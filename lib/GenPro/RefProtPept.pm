use 5.10.0;
use strict;
use warnings;

package GenPro::RefProtPept;

use MCE::Loop;
use DDP;

use Mouse 2;
use namespace::autoclean;

use GenPro::DBManager;

use Seq::Tracks;
use Seq::DBManager;
use Seq::Tracks::Gene::Site;
use Seq::Tracks::Gene::Site::CodonMap;
use Test::More;

# TODO: Could inject the Seq::Tracks object
# Trying to go for loose coupling, but this may introduce hidden deps

# Private variables
my %metaKeys = (
  refTxsWritten => 'refTxsWritten',
  refPeptidesWritten => 'refPeptidesWritten',
);

my $metaEnv = 'genpro_ref_prot_peptide_meta';
my $metaConfig = {stringKeys => 1};

# Required variables, and exports
has chromosomes => (is => 'ro', isa => 'ArrayRef', required => 1);
has digest => (is => 'ro', isa => 'GenPro::Digest', required => 1);
has geneTrack => ( is => 'ro', isa => 'Str', required => 1 );

has refProtEnv => (is => 'ro', isa => 'Str', init_arg => undef, default => 'referenceProteins');
has refPeptEnv => (is => 'ro', isa => 'Str', init_arg => undef, default => 'referencePeptides');

has refProtDbConfig => (is => 'ro', isa => 'HashRef', init_arg => undef, lazy => 1,
default => sub {
  my $self = shift;

  return {
    maxDbs => scalar @{$self->chromosomes}
  }
});

has refPeptDbConfig => (is => 'ro', isa => 'HashRef', init_arg => undef, lazy => 1, default => sub {
  my $self = shift;

  return {
    maxDbs => scalar @{$self->digest->enzymes},
    stringKeys => 1,
  }
});

# Optional
has maxThreads => (is => 'ro', isa => 'Int', default => sub { Sys::CpuAffinity::getNumCpus() });

sub go {
  my $self = shift;
  my $wantedTxNumHref = shift;

  if(!ref $wantedTxNumHref || ref $wantedTxNumHref ne 'HASH') {
    return "RefProtPept->go needs a hash reference to ( chr => txNums <HashRef> )";
  }

  my $err = $self->makeReferenceProtDb($wantedTxNumHref);

  if($err) {
    return $err;
  }

  $err = $self->makeReferenceUniquePeptides($wantedTxNumHref);

  return $err;
}

# Generates everything that GenPro_make_refprotdb does
# Also pre-calculates all digested peptides, and stores those
sub makeReferenceProtDb {
  my ($self, $wantedTxNumHref) = @_;

  my (undef, $writtenTxNumsHref, $txNumsAref)
    = $self->_checkHasTxNumToWrite($wantedTxNumHref, $metaEnv, $metaConfig, $metaKeys{refTxsWritten});

  if(@$txNumsAref == 0) {
    return;
  }

  my $geneTrack = $self->geneTrack;

  my $tracks = Seq::Tracks->new();
  my $geneTrackGetter = $tracks->getTrackGetterByName($geneTrack);
  my $geneTrackGetterDbName = $geneTrackGetter->dbName;

  # Ensure that we crate the necessary databases before entering threads
  # It *seems* I don't need to pre-make databases
  # TODO: write tests

  # Pre-load all transcript info
  my %regionDb;
  {
    # Control scope (Rust)
    my $db = Seq::DBManager->new();

    for my $chr ( keys %$wantedTxNumHref ) {
      $regionDb{$chr} //= $db->dbReadAll( $geneTrackGetter->regionTrackPath($chr) );
    }
  }

  my $siteUnpacker = Seq::Tracks::Gene::Site->new();
  my $siteTypeMap  = Seq::Tracks::Gene::Site::SiteTypeMap->new();
  my $codonMap     = Seq::Tracks::Gene::Site::CodonMap->new();

  my $strandSiteIdx        = $siteUnpacker->strandIdx;
  my $siteTypeSiteIdx      = $siteUnpacker->siteTypeIdx;
  my $codonSequenceSiteIdx = $siteUnpacker->codonSequenceIdx;
  my $codonPositionSiteIdx = $siteUnpacker->codonPositionIdx;
  my $codonNumberSiteIdx   = $siteUnpacker->codonNumberIdx;

  # These are always relative to the region database
  # The computed features are totally separate
  my $nameFeatIdx = $geneTrackGetter->getFieldDbName('name');
  my $name2FeatIdx = $geneTrackGetter->getFieldDbName('name2');
  my $strandFeatIdx = $geneTrackGetter->getFieldDbName('strand');
  my $txErrorFeatIdx = $geneTrackGetter->getFieldDbName('txError');

  my $exonStartsFeatIdx = $geneTrackGetter->getFieldDbName('exonStarts');
  my $exonEndsFeatIdx = $geneTrackGetter->getFieldDbName('exonEnds');

  my $cdsStartFeatIdx = $geneTrackGetter->getFieldDbName('cdsStart');
  my $cdsEndFeatIdx = $geneTrackGetter->getFieldDbName('cdsEnd');

  MCE::Loop::init {
    max_workers => $self->maxThreads || 8,

    # bystro-vcf outputs a very small row; fully annotated through the alt column (-ref -discordant)
    # so accumulate less than we would if processing full .snp
    chunk_size => 'auto',
    gather     => $self->_makeRefProgressFunc($writtenTxNumsHref)
  };

  mce_loop {
    my ($mce, $chunk_ref, $chunk_id) = @_;

    my $db = Seq::DBManager->new();
    my $personalDb = GenPro::DBManager->new();

    my %cursors;
    my %seenTxNums;
    my %dbs;
    for my $txNumInfo (@{ $chunk_ref }) {
      my $chr = $txNumInfo->[0];
      my $txNumber = $txNumInfo->[1];

      # TODO: refactor into function so that same options always passed
      # or at least hash of options for this environment
      # Must pass nChrs, or get -30791 MDB_DBS_FULL (no more dbs)
      # TODO: add more error checks, such that if we don't have 1 db per chr
      # that we fail...
      $dbs{$chr} //= $personalDb->_getDbi( $self->refProtEnv, $chr, $self->refProtDbConfig );

      # To save
      my $tx = $regionDb{$chr}->[$txNumber];
      my $txName = $tx->{$nameFeatIdx};
      my $strand = $tx->{$strandFeatIdx};

      if(!defined $strand) {
        die "No strand, for tx #$txNumber";
      }

      my ($cdsStart, $cdsEnd, @data);

      # say "CDS START: $cdsStart";
      # Should also be half-open, just like exonEnds
      # http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/
      # On the negative strand, $cdsStart = $cdsEnd
      # But weirdly, cdsEnd I believe is still always half-open
      $cdsStart = $tx->{$cdsStartFeatIdx};
      $cdsEnd = $tx->{$cdsEndFeatIdx};
      @data = ($cdsStart .. $cdsEnd - 1);

      if($strand eq '-') {
        @data = reverse @data;
      }

      $cursors{$chr} //= $db->dbStartCursorTxn($chr);
      $db->dbReadCursorUnsafe($cursors{$chr}, \@data);

      my ( $fetchedTxNumbers, $siteData, $wantedSiteData );
      my (@codonSequence, @aaSequence, @codonPos, @codonNums, @txNumOut);

      my $lastCodonNumber;
      for my $site (@data) {
        if(!defined) {
          die 'Err: no gene site';
        }

        ( $fetchedTxNumbers, $siteData ) = $siteUnpacker->unpack( $site->[ $geneTrackGetterDbName ] );

        if(!ref $fetchedTxNumbers) {
          if ($fetchedTxNumbers != $txNumber) {
            die "Err, expected all positions cdsStart ... cdsEnd - 1 to have an entry for tx $txNumber";
          }

          $wantedSiteData = $siteData;
        } else {
          my $i = 0;
          my $found;
          for my $num (@$fetchedTxNumbers) {
            if($num == $txNumber) {
              $found = 1;
              last;
            }

            $i++;
          }

          if(!$found) {
            die "Err: Couldn't find txNumber $txNumber at Bystro site";
          }

          $wantedSiteData = $siteData->[$i];
        }

        my $codonNumber = $wantedSiteData->[$codonNumberSiteIdx];

        if(!defined $codonNumber || (defined $lastCodonNumber && $lastCodonNumber == $codonNumber)) {
          next;
        }

        $lastCodonNumber = $codonNumber;

        my $codonSequence = $wantedSiteData->[$codonSequenceSiteIdx];

        my $aa = $codonMap->codon2aa( $codonSequence );

        if(!defined $aa) {
          p $codonSequence;
          p $aa;
          p $siteData;
          die "Couldn't make amino acid for $txNumber";
        }

        # my $codonPosition = $wantedSiteData->[$codonPositionSiteIdx];

        push @codonSequence, $codonSequence;
        push @aaSequence, $aa;

        # Thomas did this...figure out exactly why
        # But apparently on occassion there are multiple stop codons...
        # Why is this?
        # What import do they have?
        # Related to NMD
        # TODO: Write tests for such proteins
        # TODO: Maybe report errors if this is not at the end of the transcript
        # Convo w/ Thomas @ 11/13/18: Return the fragment, warn the user
        if($aa eq '*') {
          last;
        }
      }

      if(!@aaSequence) {
        p $tx;
        p @codonSequence;
        die "Err: ouldn't make for tx $txNumber";
      }

      if(@aaSequence == 1) {
        p $tx;
        # p $userSubstitutions;
        p $txNumber;
        p @aaSequence;
        p @codonSequence;
        die "Err";
      }

      # TODO: Move to testing package
      if($regionDb{$chr}[$txNumber]{$nameFeatIdx} eq 'NM_033467') {
        # https://www.ncbi.nlm.nih.gov/nuccore/NM_033467
        ok(join('', @aaSequence) eq 'MGKSEGPVGMVESAGRAGQKRPGFLEGGLLLLLLLVTAALVALGVLYADRRGKQLPRLASRLCFLQEERTFVKRKPRGIPEAQEVSEVCTTPGCVIAAARILQNMDPTTEPCDDFYQFACGGWLRRHVIPETNSRYSIFDVLRDELEVILKAVLENSTAKDRPAVEKARTLYRSCMNQSVIEKRGSQPLLDILEVVGGWPVAMDRWNETVGLEWELERQLALMNSQFNRRVLIDLFIWNDDQNSSRHIIYIDQPTLGMPSREYYFNGGSNRKVREAYLQFMVSVATLLREDANLPRDSCLVQEDMVQVLELETQLAKATVPQEERHDVIALYHRMGLEELQSQFGLKGFNWTLFIQTVLSSVKIKLLPDEEVVVYGIPYLQNLENIIDTYSARTIQNYLVWRLVLDRIGSLSQRFKDTRVNYRKALFGTMVEEVRWRECVGYVNSNMENAVGSLYVREAFPGDSKSMVRELIDKVRTVFVETLDELGWMDEESKKKAQEKAMSIREQIGHPDYILEEMNRRLDEEYSNLNFSEDLYFENSLQNLKVGAQRSLRKLREKVDPNLWIIGAAVVNAFYSPNRNQIVFPAGILQPPFFSKEQPQALNFGGIGMVIGHEITHGFDDNGRNFDKNGNMMDWWSNFSTQHFREQSECMIYQYGNYSWDLADEQNVNGFNTLGENIADNGGVRQAYKAYLKWMAEGGKDQQLPGLDLTHEQLFFINYAQVWCGSYRPEFAIQSIKTDVHSPLKYRVLGSLQNLAAFADTFHCARGTPMHPKERCRVW*');
      }

      # store the primary key as well as denormalized info
      # since we won't be using this db as a traditional, relational db
      # although we could just ignore the de-normalized data
      # my $txn = $dbs{$chr}{env}->BeginTxn();
      $personalDb->dbPut($dbs{$chr}, $txNumber, [\@aaSequence, \@codonSequence, [$txName, $cdsStart, $cdsEnd]]);

      $seenTxNums{$chr} //= {};
      $seenTxNums{$chr}{$txNumber} = 1;
    }

    MCE->gather(\%seenTxNums);

    # Doing this will lead to EINVAL during creation...race condition?
    # $personalDb->cleanUp();
  } @$txNumsAref;

  MCE::Loop::finish();

  $self->_recordTxNumsWritten($writtenTxNumsHref, $metaEnv, $metaConfig, $metaKeys{refTxsWritten});
}

sub _recordTxNumsWritten {
  my ($self, $writtenTxNumsHref, $metaEnv, $metaConfig, $metaKey) = @_;

  my $personalDb = GenPro::DBManager->new();

  my $metaDb = $personalDb->_getDbi( $metaEnv, undef, $metaConfig );
    $personalDb->dbPut( $metaDb, $metaKey, $writtenTxNumsHref );
  undef $metaDb;

  $personalDb->cleanUp();
}

sub _checkHasTxNumToWrite {
  my ($self, $wantedTxNumHref, $metaEnv, $metaConfig, $metaKey) = @_;

  my $personalDb = GenPro::DBManager->new();

  my $metaDb = $personalDb->_getDbi( $metaEnv, undef,$metaConfig);
    my $previouslyWritten = $personalDb->dbReadOne( $metaDb, $metaKey );
  undef $metaDb;

  my %writtenChrs = $previouslyWritten ? %$previouslyWritten : ();

  my @txNums;
  my %wantedChrs;
  for my $chr (sort { $a cmp $b } keys %$wantedTxNumHref) {
    for my $txNum (sort { $a <=> $b } keys %{$wantedTxNumHref->{$chr}}) {
      if(!($writtenChrs{$chr} && $writtenChrs{$chr}{$txNum})) {
        push @txNums, [$chr, $txNum];
      }
    }

    $wantedChrs{$chr} //= 1;
  }

  # Clear the singleton instance, ensure that threads don't copy any memory
  $personalDb->cleanUp();
  undef $personalDb;

  return (\%wantedChrs, \%writtenChrs, \@txNums);
}

sub _makeRefProgressFunc {
  my ($self, $writtenTxNumsHref) = @_;

  return sub {
    my $seenTxNums = shift;

    for my $chr (keys %$seenTxNums) {
      $writtenTxNumsHref->{$chr} //= {};

      for my $txNum (keys %{$seenTxNums->{$chr}}) {
        $writtenTxNumsHref->{$chr}{$txNum} = 1;
      }
    }
  };
}

# Generates everything that GenPro_make_refprotdb does
# Also pre-calculates all digested peptides, and stores those
sub makeReferenceUniquePeptides {
  my ($self, $wantedTxNumHref) = @_;

  my $wantedEnzymesAref = $self->digest->enzymes;

  my %digestFuncs;
  for my $enzyme (@$wantedEnzymesAref) {
    $digestFuncs{$enzyme} = $self->digest->makeDigestFunc($enzyme);
  }

  # for now assume we'll have trypsin, chymotrypsin, lyse-c
  my $nEnzymeTables = 3;

  my (undef, $writtenTxNumsHref, $txNumsAref) =
    $self->_checkHasTxNumToWrite($wantedTxNumHref, $metaEnv, $metaConfig, $metaKeys{refPeptidesWritten});

  # TODO: size the database based on the number of possible enzymes
  # TODO: think about combining refProtEnv and refPeptEnv
  # It doesn't seem like this is necessary
  # https://www.openldap.org/lists/openldap-technical/201608/msg00067.html
  # Not 100% sure that we don't have a race condition during database creation
  # for my $enzyme (@$wantedEnzymesAref) {
  #   $personalDb->_getDbi( $self->refPeptEnv, $enzyme, $self->refPeptDbConfig );
  # }

  if(@$txNumsAref == 0) {
    return;
  }

  my %refProtRdOnlyConfig = (%{$self->refProtDbConfig}, (readOnly => 1));

  MCE::Loop::init {
    max_workers => $self->maxThreads || 8,

    # bystro-vcf outputs a very small row; fully annotated through the alt column (-ref -discordant)
    # so accumulate less than we would if processing full .snp
    chunk_size => 'auto',
    gather     => $self->_makeRefProgressFunc($writtenTxNumsHref),
  };

  mce_loop {
    my ($mce, $chunk_ref, $chunk_id) = @_;

    # Thread-local instance
    my $personalDb = GenPro::DBManager->new();

    my %dbs;
    my %refProtDbs;
    my $existing;
    my @peptides;
    my $err;
    my $data;
    my $dbi;
    my %seenTxNums;
    for my $txNumInfo (@{ $chunk_ref }) {
      my $chr = $txNumInfo->[0];
      my $txNumber = $txNumInfo->[1];

      if($seenTxNums{$chr}) {
        $seenTxNums{$chr}{$txNumber} //= 1;
      } else {
        $seenTxNums{$chr} = {$txNumber => 1};
        # TODO: refactor this into separate function so that same options passed each time
        $refProtDbs{$chr} //= $personalDb->_getDbi( $self->refProtEnv, $chr, \%refProtRdOnlyConfig );
      }

      my $seqInfo = $personalDb->dbReadOne($refProtDbs{$chr}, $txNumber);

      my $aaSeqAref = $seqInfo->[0];
      my $txInfo = $seqInfo->[2];
      for my $enzyme (@$wantedEnzymesAref) {
        $dbs{$enzyme} //= $personalDb->_getDbi( $self->refPeptEnv, $enzyme, $self->refPeptDbConfig );
        $dbi = $dbs{$enzyme}{dbi};

        for my $peptide ($digestFuncs{$enzyme}->($aaSeqAref)) {
          # say STDERR $peptide . "\n";
          my $txn = $dbs{$enzyme}{env}->BeginTxn();

          # TODO: Check if one transcript produces the same peptide multiple times??
          # We will see it here...but will need to do linear search O(n) to check
          $data = $personalDb->dbReadOneRaw($txn, $dbi, $peptide);

          if(!defined $data) {
            $personalDb->dbPutRaw($txn, $dbi, $peptide, [$txInfo]);
          } else {
            push @{$data}, $txInfo;

            $personalDb->dbPutRaw($txn, $dbi, $peptide, $data);
          }

          $err = $txn->commit();

          if($err) {
            die "Error on commit: $err";
          }

          # if ($LMDB_File::last_err) {
          #   if ( $LMDB_File::last_err != MDB_KEYEXIST ) {

          #     # TODO: return error from thread, synchronize with other threads $$abortErr
          #     die "dbPatch put or commit LMDB error $LMDB_File::last_err";
          #   }

          #   #reset the class error variable, to avoid crazy error reporting later
          #   $LMDB_File::last_err = 0;
          # }
        }
      }
    }

    MCE->gather(\%seenTxNums);

    # Our cleanup function seems to cause issues
    # undef %refProtDbs;
    # undef %dbs;
    # $personalDb->cleanUp()
  } @$txNumsAref;

  MCE::Loop::finish();

  $self->_recordTxNumsWritten($writtenTxNumsHref, $metaEnv, $metaConfig, $metaKeys{refPeptidesWritten});
}

__PACKAGE__->meta->make_immutable;
1;