# TODO: Put each database in a separate environemnet, aka folder
# TODO: Check that we don't want stopGain, stopLoss, startGain, or spliceD/A

# or at least put the meta information as a subfolder, say refProteins/meta

use 5.10.0;
use strict;
use warnings;

package GenProSeq;

our $VERSION = '0.001';

use Mouse 2;
use lib './lib';

use List::Util qw/sum/;
use Types::Path::Tiny qw/AbsFile/;

use namespace::autoclean;

use MCE::Loop;

use Seq::InputFile;
use Seq::Output;
use Seq::Output::Delimiters;
use Seq::Tracks::Build::CompletionMeta;
use Seq::Headers;


use GenPro::DBManager;
use GenPro::Digest;
use GenPro::SampleBuilder;

use Seq::Tracks::Base::MapFieldNames;
use Seq::Tracks::Gene::Site;
use Seq::Tracks::Gene::Site::CodonMap;
use Seq::DBManager;
use Path::Tiny;

use DDP;
use Test::More;

extends 'Seq::Base';

# Bystro db is read-only
has '+readOnly' => (init_arg => undef, default => 1);

has fileProcessors => ( is => 'ro', isa => 'HashRef', required => 1);

has chromosomes => ( is => 'ro', isa => 'ArrayRef', required => 1);

has config => (is => 'ro', required => 1, isa => 'Str');

# Defines most of the properties that can be configured at run time
# Needed because there are variations of Seq.pm, ilke SeqFromQuery.pm
# Requires logPath to be provided (currently found in Seq::Base)
with 'Seq::Definition', 'Seq::Role::Validator';

has input_file => (is => 'ro' );

my $metaEnv = 'stats';
my %metaKeys = (
  samplesWritten => 'samplesWritten',
  wantedTxs => 'wantedTxs',
  refTxsWritten => 'refTxsWritten',
  refPeptidesWritten => 'refPeptidesWritten',
  samplePeptidesWritten => 'samplePeptidesWritten',
);

my $metaConfig = {stringKeys => 1};
my $refPeptideConfig = {stringKeys => 1, maxDbs => 3} ;

my $refProtEnv = 'referenceProteins';
my $refPeptideEnv = 'referencePeptides';

my $minPeptideLength = 6;
my $maxPeptideLength = 40;

my $trim_end_regex     = qr{\*[\*\w]*\z};

my @allDigestFuncs = ('trypsin');

has refProtConfig => (is => 'ro', isa => 'HashRef', init_arg => undef, lazy => 1,
default => sub {
  my $self = shift;

  return {
    maxDbs => scalar @{$self->chromosomes}
  }
});

# TODO: Move sample stuff into separate package
# my $sampleFieldNamer = Seq::Tracks::Base::MapFieldNames->new({
#   name => 'sampleVariant',
# });


# my $mp = Data::MessagePack->new()->prefer_integer()->prefer_float32();

# TODO: further reduce complexity
sub BUILD {
  my $self = shift;

  ################## Make the full output path ######################
  # The output path always respects the $self->output_file_base attribute path;
  $self->{_outPath} =
  $self->_workingDir->child( $self->outputFilesInfo->{annotation} );

  # Must come before statistics, which relies on a configured Seq::Tracks
  #Expects DBManager to have been given a database_dir
  $self->{_db} = Seq::DBManager->new();

  # Initializes the tracks
  # This ensures that the tracks are configured, and headers set
  $self->{_tracks} = $self->tracksObj;

  # TODO: Convert genpro to Bystro track, give it its own db folder
  GenPro::DBManager::initialize(
    {
      databaseDir => path( $self->database_dir )->parent()->child("genpro")->stringify()
    }
  );
}

sub annotate {
  my $self = shift;

  # TODO: Handle Sig Int (Ctrl + C) to close db, clean up temp dir
  # local $SIG{INT} = sub {
  #   my $message = shift;
  # };

  my ($err, $sampleList, $wantedTxNums);

  # TODO: Allow configuration of track, features
  # TODO: extract config needed by makeSampleUniquePeptides?
  my $sampleBuilder = GenPro::SampleBuilder->new({
    inputFile => $self->input_file,
    geneTrack => 'refSeq',
    fileProcessors => $self->fileProcessors,
    # TODO: this is clunky (from Seq::Definition)
    output_file_base => $self->output_file_base,
    temp_dir => $self->temp_dir,
    chromosomes => $self->chromosomes,
    # assembly => $self->assembly,
  });

  ( $err, $sampleList, $wantedTxNums ) = $sampleBuilder->go();

  say STDERR "FINISHED STEP 1 (make personal replacement db)";

  ( $err ) = $self->makeReferenceProtDb($wantedTxNums);

  say STDERR "FINISHED STEP 2 (make reference protein db for requested txNumbers)";

  ( $err ) = $self->makeReferenceUniquePeptides($sampleList, $wantedTxNums);
  say STDERR "FINISHED STEP 3 (make reference peptide database for trypsin for requested txNumbers)";
  
  ( $err ) = $self->makeSampleUniquePeptides($sampleList, $wantedTxNums, {
    sample => {
      dbConfig => $sampleBuilder->dbConfig,
      featureDbIdx => $sampleBuilder->featureDbIdx
    },
    refProtein => {
      dbConfig => $self->refProtConfig,
      env => $refProtEnv,
    },
    refPeptide => {
      dbConfig => $refPeptideConfig,
      env => $refPeptideEnv,
    }
  });
  say STDERR "FINISHED STEP 4 (make per-sample database of unique peptides digested by trypsin)";
  # # We may get rid of this step
  # # ( $err ) = $self->createPersProtPermutations($sampleList, $wantedTxNums);

  # TODO: Inspect vcf header

  # TODO: support any other file, by checking the extension

  # TODO: we don't really check for valid vcf, just assume it is
  # So this message is never reached
  # $self->_errorWithCleanup("File type isn\'t vcf or snp. Please use one of these files");
  # return ( "File type isn\'t vcf or snp. Please use one of these files", undef );
}

# Generates everything that GenPro_make_refprotdb does
# Also pre-calculates all digested peptides, and stores those
sub makeReferenceProtDb {
  my ($self, $wantedTxNumHref) = @_;

  my $personalDb = GenPro::DBManager->new();
  
  my $metaDb = $personalDb->_getDbi( $metaEnv, undef,$metaConfig);
    my $previouslyWritten = $personalDb->dbReadOne( $metaDb, $metaKeys{refTxsWritten} );
  undef $metaDb;

  my %writtenChrs = $previouslyWritten ? %$previouslyWritten : ();

  my @txNums;
  my %wantedChrs;
  for my $chr (sort { $a cmp $b } keys %$wantedTxNumHref) {
    for my $txNum (sort { $a <=> $b } keys %{$wantedTxNumHref->{$chr}}) {
      push @txNums, [$chr, $txNum];
    }

    $wantedChrs{$chr} //= 1;
  }

  # TODO: get the name string in some other way, maybe from YAML
  # or as an argument
  my $db              = $self->{_db};
  my $geneTrackGetter = $self->tracksObj->getTrackGetterByName('refSeq');
  my $geneTrackGetterDbName = $geneTrackGetter->dbName;

  my %regionDb;
  my $nChrs = %wantedChrs;
  my $haveWritten = 0;

  # Ensure that we crate the necessary databases before entering threads
  # It *seems* I don't need to pre-make databases
  # TODO: write tests
  for my $chr ( keys %wantedChrs ) {
    $regionDb{$chr} //= $db->dbReadAll( $geneTrackGetter->regionTrackPath($chr) );

    # Run here, in case we have yet to create this environment, before
    # we get to multi-threaded code.
    # Doesn't seem to be necessary
    # $personalDb->_getDbi( $refProtEnv, 0, 0, $chr, $nChrs );

    if(defined $writtenChrs{$chr}) {
      $haveWritten++;
    }
  }

  # Clear the singleton instance, ensure that threads don't copy any memory
  $personalDb->cleanUp();
  undef $personalDb;

  if($haveWritten == keys %wantedChrs) {
    return;
  }

  my $progressFunc = sub {
    my $seenChrs = shift;

    for my $chr (keys %$seenChrs) {
      $writtenChrs{$chr} //= 1;
    }
  };

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
    gather     => $progressFunc,
  };

  mce_loop {
    my ($mce, $chunk_ref, $chunk_id) = @_;

    $personalDb = GenPro::DBManager->new();
  
    my %seenCodonRanges;
    my %cursors;
    my %seenChrs;
    my %dbs;
    for my $txNumInfo (@{ $chunk_ref }) {
      my $chr = $txNumInfo->[0];
      my $txNumber = $txNumInfo->[1];

      $seenChrs{$chr} //= 1;

      # TODO: refactor into function so that same options always passed
      # or at least hash of options for this environment
      # Must pass nChrs, or get -30791 MDB_DBS_FULL (no more dbs)
      # TODO: add more error checks, such that if we don't have 1 db per chr
      # that we fail...
      $dbs{$chr} //= $personalDb->_getDbi( $refProtEnv, $chr, $self->refProtConfig );

      # To save
      my $tx = $regionDb{$chr}->[$txNumber];
      my $txName = $tx->{$nameFeatIdx};
      my $strand = $tx->{$strandFeatIdx};
      # p $tx;

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

      $seenCodonRanges{$chr} //= {};
      $seenCodonRanges{$chr}{$cdsStart} //= {};

      # Don't do this!!!
      # if($seenCodonRanges{$chr}{$cdsStart}{$cdsEnd}) {
      #   next;
      # }

      $seenCodonRanges{$chr}{$cdsStart}{$cdsEnd} = 1;

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
      
      # if($txNumber == 438) {
      #   say STDERR "Put in txNumber 438";

      #   my $val =  $personalDb->dbReadOne($dbs{$chr}, $txNumber);

      #   p $val;
      #   sleep(5);
      # }
      
      # $txn->commit();
    }

    MCE->gather(\%seenChrs);

    # Remove our pointers so LMDB can clean up properly
    # undef %dbs;
    # $personalDb->cleanUp();
  } @txNums;

  MCE::Loop::finish();

  $personalDb = GenPro::DBManager->new();

  $metaDb = $personalDb->_getDbi( $metaEnv, undef, $metaConfig );
    $personalDb->dbPut( $metaDb, $metaKeys{refTxsWritten}, \%writtenChrs );
  undef $metaDb;

  $personalDb->cleanUp();
  $db->cleanUp();
}

# Generates everything that GenPro_make_refprotdb does
# Also pre-calculates all digested peptides, and stores those
sub makeReferenceUniquePeptides {
  my ($self, $wantedSamples, $wantedTxNumHref, $wantedEnzymesAref) = @_;

  $wantedEnzymesAref //= ['trypsin'];

  my %digestFuncs;
  my $digest = GenPro::Digest->new();

  for my $enzyme (@$wantedEnzymesAref) {
    $digestFuncs{$enzyme} = $digest->makeDigestFunc($enzyme);
  }

  # for now assume we'll have trypsin, chymotrypsin, lyse-c
  my $nEnzymeTables = 3;

  my $personalDb = GenPro::DBManager->new();

  # TODO: size the database based on the number of possible enzymes
  # TODO: think about combining refProtEnv and refPeptideEnv
  for my $enzyme (@$wantedEnzymesAref) {
    $personalDb->_getDbi( $refPeptideEnv, $enzyme, $refPeptideConfig );
  }

  my $metaDb = $personalDb->_getDbi( $metaEnv, undef, $metaConfig );
    my $previouslyWritten = $personalDb->dbReadOne( $metaDb, $metaKeys{refPeptidesWritten} );
  undef $metaDb;

  # Clear the singleton instance, ensure that threads don't copy any memory
  $personalDb->cleanUp();
  undef $personalDb;

  my %writtenTxNums = $previouslyWritten ? %$previouslyWritten : ();

  my @txNums;
  my %wantedChrs;
  for my $chr (sort { $a cmp $b } keys %$wantedTxNumHref) {
    $wantedChrs{$chr} //= 1;

    for my $txNum (sort { $a <=> $b } keys %{$wantedTxNumHref->{$chr}}) {
      if($writtenTxNums{$chr} && $writtenTxNums{$chr}{$txNum}) {
        next;
      }

      push @txNums, [$chr, $txNum];
    }
  }

  if(@txNums == 0) {
    return;
  }

  # TODO: when working with named databases, must open each time with same flags
  # and I think, same number of maxDbs (or at least, no fewer max dbs?)
  my $nChrs = %wantedChrs;
  
  my %refProtRdOnlyConfig = (%{$self->refProtConfig}, (readOnly => 1));
  my $progressFunc = sub {
    my $seenTxNums = shift;

    for my $chr (keys %$seenTxNums) {
      $writtenTxNums{$chr} //= {};

      for my $txNum (keys %{$seenTxNums->{$chr}}) {
        $writtenTxNums{$chr}{$txNum} = 1;
      }
    }
  };

  MCE::Loop::init {
    max_workers => $self->maxThreads || 8,

    # bystro-vcf outputs a very small row; fully annotated through the alt column (-ref -discordant)
    # so accumulate less than we would if processing full .snp
    chunk_size => 'auto',
    gather     => $progressFunc,
  };

  mce_loop {
    my ($mce, $chunk_ref, $chunk_id) = @_;

    # Thread-local instance
    $personalDb = GenPro::DBManager->new();

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
        $refProtDbs{$chr} //= $personalDb->_getDbi( $refProtEnv, $chr, \%refProtRdOnlyConfig );
      }

      my $seqInfo = $personalDb->dbReadOne($refProtDbs{$chr}, $txNumber);

      my $aaSeqAref = $seqInfo->[0];
      my $txInfo = $seqInfo->[2];
      for my $enzyme (@$wantedEnzymesAref) {
        $dbs{$enzyme} //= $personalDb->_getDbi( $refPeptideEnv, $enzyme, $refPeptideConfig );
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
          # p $err;

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
  } @txNums;

  MCE::Loop::finish();

  $personalDb = GenPro::DBManager->new();

  $metaDb = $personalDb->_getDbi($metaEnv, undef, $metaConfig);
    $personalDb->dbPut( $metaDb, $metaKeys{refPeptidesWritten}, \%writtenTxNums );
  undef $metaDb;

  $personalDb->cleanUp();
}

# TODO: check if number wantedTxNumHref has changed
# and if so, don't go by the sample meta
# TODO: pass in 
sub makeSampleUniquePeptides {
  my ($self, $wantedSamplesHref, $wantedTxNumHref, $configs) = @_;

  my $wantedEnzymesAref = ['trypsin'];

  my $sampleConfig = $configs->{sample};
  my $sampleFeatureDbIdx = $sampleConfig->{featureDbIdx};
  my %sampleReadOnlyConfig = (%{$sampleConfig->{dbConfig}}, (readOnly => 1));

  my $refPeptideConfig = $configs->{refPeptide};

  my $refProteinConfig = $configs->{refProtein};
  my %refProtRdOnlyConfig = (%{$refProteinConfig->{dbConfig}}, (readOnly => 1));
  my $refProtEnv = $refProteinConfig->{env};

  # When working with named databases, need to pass same configuration each time
  # This is to know the maxDbs config for the reference peptide database
  # which stores n peptides ('for now just trypsin');
  # and maybe needs to be at least no less than the largest number it was last
  # opened with (or maybe needs to be constant)
  # TODO: Make this a configuration owned by the reference peptide database
  # my $nEnzymeTables = 3;

  my $personalDb = GenPro::DBManager->new();

  my $metaDb = $personalDb->_getDbi( $metaEnv, undef, $metaConfig );
    my $previouslyWritten = $personalDb->dbReadOne( $metaDb, $metaKeys{samplePeptidesWritten} );
  undef $metaDb;

  # Clear the singleton instance, ensure that threads don't copy any memory
  $personalDb->cleanUp();
  undef $personalDb;

  my %writtenSamples = $previouslyWritten ? %$previouslyWritten : ();

  my @toGetSamples;
  for my $sample (keys %$wantedSamplesHref) {
    if(defined $writtenSamples{$sample}) {
      next;
    }

    push @toGetSamples, $sample;
  }

  if(@toGetSamples == 0) {
    return;
  }

  my $txNumberIdx = $sampleFeatureDbIdx->{txNumber};

  my $progressFunc = sub {
    my $seenSamples = shift;

    for my $sample (keys %$seenSamples) {
      for my $chr (keys %{$seenSamples->{$sample}}) {
        $writtenSamples{$sample}{$chr} = 1;
      }
    }
  };

  MCE::Loop::init {
    max_workers => $self->maxThreads || 8,

    # bystro-vcf outputs a very small row; fully annotated through the alt column (-ref -discordant)
    # so accumulate less than we would if processing full .snp
    chunk_size => 'auto',
    gather     => $progressFunc,
  };

  my @wantedChrs = keys %$wantedTxNumHref;
  mce_loop {
    my ($mce, $chunk_ref, $chunk_id) = @_;

    # Thread-local instance
    $personalDb = GenPro::DBManager->new();

  
    my %refProteinDbs;
    
    my ($txHasUniquePeptidesFn, $cleanUp) = $self->_configureVarProtFn($sampleFeatureDbIdx, $refPeptideConfig);

    for my $sample (@{ $chunk_ref }) {
      my @finalRecords;

      my $samplePeptideDb = $personalDb->_getDbi( "$sample/peptide", 'trypsin', {maxDbs => 3, stringKeys => 1} );
      
      for my $chr (@wantedChrs) {
        $refProteinDbs{$chr} //= $personalDb->_getDbi( $refProtEnv, $chr, \%refProtRdOnlyConfig );

        my $sampleDb = $personalDb->_getDbi( $sample, $chr, \%sampleReadOnlyConfig );

        my $variantsAref = $personalDb->dbReadAll($sampleDb);

        if(!defined $variantsAref) {
          say STDERR "Found no txs for $sample on chr $chr";
          next;
        }
        
        # Place all data into
        # txName => [first modification, 2nd modification, etc]
        # Where each modification = contains a reference to the record
        # from the 
        my $userTxHref = _orderByDesired($variantsAref, $txNumberIdx);

        for my $txNum (sort {$a <=> $b} keys %$userTxHref) {
          my $userVars = $userTxHref->{$txNum};

          # p $userVars;
          if(@$userVars > 21) {
            say STDERR "Tx number $txNum on chr $chr has too many variants ( " . @$userVars . ')';
            next;
          } 

          my $refSequenceInfo = $personalDb->dbReadOne($refProteinDbs{$chr}, $txNum);

          if(!$refSequenceInfo) {
            die "Couldn't find transcript number $txNum for chr $chr";
          }
  
          my $uniquePeptidesAref = $txHasUniquePeptidesFn->($userVars, $refSequenceInfo->[0]);
          # say STDERR "UNIQUE PEPTIDES";
          # p $uniquePeptidesAref;
          # sleep(1);
          # exit;
          if(!@$uniquePeptidesAref) {
            # say STDERR "Couldn't find any peptides";
            next;
          }

          # This replicates select_entries as well as create_per_prot_rec
          # Since create_per_prot_rec is no longer necessary,
          # as $txHasUniquePeptidesFn stores 
          # \@{all affected txs}, \@{alt aa sequence}, \@uniquePeptides
          # where uniquePeptides is also modified by sortByUniquePeptides 
          my @sortedRecs = sortUniquePeptides($uniquePeptidesAref, $txNum, $personalDb, $samplePeptideDb);

          # p @sortedRecs;
          # sleep(1);

          while(@sortedRecs) {
            my $rec = shift @sortedRecs;
            push @finalRecords, $rec;

            # Note: This is ok, because sortByUnique peptides
            # will remove all peptides previously seen in this sample
            for my $peptide (@{$rec->[2]}) {
              $personalDb->dbPut($samplePeptideDb, $peptide, [1, [$txNum]]);
              # say STDERR "PUTTING: " . $peptide;
            }

            # Checks all of the pre-generated peptides
            # And:
            # 1) returns only those records with at least 1 unique peptide
            # 2) for each of those records, index 2 points to an array
            # of the unique peptides pertaining to that record
            # Should replicate Thomas' code, except avoiding 2nd round of digestion
            # and 2nd round of checking against the reference db, which is invariant
            # (so this function checks only against the sample db, which is
            # updated in the line above thise)
            @sortedRecs = sortUniquePeptides(\@sortedRecs, $txNum, $personalDb, $samplePeptideDb);  
          }
          # say STDERR "Final records are";
          # p @finalRecords;
        }
      }
    }

    $cleanUp->();

  } @toGetSamples;

  MCE::Loop::finish();

  $personalDb = GenPro::DBManager->new();

  # $metaDb = $personalDb->_getDbi($metaEnv, undef, 1);
  #   $personalDb->dbPut( $metaDb, $metaKeys{refPeptidesWritten}, \%writtenSamples );
  # undef $metaDb;

  $personalDb->cleanUp();
}



# Takes a variant record, with N affected transcripts,
# de-convolutes those,
# returning per-transcript data, indexed on some desired transcript value
# ex: _orderBy($dataAref, txNumberIdx);
# @returns <HashRef> : per-desired tx value variant records;
sub _orderByDesired {
  my ($dataAref, $desiredIdx) = @_;

  my %ret;
  for my $data (@$dataAref) {
    my $desired = $data->[$desiredIdx];

    if(ref $desired) {
      # p $data;
      for (my $i = 0; $i < @$desired; $i++) {
        $ret{$desired->[$i]} //= [];

        # A few of the properties are not arrays (chr, pos, ref, alt, type)
        # and all others are
        # TODO: be more precise about this...potentially wait for Rust
        push @{$ret{$desired->[$i]}}, [ map { ref $_ ? $_->[$i] : $_ } @$data ];
      }
    } else {
      $ret{$desired} //= [];

      push @{$ret{$desired}}, $data;
    }
  }

  return \%ret;
}

# Only sorts peptides, and updates counts, txNumbers of non-unique peptides
sub sortUniquePeptides {
  my ($recsAref, $txNum, $dbManager, $db) = @_;

  my @finalRecords;
  # say STDERR "HELLLLO";
  # p $recsAref;
  # say STDERR "PAST";

  for my $rec (sort { @{$b->[2]} <=> @{$a->[2]} } @$recsAref) {
    # say "Checking ref";
    # p $rec;
    # p $rec->[2];

    # Check each peptide against 
    my @uniquePeptides;

    for my $peptide (@{$rec->[2]}) {
      my $seenAref = $dbManager->dbReadOne($db, $peptide);

      if(defined $seenAref) {
        # p $seenAref;
        $seenAref->[0] += 1;
        push @{$seenAref->[1]}, $txNum;

        $dbManager->dbPut($db, $peptide, $seenAref);

        # say STDERR "defined";
        # p $seenAref;

        next;
      }

      push @uniquePeptides, $peptide;
    }

    if(@uniquePeptides > 0) {
      $rec->[2] = \@uniquePeptides;

      push @finalRecords, $rec;
    }
  }

  # say STDERR "FINAL";
  # p @finalRecords;
  return @finalRecords;
}

sub storeUniquePeptides {
  my ($rec, $txNum, $dbManager, $db);

  for my $peptide (@${$rec->[2]}) {
    $dbManager->dbPut($db, $peptide, [1, $txNum]);
  }

  return;
}

sub sortAndStoreUniquePeptides {
  my ($recsAref, $txNum, $dbManager, $db);

  my @finalRecords;

  for my $rec (sort { @{$b->[2]} <=> @{$a->[2]} } @$recsAref) {
    # Check each peptide against 
    my @uniquePeptides;

    for my $peptide (@${$rec->[2]}) {
      my $seenAref = $dbManager->dbReadOne($db, $peptide, 1);

      
      if(!defined $seenAref) {
        push @uniquePeptides, $peptide;

        $seenAref = [ 1, [$txNum] ];
      } else {
        $seenAref->[0] += 1;
        push @{$seenAref->[1]}, $txNum;
      }

      $dbManager->dbPut($db, $peptide, $seenAref);
    }

    if(@uniquePeptides > 0) {
      $rec->[2] = \@uniquePeptides;

      push @finalRecords, $rec;
    }
  }

  return @finalRecords;
}

sub makeVarProtForTxFn {
  my $config = shift;
  
  my $altAaIdx = $config->{altAaIdx};
  my $refAaIdx = $config->{refAaIdx};
  my $codonNumIdx = $config->{codonNumIdx};
  my $cutsHref = $config->{cutsHref};
  my $globalUniquePeptideFn = $config->{globalUniquePeptideFn};
  my $maxPeptideLength = $config->{maxPeptideLength};

  if(!%$cutsHref || !$globalUniquePeptideFn || !defined $altAaIdx || !defined $refAaIdx ||
  !defined $codonNumIdx || !defined $maxPeptideLength) {
    die "makeVarProtForTxFn requires all config properties";
  }

  return sub {
    my ($varTxAref, $refAaSeq) = @_;

    # p $varTxAref;
    # p $refAaSeq;
    # exit;
    
    my (@uniqueRecords, %seqWithUnique);

    my $lastCodonNum;

    if(!$refAaSeq->[-1] eq '*') {
      say STDERR "Got a transcript without a trailing stop";
      $seqWithUnique{-9} = [ undef, [@$refAaSeq], undef ];
    } else {
      $seqWithUnique{-9} = [ undef, [ @$refAaSeq[ 0 .. $#$refAaSeq - 1 ] ], undef ];
    }
  
    my $max = -9 + sum map { $_->[$codonNumIdx] } @$varTxAref;
    # my $max = "-9" . join(" ", @$varTxAref);

    for my $var (sort { $a->[$codonNumIdx] <=> $b->[$codonNumIdx] } @$varTxAref) {
      my $codonNum = $var->[$codonNumIdx];

      # say STDERR "VAR IS for codon $codonNum";
      # p $var;
      # sleep(5);
      # TODO: CHECK ON THIS: THOMAS HAS THIS INITIALIZED TO 1, WHICH I THINK
      # MAY CAUSE SOME COMBINATIONS NOT TO BE SEEN (WHEN THE FIRST CODON)
      # IS > maxPeptideLength away
      $lastCodonNum //= $codonNum;

      # We can avoid running this condition if we 
      # Both of these are 1-based, so add 1 to be consistent with other 
      # parts of the application, or more succinctly choose > instead of >=
      # say STDERR "DIFF : " . ($codonNum - $lastCodonNum);
      # sleep(1);

      # This exits too early
      # if($codonNum - $lastCodonNum > $maxPeptideLength) {
      #   last;
      # }

      $lastCodonNum = $codonNum;

      # say STDERR "PAST . does k Cut" . ($cutsHref->{'K'} ? 'yes' : ' no');

      # p $cutsHref;
      # sleep(1);

      my $altAa = $var->[$altAaIdx];
# say STDERR "alt is $altAa";
# sleep(1);      # Save cycles
      if(!$cutsHref->{$altAa}) {
        next;
      }

      my $refAa = $var->[$refAaIdx];
# say STDERR "altAa is $altAa and it is cuttable, ref is $refAa";
# sleep(1);
      if($refAa ne $seqWithUnique{-9}[1][$codonNum - 1]) {
        die "NOT EQUAL at codon $codonNum of $seqWithUnique{-9}";
      }

      $lastCodonNum = $codonNum;
      # say STDERR "seqWithUnique has";
      # p %seqWithUnique;

      for my $i (keys %seqWithUnique) {
        # say STDERR "Checking seqWithUnique key $i";
        # Copy this mutated allele (or reference if first iteration)
        # Mutate that 
        my @newSeq = @{$seqWithUnique{$i}[1]};
      # say STDERR "NewSeq is ";
      # p @newSeq;
     
        # Thomas did this in the seq_of_per_prot, aka the seqWithUnique loop
        # but that seems unnecessary, since we should by definition
        # never have a replacement site at a stop codon
        # and our stop codon should be at the end of the reference
        # so to save time on if statements and indexing
        # we should chomp the array
        # if($newSeq[$codonNum - 1] eq '*') {
        #   next;
        # }

        $newSeq[$codonNum - 1] = $altAa;
  # p @newSeq;
        my @uniquePeptides = $globalUniquePeptideFn->(\@newSeq);
#         say "Unique peptides are";
#         p @uniquePeptides;
#  sleep(1);
        if(@uniquePeptides == 0) {
          if($i != -9) {
            # always keep the reference around
            delete $seqWithUnique{$i};
          }

          next;
        }

        my $varRecAref;
        
        if($i == -9) {
          $varRecAref = [$var];
        } else {
          $varRecAref = [@{$seqWithUnique{$i}[0]}, $var];
        }

        my $recAref = [$varRecAref, \@newSeq, \@uniquePeptides];
      # say STDERR "recAref for $i is ";
      # p $recAref;
        # We keep a reference to the unique peptides because this allows
        # future steps that need to find the unique peptides conditioned
        # on not being in the set of the previous most unique sequences' unique peptides
        # which allows us to dramatically reduce our search space in that step
        push @uniqueRecords, $recAref;
# say STDERR "now uniqueRecords has";
# p @uniqueRecords;
# sleep(1);
        # The value of the key is only import in that we want to know
        # whether we've made a record with all substitutions
        $seqWithUnique{$i + $codonNum} = $recAref;
      }
    }

    if(!$seqWithUnique{$max}) {
      # Mutates the reference sequence copy we made, but we'll delete that
      my $newSeq = $seqWithUnique{-9}[1];

      for my $var (@$varTxAref) {
        $newSeq->[ $var->[$codonNumIdx] - 1] = $var->[$altAaIdx];
      }

      my @uniquePeptides = $globalUniquePeptideFn->($newSeq);

        # say "Unique ppetides";
        # p @uniquePeptides;
        # sleep(1);
      if(@uniquePeptides) {
        push @uniqueRecords, [$varTxAref, $newSeq, \@uniquePeptides];
      }
    }


    # delete $seqWithUnique{-9};
    # say STDERR "Unique records";
    # p %seqWithUnique;
    # p @uniqueRecords;
    # sleep(1);
    return \@uniqueRecords;
  }
}

# TODO: add other enzymes
sub _configureVarProtFn {
  my $self = shift;
  # The sample variant features
  my $featureDbIdx = shift;
 
  my $refPeptideConfig = shift;

  my $altAaIdx = $featureDbIdx->{altAminoAcid};
  my $refAaIdx = $featureDbIdx->{refAminoAcid};
  my $codonNumIdx = $featureDbIdx->{codonNumber};

  my $dbManager = GenPro::DBManager->new();

  my %refPeptideRdOnlyConfig = (%{$refPeptideConfig->{dbConfig}}, (readOnly => 1));
  my $refPeptideDb = $dbManager->_getDbi( $refPeptideConfig->{env}, 'trypsin', \%refPeptideRdOnlyConfig );
  my $refPeptideDbi = $refPeptideDb->{dbi};
  my $refPeptideTxn = $refPeptideDb->{env}->BeginTxn();

  my $digest = GenPro::Digest->new();
  my $digestFunc = $digest->makeDigestFunc('trypsin');
  my $cutsHref = $digest->digestLookups->{trypsin}{cut};

  if(!$digestFunc) {
    die 'No digest func for trypsin';
  }

  my $globalUniquePeptideFn = sub {
    my $aaAref = shift;

    my @unique;
    for my $peptide ($digestFunc->($aaAref)) {
      # say STDERR "PEPTIDE: $peptide";
      if(!defined $dbManager->dbReadOneRaw($refPeptideTxn, $refPeptideDbi, $peptide)) {
        # say STDERR "NOT DEFINED $peptide";
        push @unique, $peptide;
      }
    }

    return @unique;
  };
  

  my %txHasUniqueConfig= (
    altAaIdx => $altAaIdx,
    refAaIdx => $refAaIdx,
    codonNumIdx => $codonNumIdx,
    cutsHref => $cutsHref,
    maxPeptideLength => $digest->maxPeptideLength,
    globalUniquePeptideFn => $globalUniquePeptideFn
  );

  my $txHasUniquePeptidesFn = makeVarProtForTxFn(\%txHasUniqueConfig);

  my $cleanUp = sub {
    $refPeptideTxn->abort();
  };

  return ($txHasUniquePeptidesFn, $cleanUp);
}

__PACKAGE__->meta->make_immutable;

1;
