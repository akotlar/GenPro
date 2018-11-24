# TODO: Put each database in a separate environemnet, aka folder
# or at least put the meta information as a subfolder, say refProteins/meta

use 5.10.0;
use strict;
use warnings;

package GenProSeq;

our $VERSION = '0.001';

use Mouse 2;
use Types::Path::Tiny qw/AbsFile/;

use namespace::autoclean;

use MCE::Loop;

use Seq::InputFile;
use Seq::Output;
use Seq::Output::Delimiters;
use Seq::Tracks::Build::CompletionMeta;
use Seq::Headers;

use lib './lib';
use GenPro::DBManager;
use GenPro::Digest;

use Seq::Tracks::Gene::Site;
use Seq::Tracks::Gene::Site::CodonMap;
use Seq::DBManager;
use Path::Tiny;
use Scalar::Util qw/looks_like_number/;
use DDP;
use Test::More;

use List::Util qw/sum/;

use LMDB_File qw(:all);
$LMDB_File::die_on_err = 0;

extends 'Seq::Base';

# We  add a few of our own annotation attributes
# These will be re-used in the body of the annotation processor below
# Users may configure these
has input_file => (
  is       => 'rw',
  isa      => AbsFile,
  coerce   => 1,
  required => 1,
  handles  => { inputFilePath => 'stringify' },
  writer   => 'setInputFile'
);

has chromosomes => ( is => 'ro', required => 1, isa => 'ArrayRef' );

# Maximum (signed) size of del allele
has maxDel => ( is => 'ro', isa => 'Int', default => -32, writer => 'setMaxDel' );

# TODO: formalize: check that they have name and args properties
has fileProcessors => ( is => 'ro', isa => 'HashRef', default => 'bystro-vcf' );

# Defines most of the properties that can be configured at run time
# Needed because there are variations of Seq.pm, ilke SeqFromQuery.pm
# Requires logPath to be provided (currently found in Seq::Base)
with 'Seq::Definition', 'Seq::Role::Validator';

has '+readOnly' => ( init_arg => undef, default => 1 );

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

has sampleDbConfig => (is => 'ro', isa => 'HashRef', init_arg => undef, lazy => 1,
default => sub {
  my $self = shift;

  my $nChrs = @{$self->chromosomes};

  return {
    # + 5 because we need 1 global peptide database, and want to leave some headroom,
    maxDbs => $nChrs,
  }
});

has refProtConfig => (is => 'ro', isa => 'HashRef', init_arg => undef, lazy => 1,
default => sub {
  my $self = shift;

  return $self->sampleDbConfig;
});

# my $mp = Data::MessagePack->new()->prefer_integer()->prefer_float32();

# TODO: further reduce complexity
sub BUILD {
  my $self = shift;

  if ( $self->maxDel > 0 ) {
    $self->setMaxDel( -$self->maxDel );
  }

  ################## Make the full output path ######################
  # The output path always respects the $self->output_file_base attribute path;
  $self->{_outPath} =
  $self->_workingDir->child( $self->outputFilesInfo->{annotation} );

  # Must come before statistics, which relies on a configured Seq::Tracks
  #Expects DBManager to have been given a database_dir
  $self->{_db} = Seq::DBManager->new();

  # # TODO: Convert genpro to Bystro track, give it its own db folder
  # TODO: Convert genpro to Bystro track, give it its own db folder
  GenPro::DBManager::initialize(
    {
      databaseDir => path( $self->database_dir )->parent()->child("genpro")->stringify()
    }
  );

  # Initializes the tracks
  # This ensures that the tracks are configured, and headers set
  $self->{_tracks} = $self->tracksObj;
}

sub annotate {
  my $self = shift;

  $self->log( 'info', 'Checking input file format' );

  my $err;
  ( $err, $self->{_chunkSize} ) =
  $self->getChunkSize( $self->input_file, $self->maxThreads, 512, 5000 );

  if ($err) {
    $self->_errorWithCleanup($err);
    return ( $err, undef );
  }

  $self->log( 'debug', "chunk size is $self->{_chunkSize}" );

  #################### Validate the input file ################################
  # Converts the input file if necessary
  ( $err, my $fileType ) = $self->validateInputFile( $self->input_file );

  if ($err) {
    $self->_errorWithCleanup($err);
    return ( $err, undef );
  }

  # TODO: Handle Sig Int (Ctrl + C) to close db, clean up temp dir
  # local $SIG{INT} = sub {
  #   my $message = shift;
  # };

  my ($sampleList, $wantedTxNums);

  # TODO: For efficiency sake, return the tx names (name field) as well
  ( $err, $sampleList, $wantedTxNums ) = $self->makePersonalProtDb($fileType);

  say STDERR "FINISHED STEP 1 (make personal replacement db)";

  ( $err ) = $self->makeReferenceProtDb($wantedTxNums);

  say STDERR "FINISHED STEP 2 (make reference protein db for requested txNumbers)";

  
  ( $err ) = $self->makeReferenceUniquePeptides($sampleList, $wantedTxNums);
  say STDERR "FINISHED STEP 3 (make reference peptide database for trypsin for requested txNumbers)";
  
  ( $err ) = $self->makeSampleUniquePeptides($sampleList, $wantedTxNums);
  say STDERR "FINISHED STEP 4 (make per-sample database of unique peptides digested by trypsin)";
  # We may get rid of this step
  # ( $err ) = $self->createPersProtPermutations($sampleList, $wantedTxNums);

  # TODO: Inspect vcf header

  # TODO: support any other file, by checking the extension

  # TODO: we don't really check for valid vcf, just assume it is
  # So this message is never reached
  $self->_errorWithCleanup("File type isn\'t vcf or snp. Please use one of these files");
  return ( "File type isn\'t vcf or snp. Please use one of these files", undef );
}

sub makePersonalProtDb {

  #Inspired by T.S Wingo: https://github.com/wingolab-org/GenPro/blob/master/bin/vcfToSnp
  my $self = shift;
  my $type = shift;

  # This replicates sub ProcessIds {
  # which just returns an array of ids
  my ( $err, $fh, $outFh, $statsFh, $sampleListAref ) =
  $self->_getFileHandlesAndSampleList($type);

  if ($err) {
    $self->_errorWithCleanup($err);
    return ( $err, undef );
  }

  if ( !@$sampleListAref ) {
    $err =
    "Need samples; did you pass \%sampleList% as an argument to the fileProcessor?";

    $self->_errorWithCleanup($err);
    return ( $err, undef );
  }

  my $personalDb = GenPro::DBManager->new();

  my $metaDb = $personalDb->_getDbi( $metaEnv, undef, $metaConfig );
    my $c = $personalDb->dbReadOne( $metaDb, $metaKeys{samplesWritten} );
    my $previousTxNums = $personalDb->dbReadOne( $metaDb, $metaKeys{wantedTxs} );

  # Need to clean up after ourselves
  # Otherwise we'll have a dangling reference to a free'd environment in LMDB
  undef $metaDb;

  # Create the databases, 1 database per sample,
  # and each sample database was 25 tables, 1 per chromosome
  
  my %completed = $c ? %$c : ();
  my %wantedTxNumbers = $previousTxNums ? %$previousTxNums : ();

  if ( keys %completed == @$sampleListAref ) {
    return ( undef, \%completed, \%wantedTxNumbers);
  }

  # # Ensure that the database is opened
  my $nChr = @{ $self->chromosomes };

  # It *SEEMS* I don't have to pre-make the databases....
  # TODO: Write test for this...
  # for my $sample (@$sampleListAref) {
  #   if ( $completed{$sample} ) {
  #     next;
  #   }

  #   for my $chr ( @{ $self->chromosomes } ) {
  #     $personalDb->_getDbi( $sample, 0, 0, $chr, $nChr );
  #   }
  # }
  
  # Clear the singleton instance, ensure that threads don't copy any memory
  $personalDb->cleanUp();
  undef $personalDb;

  ######################## Build the fork pool #################################

  # Report every 1e4 lines, to avoid thrashing receiver
  my $progressFunc = $self->makeDbWriter( \%completed, \%wantedTxNumbers );

  ########################## Write the header ##################################
  my $header = <$fh>;
  $self->setLineEndings($header);

  # Register the output of the intermediate fileProcessor output with GenPro
  my $finalHeader = $self->_setFinalHeader($header);

  MCE::Loop::init {
    max_workers => $self->maxThreads || 8,
    use_slurpio => 1,

    # bystro-vcf outputs a very small row; fully annotated through the alt column (-ref -discordant)
    # so accumulate less than we would if processing full .snp
    chunk_size => $self->{_chunkSize} . "K",
    gather     => $progressFunc,
  };

  # We separate out the reference track getter so that we can check for discordant
  # bases, and pass the true reference base to other getters that may want it (like CADD)

  # To avoid the Moose/Mouse accessor penalty, store reference to underlying data

  my $refTrackGetter = $self->tracksObj->getRefTrackGetter();

  # TODO: specify the refSeq track name in the YAML config
  # A property on the GenPro track
  my $geneTrackGetter = $self->tracksObj->getTrackGetterByName('refSeq');

  my %wantedChromosomes = %{ $refTrackGetter->chromosomes };
  my $maxDel            = $self->maxDel;

  #Accessors are amazingly slow; it takes as long to call ->name as track->get
  #after accounting for the nubmer of calls to ->name
  # my @trackNames = map { $_->name } @{$self->_tracks};

  # TODO: don't annotate MT (GRCh37) if MT not explicitly specified
  # to avoid issues between GRCh37 and hg19 chrM vs MT
  # my %normalizedNames = %{$self->normalizedWantedChrs};

  # TODO: Grab it from the intermediate output header ($header) and complain if not found
  my $headers = Seq::Headers->new();

  # my $headerIndices = $headers->getParentIndices();

  my $hetIdx = $finalHeader->getFeatureIdx(undef, 'heterozygotes');
  my $homIdx = $finalHeader->getFeatureIdx(undef, 'homozygotes');

  if ( !( defined $homIdx && defined $hetIdx ) ) {
    my $err = "Require 'heterozygotes' and 'homozygotes' output for $type fileProcessor";
    $self->_errorWithCleanup($err);
    die $err;
    return ( $err, undef );
  }

  my $replacement = $geneTrackGetter->siteTypes->{replacement};

  # TODO: Get all these by reading the header, not by using the private
  # interface

  # from the db raw input src (ucsc records)
  my ( $nameIdx, $name2Idx, $txErrorIdx);

  # computed features
  my ($exonicAlleleFuncIdx,  $codonNumIdx , $refCodonIdx, $strandIdx,
    $altCodonIdx , $refAaIdx, $altAaIdx, $codonPosIdx, $txNumberIdx);

  my $trackName = $geneTrackGetter->name;

  $strandIdx = $headers->getFeatureIdx($trackName, 'strand');
  $nameIdx = $headers->getFeatureIdx($trackName, 'name');
  $name2Idx = $headers->getFeatureIdx($trackName, 'name2');
  $txErrorIdx = $headers->getFeatureIdx($trackName, 'txError');

  $exonicAlleleFuncIdx = $headers->getFeatureIdx($trackName, $geneTrackGetter->exonicAlleleFunctionKey);
  $codonNumIdx = $headers->getFeatureIdx($trackName, $geneTrackGetter->codonNumberKey);
  $codonPosIdx = $headers->getFeatureIdx($trackName, $geneTrackGetter->codonPositionKey);
  $refCodonIdx = $headers->getFeatureIdx($trackName, $geneTrackGetter->codonSequenceKey);
  $altCodonIdx = $headers->getFeatureIdx($trackName, $geneTrackGetter->newCodonKey);
  $refAaIdx = $headers->getFeatureIdx($trackName, $geneTrackGetter->refAminoAcidKey);
  $altAaIdx = $headers->getFeatureIdx($trackName, $geneTrackGetter->newAminoAcidKey);
  $txNumberIdx = $headers->getFeatureIdx($trackName, $geneTrackGetter->txNumberKey);

  # Defines feature insertion order
  # This defines what our per-user record looks like at any given position
  # Except we add 5 fields to beginning (chr, pos, type, ref, alt)
  # And 1 field to end (het/hom status
  my @wanted = (
    $exonicAlleleFuncIdx, $strandIdx, $codonNumIdx, $refCodonIdx, $altCodonIdx, $refAaIdx, $altAaIdx,
    $nameIdx, $name2Idx, $codonPosIdx, $txNumberIdx
  );

  my $db = $self->{_db};

  my %sConfig = %{$self->sampleDbConfig};
  mce_loop_f {

    my ( $mce, $slurp_ref, $chunk_id ) = @_;
    #    $_[0], $_[1],     $_[2]
    # open my $MEM_FH, '<', $slurp_ref;
    # binmode $MEM_FH, ':raw';
    open my $MEM_FH, '<', $_[1];
    binmode $MEM_FH, ':raw';

    my $total = 0;

    my $dataFromDbAref;
    my $zeroPos;
    my $chr;
    my @fields;

    #my @results;
    my %seen;
    my %cursors  = ();
    my %pCursors = ();

    my %txNumbers;
    my %dbs;

    # since we cleared db singleton using cleanUp, we now have a new instance
    # local to this thread
    $personalDb  =  GenPro::DBManager->new();
   
    ####### Read #######
    # The ouptut of the intermediate fileProcessor is expected to be
    # chrom \t pos \t type \t inputRef \t alt \t hets \t heterozygosity \t homs ...

    # For GenPro, we store all replacement sites for a given samples, with
    # some informaiton about the site and the affected transcripts

    READ_LOOP: while ( my $line = $MEM_FH->getline() ) {
      chomp $line;

      @fields = split '\t', $line;
      $total++;

      if ( !$wantedChromosomes{ $fields[0] } ) {
        next;
      }

      $chr     = $fields[0];
      $zeroPos = $fields[1] - 1;

      # Caveat: It seems that, per database ($chr), we can have only one
      # read-only transaction; so ... yeah can't combine with dbRead, dbReadOne
      if ( !$cursors{$chr} ) {
        $cursors{$chr} = $db->dbStartCursorTxn($chr);
      }

      $dataFromDbAref = $db->dbReadOneCursorUnsafe( $cursors{$chr}, $zeroPos );

      if ( !defined $dataFromDbAref ) {
        $self->_errorWithCleanup("Wrong assembly? $chr\: $fields[1] not found.");

        # Store a reference to the error, allowing us to exit with a useful fail message
        MCE->gather( 0, 0, "Wrong assembly? $chr\: $fields[1] not found." );
        $_[0]->abort();
        return;
      }

      # For now we don't handle indesl
      if ( length( $fields[4] ) > 1 ) {
        next;
      }

      # for my $trackIndex (@trackIndicesExceptReference) {
      my @out;

      $geneTrackGetter->get(
        $dataFromDbAref,
        $chr,          # chr
        $fields[3],    #ref
        $fields[4],    #alt
        0,
        \@out,
        $zeroPos
      );

      if ( $refTrackGetter->get($dataFromDbAref) ne $fields[3] ) {
        $self->log( 'info', "Discordant: $chr: $fields[1]" );
        next;
      }

      ####### Some or all of our transcripts at this location will not be replacement ####
      # When we have a single transcript, we will have a single value for each
      # feature, unless that feature is nested for that single transcript
      # (ex: sometimes 2 spDisplayID entries for 1 transcript, they will be [1,2])
      # Since exonicAlleleFunction is always 1:1 transcript, a scalar
      # Remember all indices that are replacement, and build up a result record
      # of just those indices
      # All features retain order relative to the transcript,
      # even features that have a 1 tx : M subfeatures (like spDisplayID) mapping
      # So all we need is the indices of the exonicAlleleFunction that are replacement
      my $func = $out[$exonicAlleleFuncIdx][0];

      if ( !defined $func ) {
        next;
      }

      # Check that we don't want stopGain, stopLoss, startGain, or spliceD/A

      my @ok;
      if ( !ref $func ) {
        # When scalar $out[$txErrorIdx][0] will be an array of errors, or undefined
        if( $func ne $replacement || defined $out[$txErrorIdx][0]) {
          next;
        }
      } else {
        my $found = 0;

        for ( my $idx = 0 ; $idx < @$func ; $idx++ ) {
          # Last condition is to ensure we don't have any tx errors entering
          # Any stored sites in the personal protein db will then
          # be directly useful, reliable
          if (
            defined $func->[$idx] &&
            $func->[$idx] eq $replacement &&
            !defined $out[$txErrorIdx][0][$idx]
          ) {
            push @ok, $idx;
            $found ||= 1;
          }
        }

        if ( !$found ) {
          next;
        }
      }

      ################### What we aim to write ########################
      # TODO: Can optimize away insertion of exonicAlleleFunction is only
      # interested in replacement sites
      # We're going to do this as an array, following @wanted, defined above
      # my $record_href = {
      #   aa_residue => $aa_residue,
      #   old_aa     => $codon_2_aa{$codon},
      #   new_aa     => $codon_2_aa{$new_codon},
      #   chr        => $chr,
      #   chr_pos    => $pos,
      #   codon_pos  => $codon_pos,
      #   ref_allele => $ref_allele,
      #   min_allele => $min_allele,
      # };

      # TODO: Optimize allocations
      # flexible features (such as description, not necessary for the Gene class)
      my @record;

      # First 5 records are chr, pos, type, ref, alt; last record is whether it's a het or a hom
      $#record = $#wanted + 6;

      $record[0] = $fields[0];    #chr
      $record[1] = $fields[1];    #pos
      $record[2] = $fields[2];    #type (SNP, INS, DEL, etc), could remove
      $record[3] = $fields[3];    #ref
      $record[4] = $fields[4];    #alt

      # could be more elegant
      # build up the array of results
      if ( @ok == 0 ) {
        my $outIdx = 5;
        for my $fIdx (@wanted) {
          # 0 because non-indels have only 1 position index
          $record[$outIdx] = $out[$fIdx][0];
          $outIdx++;
        }
      } elsif ( @ok == 1 ) {
        my $txIdx  = $ok[0];
        my $outIdx = 5;

        for my $fIdx (@wanted) {
          # 0 because non-indels have only 1 position index
          $record[$outIdx] = $out[$fIdx][0][$txIdx];
          $outIdx++;
        }
      } else {
        my $seenFirst;

        for my $txIdx (@ok) {

          if ( !$seenFirst ) {
            my $outIdx = 5;
            for my $fIdx (@wanted) {
              # 0 because non-indels have only 1 position index
              $record[$outIdx] = [ $out[$fIdx][0][$txIdx] ];
              $outIdx++;
            }

            $seenFirst = 1;
            next;
          }

          my $outIdx = 5;
          for my $fIdx (@wanted) {
            # 0 because non-indels have only 1 position index
            push @{ $record[$outIdx] }, $out[$fIdx][0][$txIdx];
            $outIdx++;
          }
        }
      }

      # Todo use the configured delimiters from Seq::Ouput->delimiters
      my %hets;
      if ( $fields[$hetIdx] ne '!' ) {
        %hets = map { $_ => 1 } split( ';', $fields[$hetIdx] );
      }

      my %homs;
      if ( $fields[$homIdx] ne '!' ) {
        %homs = map { $_ => 1 } split( ';', $fields[$homIdx] );
      }

      my $cnt = 0;
      # Can no longer use keys % in scalar context
      # This seems to work
      my $max = %hets + %homs;

      my $sIdx;
      my @samples;
      my $het = 1;

      $txNumbers{$chr} //= {};
      my @nums = ref $record[-2] ? @{$record[-2]} : $record[-2];
      for my $num (@nums) {
        $txNumbers{$chr}{$num} //= 1;
      }

      # my @wantedSamples = 

      # Optimize, only read over hets + homs
      for my $sample (@$sampleListAref) {
        if ( $cnt == $max ) {
          # say STDERR "count is $cnt and max is $max: $line";
          last;
        }

        # Store whether het or not
        if ( $hets{$sample} ) {
          $record[-1] = 1;
        } elsif ( $homs{$sample} ) {
          $record[-1] = 0;
        } else {
          next;
        }

        $seen{$sample} //= 1;

        $cnt++;
        # Doesn't work yet
        #my (  $cursor, $chr, $trackKey, $pos, $newValue, $mergeFunc) = @_;
        # $personalDb->dbPatchCursorUnsafe( $pCursors{$sample}{$chr},
        #   $chr, 0, $zeroPos, \@record );

        # dbPatch with per-op commit seems to work fine

        # Here 0 is the index of the "Genpro" database
        # This is of course not really necessary, it's just what Bystro does
        # If we ended up adding a GenPro track to Bystro, this could be the
        # sample index
        # TODO: Use cursor; don't add unnecessary track key, add merge function
        # in case of multiallelics
        if(!$dbs{$sample}) {
          $dbs{$sample} = { $chr => $personalDb->_getDbi( $sample, $chr, \%sConfig ) };
        } elsif(!$dbs{$sample}{$chr}) {
          $dbs{$sample}{$chr} = $personalDb->_getDbi( $sample, $chr, \%sConfig );
        }

        # We don't currently support multiallelics
        $personalDb->dbPut( $dbs{$sample}{$chr}, $zeroPos, \@record );
      }
    }

    MCE->gather( \%seen, \%txNumbers );

    # Get in a good habit; need to free environment pointers for 
    # LMDBB to clean up properly
    # undef %dbs;
    # $personalDb->cleanUp();

    close $MEM_FH;
  }
  $fh;

  # Force flush of output
  # $progressFunc->(undef);

  MCE::Loop::finish();

  $personalDb = GenPro::DBManager->new();

  $metaDb = $personalDb->_getDbi( $metaEnv, undef, $metaConfig);
    $personalDb->dbPut($metaDb, $metaKeys{samplesWritten}, \%completed );
    $personalDb->dbPut($metaDb, $metaKeys{wantedTxs}, \%wantedTxNumbers );
  undef $metaDb;

  $db->cleanUp();
  $personalDb->cleanUp();

  return ( $err, \%completed, \%wantedTxNumbers );
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
  my $nChrs = @{$self->chromosomes};
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

      if($txNumber == 438) {
        say STDERR "Working on txNumber 438";
      }
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
    if(defined $writtenTxNums{$chr}) {
      next;
    }

    for my $txNum (sort { $a <=> $b } keys %{$wantedTxNumHref->{$chr}}) {
      if($writtenTxNums{$chr}{$txNum}) {
        next;
      }

      push @txNums, [$chr, $txNum];
      $wantedChrs{$chr} //= 1;
    }
  }

  if(@txNums == 0) {
    return;
  }
  
  # TODO: when working with named databases, must open each time with same flags
  # and I think, same number of maxDbs (or at least, no fewer max dbs?)
  my $nChrs = @{$self->chromosomes};
  
  my %refProtRdOnlyConfig = (%{$self->refProtConfig}, (readOnly => 1));
  my $progressFunc = sub {
    my $seenTxNums = shift;

    for my $chr (keys %$seenTxNums) {
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

        for my $peptide (trypsin($aaSeqAref)) {
          # say STDERR $peptide . "\n";
          my $txn = $dbs{$enzyme}{env}->BeginTxn();

          $data = $personalDb->dbReadOneRaw($txn, $dbi, $peptide);
          
          if(!defined $data) {
            $personalDb->dbPutRaw($txn, $dbi, $peptide, [1, [$txInfo]]);
            # say STDERR "Count: 0";
          } else {
            $data->[0] = $data->[0] + 1;
            push @{$data->[1]}, $txInfo;

            $personalDb->dbPutRaw($txn, $dbi, $peptide, $data);
            # say STDERR "Count: ref: $data->[0]";
            # p $data;
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
sub makeSampleUniquePeptides {
  my ($self, $wantedSamplesHref, $wantedTxNumHref, $wantedEnzymesAref) = @_;

  $wantedEnzymesAref //= ['trypsin'];

  

  # When working with named databases, need to pass same configuration each time
  # This is to know the maxDbs config for the reference peptide database
  # which stores n peptides ('for now just trypsin');
  # and maybe needs to be at least no less than the largest number it was last
  # opened with (or maybe needs to be constant)
  # TODO: Make this a configuration owned by the reference peptide database
  my $nEnzymeTables = 3;

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

  # When working with named databases, need to pass same configuration each time
  # This is to know the maxDbs config for the reference protein database
  # and maybe needs to be at least no less than the largest number it was last
  # opened with (or maybe needs to be constant)
  # TODO: Make this a configuration owned by the reference protein database
  my $nChrs = @{$self->chromosomes};

  my @wanted = (
    'chr', 'pos', 'type', 'ref', 'alt',
    'exonicAlleleFunction', 'strand', 'codonNum', 'refCodon', 'altCodon', 'refAa', 'altAa',
    'name', 'name2', 'codonPos', 'txNumber', 'het'
  );

  my %wanted;
  for (my $i = 0; $i < @wanted; $i++) {
    $wanted{$wanted[$i]} = $i;
  }

  # TODO: add other enzymes
  my $digest = GenPro::Digest->new();
  my $digestFunc = $digest->makeDigestFunc('trypsin');
  my $cutsHref = $digest->digestLookups->{'trypsin'}{'cut'};

  if(!$digestFunc) {
    die 'No digest func for trypsin';
  }

  my $txNumberIdx = $wanted{'txNumber'};
  my $altAaIdx = $wanted{'altAa'};
  my $refAaIdx = $wanted{'refAa'};
  my $codonNumIdx = $wanted{'codonNum'};
  
  my %txHasUniqueConfigTemplate = (
    altAaIdx => $altAaIdx,
    refAaIdx => $refAaIdx,
    codonNumIdx => $codonNumIdx,
    cutsHref => $cutsHref,
    maxPeptideLength => $digest->maxPeptideLength,
  );
  
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

    my $existing;
    my @peptides;
    my %seenSamples;

    my %refPeptideRdOnlyConfig = (%$refPeptideConfig, (readOnly => 1));
    my %refProtRdOnlyConfig = (%{$self->refProtConfig}, (readOnly => 1));
    my %sampleReadOnlyConfig = (%{$self->sampleDbConfig}, (readOnly => 1));

    my $refPeptideDb = $personalDb->_getDbi( $refPeptideEnv, 'trypsin', \%refPeptideRdOnlyConfig );
    my $refPeptideDbi = $refPeptideDb->{dbi};
    my $refPeptideTxn = $refPeptideDb->{env}->BeginTxn();

    my %refProteinDbs;

    my $globalUniquePeptideFn = sub {
      my $aaAref = shift;

      my @unique;
      for my $peptide (trypsin($aaAref)) {
        # say STDERR "PEPTIDE: $peptide";
        if(!defined $personalDb->dbReadOneRaw($refPeptideTxn, $refPeptideDbi, $peptide)) {
          # say STDERR "NOT DEFINED $peptide";
          push @unique, $peptide;
        }
      }

      return @unique;
    };

    my %txHasUniqueConfig = (%txHasUniqueConfigTemplate, (globalUniquePeptideFn => $globalUniquePeptideFn));
    my $txHasUniquePeptidesFn = makeVarProtForTxFn(\%txHasUniqueConfig);

    # p %txHasUniqueConfig;

    for my $sample (@{ $chunk_ref }) {
      my @finalRecords;

      my $samplePeptideDb = $personalDb->_getDbi( "$sample/peptide", undef, $refPeptideConfig );

      # TODO: take chromosomes from those written
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

sub trypsin {
  my $aaAref = shift;

  my @cutSites;
  my $lastCutSite = 0;

  my $i = -1;
  my $start;
  my $seq;
  my @peptides;
  for my $aa (@$aaAref) {
    $i++;

    if($aa eq '*') {
      last;
    }

    if($aa eq 'K' || $aa eq 'R') {
      push @cutSites, $i + 1;
      $lastCutSite = $i + 1;
    }
  }

  # if a K or R at the end, don't include it twice but also include the
  # end of the protein
  if( $lastCutSite < @$aaAref) {
    push @cutSites, scalar @$aaAref;
  }

  # reset
  $lastCutSite = 0;
  
  $i = -1;
  for my $end (@cutSites) {
    $i++;

    # Don't add + 1 because we take lastCutSite .. $end - 1, so this is already +1 length
    # and because end in the above is $i + 1
    if($end - $lastCutSite >= $minPeptideLength && $end - $lastCutSite <= $maxPeptideLength) {
      # No need to join into a string...that's a byte array anyway
      # Just keep the reference to an array reference
      push @peptides, join('', @$aaAref[ $lastCutSite .. $end - 1 ]);
    }

    if ( $end == @$aaAref ) {
      last;
    }

    # If the next base is a proline that blocks
    if( $aaAref->[$end] eq 'P' ) {
      # any more cut sites?
      if( $i + 1 < @cutSites) {
        $end = $cutSites[ $i + 1 ];

        if($end - $lastCutSite >= $minPeptideLength && $end - $lastCutSite <= $maxPeptideLength) {
          # No need to join into a string...that's a byte array anyway
          # Just keep the reference to an array reference
          push @peptides, join('', @$aaAref[ $lastCutSite .. $end - 1 ] );
        }
      }
    }

    $lastCutSite = $end;
  }

  return @peptides;
}

sub makeDbWriter {
  my ( $self, $samplesSeenHref, $txNumbersSeenHref ) = @_;

  return sub {
    my ($newSeen, $txNumbersHref) = @_;

    for my $sample ( keys %$newSeen ) {
      $samplesSeenHref->{$sample} //= 1;
    }

    for my $chr (keys %$txNumbersHref) {
      $txNumbersSeenHref->{$chr} //= {};

      for my $txNum (keys %{$txNumbersHref->{$chr}}) {
        $txNumbersSeenHref->{$chr}{$txNum} //= 1;
      }
    }
  }
}

################################################################################

sub _getFileHandles {
  my ( $self, $type ) = @_;

  my ( $outFh, $statsFh, $inFh );
  my $err;

  ( $err, $inFh ) = $self->_openAnnotationPipe($type);

  if ($err) {
    return ( $err, undef, undef, undef );
  }

  if ( $self->run_statistics ) {
    ########################## Tell stats program about our annotation ##############
    # TODO: error handling if fh fails to open
    my $statArgs = $self->_statisticsRunner->getStatsArguments();

    $err = $self->safeOpen( $statsFh, "|-", $statArgs );

    if ($err) {
      return ( $err, undef, undef, undef );
    }
  }

  # $fhs{stats} = $$statsFh;
  ( $err, $outFh ) = $self->getWriteFh( $self->{_outPath} );

  if ($err) {
    return ( $err, undef, undef, undef );
  }

  return ( undef, $inFh, $outFh, $statsFh );
}

# Open the pipe to the annotation pre-processor, then close it
# All file setup, including generating sample list happens
# Regardless of how many lines of the output we read;
# The sample list will get generated, if requested,
# Even without reading from the stdout filehandle
sub _getFileHandlesAndSampleList {
  my ( $self, $type ) = @_;

  my ( $err, $fh, $outFh, $statsFh ) = $self->_getFileHandles($type);

  if ($err) {
    return ( $err, undef );
  }

  my $sampleListPath =
  $self->_workingDir->child( $self->outputFilesInfo->{sampleList} )->stringify();

  my $sampleListFh;
  $err = $self->safeOpen( $sampleListFh, '<', $sampleListPath );

  if ($err) {
    return ( $err, undef );
  }

  my @samples = <$sampleListFh>;

  $self->safeClose($sampleListFh);

  chomp @samples;

  return ( undef, $fh, $outFh, $statsFh, \@samples );
}

sub _openAnnotationPipe {
  my ( $self, $type ) = @_;

  my ( $errPath, $inPath, $echoProg ) = $self->_getFileInfo();

  my $fh;

  my $args = $self->_getFileProcessor($type);

  # TODO:  add support for GQ filtering in vcf
  my $err = $self->safeOpen( $fh, '-|', "$echoProg $inPath | $args 2> $errPath" );

  return ( $err, $fh );
}

sub _getFileInfo {
  my $self = shift;

  my $errPath =
  $self->_workingDir->child( $self->input_file->basename . '.file-log.log' );

  my $inPath = $self->inputFilePath;
  my $echoProg =
    $self->isCompressedSingle($inPath)
  ? $self->gzip . ' ' . $self->decompressArgs
  : 'cat';

  return ( $errPath, $inPath, $echoProg );
}

sub _getFileProcessor {
  my ( $self, $type ) = @_;

  if ( !$self->fileProcessors->{$type} ) {
    $self->_errorWithCleanup("No fileProcessors defined for $type file type");
  }

  my $fp   = $self->fileProcessors->{$type};
  my $args = $fp->{program} . " " . $fp->{args};

  for my $type ( keys %{ $self->outputFilesInfo } ) {
    if ( index( $args, "\%$type\%" ) > -1 ) {
      substr( $args, index( $args, "\%$type\%" ), length("\%$type\%") ) =
      $self->_workingDir->child( $self->outputFilesInfo->{$type} );
    }
  }

  return $args;
}

sub _setFinalHeader {
  my ( $self, $header ) = @_;
  ######### Build the header, and write it as the first line #############
  my $finalHeader = Seq::Headers->new();

  chomp $header;

  my @headerFields = split '\t', $header;

  # Our header class checks the name of each feature
  # It may be, more than likely, that the pre-processor names the 4th column 'ref'
  # We replace this column with trTv
  # This not only now reflects its actual function
  # but prevents name collision issues resulting in the wrong header idx
  # being generated for the ref track
  $headerFields[3] = $self->discordantField;

  # Bystro takes data from a file pre-processor, which spits out a common
  # intermediate format
  # This format is very flexible, in fact Bystro doesn't care about the output
  # of the pre-processor, provided that the following is found in the corresponding
  # indices:
  # idx 0: chromosome,
  # idx 1: position
  # idx 3: the reference (we rename this to discordant)
  # idx 4: the alternate allele

  # Prepend all of the headers created by the pre-processor
  $finalHeader->addFeaturesToHeader( \@headerFields, undef, 1 );

  return $finalHeader;
}

sub _errorWithCleanup {
  my ( $self, $msg ) = @_;

  $self->log( 'error', $msg );

  $self->{_db}->cleanUp();

  return $msg;
}

__PACKAGE__->meta->make_immutable;

1;
