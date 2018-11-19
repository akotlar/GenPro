use 5.10.0;
use strict;
use warnings;

package GenProSeq;

our $VERSION = '0.001';

# ABSTRACT: Annotate a snp file

# TODO: make temp_dir handling more transparent
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
use Seq::Tracks::Gene::Site;
use Seq::Tracks::Gene::Site::CodonMap;
use Seq::DBManager;
use Path::Tiny;
use Scalar::Util qw/looks_like_number/;
use DDP;
use Test::More;

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

my $refProtEnv = 'referenceProteins';
my $refPeptideEnv = 'referencePeptides';


my $min_peptide_length = 6;
my $max_peptide_length = 40;
my $trim_end_regex     = qr{\*[\*\w]*\z};

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

  my $db = $self->{_db};

  my $personalDb = GenPro::DBManager->new();

  $personalDb->_getDbi( $metaEnv, undef, 1);

  my %okToBuild;
  my %seen;
  # Create the databases, 1 database per sample,
  # and each sample database was 25 tables, 1 per chromosome
  my $nChr = @{ $self->chromosomes };

  my $toBuild = 0;
  my %completionMeta;
  # my $nSamples = @$sampleListAref;
  my $c = $personalDb->dbReadOne( $metaEnv, undef, $metaKeys{samplesWritten} );
  my $previousTxNums = $personalDb->dbReadOne( $metaEnv, undef, $metaKeys{wantedTxs} );

  my %completed = $c ? %$c : ();
  my %wantedTxNumbers = $previousTxNums ? %$previousTxNums : ();

  if ( keys %completed == @$sampleListAref ) {
    return ( undef, \%completed, \%wantedTxNumbers);
  }

  for my $sample (@$sampleListAref) {
    if ( $completed{$sample} ) {
      next;
    }

    for my $chr ( @{ $self->chromosomes } ) {
      $personalDb->_getDbi( $sample, 0, 0, $chr, $nChr );
    }
  }

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

  mce_loop_f {

    my ( $mce, $slurp_ref, $chunk_id ) = @_;
    #    $_[0], $_[1],     $_[2]
    # open my $MEM_FH, '<', $slurp_ref;
    # binmode $MEM_FH, ':raw';
    open my $MEM_FH, '<', $_[1];
    binmode $MEM_FH, ':raw';

    my $total = 0;

    # my @indelDbData;
    # my @indelRef;
    # my @lines;
    my $dataFromDbAref;
    my $zeroPos;
    my $chr;
    my @fields;

    #my @results;
    my %seen;
    my %cursors  = ();
    my %pCursors = ();

    my %txNumbers;
    # GenPro::DBManager::initialize(
    #   {
    #     databaseDir => path( $self->database_dir )->parent()->child("genpro")->stringify()
    #   }
    # );

    # my $personalDb = GenPro::DBManager->new();

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
        $personalDb->dbPatch( $sample, $chr, 0, $zeroPos, \@record );

        # my $name = $record[$nameIdx];
        # my $altAa = $record[$altAaIdx];
        # my $altCodon = $record[$altCodonIdx];
        # my $codonPos = $record[$codonPosIdx];
        # my $codonNum = $record[$codonNumIdx];
        # my $pos = $zeroPos;

        # my $txNumber = $record[$txNumberIdx];
        # my $strand = $record[$strandIdx];

        # if(ref $name) {
        #   # p $data;
        #   for (my $i = 0; $i < @$name; $i++) {
        #     $userTx{$txNumber->[$i]} //= [];

        #     my @details = (
        #       $chr, $pos, $txNumber->[$i], $name->[$i], $codonNum->[$i],
        #       $codonPos->[$i], $altCodon->[$i], $altAa->[$i], $strand->[$i]
        #     );

        #     push @{$userTx{$txNumber->[$i]}}, \@details;
        #   }
        # } else {
        #   $userTx{$txNumber} //= [];

        #   # If any of the inner stuff is an array, it is guaranteed to
        #   # correspond to the $name
        #   # for instance, if multiple spId's, they all belong to the 1 $name
        #   my @details = (
        #     $chr, $pos, $txNumber, $name, $codonNum, $codonPos,
        #     $altCodon, $altAa, $strand
        #   );

        #   push @{$userTx{$txNumber}}, \@details;
        # }

        # # push @txNumbers, $record[$txNumberIdx];
        # push @{$userTx{$txNumber->[$i]}}, \@details;
        # Uncomment to test
        # my $val = $personalDb->dbReadOne( $sample, $chr, $zeroPos );
      }

      # If it becomes unsafe to use NOTLS and write from multiple processes
      # push @results, [ \@samples, $chr, $zeroPos, \@record ];

    }

    # say STDERR "Done";

    # for my $sample (@$sampleListAref) {
    #   $personalDb->dbEndCursorTxn($sample);
    # }
    # $personalDb->dbEndCursorTxn()

    MCE->gather( \%seen, \%txNumbers );

    close $MEM_FH;
  }
  $fh;

  # Force flush of output
  # $progressFunc->(undef);

  MCE::Loop::finish();

  $personalDb->dbPut( $metaEnv, undef, $metaKeys{samplesWritten}, \%completed );

  $personalDb->dbPut( $metaEnv, undef, $metaKeys{wantedTxs}, \%wantedTxNumbers );

  $db->cleanUp();

  # Unfortunately, MCE::Shared::stop() removes the value of $abortErr
  # according to documentation, and I did not see mention of a way
  # to copy the data from a scalar, and don't want to use a hash for this alone
  # So, using a scalar ref to abortErr in the gather function.
  # if ($abortErr) {
  #   say "Aborted job due to $abortErr";
  #   return ( 'Job aborted due to error', undef );
  # }

  ################ Finished writing file. If statistics, print those ##########
  # Sync to ensure all files written
  # This simply tries each close/sync/move operation in order
  # And returns an error early, or proceeds to next operation
  # $err =
  #    $self->safeClose($outFh)
  # || ( $statsFh && $self->safeClose($statsFh) )
  # || $self->safeSystem('sync')
  # || $self->_moveFilesToOutputDir();

  return ( $err, \%completed, \%wantedTxNumbers );
}

# Generates everything that GenPro_make_refprotdb does
# Also pre-calculates all digested peptides, and stores those
sub makeReferenceProtDb {
  my ($self, $wantedTxNumHref) = @_;

  my $personalDb = GenPro::DBManager->new();
  $personalDb->_getDbi( $metaEnv, 'completed', 1 );

  my $previouslyWritten = $personalDb->dbReadOne( $metaEnv, undef, $metaKeys{refTxsWritten} );
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
  for my $chr ( keys %wantedChrs ) {
    $regionDb{$chr} //= $db->dbReadAll( $geneTrackGetter->regionTrackPath($chr) );

    $personalDb->_getDbi( $refProtEnv, 0, 0, $chr, $nChrs );

    if(defined $writtenChrs{$chr}) {
      $haveWritten++;
    }
  }

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

    my %seenCodonRanges;
    my %cursors;
    my %seenChrs;
    for my $txNumInfo (@{ $chunk_ref }) {
      my $chr = $txNumInfo->[0];
      my $txNumber = $txNumInfo->[1];

      $seenChrs{$chr} //= 1;

      # To save
      my $tx = $regionDb{$chr}->[$txNumber];
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

      $seenCodonRanges{$chr} //= {};
      $seenCodonRanges{$chr}{$cdsStart} //= {};

      if($seenCodonRanges{$chr}{$cdsStart}{$cdsEnd}) {
        next;
      }

      $seenCodonRanges{$chr}{$cdsStart}{$cdsEnd} = 1;

      if($strand eq '-') {
        @data = reverse @data;
      }

      $cursors{$chr} //= $db->dbStartCursorTxn($chr);
      $db->dbReadCursorUnsafe($cursors{$chr}, \@data);

      my ( $fetchedTxNumbers, $siteData, $wantedSiteData );
      my (@codonSequence, @aaSequence, @codonPos, @codonNums);

      my $lastCodonNumber;
      for my $site (@data) {
        if(!defined) {
          die 'Err: no gene site';
        }

        ( $fetchedTxNumbers, $siteData ) = $siteUnpacker->unpack( $site->[ $geneTrackGetterDbName ] );

        # if($txNumber == 2211) {
        #   p $fetchedTxNumbers;
        #   p $siteData;
        # }
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
        ok(join('', @aaSequence eq 'MGKSEGPVGMVESAGRAGQKRPGFLEGGLLLLLLLVTAALVALGVLYADRRGKQLPRLASRLCFLQEERTFVKRKPRGIPEAQEVSEVCTTPGCVIAAARILQNMDPTTEPCDDFYQFACGGWLRRHVIPETNSRYSIFDVLRDELEVILKAVLENSTAKDRPAVEKARTLYRSCMNQSVIEKRGSQPLLDILEVVGGWPVAMDRWNETVGLEWELERQLALMNSQFNRRVLIDLFIWNDDQNSSRHIIYIDQPTLGMPSREYYFNGGSNRKVREAYLQFMVSVATLLREDANLPRDSCLVQEDMVQVLELETQLAKATVPQEERHDVIALYHRMGLEELQSQFGLKGFNWTLFIQTVLSSVKIKLLPDEEVVVYGIPYLQNLENIIDTYSARTIQNYLVWRLVLDRIGSLSQRFKDTRVNYRKALFGTMVEEVRWRECVGYVNSNMENAVGSLYVREAFPGDSKSMVRELIDKVRTVFVETLDELGWMDEESKKKAQEKAMSIREQIGHPDYILEEMNRRLDEEYSNLNFSEDLYFENSLQNLKVGAQRSLRKLREKVDPNLWIIGAAVVNAFYSPNRNQIVFPAGILQPPFFSKEQPQALNFGGIGMVIGHEITHGFDDNGRNFDKNGNMMDWWSNFSTQHFREQSECMIYQYGNYSWDLADEQNVNGFNTLGENIADNGGVRQAYKAYLKWMAEGGKDQQLPGLDLTHEQLFFINYAQVWCGSYRPEFAIQSIKTDVHSPLKYRVLGSLQNLAAFADTFHCARGTPMHPKERCRVW'));
      }

      $personalDb->dbPut($refPeptideEnv, $chr, $txNumber, [\@aaSequence, \@codonSequence]);
    }

    MCE->gather(\%seenChrs);
  } @txNums;

  MCE::Loop::finish();

  $personalDb->dbPut( $metaEnv, undef, $metaKeys{refTxsWritten}, \%writtenChrs );

  $personalDb->cleanUp();
}

# Generates everything that GenPro_make_refprotdb does
# Also pre-calculates all digested peptides, and stores those
sub makeReferenceUniquePeptides {
  my ($self, $wantedSamples, $wantedTxNumHref) = @_;
# p $wantedTxNumHref;
  my $personalDb = GenPro::DBManager->new();
  $personalDb->_getDbi( $refPeptideEnv, undef, 1 );

  $personalDb->_getDbi( $metaEnv, undef, 1 );

  my $previouslyWritten = $personalDb->dbReadOne( $metaEnv, undef, $metaKeys{refPeptidesWritten} );

  my %writtenTxNums = $previouslyWritten ? %$previouslyWritten : ();
# p %writtenTxNums;
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
# p @txNums;
  if(@txNums == 0) {
    return;
  }

  my $progressFunc = sub {
    my $seenChrs = shift;

    for my $chr (keys %$seenChrs) {
      for my $txNum (keys %{$seenChrs->{$chr}}) {
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

    my %seenChrs;
    for my $txNumInfo (@{ $chunk_ref }) {
      my $chr = $txNumInfo->[0];
      my $txNumber = $txNumInfo->[1];

      $seenChrs{$chr} //= 1;

      my $seqInfo = $personalDb->dbReadOne($refProtEnv, $chr, $txNumber);

      # p $seqInfo;

    }

    MCE->gather(\%seenChrs);
  } @txNums;

  MCE::Loop::finish();

  # $personalDb->dbPut( $metaEnv, undef, $metaKeys{refTxsWritten}, \%writtenChrs );
}

# Replicates:
# 1) sub var_prot_for_tx, which:
## # var_prot_for_tx takes a variant record and returns an array of variant proteins;
#### the caller may select the most parsimonious elements
sub createPersProtPermutations {
  my $self           = shift;
  my $sampleListHref = shift;
  my $wantedTxNumHref = shift;

  my $err;

  my $personalDb = GenPro::DBManager->new();

  my $db              = $self->{_db};
  my $geneTrackGetter = $self->tracksObj->getTrackGetterByName('refSeq');

  my $completed = $personalDb->dbReadOne( $metaEnv, undef, $metaKeys{samplePeptidesWritten} );
 
  my $found = 0;
  for my $sample (keys %$sampleListHref) {
    if(defined $completed->{$sample}) {
      $found++;
    }
  }

  if($found == %$sampleListHref) {
    return;
  }

  # my $start = time();
  my %refProtSeqData;
  for my $chr (sort { $a cmp $b } keys %$wantedTxNumHref) {
    $refProtSeqData{$chr} //= {};
    for my $txNum (sort { $a <=> $b } keys %{$wantedTxNumHref->{$chr}}) {
      # push @txNums, [$chr, $txNum];
      $refProtSeqData{$chr}{$txNum} =  $personalDb->dbReadOne($refProtEnv, $chr, $txNum);
    }
  }
  
  # say STDERR "PAST: " . (time() - $start) . " : num stuff: " . (keys %{$refProtSeqData{chr1}});
  # sleep(5);
  my %regionDb;
  my %byTxRegion;

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


  # TODO: store the order of these in meta table
  # since that is effectively the feature index
  # my @wanted = (
  #   $exonicAlleleFuncIdx, $strandIdx, $codonNumIdx, $refCodonIdx, $altCodonIdx, $refAaIdx, $altAaIdx,
  #   $nameIdx, $name2idx, $codonPosIdx
  # );
  my @wanted = (
    'chr', 'pos', 'type', 'ref', 'alt',
    'exonicAlleleFunction', 'strand', 'codonNum', 'refCodon', 'altCodon', 'refAa', 'altAa',
    'name', 'name2', 'codonPos', 'txNumber', 'het'
  );

  my %wanted;
  for (my $i = 0; $i < @wanted; $i++) {
    $wanted{$wanted[$i]} = $i;
  }

  # These are never arrays
  my $chrIdx = $wanted{'chr'};
  my $posIdx = $wanted{'pos'};
  my $typeIdx = $wanted{'type'};
  my $refIdx = $wanted{'ref'};
  my $altIdx = $wanted{'alt'};
  my $hetIdx = $wanted{'het'};

  # These can be arrays, if there are multiple transcripts
  my $nameIdx = $wanted{'name'};
  my $altAaIdx = $wanted{'altAa'};
  my $altCodonIdx = $wanted{'altCodon'};
  my $codonPosIdx = $wanted{'codonPos'};
  my $codonNumIdx = $wanted{'codonNum'};
  my $txNumberIdx = $wanted{'txNumber'};
  my $strandIdx = $wanted{'strand'};

  my $siteUnpacker = Seq::Tracks::Gene::Site->new();
  my $siteTypeMap  = Seq::Tracks::Gene::Site::SiteTypeMap->new();
  my $codonMap     = Seq::Tracks::Gene::Site::CodonMap->new();

  my $strandSiteIdx        = $siteUnpacker->strandIdx;
  my $siteTypeSiteIdx      = $siteUnpacker->siteTypeIdx;
  my $codonSequenceSiteIdx = $siteUnpacker->codonSequenceIdx;
  my $codonPositionSiteIdx = $siteUnpacker->codonPositionIdx;
  my $codonNumberSiteIdx   = $siteUnpacker->codonNumberIdx;

  my $geneTrackGetterDbName = $geneTrackGetter->dbName;

  # Report every 1e4 lines, to avoid thrashing receiver
  my %wantedTranscripts;
  my %sampleData;
  my $progressFunc = sub {
    my ($chr, $sample, $userTxHref);

    my @txNums = keys %$userTxHref;

    $wantedTranscripts{$chr} //= {};


    for my $txNum (@txNums) {
      $wantedTranscripts{$chr}{$txNum} //= 1
    }

    $sampleData{$sample} //= {};
    $sampleData{$sample} = $userTxHref;
  };

  MCE::Loop::init {
    max_workers => $self->maxThreads || 8,

    # bystro-vcf outputs a very small row; fully annotated through the alt column (-ref -discordant)
    # so accumulate less than we would if processing full .snp
    chunk_size => 1,
    gather     => $progressFunc,
  };

  mce_loop {
    my ( $mce, $chunk_ref, $chunk_id ) = @_;

    my $sample = $_;
    # Each thread gets its own cursor
    my  %cursors = ();

    # TODO: optimize wanted chromosomes
    for my $chr ( @{ $self->chromosomes } ) {
      my $dataAref = $personalDb->dbReadAll( $sample, $chr );

      p $dataAref;

      if(!defined $dataAref) {
        next;
      }

      # Place all data into
      # txName => [first modification, 2nd modification, etc]
      # Where each modification = [chr, pos, name, codonNum, codonPos, altCodon, altAa]
      my %userTx;

      # Caveat: It seems that, per database ($chr), we can have only one
      # read-only transaction; so ... yeah can't combine with dbRead, dbReadOne
      if ( !$cursors{$chr} ) {
        $cursors{$chr} = $db->dbStartCursorTxn($chr);
      }

      for my $data (@$dataAref) {
        my $name = $data->[0][$nameIdx];
        my $altAa = $data->[0][$altAaIdx];
        my $altCodon = $data->[0][$altCodonIdx];
        # my $codonPos = $data->[0][$codonPosIdx];
        my $codonNum = $data->[0][$codonNumIdx];
        # my $pos = $data->[0][$posIdx];

        my $txNumber = $data->[0][$txNumberIdx];
        my $strand = $data->[0][$strandIdx];

        if(ref $name) {
          # p $data;
          for (my $i = 0; $i < @$name; $i++) {
            $userTx{$txNumber->[$i]} //= [];

            my @details = (
              $txNumber->[$i], $name->[$i], $altCodon->[$i], $altAa->[$i], $strand->[$i]
            );

            push @{$userTx{$txNumber->[$i]}}, \@details;
          }
        } else {
          $userTx{$txNumber} //= [];

          # If any of the inner stuff is an array, it is guaranteed to 
          # correspond to the $name
          # for instance, if multiple spId's, they all belong to the 1 $name
          my @details = (
            $txNumber, $name, $altCodon, $altAa, $strand
          );

          push @{$userTx{$txNumber}}, \@details;
        }
      }

      # TODO: it may be more efficient to build all transcripts at
      # once, at cost of memory
      # To not discard info / waste iterations when transcripts overlap
      my @txNums = keys %userTx;
      
      for my $txNum (@txNums) {
        my $refSeqData = $refProtSeqData{$chr}{$txNum};

        p $refSeqData;
        sleep(10);
      }

      # MCE->gather($sample, $chr, \%userTx);
    }
  }
  sort { $a cmp $b } keys %$sampleListHref;

  MCE::Loop::finish();

  return ($err);
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

################ Borrowed from Thomas, with minor modifications ################
# Choose informativ
sub select_entries {
  my $recs_aref = shift;

  my @final_records;

  my @s_records = sort_by_uniq_peptides(@$recs_aref);

  while (@s_records) {
    my $rec = shift @s_records;
    push @final_records, $rec;

    # add the tryptic peptides to the global db, which will be used when
    # creating the next round of counts
    add_seq_to_trp_db( $rec->{seq} );
    @s_records = sort_by_uniq_peptides(@s_records);
  }
  return \@final_records;
}

sub sort_by_uniq_peptides {
  my (@records) = @_;

  my @s_records;
  my @s_records_with_count = sort { $b->[0] <=> $a->[0] }
  map { [ new_global_peptide( $_->{seq} ), $_ ] } @records;

  # only want records if they contribute unique peptides
  for ( my $i = 0 ; $i < @s_records_with_count ; $i++ ) {
    if ( $s_records_with_count[$i]->[0] > 0 ) {
      push @s_records, $s_records_with_count[$i]->[1];
    }
  }
  return @s_records;
}

# Will not be used
my %trp_peptides;
my %trpPepDb;
my @cut_sites;

# new_global_peptide takes a string and returns a count of the number of globally
# unique peptides the string possessesk
sub new_global_peptide {
  my $seq = shift;

  my $existing_peptide_count = 0;

  my %trp_peptides = Trypsin($seq);
  if ( !%trp_peptides ) {
    return $existing_peptide_count;
  } else {
    for my $pep ( values %trp_peptides ) {
      if ( !exists $trpPepDb{$pep} ) {
        $existing_peptide_count++;
      }
    }
  }
  return $existing_peptide_count;
}

# in-silico trypsin digestion

# Trypsin takes a string and cuts it into peptides as fully tryptic peptides
# while allowing for blockage by protline, results returned as either hash or
# hash reference depending on calling context
sub Trypsin {
  # An array reference where each entry is an amino acid
  my $aaAref = shift;

  my %digest;

  # Dont' check for "*", we'll just skip the "*" which should be last entry
  # May be safe to just skip last entry, but let's check in case
  # # always trim off anything after a '*'; the if ( index ...) > -1) is (probably)
  # # speeding things up since most of the time there's no '*' and using the
  # # regex is comparitively slow versus just checking whether the string exists
  # if ( index( $peptide, '*' ) > -1 ) {
  #   $peptide =~ s/$trim_end_regex//xm;
  # }
  
  # replaced by $aaAref my @peptide = split( //, $peptide );
  
  my @peptides;
  my $last_cut_site = 0;
  
  # The first base in the peptide we're building
  my $nextStart = 0;

  my $i = -1;
  my @peptide;
  for my $aa (@$aaAref) {
    $i++;

    if($aa eq '*') {
      last;
    }

    push @peptide, $aa;

    if($aa eq 'K' || $aa eq 'R') {
      push @peptides, \@peptide;
      @peptide = ();

      $nextStart = $i + 1;
    
    }
  }
  for ( my $i = 0 ; $i < @$aaAref ; $i++ ) {
    if ( $aaAref ) {
      push @cut_sites, $i + 1;
      $last_cut_site = $i + 1;
    }
  }

  # if a K or R at the end, don't include it twice but also include the
  # end of the protein
  if ( $last_cut_site < @peptide ) {
    push @cut_sites, scalar @peptide;
  }

  # reset
  $last_cut_site = 0;

  for ( my $i = 0 ; $i < @cut_sites ; $i++ ) {
    my $start = $last_cut_site;
    my $end   = $cut_sites[$i];
    my $seq   = join "", @peptide[ $start .. $end - 1 ];

    if ( length $seq >= $min_peptide_length && length $seq <= $max_peptide_length ) {
      my $coord = sprintf( "%d:%d", $start + 1, $end );
      $digest{$coord} = $seq;
    }

    if ( @peptide == $end ) {
      last;
    }

    if ( BlockCutAa( $peptide[$end] ) ) {

      # any more cut sites?
      if ( $i + 1 <= @cut_sites ) {
        my $end = $cut_sites[ $i + 1 ];
        my $seq = join "", @peptide[ $start .. $end - 1 ];
        if ( length $seq >= $min_peptide_length && length $seq <= $max_peptide_length ) {
          my $coord = sprintf( "%d:%d", $start + 1, $end );
          $digest{$coord} = $seq;
        }
      }
    }
    $last_cut_site = $end;
  }
  if (wantarray) {
    return %digest;
  } elsif ( defined wantarray ) {
    return \%digest;
  } else {
    die "Trypsin() should be called in the list or scalar context";
  }
}

sub CutAa {
  my $aa = shift;

  if ( $aa eq 'K' || $aa eq 'R' ) {
    return 1;
  }
  return;
}

sub BlockCutAa {
  my $aa = shift;

  if ( $aa eq 'P' ) {
    return 1;
  }
  return;
}

sub MakeVarProt {
  my $path_db    = shift;
  my $path_out   = shift;
  my $idListAref = shift;
  my $outDbName  = shift;
  my $okChrsAref = shift;

  my @records;

  my $field_sep_char = ";";
  my $rec_sep_char   = "|";

  # defined acceptable chromosomes, which is really dependent on the binary db
  if ( !defined $okChrsAref ) {
    my @chrs = map { "chr$_" } ( 1 .. 22, 'X', 'Y', 'M' );
    $okChrsAref = \@chrs;
  }

  for my $id (@$idListAref) {
    # read data for each chromosome
    for my $chr (@$okChrsAref) {
      Log("Working on Chr: $chr");
      my $per_rec_aref = ReadPerDb( $path_db, $id, $chr );
      Log( "Read Tx from db: ", scalar @$per_rec_aref );

      for my $per_rec_href (@$per_rec_aref) {
        Log( "Working on", $per_rec_href->{name} );

        # generate all permutations
        my $var_prot_aref = var_prot_for_tx($per_rec_href);
        #say dump( { before => $var_prot_aref } );

        # select the most informative entries
        $var_prot_aref = select_entries($var_prot_aref);
        #say dump({ after => $var_prot_aref });
        for my $r_href (@$var_prot_aref) {
          my $final_href = create_per_prot_rec( $r_href, $per_rec_href, $field_sep_char );
          push @records, $final_href;
        }
      }
    }
  }

  Log("Started writing personal protein entries");
  WritePerProt( $outDbName, $path_out, \@records, $rec_sep_char );
  Log("Finished writing personal protein entries");
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
