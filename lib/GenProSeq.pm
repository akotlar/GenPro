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
use Seq::Headers;

use lib './lib';
use GenPro::DBManager;
use Seq::DBManager;
use Path::Tiny;
use Scalar::Util qw/looks_like_number/;
use DDP;
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
  p $self->{_chunkSize};
  sleep(1);
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

  if ( $fileType eq 'snp' ) {
    $self->log( 'info', 'Beginning annotation' );
    return $self->annotateFile('snp');
  }

  if ( $fileType eq 'vcf' ) {
    $self->log( 'info', 'Beginning annotation' );
    return $self->annotateFile('vcf');
  }

  # TODO: Inspect vcf header

  # TODO: support any other file, by checking the extension

  # TODO: we don't really check for valid vcf, just assume it is
  # So this message is never reached
  $self->_errorWithCleanup("File type isn\'t vcf or snp. Please use one of these files");
  return ( "File type isn\'t vcf or snp. Please use one of these files", undef );
}

sub annotateFile {

  #Inspired by T.S Wingo: https://github.com/wingolab-org/GenPro/blob/master/bin/vcfToSnp
  my $self = shift;
  my $type = shift;

  # TODO: Convert genpro to Bystro track, give it its own db folder
  GenPro::DBManager::initialize(
    {
      databaseDir => path( $self->database_dir )->parent()->child("genpro")->stringify()
    }
  );

  my $personalDb = GenPro::DBManager->new();

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

  # Create the databases, 1 database per sample,
  # and each sample database was 25 tables, 1 per chromosome
  my $numChr = @{ $self->chromosomes };
  for my $sample (@$sampleListAref) {
    for my $chr ( @{ $self->chromosomes } ) {
      $personalDb->_getDbi( $sample, 0, 0, $chr, $numChr );
    }
  }

  $personalDb->cleanUp();

  ########################## Write the header ##################################
  my $header = <$fh>;
  $self->setLineEndings($header);

  # Register the output of the intermediate fileProcessor output with GenPro
  $self->_setFinalHeader($header);
  ######################## Build the fork pool #################################
  my $abortErr;

  my $messageFreq = ( 2e4 / 4 ) * $self->maxThreads;

  # Report every 1e4 lines, to avoid thrashing receiver
  my $progressFunc =
  $self->makeLogProgressAndPrint( \$abortErr, $outFh, $statsFh, $messageFreq );
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
  my $headerIndices = Seq::Headers::getParentIndices();

  my $hetIdx = $headerIndices->{heterozygotes};
  my $homIdx = $headerIndices->{homozygotes};

  if ( !( defined $homIdx && defined $hetIdx ) ) {
    my $err = "Require 'heterozygotes' and 'homozygotes' output for $type fileProcessor";
    $self->_errorWithCleanup($err);
    return ( $err, undef );
  }

  my $replacement = $geneTrackGetter->siteTypes->{replacement};

  my $exonicAlleleFuncIdx = $geneTrackGetter->{_alleleFuncFidx};
  my $codonNumIdx         = $geneTrackGetter->{_codonNumFidx};
  my $refAaIdx            = $geneTrackGetter->{_refAaFidx};
  my $altAaIdx            = $geneTrackGetter->{_altAaFidx};
  my $descIdx             = $geneTrackGetter->{_featureIdxMap}{description};
  my $nameIdx             = $geneTrackGetter->{_featureIdxMap}{name};
  my $name2idx            = $geneTrackGetter->{_featureIdxMap}{name2};
  my $codonPosIdx         = $geneTrackGetter->{_codonPosFidx};

  # Defines feature insertion order
  # This defines what our per-user record looks like at any given position
  # Except we add 5 fields to beginning (chr, pos, type, ref, alt)
  # And 1 field to end (het/hom status
  my @wanted = (
    $exonicAlleleFuncIdx, $codonNumIdx, $refAaIdx, $altAaIdx,
    $descIdx, $nameIdx, $name2idx, $codonPosIdx
  );

  mce_loop_f {

    #my ($mce, $slurp_ref, $chunk_id) = @_;
    #    $_[0], $_[1],     $_[2]
    #open my $MEM_FH, '<', $slurp_ref; binmode $MEM_FH, ':raw';
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

    my %cursors = ();

    ####### Read #######
    # The ouptut of the intermediate fileProcessor is expected to be
    # chrom \t pos \t type \t inputRef \t alt \t hets \t heterozygosity \t homs ...

    # For GenPro, we store all replacement sites for a given samples, with
    # some informaiton about the site and the affected transcripts

    while ( my $line = $MEM_FH->getline() ) {
      chomp $line;

      @fields = split '\t', $line;
      $total++;

      if ( !$wantedChromosomes{ $fields[0] } ) {
        next;
      }

      # if ( $chr ne $fields[0] ) {
      #   $personalDb->dbEndCursorTxn($chr);
      #   delete $wCursors{$chr};
      # }

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

      # Check that we don't want stopGain, stopLoss, startGain, or spliceD/A
      if ( !defined $func || ( !ref $func && $func ne $replacement ) ) {
        next;
      }

      my $found = 0;
      my @ok;
      for ( my $idx = 0 ; $idx < @$func ; $idx++ ) {
        if ( !defined $func->[$idx] || $func->[$idx] ne $replacement ) {
          push @ok, $idx;
          next;
        }

        $found ||= 1;

        # $idx++;
      }

      if ( !$found ) {
        next;
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

      # WriteToDb( $id, $chr, $transcript_name, $record_href );
      #                   }
      #                 }
      #               }
      #             }
      #             $dat_record_len -= $dat_typedef_size;
      #           }
      #         }
      #       }
      #     }
      #   }

      # TODO: Optimize allocations
      # flexible features (such as description, not necessary for the Gene class)
      my @record;

      # First 4 record; last record is whether it's a het or a hom
      $#record = $#wanted + 5;

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
        my $txIdx = @ok;

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
      for my $sample (@$sampleListAref) {
        if ( $cnt == $max ) {
          say STDERR "count is $cnt and max is $max: $line";
          last;
        }

        # Store whether het or not
        if ( $hets{$sample} ) {
          $record[-1] = 0;
        } elsif ( $homs{$sample} ) {
          $record[-1] = 1;
        } else {
          next;
        }

        $cnt++;

        $personalDb->dbPatch( $sample, $chr, 0, $zeroPos, \@record );

        # Uncomment to test
        # my $val = $personalDb->dbReadOne( $sample, $chr, $zeroPos );
        # p $val;
      }

      $overallCnt += $cnt;
    }

    $personalDb->cleanUp();
    close $MEM_FH;


  }
  $fh;

  # Force flush of output
  $progressFunc->( 0, 0, undef, undef, 1 );

  MCE::Loop::finish();

  # Unfortunately, MCE::Shared::stop() removes the value of $abortErr
  # according to documentation, and I did not see mention of a way
  # to copy the data from a scalar, and don't want to use a hash for this alone
  # So, using a scalar ref to abortErr in the gather function.
  if ($abortErr) {
    say "Aborted job due to $abortErr";

    # Database & tx need to be closed
    $db->cleanUp();

    return ( 'Job aborted due to error', undef );
  }

  ################ Finished writing file. If statistics, print those ##########
  # Sync to ensure all files written
  # This simply tries each close/sync/move operation in order
  # And returns an error early, or proceeds to next operation
  $err =
     $self->safeClose($outFh)
  || ( $statsFh && $self->safeClose($statsFh) )
  || $self->safeSystem('sync')
  || $self->_moveFilesToOutputDir();

  if ($err) {
    $self->_errorWithCleanup($err);
    return ( $err, undef );
  }

  $db->cleanUp();

  return ( $err, $self->outputFilesInfo );
}

sub makeLogProgressAndPrint {
  my ( $self, $abortErrRef, $outFh, $statsFh, $throttleThreshold ) = @_;

  my $totalAnnotated = 0;
  my $totalSkipped   = 0;

  my $publish = $self->hasPublisher;

  my $thresholdAnn     = 0;
  my $thresholdSkipped = 0;

  if ( !$throttleThreshold ) {
    $throttleThreshold = 1e4;
  }
  return sub {

    #<Int>$annotatedCount, <Int>$skipCount, <Str>$err, <Str>$outputLines, <Bool> $forcePublish = @_;
    ##    $_[0],          $_[1]           , $_[2],     $_[3].           , $_[4]
    if ( $_[2] ) {
      $$abortErrRef = $_[2];
      return;
    }

    if ($publish) {
      $thresholdAnn     += $_[0];
      $thresholdSkipped += $_[1];

      if ( $_[4]
        || $thresholdAnn + $thresholdSkipped >= $throttleThreshold )
      {
        $totalAnnotated += $thresholdAnn;
        $totalSkipped   += $thresholdSkipped;

        $self->publishProgress( $totalAnnotated, $totalSkipped );

        $thresholdAnn     = 0;
        $thresholdSkipped = 0;
      }
    }

    if ( $_[3] ) {
      if ($statsFh) {
        print $statsFh $_[3];
      }

      print $outFh $_[3];
    }
  }
}

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
