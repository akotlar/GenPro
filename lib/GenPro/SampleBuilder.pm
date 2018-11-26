use 5.10.0;
use strict;
use warnings;

package GenPro::SampleBuilder;

use Mouse 2;
use lib './lib';

use Types::Path::Tiny qw/AbsFile/;

use namespace::autoclean;

use MCE::Loop;

use Seq::Headers;
use Seq::DBManager;
use Seq::Tracks;
use Seq::Tracks::Gene::Definition;

use GenPro::DBManager;

use Path::Tiny;
use DDP;

my $metaEnv = 'samples_meta';
my %metaKeys = (
  samplesWritten => 'samplesWritten',
  wantedTxs => 'wantedTxs',
);

my $metaConfig = {stringKeys => 1};
# We  add a few of our own annotation attributes
# These will be re-used in the body of the annotation processor below
# Users may configure these
has inputFile => (
  is       => 'ro',
  isa      => AbsFile,
  coerce   => 1,
  required => 1,
  handles  => { inputFilePath => 'stringify' },
);

has geneTrack => ( is => 'ro', isa => 'Str', required => 1 );

has chromosomes => ( is => 'ro', isa => 'ArrayRef', required => 1);
# TODO: formalize: check that they have name and args properties
has fileProcessors => ( is => 'ro', isa => 'HashRef', default => 'bystro-vcf' );

# TODO: set the feature names based on the tx track
has features => (is => 'ro', isa => 'ArrayRef', init_arg => undef, lazy => 1, default => sub {
  my $self = shift;

  my $locationFeatures = $self->locationFeatures;
  my $geneTrackFeatures = $self->geneTrackFeatures;
  my $computedFeatures = $self->computedFeatures;

  return [@$locationFeatures, @$computedFeatures, @$geneTrackFeatures];
});

has dbConfig => (is => 'ro', isa => 'HashRef', init_arg => undef, lazy => 1,
default => sub {
  my $self = shift;

  my $nChrs = @{$self->chromosomes};

  return {
    # + 5 because we need 1 global peptide database, and want to leave some headroom,
    maxDbs => $nChrs,
  }
});

has locationFeatures => (is => 'ro', isa => 'ArrayRef', lazy => 1, default => sub {
  return ['chr', 'pos', 'type', 'ref', 'alt'];
});

has computedFeatures => (is => 'ro', isa => 'ArrayRef', lazy => 1, default => sub {
  return ['isHet'];
});

# TODO: allow these to be use the gene track getter defaults
has geneTrackFeatures => (is => 'ro', isa => 'ArrayRef', lazy => 1, default => sub {
  return [
    "txNumber", "name", "name2", "protAcc", "spID", "spDisplayID", "mRNA", "kgID", "ensemblID",
    "exonicAlleleFunction", "strand",
    "codonNumber", "codonPosition", "refCodon", "altCodon", "refAminoAcid", "altAminoAcid"
  ];
});

# A simpler version that doesn't rely on the meta table
# Could use the getFieldDbName instead of the index
has featureDbIdx => (is => 'ro', isa => 'HashRef', init_arg => undef, lazy => 1, default => sub {
  my $self = shift;

  my $trackName = $self->geneTrack;

  my %features;

  my $i = -1;
  for my $feat (@{$self->features}) {
    $i++;

    $features{$feat} = $i;

    my $parentIdx = index($feat, "$trackName.");

    if($parentIdx > -1) {
      $features{substr($feat, $parentIdx + 1)} = $i
    }
  }

  return \%features;
});

# Defines most of the properties that can be configured at run time
# Needed because there are variations of Seq.pm, ilke SeqFromQuery.pm
# Requires logPath to be provided (currently found in Seq::Base)
with 'Seq::Definition', 'Seq::Role::Validator';

# TODO: further reduce complexity
sub BUILD {
  my $self = shift;

  # Pre-build/cache features
  # If we need db access, do before child processes fork
  $self->featureDbIdx;

  # Make our use of Seq::Tracks clear
  my $tracks = Seq::Tracks->new({gettersOnly => 1});
  $self->{_geneGetter} = $tracks->getTrackGetterByName($self->geneTrack);
  $self->{_refGetter} = $tracks->getRefTrackGetter();

  if(!$self->{_geneGetter}) {
    die "Unrecognized track " . $self->geneTrack;
  }

  if(!$self->{_refGetter}) {
    die "Couldn't find reference track ";
  }
}

sub go {
  my $self = shift;

  $self->log( 'info', 'Checking input file format' );

  my $err;
  ( $err, $self->{_chunkSize} ) =
  $self->getChunkSize( $self->inputFile, $self->maxThreads, 512, 5000 );

  if ($err) {
    return ( $err, undef, undef );
  }

  $self->log( 'debug', "chunk size is $self->{_chunkSize}" );

  #################### Validate the input file ################################
  # Converts the input file if necessary
  ( $err, my $fileType ) = $self->validateInputFile( $self->inputFile );

  if ($err) {
    return ( $err, undef, undef );
  }

  # This replicates sub ProcessIds {
  # which just returns an array of ids
  my ($fh, $sampleListAref);
  ($err, $fh, $sampleListAref) = $self->_getFhAndSampleList($fileType);

  if ($err) {
    return ( $err, undef, undef );
  }

  if ( !@$sampleListAref ) {
    $err =
    "Need samples; did you pass \%sampleList% as an argument to the fileProcessor?";

    return ( $err, undef, undef );
  }

  return $self->makePersonalProtDb($fh, $sampleListAref);
}

sub makePersonalProtDb {
  #Inspired by T.S Wingo: https://github.com/wingolab-org/GenPro/blob/master/bin/vcfToSnp
  my ($self, $fh, $sampleListAref) = @_;

  ########################## Write the header ##################################
  my $header = <$fh>;
  $self->setLineEndings($header);

  # Register the output of the intermediate fileProcessor output with GenPro
  my $headers = $self->_setFinalHeader($header);

  my ($completedHref, $recordedTxNumsHref) = $self->_getSampleTodo();

  if ( keys %$completedHref == @$sampleListAref ) {
    return ( undef, $completedHref, $recordedTxNumsHref);
  }

  # What we could possibly annotate
  my %possibleChromosomes = %{ $self->{_refGetter}->normalizedWantedChr };

  my $replacement = $self->{_geneGetter}->siteTypes->{replacement};

  # TODO: make this less clumsy
  my $geneDef = Seq::Tracks::Gene::Definition->new();

  my $geneTrackName = $self->{_geneGetter}->name;

  my $txErrorInIdx = $headers->getFeatureIdx($geneTrackName, $geneDef->txErrorName);
  my $exonicAlleleFuncInIdx = $headers->getFeatureIdx($geneTrackName, $self->{_geneGetter}->exonicAlleleFunctionKey);
  my $hetInIdx = $headers->getFeatureIdx(undef, 'heterozygotes');
  my $homInIdx = $headers->getFeatureIdx(undef, 'homozygotes');

  if ( !( defined $txErrorInIdx && defined $homInIdx && defined $hetInIdx && defined $exonicAlleleFuncInIdx ) ) {
    my $err = "Require 'heterozygotes' and 'homozygotes' and "
      . $self->{_geneGetter}->exonicAlleleFunctionKey . " output ";

    return ( $err, undef, undef );
  }

  my %featureDbIdx = %{$self->featureDbIdx};

  my @geneIdx;
  for my $feat (@{$self->geneTrackFeatures}) {
    # Where the data will be stored once we call the gene track getter
    my $inIdx = $headers->getFeatureIdx($geneTrackName, $feat);
    # Where the data will be stored in our personal variant db
    my $outIdx = $featureDbIdx{$feat};

    if(!(defined $inIdx && defined $outIdx)) {
      die "Expected feat $feat to be in the header: " . $headers->getString();
    }

    push @geneIdx, [$inIdx, $outIdx];
  }

  my $chrIdx = $featureDbIdx{chr};
  my $posIdx = $featureDbIdx{pos};
  my $typeIdx = $featureDbIdx{type};
  my $refIdx =  $featureDbIdx{ref};
  my $altIdx = $featureDbIdx{alt};
  my $txNumberIdx = $featureDbIdx{txNumber};
  my $isHetIdx = $featureDbIdx{isHet};

  # TODO: We don't actually require the names of these to be fixed
  # but maybe we should
  # order is however fixed, so we could use 0 - 4 here instead
  my $chrInIdx  = $headers->getFeatureIdx(undef, 'chrom');
  my $posInIdx  = $headers->getFeatureIdx(undef, 'pos');
  my $typeInIdx = $headers->getFeatureIdx(undef, 'type');
  # We overwrite the reference with 'discordant';
  # So in the header the supplied reference field is replaced with 'discordant';
  my $refInIdx  = $headers->getFeatureIdx(undef, 'discordant');
  my $altInIdx  = $headers->getFeatureIdx(undef, 'alt');

  ######################## Define variables we want COW ########################
  my $sConfig = $self->dbConfig;
  my $nFeatures = @{ $self->features };
  my $geneTrackGetter = $self->{_geneGetter};
  my $refTrackGetter = $self->{_refGetter};

  ######################## Build the fork pool #################################

  # Report every 1e4 lines, to avoid thrashing receiver
  my $abortErr;
  my $progressFunc = $self->_makeDbWriter( $completedHref, $recordedTxNumsHref, \$abortErr );

  MCE::Loop::init {
    max_workers => $self->maxThreads || 8,
    use_slurpio => 1,

    # bystro-vcf outputs a very small row; fully annotated through the alt column (-ref -discordant)
    # so accumulate less than we would if processing full .snp
    chunk_size => $self->{_chunkSize} . "K",
    gather     => $progressFunc,
  };

  mce_loop_f {

    my ( $mce, $slurp_ref, $chunk_id ) = @_;
    #    $_[0], $_[1],     $_[2]
    # open my $MEM_FH, '<', $slurp_ref;
    # binmode $MEM_FH, ':raw';
    open my $MEM_FH, '<', $_[1];
    binmode $MEM_FH, ':raw';

    my $total = 0;

    my ($dataFromDbAref, $zeroPos, $chr, @fields);

    my (%cursors, %dbs);

    my (%seen, %txNumbers);

    # After fork, so that each process opens its environments only once
    # and doesnt share a database handle across threads
    my $personalDb  =  GenPro::DBManager->new();
    my $db = Seq::DBManager->new();

    ####### Read #######
    # The ouptut of the intermediate fileProcessor is expected to be
    # chrom \t pos \t type \t inputRef \t alt \t hets \t heterozygosity \t homs ...

    # For GenPro, we store all replacement sites for a given samples, with
    # some informaiton about the site and the affected transcripts

    READ_LOOP: while ( my $line = $MEM_FH->getline() ) {
      chomp $line;

      @fields = split '\t', $line;
      $total++;

      if ( !$possibleChromosomes{ $fields[0] } ) {
        next;
      }

      $chr     = $fields[$chrInIdx];
      $zeroPos = $fields[$posInIdx] - 1;

      # Caveat: It seems that, per database ($chr), we can have only one
      # read-only transaction; so ... yeah can't combine with dbRead, dbReadOne
      if ( !$cursors{$chr} ) {
        $cursors{$chr} = $db->dbStartCursorTxn($chr);
      }

      $dataFromDbAref = $db->dbReadOneCursorUnsafe( $cursors{$chr}, $zeroPos );

      if ( !defined $dataFromDbAref ) {
        # $self->_errorWithCleanup("Wrong assembly? $chr\: $fields[$posInIdx] not found.");

        # Store a reference to the error, allowing us to exit with a useful fail message
        MCE->gather( 0, 0, "Wrong assembly? $chr\: $fields[$posInIdx] not found." );
        $_[0]->abort();
        return;
      }

      # For now we don't handle indesl
      if ( length( $fields[$altInIdx] ) > 1 ) {
        next;
      }

      # for my $trackIndex (@trackIndicesExceptReference) {
      my @out;

      $geneTrackGetter->get(
        $dataFromDbAref,
        $chr,          # chr
        $fields[$refInIdx],    #ref
        $fields[$altInIdx],    #alt
        0,
        \@out,
        $zeroPos
      );

      if ( $refTrackGetter->get($dataFromDbAref) ne $fields[$refInIdx] ) {
        $self->log( 'info', "Discordant: $chr: $fields[$posInIdx]" );
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
      my $func = $out[$exonicAlleleFuncInIdx][0];

      if ( !defined $func ) {
        next;
      }

      my @ok;
      if ( !ref $func ) {
        # When scalar $out[$txErrorIdx][0] will be an array of errors, or undefined
        if( $func ne $replacement || defined $out[$txErrorInIdx][0]) {
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
            !defined $out[$txErrorInIdx][0][$idx]
          ) {
            push @ok, $idx;
            $found ||= 1;
          }
        }

        if ( !$found ) {
          next;
        }
      }

      # TODO: Optimize allocations
      # flexible features (such as description, not necessary for the Gene class)
      my @record;

      # First 5 records are chr, pos, type, ref, alt; last record is whether it's a het or a hom
      $#record = $nFeatures - 1;

      # TODO: Use the index from the intermediate output
      $record[$chrIdx]  = $fields[0];    #chr
      $record[$posIdx]  = $fields[1];    #pos
      $record[$typeIdx] = $fields[2];    #type (SNP, INS, DEL, etc), could remove
      $record[$refIdx]  = $fields[3];    #ref
      $record[$altIdx]  = $fields[4];    #alt

      # could be more elegant
      # build up the array of results
      if ( @ok == 0 ) {
        for my $fIdx (@geneIdx) {
          # 0 because non-indels have only 1 position index
          $record[$fIdx->[1]] = $out[ $fIdx->[0] ][0];
        }
      } elsif ( @ok == 1 ) {
        my $txIdx  = $ok[0];

        for my $fIdx (@geneIdx) {
          # 0 because non-indels have only 1 position index
          $record[$fIdx->[1]] = $out[ $fIdx->[0] ][0][$txIdx];
        }
      } else {
        my $seenFirst;

        for my $txIdx (@ok) {
          if ( !$seenFirst ) {
            for my $fIdx (@geneIdx) {
              # 0 because non-indels have only 1 position index
              $record[$fIdx->[1]] = [ $out[ $fIdx->[0] ][0][$txIdx] ];
            }

            $seenFirst = 1;
            next;
          }

          for my $fIdx (@geneIdx) {
            # 0 because non-indels have only 1 position index
            push @{ $record[$fIdx->[1]] }, $out[ $fIdx->[0] ][0][$txIdx];
          }
        }
      }

      # Todo use the configured delimiters from Seq::Ouput->delimiters
      my %hets;
      if ( $fields[$hetInIdx] ne '!' ) {
        %hets = map { $_ => 1 } split( ';', $fields[$hetInIdx] );
      }

      my %homs;
      if ( $fields[$homInIdx] ne '!' ) {
        %homs = map { $_ => 1 } split( ';', $fields[$homInIdx] );
      }

      my $cnt = 0;

      # Can no longer use keys % in scalar context
      my $max = %hets + %homs;

      my $het = 1;

      $txNumbers{$chr} //= {};
      my @nums = ref $record[$txNumberIdx] ? @{$record[$txNumberIdx]} : $record[$txNumberIdx];
      for my $num (@nums) {
        $txNumbers{$chr}{$num} //= 1;
      }

      # Optimize, only read over hets + homs
      for my $sample (@$sampleListAref) {
        if ( $cnt == $max ) {
          # say STDERR "count is $cnt and max is $max: $line";
          last;
        }

        # Store whether het or not
        if ( $hets{$sample} ) {
          $record[$isHetIdx] = 1;
        } elsif ( $homs{$sample} ) {
          $record[$isHetIdx] = 0;
        } else {
          next;
        }

        $seen{$sample} //= 1;

        $cnt++;

        # dbPatch with per-op commit seems to work fine

        if(!$dbs{$sample}) {
          $dbs{$sample} = { $chr => $personalDb->_getDbi( $sample, $chr, $sConfig ) };
        } elsif(!$dbs{$sample}{$chr}) {
          $dbs{$sample}{$chr} = $personalDb->_getDbi( $sample, $chr, $sConfig );
        }

        # We don't currently support multiallelics
        $personalDb->dbPut( $dbs{$sample}{$chr}, $zeroPos, \@record );
      }
    }

    MCE->gather( \%seen, \%txNumbers, undef );

    # If we issue personalDb->cleanUp() here, for some reason
    # we can eventually get an 22 EINVAL error
    # maybe too many concurrent openings, or a race condition

    close $MEM_FH;
  }
  $fh;

  MCE::Loop::finish();

  if($abortErr) {
    say STDERR "Abort err";
    p $abortErr;
    return ($abortErr, undef, undef);
  }

  $self->_writeSampleCompleted($completedHref, $recordedTxNumsHref);

  return ( undef, $completedHref, $recordedTxNumsHref );
}

sub _writeSampleCompleted {
  my ($self, $completedHref, $recordedTxNumsHref) = @_;

  my $personalDb = GenPro::DBManager->new();

  my $metaDb = $personalDb->_getDbi( $metaEnv, undef, $metaConfig);
    $personalDb->dbPut($metaDb, $metaKeys{samplesWritten}, $completedHref );
    $personalDb->dbPut($metaDb, $metaKeys{wantedTxs}, $recordedTxNumsHref );
  undef $metaDb;

  $personalDb->cleanUp();
}

sub _getSampleTodo {
  my $personalDb = GenPro::DBManager->new();

  my $metaDb = $personalDb->_getDbi( $metaEnv, undef, $metaConfig );
    my $c = $personalDb->dbReadOne( $metaDb, $metaKeys{samplesWritten} ) || {};
    my $previousTxNums = $personalDb->dbReadOne( $metaDb, $metaKeys{wantedTxs} ) || {};

  # Need to clean up after ourselves
  # Otherwise we'll have a dangling reference to a free'd environment in LMDB
  undef $metaDb;

  # Clear the singleton instance, ensure that threads don't copy any memory
  $personalDb->cleanUp();

  return ($c, $previousTxNums);
}

sub _makeDbWriter {
  my ( $self, $samplesSeenHref, $txNumbersSeenHref, $abortErrRef ) = @_;

  return sub {
    my ($newSeen, $txNumbersHref, $err) = @_;

    if($err) {
      $$abortErrRef = $err;
      return;
    }

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

sub _getFh {
  my ( $self, $type ) = @_;

  my ( $err, $inFh ) = $self->_openAnnotationPipe($type);

  if ($err) {
    return ( $err, undef );
  }

  return ( undef, $inFh);
}

# Open the pipe to the annotation pre-processor, then close it
# All file setup, including generating sample list happens
# Regardless of how many lines of the output we read;
# The sample list will get generated, if requested,
# Even without reading from the stdout filehandle
sub _getFhAndSampleList {
  my ( $self, $type ) = @_;

  my ( $err, $fh ) = $self->_getFh($type);

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

  return ( undef, $fh, \@samples );
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
  $self->_workingDir->child( $self->inputFile->basename . '.file-log.log' );

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
  # idx 2: type
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