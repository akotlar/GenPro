# TODO: Put each database in a separate environemnet, aka folder
# TODO: Check that we don't want stopGain, stopLoss, startGain, or spliceD/A

# or at least put the meta information as a subfolder, say refProteins/meta

use 5.10.0;
use strict;
use warnings;

package GenPro;

our $VERSION = '0.001';

use Mouse 2;
use lib './lib';

use List::Util qw/sum/;
use Types::Path::Tiny qw/AbsFile/;

use namespace::autoclean;

use MCE::Loop;

use GenPro::DBManager;

use GenPro::Digest;
use GenPro::SampleBuilder;
use GenPro::RefProtPept;

use Seq::DBManager;

use Path::Tiny;
use DDP;

use Types::Path::Tiny qw/AbsPath AbsFile AbsDir/;
use Sys::CpuAffinity;

extends 'Seq::Base';

# Bystro db is read-only
has '+readOnly' => (init_arg => undef, default => 1);

# Private
my $metaEnv = 'sample_peptide_meta';
my %metaKeys = (
  samplePeptidesWritten => 'samplePeptidesWritten',
);

my $metaConfig = {stringKeys => 1};

# Taken as arguments / config 
has fileProcessors => ( is => 'ro', isa => 'HashRef', required => 1);
has chromosomes => ( is => 'ro', isa => 'ArrayRef', required => 1);

has input_file => (is => 'ro', required => 1);
has output_file_base => ( is => 'ro', isa => AbsPath, coerce => 1, required => 1,
  handles => { outDir => 'parent', outBaseName => 'basename' });
has temp_dir => (is => 'ro', isa => 'Maybe[Str]');

has geneTrack => ( is => 'ro', isa => 'Str', default => 'refSeq');

has minPeptideLength => (is => 'ro', isa => 'Int', default => 6);
has maxPeptideLength => (is => 'ro', isa => 'Int', default => 40);

has maxThreads => (is => 'ro', isa => 'Int', lazy => 1, default => sub {
  return Sys::CpuAffinity::getNumCpus();
});

# TODO: further reduce complexity
sub BUILD {
  my $self = shift;

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

  # Could allow configuration using a digest key
  $self->{_configuredDigest} = GenPro::Digest->new({
    minPeptideLength => $self->minPeptideLength,
    maxPeptideLength => $self->maxPeptideLength,
  });
}

sub annotate {
  my $self = shift;

  my ($err, $sampleList, $wantedTxNums);

  # TODO: Allow configuration of track, features
  # TODO: extract config needed by makeSampleUniquePeptides?
  my $sampleBuilder = GenPro::SampleBuilder->new({
    inputFile => $self->input_file,
    geneTrack => $self->geneTrack,
    fileProcessors => $self->fileProcessors,
    # TODO: this is clunky (from Seq::Definition)
    output_file_base => $self->output_file_base,
    temp_dir => $self->temp_dir,
    chromosomes => $self->chromosomes,
    # assembly => $self->assembly,
  });

  ( $err, $sampleList, $wantedTxNums ) = $sampleBuilder->go();

  say STDERR "FINISHED STEP 1 (make personal replacement db)";

  my $refBuilder = GenPro::RefProtPept->new({
    digest => $self->{_configuredDigest},
    chromosomes => $self->chromosomes,
    geneTrack => $self->geneTrack,
  });

  $err = $refBuilder->go($wantedTxNums);

  say STDERR "FINISHED STEP 2 (make reference protein and peptide db for requested txNumbers)";

  $err = $self->makeSampleUniquePeptides($sampleList, $wantedTxNums, {
    sample => {
      dbConfig => $sampleBuilder->dbConfig,
      featureDbIdx => $sampleBuilder->featureDbIdx
    },
    refProtein => {
      dbConfig => $refBuilder->refProtDbConfig,
      env => $refBuilder->refProtEnv,
    },
    refPeptide => {
      dbConfig => $refBuilder->refPeptDbConfig,
      env => $refBuilder->refPeptEnv,
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

# TODO: check if number wantedTxNumHref has changed
# and if so, don't go by the sample meta
# TODO: pass in
sub makeSampleUniquePeptides {
  my ($self, $wantedSamplesHref, $wantedTxNumHref, $configs) = @_;

  my $digest = $self->{_configuredDigest};
  my $wantedEnzymesAref = $digest->enzymes;

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

          if(@$userVars > 21) {
            say STDERR "Tx number $txNum on chr $chr has too many variants ( " . @$userVars . ')';
            next;
          }

          my $refSequenceInfo = $personalDb->dbReadOne($refProteinDbs{$chr}, $txNum);

          if(!$refSequenceInfo) {
            die "Couldn't find transcript number $txNum for chr $chr";
          }

          my $uniquePeptidesAref = $txHasUniquePeptidesFn->($userVars, $refSequenceInfo->[0]);

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

  return;
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
  for my $rec (sort { @{$b->[2]} <=> @{$a->[2]} } @$recsAref) {
    # Check each peptide against
    my @uniquePeptides;

    for my $peptide (@{$rec->[2]}) {
      my $seenAref = $dbManager->dbReadOne($db, $peptide);

      if(defined $seenAref) {
        $seenAref->[0] += 1;
        push @{$seenAref->[1]}, $txNum;

        $dbManager->dbPut($db, $peptide, $seenAref);

        next;
      }

      push @uniquePeptides, $peptide;
    }

    if(@uniquePeptides > 0) {
      $rec->[2] = \@uniquePeptides;

      push @finalRecords, $rec;
    }
  }

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

      # TODO: CHECK ON THIS: THOMAS HAS THIS INITIALIZED TO 1, WHICH I THINK
      # MAY CAUSE SOME COMBINATIONS NOT TO BE SEEN (WHEN THE FIRST CODON)
      # IS > maxPeptideLength away
      $lastCodonNum //= $codonNum;

      # We can avoid running this condition if we
      # Both of these are 1-based, so add 1 to be consistent with other
      # parts of the application, or more succinctly choose > instead of >=

      # This exits too early
      # if($codonNum - $lastCodonNum > $maxPeptideLength) {
      #   last;
      # }

      $lastCodonNum = $codonNum;

      my $altAa = $var->[$altAaIdx];

      if(!$cutsHref->{$altAa}) {
        next;
      }

      my $refAa = $var->[$refAaIdx];

      if($refAa ne $seqWithUnique{-9}[1][$codonNum - 1]) {
        die "NOT EQUAL at codon $codonNum of $seqWithUnique{-9}";
      }

      $lastCodonNum = $codonNum;

      for my $i (keys %seqWithUnique) {
        # say STDERR "Checking seqWithUnique key $i";
        # Copy this mutated allele (or reference if first iteration)
        # Mutate that
        my @newSeq = @{$seqWithUnique{$i}[1]};

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
        my @uniquePeptides = $globalUniquePeptideFn->(\@newSeq);

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

        # We keep a reference to the unique peptides because this allows
        # future steps that need to find the unique peptides conditioned
        # on not being in the set of the previous most unique sequences' unique peptides
        # which allows us to dramatically reduce our search space in that step
        push @uniqueRecords, $recAref;

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

      if(@uniquePeptides) {
        push @uniqueRecords, [$varTxAref, $newSeq, \@uniquePeptides];
      }
    }

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

  my $digest = $self->{_configuredDigest};
  my $digestFunc = $digest->makeDigestFunc('trypsin');
  my $cutsHref = $digest->digestLookups->{trypsin}{cut};

  if(!$digestFunc) {
    die 'No digest func for trypsin';
  }

  my $globalUniquePeptideFn = sub {
    my $aaAref = shift;

    my @unique;
    for my $peptide ($digestFunc->($aaAref)) {
       if(!defined $dbManager->dbReadOneRaw($refPeptideTxn, $refPeptideDbi, $peptide)) {
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
