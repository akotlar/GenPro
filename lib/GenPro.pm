# TODO: Put each database in a separate environemnet, aka folder
# TODO: Check that we don't want stopGain, stopLoss, startGain, or spliceD/A
# TODO: Clean up the sample peptide database after done? So that we can always
# Get unique set... or keep around and assume we'll want to re-check what was previously shared
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

  my $p = $self->output_file_base->stringify();
  $self->{_tsvOut} = "$p.tsv";
  $self->{_fastaOut} = "$p.fasta";

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

  say STDERR "FINISHED STEP 3 (make per-sample database of unique peptides digested by trypsin)";
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
  my %samplePeptideConfig  = (maxDbs => scalar @{$digest->enzymes}, stringKeys => 1);

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

  my ($fastaFn, $tsvHeader) = _makeFastaPrintFn($sampleFeatureDbIdx);

  my $txNumberIdx = $sampleFeatureDbIdx->{txNumber};

  open(my $fastaFh, '>', $self->{_fastaOut});
  open(my $tsvFh, '>', $self->{_tsvOut});

  say $tsvFh $tsvHeader;

  my $progressFunc = sub {
    # my ($fastaStr, $tsvStr) = @_;
    #     $_[0],     $_[1]

    say $fastaFh $_[0];
    say $tsvFh $_[1];
  };

  MCE::Loop::init {
    max_workers => $self->maxThreads || 8,

    # bystro-vcf outputs a very small row; fully annotated through the alt column (-ref -discordant)
    # so accumulate less than we would if processing full .snp
    chunk_size => 'auto',
    gather     => $progressFunc,
  };

  my @wantedChrs = keys %$wantedTxNumHref;
  my $sampleDbConfig = $self->{_sampleDbConfig};

  mce_loop {
    my ($mce, $chunk_ref, $chunk_id) = @_;

    # Thread-local instance
    $personalDb = GenPro::DBManager->new();

    my %refProteinDbs;

    my @finalRecords;
    # my $sample = $_;
    for my $sample (@{ $chunk_ref }) {
      my $samplePeptideDb = $personalDb->_getDbi( "$sample/peptide", 'trypsin', \%samplePeptideConfig );

      my $config = {
        refPeptideConfig => $refPeptideConfig,
        samplePeptideDb => $samplePeptideDb,
      };

      # Something like  20% of the peptides can be  skipped by configuring sampleDb
      # (to check if sample's peptide seen before)
      my ($txHasUniquePeptidesFn, $cleanUp) = $self->_configureVarProtFn($sampleFeatureDbIdx, $config);

      for my $chr (@wantedChrs) {
        $refProteinDbs{$chr} //= $personalDb->_getDbi( $refProtEnv, $chr, \%refProtRdOnlyConfig );

        my $sampleDb = $personalDb->_getDbi( $sample, $chr, \%sampleReadOnlyConfig );

        # At this point, memory usage is 8MB per thread

        # This step eats a few mb for 1M variants on chr1, negligible.
        my $variantsAref = $personalDb->dbReadAll($sampleDb);

        # Go out of scope early,
        undef $sampleDb;

        if(!defined $variantsAref) {
          say STDERR "Found no txs for $sample on chr $chr";
          next;
        }

        # Place all data into
        # txName => [first modification, 2nd modification, etc]
        # Where each modification = contains a reference to the record
        # from the
        # This step seems to eat about 1 more megabyte
        my $userTxHref = _orderByDesired($variantsAref, $txNumberIdx);

        for my $txNum (sort {$a <=> $b} keys %$userTxHref) {
          my $userVars = $userTxHref->{$txNum};

          if(@$userVars > 21) {
            say STDERR "tooManyVariants;txNum:$txNum;chr:$chr;n: " . @$userVars;
            next;
          }

          # This seems to also eat no meaningful memory; memory is still 8MB / thread
          my $refSequenceInfo = $personalDb->dbReadOne($refProteinDbs{$chr}, $txNum);

          if(!$refSequenceInfo) {
            die "Couldn't find transcript number $txNum for chr $chr";
          }

          # Tracking unique peptides at the tx level reduces memory usage
          # by roughly 50-75%; max real memory of ~100MB vs 400 without per thread
          my %seen;
          # Memory usage is essentially flat here too. + 5MB
          my $uniquePeptidesAref = $txHasUniquePeptidesFn->($userVars, $refSequenceInfo->[0], \%seen);

          if(!@$uniquePeptidesAref) {
            next;
          }

          # This replicates select_entries as well as create_per_prot_rec
          # Since create_per_prot_rec is no longer necessary,
          # as $txHasUniquePeptidesFn stores
          # \@{all affected txs}, \@{alt aa sequence}, \@uniquePeptides
          # where uniquePeptides is also modified by sortByUniquePeptides

          # This  also seems to eat no meaningful memory;
          # 1 thread jumps to 22MB, the others between 9 and 16MB; usage grows slowly
          # Then spikes to 63MB; then main process goes up to 800MB..... then lots more
          # Not sure what's happening. Maybe a lot of garbage being generated.
          # Still only 2 threads have a lot of usage;  The first one (624MB) and the 2nd 435MB
          # My first guess is that this is just the memory being read by LMDB, not actual usage
          my @sortedRecs = _sortUniquePeptides($uniquePeptidesAref, $txNum, $personalDb, $samplePeptideDb);

          @finalRecords =  ();
          while(@sortedRecs) {
            my $rec = shift @sortedRecs;
            push @finalRecords, [$rec->[0], $rec->[1]];

            # Note: This is ok, because sortByUnique peptides
            # will remove all peptides previously seen in this sample
            for my $peptide (@{$rec->[2]}) {
              $personalDb->dbPut($samplePeptideDb, $peptide, [1, [$txNum]]);
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
            @sortedRecs = _sortUniquePeptides(\@sortedRecs, $txNum, $personalDb, $samplePeptideDb);
          }

          #  TODO: figure out  if we can really get away with sort step being done per tx
          # This may save  a lot of memory
          if(@finalRecords) {
            $fastaFn->($mce, $sample, \@finalRecords);
          }
        }
      }

      # Clean up the progress function
      $cleanUp->();

      # Reduce number of open file handles
      $personalDb->closeEnv( $sample );
      $personalDb->closeEnv( "$sample/peptide" );
    }

  } @toGetSamples;

  MCE::Loop::finish();

  $personalDb = GenPro::DBManager->new();

  # $metaDb = $personalDb->_getDbi($metaEnv, undef, 1);
  #   $personalDb->dbPut( $metaDb, $metaKeys{refPeptidesWritten}, \%writtenSamples );
  # undef $metaDb;

  $personalDb->cleanUp();

  return;
}


sub _makeFastaPrintFn {
  # The sample variant features
  my $featureDbIdx = shift;

  my $chrIdx = $featureDbIdx->{chrom};
  my $codonNumIdx = $featureDbIdx->{codonNumber};
  my $codonPosIdx = $featureDbIdx->{codonPosition};
  my $posIdx = $featureDbIdx->{pos};
  my $refIdx = $featureDbIdx->{ref};
  my $altIdx = $featureDbIdx->{alt};
  my $altAaIdx = $featureDbIdx->{altAminoAcid};
  my $refAaIdx = $featureDbIdx->{refAminoAcid};

  my @names = ('name2', 'name', 'spID', 'spDisplayID', 'ensemblID', 'kgID');
  my @fIdx = map { $featureDbIdx->{$_} } @names;
  my @idx = 0 .. $#names;

  my $tsvHeader = join("\t", 'seq', 'HGVS_C', 'HGVS_P', 'HGVS_G', @names);

  my $fn = sub {
    my ($mce, $sample, $finalRecordsAref) = @_;

    my @fasta;
    my @tsv;
    for my $rec (@$finalRecordsAref) {
      # my $info = $rec->[0];

      # This assumes that the tx name, etc is same for all modificaitons
      # Which should be sound, as we are supposed to be generating these
      # per transcript number

      my (@c, @p, @g, $altAa);

      for my $r (@{$rec->[0]}) {
        my $cSite = 3 * $r->[$codonNumIdx] + $r->[$codonPosIdx];
        # my $altAa = $r->[$altAaIdx];
        # my $ref
        push @c, sprintf("c.%d%s>%s", $cSite, $r->[$refAaIdx], $r->[$altAaIdx]);
        push @p, sprintf("p.%d%s>%s", $r->[$codonNumIdx], $r->[$refAaIdx],$r->[$altAaIdx]);
        push @g, sprintf("g.%d%s>%s", $r->[$posIdx], $r->[$refIdx], $r->[$altIdx]);
      }

      my $tx = $rec->[0][0];
      my $mods = join(',', map {
        $rec->[$_]
      } 0 .. $#$rec);

      my $c = join(',', @c);
      my $p = join(',', @p);
      my $g = join(',', @g);

      my $header = ">";

      my @data;
      foreach (@idx) {
        if(!defined $tx->[$fIdx[$_]]) {
          push @data, '.';
        } elsif(ref $tx->[$fIdx[$_]]) {
          push @data, join(',', @{$tx->[$fIdx[$_]]});
        } else {
          push @data, $tx->[$fIdx[$_]];
        }
      }

      my $seq = join('', @{$rec->[1]});

      # Single concatenation for efficiency (avoid copying the string array by avoiding .=)
      push @fasta, '>' . join('|', map { $names[$_] . '=' . $data[$_] } @idx ) . "|$c|$p|$g\n$seq";

      push @tsv, join("\t", $seq, $c, $p, $g, @data)
    }

    $mce->gather(join("\n", @fasta), join("\n", @tsv));
  };

  return ($fn, $tsvHeader);
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
sub _sortUniquePeptides {
  my ($recsAref, $txNum, $dbManager, $db) = @_;

  my @finalRecords;
  my $txn = $dbManager->getTxn($db);
  for my $rec (sort { @{$b->[2]} <=> @{$a->[2]} } @$recsAref) {
    # Check each peptide against
    my @uniquePeptides;

    for my $peptide (@{$rec->[2]}) {
      my $seenAref = $dbManager->dbReadOneRaw($txn, $db->{dbi}, $peptide);

      if(defined $seenAref) {
        $seenAref->[0] += 1;
        push @{$seenAref->[1]}, $txNum;

        $dbManager->dbPutRaw($txn, $db->{dbi}, $peptide, $seenAref);
        $txn->commit();
        $txn = $dbManager->getTxn($db);
        next;
      }

      push @uniquePeptides, $peptide;
    }

    if(@uniquePeptides > 0) {
      $rec->[2] = \@uniquePeptides;

      push @finalRecords, $rec;
    }
  }

  $txn->commit();
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

  # Eats a ton of memory...not sure why
  # maybe garbage collection
  # However, it is not storing the unique peptides, or variant info that
  # takes the memory ...
  # Some of it is pushing the sequence reference to @uniqueRecords
  # But we're still using say 1GB per thread...
  # It does get collected though...
  return sub {
    my ($varTxAref, $refAaSeq, $seenHref) = @_;

    my (@uniqueRecords, %seqWithUnique);

    $seqWithUnique{-9} = [ undef, $refAaSeq, undef ];

    my $max = -9;
    my $lastCodonIdx = 0;
    my $longer = 0; my $shorter = 0;

    # TODO: do we want to sort and drop based on the lastCodon condition?
    # for my $var (sort { $a->[$codonNumIdx] <=> $b->[$codonNumIdx] } @$varTxAref) {
    for my $var (sort { $a->[$codonNumIdx] <=> $b->[$codonNumIdx] } @$varTxAref) {
      my $codonIdx = $var->[$codonNumIdx] - 1;

      # This only works if we don't break out of this loop early or skip iterations
      $max += $codonIdx;

      # This saves a large amount of memory (probably garbage) (1/3x)
      # and run time
      # if(!$cutsHref->{$var->[$altAaIdx]} && $codonIdx - $lastCodonIdx > $maxPeptideLength) {
      #   next;
      # }

      if($var->[$refAaIdx] ne $seqWithUnique{-9}[1][$codonIdx]) {
        die "Sample's ref AA != reference amino acid $seqWithUnique{-9}[1][$codonIdx] @ codonIdx: $codonIdx";
      }

      for my $i (keys %seqWithUnique) {
        # Copy this mutated allele (or reference if first iteration)
        # Mutate that
        my @newSeq = @{$seqWithUnique{$i}[1]};

        # This is needed to avoid undefined behavior if an interstitial stop
        # is identified
        if($newSeq[$codonIdx] eq '*') {
          next;
        }

        $newSeq[$codonIdx] = $var->[$altAaIdx];

        my @uniquePeptides = $globalUniquePeptideFn->(\@newSeq, $seenHref);

        if(@uniquePeptides == 0) {
          if($i != -9) {
            # always keep the reference around
            delete $seqWithUnique{$i};
          }

          next;
        }

        if($codonIdx - $lastCodonIdx > $maxPeptideLength) {
          $longer++;
        }  else {
          $shorter++;
        }
        #  say "HAD: current $codonIdx, last: $lastCodonIdx";
        # p $var;

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
        $seqWithUnique{$i + $codonIdx} = $recAref;
      }

      $lastCodonIdx = $codonIdx;
    }

    # Many times we'll never not have $max
    if(!$seqWithUnique{$max}) {
      # Mutates the reference sequence copy we made
      # This becomes our fully-mutated copy
      my @newSeq = @{$seqWithUnique{-9}[1]};

      for my $var (@$varTxAref) {
        $newSeq[ $var->[$codonNumIdx] - 1 ] = $var->[$altAaIdx];
      }

      my @uniquePeptides = $globalUniquePeptideFn->(\@newSeq);

      if(@uniquePeptides) {
        push @uniqueRecords, [$varTxAref, \@newSeq, \@uniquePeptides];
      }
    }
    undef %seqWithUnique;

    # if($longer > 0|| $shorter > 0) {
    #   say "Longer %: " . ($longer / ($longer + $shorter));
    # }

    return \@uniqueRecords;
  }
}

# TODO: add other enzymes
sub _configureVarProtFn {
  my $self = shift;
  # The sample variant features
  my $featureDbIdx = shift;

  my $config = shift;

  my $refPeptideConfig = $config->{refPeptideConfig};

  my $sampleDb = $config->{samplePeptideDb};

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

  if(!($digestFunc && $cutsHref)) {
    die 'No digest func for trypsin';
  }
  my $skipped = 0;
  my $kept = 0;
  # No effect on memory usage
  my $globalUniquePeptideFn = sub {
    my ($aaAref, $seenLocally) = @_;

    my $sampleTxn = $dbManager->getTxn($sampleDb, $dbManager->MDB_RDONLY);

    if(!$sampleTxn) {
      die "Couldn't open txn for sample";
    }

    my @unique;
    for my $peptide ($digestFunc->($aaAref)) {
       if(!defined $dbManager->dbReadOneRaw($refPeptideTxn, $refPeptideDbi, $peptide)) {
          if(($seenLocally && $seenLocally->{$peptide})
          || defined $dbManager->dbReadOneRaw($sampleTxn, $sampleDb->{dbi}, $peptide)) {
            $skipped++;
            next
          }

          push @unique, $peptide;
          $seenLocally->{$peptide} = 1;
          $kept++;
      }
    }

    $sampleTxn->abort();

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

  # my %total;
  my $cleanUp = sub {
    # $total{$sampleDb->{envName}} //= 1;

    say STDERR "SKIPPED for sample $sampleDb->{envName}: $skipped (kept $kept)";
    $refPeptideTxn->abort();
  };

  return ($txHasUniquePeptidesFn, $cleanUp);
}

__PACKAGE__->meta->make_immutable;

1;
