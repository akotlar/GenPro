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
has maxDel =>
  ( is => 'ro', isa => 'Int', default => -32, writer => 'setMaxDel' );

# TODO: formalize: check that they have name and args properties
has fileProcessors => ( is => 'ro', isa => 'HashRef', default => 'bystro-vcf' );

# Defines most of the properties that can be configured at run time
# Needed because there are variations of Seq.pm, ilke SeqFromQuery.pm
# Requires logPath to be provided (currently found in Seq::Base)
with 'Seq::Definition', 'Seq::Role::Validator';

has '+readOnly' => ( init_arg => undef, default => 0 );

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

# TODO: For now we only accept tab separated files
# We could change this, although comma separation causes may cause with our fields
# And is slower, since we cannot split on a constant
# We would also need to take care with properly escaping intra-field commas
# $err = $self->setDelimiter($firstLine);

  # if($err) {
  #   $self->_errorWithCleanup($err);
  #   return ($err, undef);
  # }

  # Calling in annotate allows us to error early
  my $err;
  ( $err, $self->{_chunkSize} ) =
    $self->getChunkSize( $self->input_file, $self->maxThreads, 512, 16384 );

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
  $self->_errorWithCleanup(
    "File type isn\'t vcf or snp. Please use one of these files");
  return ( "File type isn\'t vcf or snp. Please use one of these files",
    undef );
}

sub annotateFile {

#Inspired by T.S Wingo: https://github.com/wingolab-org/GenPro/blob/master/bin/vcfToSnp
  my $self = shift;
  my $type = shift;

  GenPro::DBManager::initialize(
    {
      databaseDir => path( $self->database_dir )->child("genpro")->stringify()
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

  my $db  = $self->{_db};
  my $idx = 0;
  my %sampleIdx;

  # How many sub-databases we need
  my $numSamples = @$sampleListAref;

  # my ( $self, $name, $dontCreate, $stringKeys, $namedDb, $maxDbs ) = @_;
  for my $sample ( @$sampleListAref[ 0 .. 100 ] ) {

    # $sampleIdx{$sample} = $idx;
    # say STDERR "Working on sample $sample";

    # $db->_getDbi("$sample");

    for my $chr ( @{ $self->chromosomes } ) {

      $personalDb->_getDbi( $sample, 0, 0, $idx, $numSamples );
      $idx++;
    }

  }

  ########################## Write the header ##################################
  my $header = <$fh>;
  $self->setLineEndings($header);

  ######################## Build the fork pool #################################
  my $abortErr;

  my $messageFreq = ( 2e4 / 4 ) * $self->maxThreads;

  # Report every 1e4 lines, to avoid thrashing receiver
  my $progressFunc =
    $self->makeLogProgressAndPrint( \$abortErr, $outFh, $statsFh,
    $messageFreq );
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

  mce_loop_f {

    #my ($mce, $slurp_ref, $chunk_id) = @_;
    #    $_[0], $_[1],     $_[2]
    #open my $MEM_FH, '<', $slurp_ref; binmode $MEM_FH, ':raw';
    open my $MEM_FH, '<', $_[1];
    binmode $MEM_FH, ':raw';

    my $total = 0;

    my @indelDbData;
    my @indelRef;
    my @lines;
    my $dataFromDbAref;
    my $zeroPos;
    my $chr;
    my @fields;

    my %cursors = ();

    # Write cursors...stored as user/chr
    my %wCursors = ();

    my $out = [];

    # Each line is expected to be
    # chrom \t pos \t type \t inputRef \t alt \t hets \t homozygotes \n
    # the chrom is always in ucsc form, chr (the golang program guarantees it)

    # For GenPro, we need the following record written to LMDB
    # my $record_href = {
    #     aa_residue => $aa_residue,
    #     old_aa     => $codon_2_aa{$codon},
    #     new_aa     => $codon_2_aa{$new_codon},
    #     chr        => $chr,
    #     chr_pos    => $pos,
    #     codon_pos  => $codon_pos,
    #     ref_allele => $ref_allele,
    #     min_allele => $min_allele,
    # };

    # WriteToDb( $id, $chr, $transcript_name, $record_href );

    # sub WriteToDb {
    #     my ( $id, $chr, $key, $href ) = @_;
    #     my $db_href = $db{$id}{$chr};

    #     if ( !defined $db_href ) {
    #         my $msg = sprintf( "Error - no db for %s %s", $id, $chr );
    #         Log( "Fatal", $msg );
    #     }
    #     my $val = $db_href->{$key};
    #     if ( !defined $val ) {
    #         $db_href->{$key} = encode_json( [$href] );
    #     }
    #     else {
    #         my $recs_aref = decode_json($val);
    #         push @$recs_aref, $href;
    #         $db_href->{$key} = encode_json($recs_aref);
    #     }
    #     }

# We then want to store all replacement sites for a given id,
# which I think is the refseq id
# # ReadPerDb takes a path tiny object of the location to the personal databases,
# an id, and a chromosome and returns all the replacement records for that id
# on that chromosome

# push @records, MergePerDbEntriesToRecord( $tx, $recs_aref ); in make_perprotdb2
# this I think is because # my $file = $path->child("$id.$chr.db");...
# each transcript id (I think this is the id), has its own database
# and MergePerDbEntriesToRecord combines that, which I think is the personal db
# with the reference db stuff

    # No! The ids are samples

    #
    # sub ProcessIds {
    # my ( $fieldsAref, $wantedIdHref, $altNameHref ) = @_;

    # my @ids;

   # # using the snpfile header line to get ids - format 'id\t\tid' and we split
   # # using '\t' in ReadSnpFile
    while ( my $line = $MEM_FH->getline() ) {
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
        $self->_errorWithCleanup(
          "Wrong assembly? $chr\: $fields[1] not found.");

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
      my $out = [];

      $geneTrackGetter->get(
        $dataFromDbAref,
        $chr,          # chr
        $fields[3],    #ref
        $fields[4],    #alt
        0,
        $out,
        $zeroPos
      );

      # TODO: Write hash of variant samples
      # TODO: Write hash of variant sample
      # Or removel imit on number of databases
      # $fields[$refTrackOutIdx][0] = $refTrackGetter->get($dataFromDbAref);

      if ( $refTrackGetter->get($dataFromDbAref) ne $fields[3] ) {
        $self->log( 'info', "Discordant: $chr: $fields[1]" );
        next;
      }

      # We need the following fields
      # my $chr  = $data{Fragment}; # chr
      # my $pos  = $data{Position}; # pos
      # my $ref  = $data{Reference}; # ref
      # my $type = $data{Type}; # type
      # $ref_for_site{$chr}{$pos} = $ref;
      # and we want to return / store
      # if ( $prob >= 0.95 ) {
      #   my $min_allele = $hIUPAC{$geno}{ $data{Reference} };
      #   if ( defined $min_allele ) {
      #     $ids{ $ids[$i] }++;
      #     push @{ $variant_sites{$chr}{$pos}{ids} },        $ids[$i];
      #     push @{ $variant_sites{$chr}{$pos}{min_allele} }, $min_allele;
      #   }
      # }
      # @ids = map { $_ } sort { $a cmp $b } ( keys %ids );
      # return ( \%variant_sites, \%ref_for_site, \@ids );

      # my ( $variant_site_href, $ref_for_site_href, $ids_aref ) =
      # ReadSnpFile( $snp_file, $wanted_chr, $wantedIdHref, $altIdForId_href );

# path_bin is path to the parent of the .idx offset files $db_name.$chr.idx
# ReplacementSites( $path_bin, $path_out, $variant_site_href, $ref_for_site_href, $ids_aref );
# Then we check, for each site, whether it's a repalcement
# If so, we store:
# my (
#   @codon,      @new_codon, $strand,     $gene_symbol, $codon_pos,
#   $codon_code, $codon,     $aa_residue, $codon_number,
# );
#

      # We drop everything that is discordant

      # We write:

#          if ( defined $new_base and defined $aa_residue and defined $strand ) {
#               @new_codon = @codon;
#               $new_codon[ $codon_pos - 1 ] = $new_base;
#               my $new_codon = join "", @new_codon;
#               my $new_aa    = $codon_2_aa{$new_codon};
#               my $old_aa    = $codon_2_aa{$codon};

#               # save the site for certain ids
#               if ( defined $new_aa and defined $old_aa and $new_aa ne $old_aa ) {
#                 # cycle through the list of IDs and minor alleles and save the transcript and
#                 # variant info for each individual with variants in the particular transcripts
#                 #                                                   AA_residue    original_AA         new_AA
#                 # e.g., push @{ r_sites{$id}{$transcript_name} }, "$aa_residue:$codon2_aa{$codon}:$codon_2_aa{$new_codon}";
#                 for ( my $i = 0; $i < @$min_alleles_aref; $i++ ) {
#                   my $id        = $these_ids_aref->[$i];
#                   my $id_allele = $min_alleles_aref->[$i];
#                   if ( $id_allele eq $min_allele ) {
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

      # Essentially, we write everything that is a repalcement

# This is what we get out of the gene track:
#  $outAccum->[$idxMap->{$feature}][$posIdx] = [map { $geneDb->[$_]{$cachedDbNames->{$feature}} } @$txNumbers];
# $outAccum->[$self->{_codonPosFidx}][$posIdx] = \@codonPos;
# $outAccum->[$self->{_codonNumFidx}][$posIdx] = \@codonNum;
# $outAccum->[$self->{_alleleFuncFidx}][$posIdx] = \@funcAccum;
# $outAccum->[$self->{_refAaFidx}][$posIdx] = \@refAA;
# $outAccum->[$self->{_altAaFidx}][$posIdx] = \@newAA;
# $outAccum->[$self->{_codonSidx}][$posIdx] = \@codonSeq;
# $outAccum->[$self->{_altCodonSidx}][$posIdx] = \@newCodon;

      # ...

# push @{ $outAccum->[$self->{_strandFidx}][$posIdx] }, $site->[$strandIdx];
#     push @{ $outAccum->[$self->{_siteFidx}][$posIdx] }, $site->[$siteTypeIdx];
# we want all of this stuff... per sample

      # and we also want chr, position, and the transcript name
      # push @lines, \@fields;

      my $dbPath;
      for my $sample (@$sampleListAref) {

        # $dbPath = "$sample/$chr";
        # p $dbPath;
        $wCursors{$dbPath} //= $db->dbStartCursorTxn($sample);

        $db->dbPutCursorUnsafe( $wCursors{$dbPath},
          $dbPath, 0, $zeroPos, $out );
      }

    }

    close $MEM_FH;

    # if (@lines) {
    #     MCE->gather(
    #       scalar @lines,
    #       $total - @lines,
    #       undef, $outputter->makeOutputString( \@lines )
    #     );
    # }
    # else {
    #     MCE->gather( 0, $total );
    # }

  }
  $fh;

  # Force flush
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
    $self->_workingDir->child( $self->outputFilesInfo->{sampleList} )
    ->stringify();

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
  p $args;

  # TODO:  add support for GQ filtering in vcf
  my $err =
    $self->safeOpen( $fh, '-|', "$echoProg $inPath | $args 2> $errPath" );

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

sub _getFinalHeader {
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
