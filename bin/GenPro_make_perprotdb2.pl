#!/usr/bin/env perl
# Name:           GenPro_make_perprotdb2
# Description:    reads the dbs created by `make_perprotdb1.pl` and creates
#                 a single json and fasta record for a specified individual
# Date Created:   Wed Aug 27 21:09:45 2014
# Date Modified:  2017-02-18
# By:             TS Wingo

use 5.10.0;
use strict;
use warnings;

use Cpanel::JSON::XS;
use Data::Dump qw/dump/;
use DB_File;
use Digest::SHA qw/sha1_hex/;
use Getopt::Long;
use Path::Tiny;
use Scalar::Util qw/reftype/;
use GenPro qw/ CleanChr /;

our $VERSION = '0.01';

my $min_peptide_length = 6;
my $max_peptide_length = 40;
my $trim_end_regex     = qr{\*[\*\w]*\z};

my (
  %refprots,     %trpPepDb,   %hIUPAC,    $id,         $dir_db,
  $dir_ref_prot, $dir_out,    $snp_file,  $verbose,    $wantedIdFile,
  $wanted_chr,   $idListFile, $outDbName, $idListAref, $debug,
);

tie %trpPepDb, 'DB_File', Path::Tiny->tempfile->stringify();

#<<< No perltidy
my %codon_2_aa = (
  "AAA" => "K", "AAC" => "N", "AAG" => "K", "AAT" => "N",
  "ACA" => "T", "ACC" => "T", "ACG" => "T", "ACT" => "T",
  "AGA" => "R", "AGC" => "S", "AGG" => "R", "AGT" => "S",
  "ATA" => "I", "ATC" => "I", "ATG" => "M", "ATT" => "I",
  "CAA" => "Q", "CAC" => "H", "CAG" => "Q", "CAT" => "H",
  "CCA" => "P", "CCC" => "P", "CCG" => "P", "CCT" => "P",
  "CGA" => "R", "CGC" => "R", "CGG" => "R", "CGT" => "R",
  "CTA" => "L", "CTC" => "L", "CTG" => "L", "CTT" => "L",
  "GAA" => "E", "GAC" => "D", "GAG" => "E", "GAT" => "D",
  "GCA" => "A", "GCC" => "A", "GCG" => "A", "GCT" => "A",
  "GGA" => "G", "GGC" => "G", "GGG" => "G", "GGT" => "G",
  "GTA" => "V", "GTC" => "V", "GTG" => "V", "GTT" => "V",
  "TAA" => "*", "TAC" => "Y", "TAG" => "*", "TAT" => "Y",
  "TCA" => "S", "TCC" => "S", "TCG" => "S", "TCT" => "S",
  "TGA" => "*", "TGC" => "C", "TGG" => "W", "TGT" => "C",
  "TTA" => "L", "TTC" => "F", "TTG" => "L", "TTT" => "F"
);

$hIUPAC{K}{G} = "T";
$hIUPAC{K}{T} = "G";
$hIUPAC{M}{A} = "C";
$hIUPAC{M}{C} = "A";
$hIUPAC{R}{A} = "G";
$hIUPAC{R}{G} = "A";
$hIUPAC{S}{C} = "G";
$hIUPAC{S}{G} = "C";
$hIUPAC{W}{A} = "T";
$hIUPAC{W}{T} = "A";
$hIUPAC{Y}{C} = "T";
$hIUPAC{Y}{T} = "C";
#>>>

die "Usage: $0
    --r|ref <reference protein db directory>
    --d|db <personal db directory>
    --i|id <sample id>
    --o|out <output directory>
    [--list <file with Ids>]
    [--name <output database name, defaults to id]
    [--max <max peptide length, default: $max_peptide_length>]
    [--min <min peptide length, default: $min_peptide_length>]\n"
unless GetOptions(
  'v|verbose'   => \$verbose,
  'd|db=s'      => \$dir_db,
  'r|ref=s'     => \$dir_ref_prot,
  'o|dir_out=s' => \$dir_out,
  'i|id=s'      => \$id,
  'l|list=s'    => \$idListFile,
  'n|name=s'    => \$outDbName,
  'max=n'       => \$max_peptide_length,
  'min=n'       => \$min_peptide_length,
  'debug'       => \$debug,
)
and $dir_db
and $dir_ref_prot
and $dir_out
and ( ( $idListFile and $outDbName ) or $id );

# check - output directory
my $path_out = path($dir_out);
if ( !$path_out->is_dir() ) {
  $path_out->mkpath();
}

# check - reference proteome location
my $path_ref_prot = path($dir_ref_prot);
if ( !$path_ref_prot->is_dir() ) {
  Log( "Error", sprintf( "'%s' does not exist", $dir_ref_prot ) );
}

# check - personal protein db
my $path_db = path($dir_db);
if ( !$path_db->is_dir() ) {
  Log( "Error", sprintf( "'%s' does not exist", $dir_db ) );
}

# read reference protein database
Log( "Started reading reference protein entries", $path_ref_prot->stringify() );
ReadRefProt($path_ref_prot);
Log( "Finished reading reference protein entries", $path_ref_prot->stringify() );

# make varaint proteins

if ( defined $idListFile ) {
  $idListAref = readIdList($idListFile);
} else {
  $idListAref = [$id];
  if ( !defined $outDbName ) {
    $outDbName = $id;
  }
}
Log(
  sprintf( "Started creating personal protein entries for %s", join ", ", @$idListAref )
);
MakeVarProt( $path_db, $path_out, $idListAref, $outDbName );
Log(
  sprintf(
    "Finished creating personal protein entries for  %s", join ", ", @$idListAref
  )
);

# MakeVarProt takes the path to the personal db and id and creates a personal
# protein database for the id; it optionally takes a list of aceptable
# chromosomes to consider
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

# readIdList takes a file and reads the contents, assuming the 1st column is
# the id of interest. A unique list of ids is made and returned as an array
# reference
sub readIdList {
  my $file = shift;
  my $p    = path($file);

  my %ids;

  if ( !$p->is_file() ) {
    Log( "Fatal", sprintf( "%s does not exist", $p->stringify() ) );
  }

  my $fh = $p->filehandle();
  while (<$fh>) {
    chomp $_;
    my @f = split /\s+/, $_;
    $ids{ $f[0] }++;
  }
  close $fh;
  return [ sort keys %ids ];
}

# create_per_prot_rec takes a personal protein and personal variatn record and
# merges them into a record for printing
sub create_per_prot_rec {
  my ( $per_prot_href, $per_rec_href, $sep_char ) = @_;

  my %m;

  # sync geno with substitutions
  my $href = variant_sites( $per_prot_href, $per_rec_href, $sep_char );
  for my $feature ( keys %$href ) {
    $m{$feature} = $href->{$feature};
  }

  # add basic stuff
  my $prot_record_href = $refprots{ $per_rec_href->{name} };
  for my $feature ( keys %$prot_record_href ) {

    # tx_seq is the reference sequence, which would be confusing in this
    # context so omit it
    if ( $feature eq "tx_seq" ) {
      next;
    }
    $m{$feature} = $prot_record_href->{$feature};
  }

  # this is a variant / personal protein, remove ref_prot key/value which chromosomes
  # from the prot_record_href
  $m{ref_prot} = 0;
  $m{prot_seq} = $per_prot_href->{seq};

  return \%m;
}

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

# var_prot_for_tx takes a variant record and returns an array of variant proteins;
# the caller may select the most parsimonious elements
sub var_prot_for_tx {
  my $variant_rec_href = shift;

  my ( @records, %seq_of_per_prot );

  my $msg = sprintf(
    "considering %d sites in %s",
    scalar @{ $variant_rec_href->{new_aa} },
    $variant_rec_href->{name}
  );
  Log($msg);

  # first key, '-9', is the reference protein
  $seq_of_per_prot{'-9'} = $variant_rec_href->{ref_seq};
  my $aa_residues_aref = $variant_rec_href->{aa_residue};
  my $old_aas_aref     = $variant_rec_href->{old_aa};
  my $new_aas_aref     = $variant_rec_href->{new_aa};

  # here, we are limiting the number of permutations to be done in a reasonable
  # time, recall how large 20!, we will still be making all substitutions on
  # the protein afterward
  if ( scalar @{ $variant_rec_href->{new_aa} } < 21 ) {
    my $last_site = 1;

    # aa_residue is the codon number
    for ( my $i = 0 ; $i < @$aa_residues_aref ; $i++ ) {
      my $site   = $aa_residues_aref->[$i];
      my $new_aa = $new_aas_aref->[$i];

      # which substitutions should be made
      # hash keys are a string list of all sites with substitutions, with -9 meaning
      # reference
      for my $geno ( keys %seq_of_per_prot ) {
        my @sites = split / /, $geno;
        my @seq   = split //,  $seq_of_per_prot{$geno};
        my $old_aa = $seq[ $site - 1 ];

        # check the site exists in the seq
        if ( !defined $old_aa ) {
          my $msg = sprintf(
            "Unable to substitute site %d in %s, expected to change %s for %s",
            $aa_residues_aref->[$i], $variant_rec_href->{name},
            $old_aas_aref->[$i],     $new_aas_aref->[$i]
          );
          Log( "Error", $msg );
        }

        # assume that nonsense mediated decay will remove mRNA that no longer
        # ends in a stop
        if ( $old_aa eq '*' ) {
          next;
        }

        # try to limit the number of permutations
        if ( !CutAa($new_aa)
          && $site - $last_site > $max_peptide_length )
        {
          Log("skipping $new_aa at $site with last site: $last_site");
          next;
        }

        # make substitution
        $seq[ $site - 1 ] = $new_aa;
        my $seq_str = join "", @seq;
        my $new_peptide_count = new_global_peptide($seq_str);

        if ( $new_peptide_count > 0 ) {
          my $new_geno = sprintf( "%s %s", $geno, $aa_residues_aref->[$i] );
          $seq_of_per_prot{$new_geno} = $seq_str;
          my $href = {
            #unique_peptide_count => $new_peptide_count,
            geno => $new_geno,
            seq  => $seq_str,
          };
          push @records, $href;
        } elsif ( $geno eq "-9" ) {
          next;
        } else {
          delete $seq_of_per_prot{$geno};
        }
      }
      $last_site = $site;
    }
  } else {
    my $msg = sprintf(
      "skipping permutations for %s due to %d sites",
      $variant_rec_href->{name},
      scalar @{ $variant_rec_href->{new_aa} }
    );
    Log($msg);
  }

  # make a protein with all substitutions if not already done
  my $all_geno = join " ", "-9", @$aa_residues_aref;
  if ( !exists $seq_of_per_prot{$all_geno} ) {
    Log("making a protein with all genotypes");
    my @seq = split //, $seq_of_per_prot{'-9'};
    for ( my $i = 0 ; $i < @$aa_residues_aref ; $i++ ) {
      my $site   = $aa_residues_aref->[$i];
      my $new_aa = $new_aas_aref->[$i];
      my $old_aa = $seq[ $site - 1 ];

      # check the site exists in the seq
      if ( !defined $old_aa ) {
        my $msg = sprintf(
          "Unable to substitute site %d in %s, expected to change %s for %s",
          $aa_residues_aref->[$i], $variant_rec_href->{name},
          $old_aas_aref->[$i],     $new_aas_aref->[$i]
        );
        Log( "Error", $msg );
      }

      # assume that nonsense mediated decay will remove mRNA that no longer
      # ends in a stop
      if ( $old_aa eq '*' ) {
        next;
      }
      $seq[ $site - 1 ] = $new_aa;
    }
    my $all_geno_seq_str = join "", @seq;
    $seq_of_per_prot{$all_geno} = $all_geno_seq_str;
    #my $new_peptide_count = new_global_peptide($all_geno_seq_str);
    my $href = {
      #unique_peptide_count => $new_peptide_count,
      geno => $all_geno,
      seq  => $all_geno_seq_str,
    };
    push @records, $href;
  }

  if (wantarray) {
    return @records;
  } elsif ( defined wantarray ) {
    return \@records;
  } else {
    die "var_prot_for_tx() should be called in the list or scalar context";
  }
}

# DigestRefProt acts on the global refprots hash and digests  reference protein
# and save the peptide fragments that are between min and max lengths (also
# global variables)
#  - note: will remove the '*' and any trailing sequence before making the digestion
sub DigestRefProt {
  for my $prot_id ( keys %refprots ) {
    my $rec = $refprots{$prot_id};
    my %d   = Trypsin( $rec->{prot_seq} );
    for my $pep ( values %d ) {
      $trpPepDb{$pep} = 1;
    }
  }
}

sub add_seq_to_trp_db {
  my $seq = shift;

  my %d = Trypsin($seq);
  for my $pep ( values %d ) {
    $trpPepDb{$pep} = 1;
  }
}

# ReadRefProt reads a directory that is expected to have all the reference
# proteins (in json format) created by make_refprotdb.pl; it populates a global
# variable, %refprots; it optionally takes a list of acceptable chromosomes
sub ReadRefProt {
  my $path       = shift;
  my $okChrsAref = shift;

  # defined acceptable chromosomes, which is really dependent on the binary db
  if ( !defined $okChrsAref ) {
    my @chrs = map { "chr$_" } ( 1 .. 22, 'X', 'Y', 'M' );
    $okChrsAref = \@chrs;
  }
  my $file_count = 0;

  # read json records and add to overall record of all proteins
  # the format is [ {record => 0}, {record => 1}, ]
  for my $chr (@$okChrsAref) {
    my $file         = $path->child("refProt.$chr.json");
    my $txt          = $file->slurp_raw();
    my $records_aref = decode_json($txt);
    foreach my $r (@$records_aref) {

      # rarely, due to mis-annotation, a protein has no actual sequence; an
      # example is transcript NM_001300891 for PRAMEF9 gene.
      if ( $r->{prot_seq} eq "" ) {
        next;
      }

      # would be nice if we could do this here - downstream consequences though.
      # E.g., malformed proteins may occur - eliminate anything after the '*'
      # the index() is just checking for the '*' in an attempt to speed the
      # processing up for the usual case of no '*' present
      #if ( index( $r->{prot_seq}, '*' ) > -1 ) {
      #  $r->{prot_seq} =~ s/$trim_end_regex//;
      #}
      $refprots{ $r->{id} } = $r;
    }
    $file_count++;
  }
  my $msg = sprintf( "Reference Protein entries: %d, Files: %d",
    scalar keys %refprots, $file_count );
  Log($msg);
}

# ReadPerDb takes a path tiny object of the location to the personal databases,
# an id, and a chromosome and returns all the replacement records for that id
# on that chromosome
sub ReadPerDb {
  my $path = shift;
  my $id   = shift;
  my $chr  = shift;

  my @records;

  my $file = $path->child("$id.$chr.db");
  if ( !$file->is_file() ) {
    Log( "Warn", "Skipping $chr - Cannot find $id.$chr.db" );
    return \@records;
  }

  my %db;
  tie %db, 'DB_File', $file->stringify();

  while ( my ( $tx, $rec_json ) = each %db ) {
    my $recs_aref = decode_json($rec_json);
    push @records, MergePerDbEntriesToRecord( $tx, $recs_aref );
  }
  untie %db;

  return \@records;
}

sub MergePerDbEntriesToRecord {
  my $tx            = shift;
  my $per_recs_aref = shift;

  my %m;
  for my $rec_href (@$per_recs_aref) {
    for my $key ( keys %$rec_href ) {
      push @{ $m{$key} }, $rec_href->{$key};
    }
  }

  # from reference db
  my $prot_record_href = $refprots{$tx};
  if ( !defined $prot_record_href ) {
    Log( "Fatal", "Could not find $tx in ref prot db - did you load refprot?" );
  }

  my $ref_seq    = $prot_record_href->{prot_seq};
  my $ref_strand = $prot_record_href->{strand};

  # get everything in natural order so reverse antisense transcript data
  if ( $ref_strand eq "-" ) {
    for my $key ( keys %m ) {
      $m{$key} = [ reverse @{ $m{$key} } ];
    }
  }

  # add some global stuff
  $m{name}       = $tx;
  $m{ref_seq}    = $ref_seq;
  $m{ref_strand} = $ref_strand;

  checkObsExpAaInTxSeq( \%m );
  return \%m;
}

sub checkObsExpAaInTxSeq {
  my $var_rec = shift;

  my $aa_residue_aref = $var_rec->{aa_residue};
  my $old_aa_aref     = $var_rec->{old_aa};
  my @ref_seq         = split //, $var_rec->{ref_seq};

  for ( my $i = 0 ; $i < @$aa_residue_aref ; $i++ ) {
    # from tx
    my $obs_aa = $ref_seq[ $aa_residue_aref->[$i] - 1 ];
    # from perdb
    my $exp_aa = $old_aa_aref->[$i];

    if ( !defined $obs_aa ) {
      my $msg = sprintf( "No site %d for %s", $var_rec->{name}, $i + 1 );
      Log( "Error", $msg, dump($var_rec), string_protein( $var_rec, $i ) );
    }

    if ( !defined $exp_aa ) {
      my $msg = sprintf( "No expected aa for site %d for %s", $i + 1, $var_rec->{name} );
      Log( "Error", $msg, dump($var_rec), string_protein( $var_rec, $i ) );
    }

    if ( $obs_aa ne $exp_aa ) {
      my $msg =
      sprintf( "Discordant expected '%s' and observed AA '%s' for site %d for %s",
        $obs_aa, $exp_aa, $i + 1, $var_rec->{name} );
      Log( "Error", $msg, dump($var_rec), string_protein( $var_rec, $i ) );
    }
  }
}

# variant_sites takes a personal protein record and personal variant record and
# returns a hash reference with the relevant variant sites enumerated for the
# personal protein record (genomic, coding mRNA, and protein)
sub variant_sites {
  my ( $per_prot_href, $per_rec_href, $sep_char ) = @_;

  my @geno = split / /, $per_prot_href->{geno};
  my %geno_present = map { $_ => 1 } @geno;

  my %m;

  my $aa_residue_aref = $per_rec_href->{aa_residue};
  my $old_aa_aref     = $per_rec_href->{old_aa};
  my $new_aa_aref     = $per_rec_href->{new_aa};
  my $chr_pos_aref    = $per_rec_href->{chr_pos};
  my $ref_allele_aref = $per_rec_href->{ref_allele};
  my $min_allele_aref = $per_rec_href->{min_allele};
  my $codon_pos_aref  = $per_rec_href->{codon_pos};

  for ( my $i = 0 ; $i < @$aa_residue_aref ; $i++ ) {
    my $aa_residue = $aa_residue_aref->[$i];
    if ( exists $geno_present{$aa_residue} ) {
      my $c_site = ( 3 * $aa_residue ) + $codon_pos_aref->[$i];
      push @{ $m{c_changes} },
      sprintf( "c.%d%s>%s", $c_site, $ref_allele_aref->[$i], $min_allele_aref->[$i] );
      push @{ $m{p_changes} },
      sprintf( "p.%d%s>%s", $aa_residue, $old_aa_aref->[$i], $new_aa_aref->[$i] );
      push @{ $m{g_changes} },
      sprintf( "g.%d%s>%s",
        $chr_pos_aref->[$i],
        $ref_allele_aref->[$i],
        $min_allele_aref->[$i] );
    }
  }

  # flatten the arrays into strings
  my %m_str;
  for my $key ( keys %m ) {
    $m_str{$key} = join $sep_char, @{ $m{$key} };
  }
  return \%m_str;
}

# printRec takes a protein entry record and turns it into a fasta an optional
# arrayref of items from the record you don't want printed is passed
#{
#    "ref_prot" : 1
#    "protAcc" : "NM_000451",
#    "spID" : "O15266",
#    "kgID" : "uc004fou.1",
#    "refseq" : "NM_000451",
#    "strand" : "+",
#    "tx_seq" : "ATG...",
#    "sha1" : "48c38981557f4138f85a83b68bfc9696432e00c1",
#    "mRNA" : "NM_000451",
#    "id" : "uc004fou.1",
#    "spDisplayID" : "SHOX_HUMAN",
#    "geneSymbol" : "SHOX",
#    "error_code" : 0,
#    "prot_seq" : "M...",
#    "desc" : "uc004fou.1;NM_000451;O15266;SHOX_HUMAN;SHOX;NM_000451;NM_000451"
# },
sub recToFastaEntry {
  my $href               = shift;
  my $ignoreFeaturesAref = shift;
  my $rec_sep            = shift;

  my ( @array, %iFeature );

  if ( defined $ignoreFeaturesAref ) {
    %iFeature = map { $_ => 1 } @$ignoreFeaturesAref;
  }

  for my $feat ( sort keys %$href ) {
    if ( exists $iFeature{$feat} ) {
      next;
    }
    my $val = $href->{$feat};
    if ( !defined $val ) {
      $val = "NA";
    }
    push @array, sprintf( "%s=%s", $feat, $val );
  }

  # trim off the '*' and beyond
  ( my $seq = $href->{prot_seq} ) =~ s/$trim_end_regex//xm;

  if ( $href->{ref_prot} == 1 ) {
    return (
      $seq,
      sprintf(
        ">%s %s\n%s", $href->{id}, join( $rec_sep, ( "ref_prot", @array ) ), $seq
      )
    );
  } elsif ( $href->{ref_prot} == 0 ) {
    return (
      $seq,
      sprintf(
        ">%s %s\n%s", $href->{id}, join( $rec_sep, ( "var_prot", @array ) ), $seq
      )
    );
  } else {
    my $msg = "Error writing fasta for: " . dump($href);
    croak $msg;
  }
}

# WritePerProt writes a json and fasta record for all personal proteins
sub WritePerProt {
  my $id        = shift;
  my $path_out  = shift;
  my $recs_aref = shift;
  my $sep_char  = shift;

  # write json
  my $json_fh = $path_out->child("$id.perProt.json")->filehandle(">");
  print {$json_fh} encode_json($recs_aref);
  close $json_fh;

  # write fasta
  my $fasta_fh = $path_out->child("$id.perProt.fasta.txt")->filehandle(">");

  my @notInHeader = qw/ id ref_prot prot_seq tx_seq desc sha1 error_code strand /;

  # hash to hold all printed sequence to ensure uniqueness of printed content
  my %prnSeq;
  tie %prnSeq, "DB_File", Path::Tiny->tempfile->stringify();

  # reference protein
  for my $protId ( sort keys %refprots ) {
    my ( $seq, $entry_str ) =
    recToFastaEntry( $refprots{$protId}, \@notInHeader, $sep_char );
    if ( exists $prnSeq{$seq} ) {
      my $msg =
      sprintf( "Exact _reference_ protein previously found in '%s' for '%s'... skipping.",
        $prnSeq{$seq}, $protId );
      Log( "Debug", $msg );
      next;
    }
    # print entry
    say {$fasta_fh} $entry_str;

    # save protein sequence as printed
    $prnSeq{$seq} = $protId;
  }

  # variant protein
  for my $r_href (@$recs_aref) {
    my ( $seq, $entry_str ) = recToFastaEntry( $r_href, \@notInHeader, $sep_char );
    if ( exists $prnSeq{$seq} ) {
      my $msg =
      sprintf( "Exact _variant_ protein previously found in '%s' for '%s'... skipping.",
        $prnSeq{$seq}, $r_href->{id} );
      Log($msg);
      next;
    }
    # print entry
    say {$fasta_fh} $entry_str;

    # save protein sequence as printed
    $prnSeq{$seq} = $r_href->{id};
  }
  untie %prnSeq;
}

# string_protein takes a variatn record hash reference and an optional site
# and returns a formatted reference protein sequence that highlights a specified
# site
sub string_protein {
  my $var_rec_href = shift;
  my $site         = shift;

  if ( !defined $site ) {
    $site = -9;
  }

  my @ref_seq = split //, $var_rec_href->{ref_seq};
  my $msg;

  $msg = sprintf( ">>> Diagnostic Info: '%s', site of interest: '%d', seq len: %d\n",
    $var_rec_href->{name}, $site, scalar @ref_seq );

  $msg .= "'";
  for ( my $j = 0 ; $j < @ref_seq ; $j++ ) {
    if ( $j % 80 == 0 && $j != 0 ) {
      $msg .= "\n";
    }
    if ( $j == $site - 1 ) {
      $msg .= ". ";
    } elsif ( $j == $site + 1 ) {
      $msg .= " .";
    }
    $msg .= $ref_seq[$j];
  }
  $msg .= "'\n\n";
  return $msg;
}

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

# Trypsin takes a string and cuts it into peptides as fully tryptic peptides
# while allowing for blockage by protline, results returned as either hash or
# hash reference depending on calling context
sub Trypsin {
  my $peptide = shift;

  my %digest;

  # always trim off anything after a '*'; the if ( index ...) > -1) is (probably)
  # speeding things up since most of the time there's no '*' and using the
  # regex is comparitively slow versus just checking whether the string exists
  if ( index( $peptide, '*' ) > -1 ) {
    $peptide =~ s/$trim_end_regex//xm;
  }
  my @peptide = split( //, $peptide );
  my @cut_sites;
  my $last_cut_site = 0;
  for ( my $i = 0 ; $i < @peptide ; $i++ ) {
    if ( CutAa( $peptide[$i] ) ) {
      push @cut_sites, $i + 1;
      $last_cut_site = $i + 1;
    }
  }

  # if a K or R at the end, don't include it twice but also include the
  # end of the protein
  if ( $last_cut_site < length $peptide ) {
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

# Log takes and array of strings and uses the first element to handle the
# message reporting - 1) Error, probably a user/input error (print and exit),
# 2) Fatal - probably an internal error, print and croak, 3) Debug is
# something that is printed if the global debug is true, and 4) Warn prints
# a message to standard out regardless of the level of verbosity. Anything else
# will be printed depending on whether the gloabal verbose value is set.
sub Log {
  my $type = shift;
  if ( $type eq 'Error' ) {
    my $msg = join ": ", $type, ( join " ", @_ );
    say STDERR $msg;
    exit(1);
  } elsif ( $type eq 'Fatal' ) {
    my $msg = join ": ", $type, ( join " ", @_ );
    croak $msg;
  } elsif ( $type eq 'Warn' ) {
    my $msg = join ": ", $type, ( join " ", @_ );
    say STDERR $msg;
  } elsif ( $type eq 'Info' ) {
    if ($verbose) {
      my $msg = join ": ", $type, ( join " ", @_ );
      say STDERR $msg;
    }
  } elsif ( $type eq 'Debug' ) {
    if ($debug) {
      my $msg = join ": ", $type, ( join " ", @_ );
      say STDERR $msg;
    }
  } else {
    if ($verbose) {
      my $msg = join " ", $type, @_;
      say STDERR "Info: ", $msg;
    }
  }
}
