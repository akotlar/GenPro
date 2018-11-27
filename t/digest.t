use 5.10.0;
use strict;
use warnings;

use Test::More;

use lib './lib';
use GenPro::Digest;
use DDP;
# TODO: Test that lengths work

my $digest = GenPro::Digest->new({
  minPeptideLength => 5,
  maxPeptideLength => 40,
});

my $fn = $digest->makeDigestFunc('trypsin');

my $prot1 =
  "MWWKQLVAGAVAGAVSRTGTAPLDRLKVFMQVHASKTNRLNILGGLRSMVLEGGIRSLWRGNGINVLKIAPESAIKFMAYEQIKRAILGQQETLHVQERFVAGSLAGATAQTIIYPMETLKNWWLQQYSHDSADPGILVLLACGTISSTCGQIASYPLALVRTRMQAQASIEGGPQLSMLGLLRHILSQEGMRGLYRGIAPNFMKVIPAVSISYVVYENMKQALGVTSR";

my @shorterFragments = sort {$a cmp $b} (
      "QLVAGAVAGAVSR",        "TGTAPLDR",
      "VFMQVHASK",            "LNILGGLR",
      "SMVLEGGIR",            "GNGINVLK",
      "IAPESAIK",             "FMAYEQIK",
      "AILGQQETLHVQER",       "FVAGSLAGATAQTIIYPMETLK",
      "MQAQASIEGGPQLSMLGLLR", "HILSQEGMR",
      "GIAPNFMK",             "VIPAVSISYVVYENMK",
      "QALGVTSR"
);

my @results = sort { $a cmp $b } $fn->([split '', $prot1]);

p @shorterFragments;
p @results;

ok(@shorterFragments == @results, "For max 40 peptide digest of prot1 found " . @results
  . " fragments (\n\tfound @results \n\texpct @shorterFragments\n) ");

for my $i (0 .. $#results ) {
  ok($results[$i] eq $shorterFragments[$i], "Fragment $i matches for prot1: (found $results[$i], expected: $shorterFragments[$i])");
}

my @longerFragments = sort {$a cmp $b} (
  "QLVAGAVAGAVSR",                             "TGTAPLDR",
  "VFMQVHASK",                                 "LNILGGLR",
  "SMVLEGGIR",                                 "GNGINVLK",
  "IAPESAIK",                                  "FMAYEQIK",
  "AILGQQETLHVQER",                            "FVAGSLAGATAQTIIYPMETLK",
  "NWWLQQYSHDSADPGILVLLACGTISSTCGQIASYPLALVR", "MQAQASIEGGPQLSMLGLLR",
  "HILSQEGMR",                                 "GIAPNFMK",
  "VIPAVSISYVVYENMK",                          "QALGVTSR"
);

my $digestMax41 = GenPro::Digest->new({
  minPeptideLength => 6,
  maxPeptideLength => 41,
});

my $fnMax41 = $digestMax41->makeDigestFunc('trypsin');


@results = sort { $a cmp $b } $fnMax41->([split '', $prot1]);

ok(@longerFragments == @results, "For max 41 peptide digest of prot1 found " . @results
  . " fragments (\n\tfound @results \n\texpct @longerFragments\n)");

for my $i (0 .. $#results ) {
  ok($results[$i] eq $longerFragments[$i], "Fragment $i matches for prot1: (found $results[$i], expected: $longerFragments[$i])");
}

my %test = (
  ''                     => { Fragments => [], },
  AAAAVVVRRRTTTTTTTTKRRT => { Fragments => [ "AAAAVVVR", "TTTTTTTTK" ], },
  AAAAARCCCCCRP          => { Fragments => [ "AAAAAR", "CCCCCR", "CCCCCRP" ], },
  AAAARPVVVR             => { Fragments => ["AAAARPVVVR"], },
  AAAARPVVVVR => { Fragments => [ "AAAARPVVVVR", "PVVVVR" ], },
  AAAARPVVVVR => { Fragments => [ "AAAARPVVVVR", "PVVVVR" ], },
);

for my $protein (keys %test) {
  my @digest = sort { $a cmp $b } $fnMax41->([split '', $protein]);
  my @expected = sort { $a cmp $b } @{$test{$protein}{Fragments}};

  ok(@digest == @expected, "For $protein found expected " . @digest . " digests (found @digest, expected @expected)");

  my $i = -1;
  for my $fragment (@digest) {
    $i++;

    ok($fragment eq $expected[$i], "Fragment $i matches for $protein: (found  $fragment, expected $expected[$i])");
  }
}


done_testing();
