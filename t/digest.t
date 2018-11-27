use 5.10.0;
use strict;
use warnings;

use Test::More;

use lib './lib';
use GenPro::Digest;
use DDP;
# TODO: Test that lengths work

my $digest = GenPro::Digest->new({
  minPeptideLength => 6,
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

ok(@shorterFragments == @results, "Found expected # of fragments (not including longer one)");

for my $i (0 .. $#results ) {
  ok($results[$i] eq $shorterFragments[$i], "Found $results[$i] ( expected: $shorterFragments[$i] )");
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

ok(@longerFragments == @results, "Found expected # of fragments (including ones longer than 40)");

for my $i (0 .. $#results ) {
  ok($results[$i] eq $longerFragments[$i], "Found $results[$i] (expected: $longerFragments[$i] )");
}


done_testing();
