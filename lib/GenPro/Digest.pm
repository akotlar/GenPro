use 5.10.0;
use strict;
use warnings;

package GenPro::Digest;

use Mouse 2;
use namespace::autoclean;
use Carp qw/croak/;
use DDP;

has minPeptideLength => (is => 'ro', isa => 'Int', default => 6);
has maxPeptideLength => (is => 'ro', isa => 'Int', default => 40);

my %digestLookups = (
  trypsin => {
    cut => {
      R => 1,
      K => 1,
    },
    block => {
      P => 1,
    }
  }
);

# Don't configure twice; derive state from our capabilities
has enzymes => (is => 'ro', isa => 'ArrayRef', init_arg => undef, default => sub {
  return [ sort {$a cmp $b} keys %digestLookups ];
});

has digestLookups => (is => 'ro', isa => 'HashRef', init_arg => undef, default => sub {
  return \%digestLookups;
});

# If a dbConfig is passed, it should have a
# dbManager and db, where db is a database configuration
# and should be callable as dbManager->dbReadOne($db, $peptide);
# Also, if provided, it will check uniqueness

sub makeDigestFunc {
  my ($self, $type) = @_;

  # Cache because Perl/Mouse accessors slow
  my $minPeptideLength = $self->minPeptideLength;
  my $maxPeptideLength = $self->maxPeptideLength;

  my %cuts = %{$digestLookups{$type}{cut}};
  my %blocks = %{$digestLookups{$type}{block}};

  if(!(%cuts && %blocks)) {
    croak "We don't know how to digest $type";
  }

  return sub {
    my $aaAref = shift;

    if(!@$aaAref) {
      return;
    }

    my @cutSites;

    # Use this to track maximum length of peptide, excluding the "*"
    # and everything after it
    my $effectiveLength = @$aaAref;

    my $i = -1;
    my $start;
    my $seq;
    my @peptides;

    for my $aa (@$aaAref) {
      $i++;

      # TODO: Should we skip stop?
      # Is this identical to $trim_end_regex qr{\*[\*\w]*\z}; $peptide =~ /$trim_end_regex//xm
      # Should be, this removes everything after an observed * (including the *)
      if($aa eq '*') {
        # This tracks the end of the non-NMD protein (effectively trims \*\w*)
        $effectiveLength = $i + 1;
        last;
      }

      if($cuts{$aa}) {
        push @cutSites, $i + 1;
      }
    }

    # if a K or R at the end, don't include it twice but also include the
    # end of the protein
    if( $cutSites[-1] < $effectiveLength) {
      push @cutSites, $effectiveLength;
    }

    # reset, and now use this to create the left bound so that we 
    # may test if our fragment is too short or too long
    my $lastCutSite = 0;
    $i = 0;
    # Iterate over $_ to avoid the confusiong that for my $end (@cutSites)
    # actually copies $end; no, instead it is a mutable reference
    # and that means mutating it will cause major issues (modify @cutSites in the middle of execution
    # This is also the fastest way
    # https://stackoverflow.com/questions/10487316/best-way-to-iterate-through-a-perl-array
    foreach (@cutSites) {
      # In this loop $_ is an element of @cutSites, is 1 past the end of the cut 

      # Don't add + 1 because we take lastCutSite .. $end - 1, so this is already +1 length
      # and because end in the above is $i + 1
      if($_ - $lastCutSite >= $minPeptideLength  && $_ - $lastCutSite <= $maxPeptideLength) {
        push @peptides, join('', @$aaAref[ $lastCutSite .. $_ - 1 ]);
      }

      if ( $_ == $effectiveLength ) {
        last;
      }

      # If the next base is a proline (or other blocker for a different cutter)
      if( $blocks{$aaAref->[$_]} ) {
        # any more cut sites?
        if( $i + 1 < @cutSites) {
          # Cannot re-assign $end, that will modify @cutSites...
          my $tailEnd = $cutSites[ $i + 1 ];

          if($tailEnd - $lastCutSite >= $minPeptideLength && $tailEnd - $lastCutSite <= $maxPeptideLength) {
            push @peptides, join('', @$aaAref[ $lastCutSite .. $tailEnd - 1 ]);
          }
        }
      }

      $lastCutSite = $_;
      $i++;
    }

    return @peptides;
  }
}

__PACKAGE__->meta->make_immutable();
1;