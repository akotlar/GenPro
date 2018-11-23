use 5.10.0;
use strict;
use warnings;

package GenPro::Digest;

use Mouse 2;
use namespace::autoclean;
use Carp qw/croak/;

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

has digestLookups => (is => 'ro', isa => 'HashRef', init_arg => undef, default => sub {
  return \%digestLookups;
});

# If a dbConfig is passed, it should have a
# dbManager and db, where db is a database configuration
# and should be callable as dbManager->dbReadOne($db, $peptide);
# Also, if provided, it will check uniqueness
sub makeDigestFunc {
  my ($self, $type, $dbConfig) = @_;

  # if($digestFuncs{$type}) {
  #   return $digestFuncs{$type};
  # }

  # Cache because Perl/Mouse accessors slow
  my $minPeptideLength = $self->minPeptideLength;
  my $maxPeptideLength = $self->maxPeptideLength;

  my %cuts = $digestLookups{$type}{cut};
  my %blocks = $digestLookups{$type}{block};

  if(!(%cuts && %blocks)) {
    croak "We don't know how to digest $type";
  }

  # my ($dbManager, $db);

  # if($dbConfig) {
  #   $dbManager = $dbConfig->{dbManager};
  #   $db = $dbConfig->{db};

  #   if(!$dbManager && $db) {
  #     die "Expected dbManager and db keys in dbConfig arg passed to makeDigestFunc";
  #   }
  # }

  return sub {
    my $aaAref = shift;
  p $aaAref;
    my @cutSites;
    my $lastCutSite = 0;

    my $i = -1;
    my $start;
    my $seq;
    my @peptides;
    for my $aa (@$aaAref) {
      $i++;

      if($aa eq '*') {
        last;
      }

      if($cuts{$aa}) {
        push @cutSites, $i + 1;
        $lastCutSite = $i + 1;
      }
    }

    # if a K or R at the end, don't include it twice but also include the
    # end of the protein
    if( $lastCutSite < @$aaAref) {
      push @cutSites, scalar @$aaAref;
    }

    # reset
    $lastCutSite = 0;
    
    $i = -1;
    for my $end (@cutSites) {
      $i++;

      # Don't add + 1 because we take lastCutSite .. $end - 1, so this is already +1 length
      # and because end in the above is $i + 1
      if($end - $lastCutSite >= $minPeptideLength && $end - $lastCutSite <= $maxPeptideLength) {
        push @peptides, join('', @$aaAref[ $lastCutSite .. $end - 1 ]);
      }

      if ( $end == @$aaAref ) {
        last;
      }

      # If the next base is a proline that blocks
      if( $blocks{$aaAref->[$end]} ) {
        # any more cut sites?
        if( $i + 1 < @cutSites) {
          $end = $cutSites[ $i + 1 ];

          if($end - $lastCutSite >= $minPeptideLength && $end - $lastCutSite <= $maxPeptideLength) {
            push @peptides, join('', @$aaAref[ $lastCutSite .. $end - 1 ]);
          }
        }
      }

      $lastCutSite = $end;
    }

    return @peptides;
  }

  # return $digestFuncs{$type};
}

__PACKAGE__->meta->make_immutable();
1;