#misc. useful code blocks

#process piped output from external command line by line
open(my $fh, '-|', 'powercfg -l') or die $!;
while (my $line = <$fh>) {
    # Do stuff with each $line.
}

#return capture directly to variable
my ($id) =  /(^\S+)/ ;



#bits from http://stackoverflow.com/questions/28443343/needleman-wunsch-algorithm-perl
# creates matrix for needleman-wunsch

use strict;
use warnings;

use List::Util 'max';

STDOUT->autoflush;

#use qw to read in multiple strings
my ($seq1, $seq2) = qw/ ACTTCAATCGGT ACTGGTCAATCGGT /;
#use of 'map'
my ($len1, $len2) = map length, $seq1, $seq2;
#read characters from string into array
my @seq1 = $seq1 =~ /./g;
my @seq2 = $seq2 =~ /./g;

#sequential bracket notation for mutlidimensional arrays (arrays of arrays)
#will work in as many dimensions as storage is available (probably some max, though)
my @matrix;
$matrix[$_][0] = -$_ for 0 .. $len2;
$matrix[0][$_] = -$_ for 0 .. $len1;

for my $x ( 1 .. $len2 ) {     # Rows (bases of sequence 2)
    for my $y ( 1 .. $len1 ) { # Columns (bases of sequence 1)

        my $match = $seq1[$y-1] eq $seq2[$x-1];

        my @scores = (
          $matrix[$x-1][$y] - 1,  # Up
          $matrix[$x][$y-1] - 1,  # Left
          $match ? $matrix[$x-1][$y-1] + 1 : $matrix[$x-1][$y-1] - 1, # Diagonal
        );

        $matrix[$x][$y] = max @scores;
    }        
}

for my $i ( 0 .. $len2 ) {
    my $row = $matrix[$i];
    print join(' ', map { sprintf '%3d', $_ } @$row), "\n";
}

