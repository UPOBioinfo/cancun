#!/usr/bin/perl
use strict;

my $F = $ARGV[0];
my $N = 81981;
my %pos;

for (my $i = 1; $i <= $N; $i++) {
  $pos{$i} = 0;
}

open in, $F;
while (<in>) {
  chomp;

  my ($p1, $p2) = (split/\t/)[5,6];
  ($p1, $p2) = ($p2, $p1) if $p1 > $p2;
  for (my $i = $p1; $i <= $p2; $i++) {
    $pos{$i}++;
  }
}
close in;

# Output
print "Position\tFrequency\n";
for (my $i = 1; $i <= $N; $i++) {
  print "$i\t$pos{$i}\n";
}

exit;
