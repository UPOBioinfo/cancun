#!/usr/bin/perl
use strict;

my $F = $ARGV[0];
my $N = 81981;
my %coo;
my %pos;

# Gather spacer types
my %t;
my $TYPE_FILE = "../spacers_ab_uniq_types2.tsv";
open in, $TYPE_FILE;
while (<in>) {
  chomp;

  my ($id, $t) = split/\t/;
  $t{$id} .= $t; # can be various types
}
close in;

# Read GFF
open in, "../gggenes/p1virus.gggenes";
while (<in>) {
  chomp;

  next if /^molecule\t/;
  my ($v, $g, $id1, $p1, $p2) = split/\t/;
  ($p1, $p2) = ($p2, $p1) if $p1 > $p2;
  push @{$coo{$g}}, "$p1\_$p2"; # array due to the 'interrupted genes'
}
close in;

# Gather clustered_proteins if cds activated
my %gn;
open in, "../../roary/ab2/clustered_proteins";
while (<in>) {
  chomp;

  my ($c1, @c) = split/\t/;
  my ($gn, $id1) = split/: /, $c1;
  push @c, $id1;
  foreach my $c (@c) {
    $gn{$c} = $gn;
  }
}
close in;

# Initialize to 0
for (my $i = 1; $i <= $N; $i++) {
  $pos{$i} = 0;
}

# Run through Blast
open in, $F;
while (<in>) {
  chomp;

  my ($sp, $id) = (split/\t/)[0, 1];

  next if $t{$sp} =~ /ifb/; # jump I-Fb spacers
  
  my $g = $gn{$id};
  foreach my $coo (@{$coo{$g}}) {
    my ($p1, $p2) = split/_/, $coo;
    for (my $i = $p1; $i <= $p2; $i++) {
      $pos{$i}++;
    }
  }
}
close in;

# Output
print "Position\tFrequency\n";
for (my $i = 1; $i <= $N; $i++) {
  
  next if $pos{$i} == 0; # remove 0's
  print "$i\t$pos{$i}\n";
}

exit;

