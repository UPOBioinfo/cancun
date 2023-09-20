#!/usr/bin/perl
use strict;

my $P = $ARGV[1] || 90;
my $C = $ARGV[2] || 90;
my $ab;

open in, $ARGV[0];
while (<in>) {
  chomp;

  if (/^ab[0-9]{5}$/) {
    $ab = $_;
  } else {
    my ($q, $s, $p, $c, $a1, $seq) = split/\t/; # v1
    
    next unless $p >= $P && $c >= $C;
    print "$ab\t$p\t$c\n";
    
    # Separate Fastas from virus
    #open FILE, ">helpers/$ab.fasta";
    #print FILE ">$ab\n$seq\n";
    #close FILE;
  }
}
close in;

exit;
