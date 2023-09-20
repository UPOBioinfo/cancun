#!/usr/bin/perl
use strict;

my $P = $ARGV[1] || 90;
my $C = $ARGV[2] || 90;
my $L = $ARGV[3] || 130000;
my $ab;
my %n;

open in, $ARGV[0];
while (<in>) {
  chomp;

  if (/^ab[0-9]{5}$/) {
    $ab = $_;
  } else {
    my ($q, $s, $p, $a1, $c, $sp1, $sp2, $slen, $seq) = split/\t/;
    
    next unless $p >= $P && $c >= $C && $slen < $L;
    print "$ab\t$s\t$p\t$c\t$slen\n";
    
    # Separate Fastas from virus
    #$n{$ab}++;
    #$seq =~ s/-//g;
    #open FILE, ">>helpers/$ab.fasta";
    #print FILE ">$ab\_$n{$ab}\n$seq\n";
    #close FILE;
  }
}
close in;

exit;
