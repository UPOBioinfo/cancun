#!/usr/bin/perl
use strict;

my $id;
my %seq;
my %as;

# Gather fasta
open in, "helper_assemblies.fasta";
while (<in>) {
  chomp;

  if (/^>(.+)/) {
    $id = $1;
  } else {
    $seq{$id} .= $_;
  }
}
close in;

# Map ab and assembky ID
open in, "helper_assemblies.tsv";
while (<in>) {
  chomp;

  my ($ab, $as) = split/\t/;
  $as{$ab} = $as;
}
close in;

# Gather coordinates and mask
open in, "helper_assemblies_coor2.tsv";
while (<in>) {
  chomp;

  my ($ab, $p1, $p2) = split/ /;
  my $as = $as{$ab};
  my $seq = $seq{$as};
  if ($p1 && $p2) {
    $p1 = $p1 - 1;
    my $le = $p2 - $p1;
    my $n = "N" x $le;
    my $subseq = substr $seq, $p1, $le;
    $seq =~ s/$subseq/$n/;
  }
  print ">$as\n$seq\n";
}
close in;

exit;
