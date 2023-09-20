#!/usr/bin/perl
use strict;

my %ref;
open in, "../roary/ab2/clustered_proteins";
while (<in>) {
  chomp;

  my ($ref, @ab) = split/\t/;
  if ($ref =~ / ([a-z][a-z]([0-9]){5}_([0-9]){5})$/) {
    $ref = $1;
  }

  push @ab, $ref;
  foreach my $ab (@ab) {
    $ref{$ab} = $ref;
  }
}
close in;

# ffn file
open in, $ARGV[0];
while (<in>) {
  chomp;

  my ($sp, $id) = split/\t/;
  print "$sp\t$ref{$id}\n";
}
close in;

exit;
