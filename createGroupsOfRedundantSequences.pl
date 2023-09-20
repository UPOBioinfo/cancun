#!/usr/bin/perl
use strict;

my $F = $ARGV[0];
my ($name, $ext) = split/\./, $F;
open FASTA, ">$name\_uniq.$ext";

my %seq;
my $id;

open in, $F;
while (<in>) {
  chomp;

  if (/^>(.+)/) {
    $id = $1;
  } else {
    next if $_ =~ /N/;
    my $s1 = $_;
    my $s2 = $s1;
    my @s = sort($s1, $s2);
    my $s12 = join "_", @s;

    $seq{$id} .= $s12;
  }
}
close in;

my %red;
foreach my $i (keys %seq) {
  push @{$red{$seq{$i}}}, $i;
}

my $n = 0;
print "Array\tNspacers\n";
foreach my $s (keys %red) {
  my $total = 0;
  my $n5 = 0;
  $n++;
  my $spid = "ab_sp_$n";

  # Check mezcla de arrays
  my %npre;
  my @out;
  my @out5;
  foreach my $seqs (@{$red{$s}}) {
    my ($pre) = split/_/, $seqs;
    $npre{$pre}++;
    $total++;
  }
  foreach my $p (sort keys %npre) {
    push @out, "$p($npre{$p})";
    if ($npre{$p} >= 5) { # minimo de spacer en un tipo de array
      push @out5, "$p";
      $n5 += $npre{$p};
    }
  }
  
  # Print 
  print FASTA ">$spid $total:";
  print FASTA join ",", @out;
  print FASTA " ";
  print FASTA join " ", @{$red{$s}};
  print FASTA "\n";
  my ($seq1) = split/_/, $s; # nos quedamos arbitrariamente con la 1ª
  print FASTA "$seq1\n";

  next unless @out5;
  print join ",", @out5;
  print "\t$n5\n";
}

close FASTA;

exit;

