#!/usr/bin/perl
use strict;

my $G = $ARGV[0] || "ssrA";
open in, "../strains.ab";
my (@AB) = <in>;
chomp @AB;
close in;

# Searching
my %g;
open in, "../genes/aba_genes_ssrA_dr.tsv";
while (<in>) {
  chomp;

  my ($ab, $id, $line) = split/\t/;
  if ($line =~ /($G)\{[a-z][a-z][0-9]{5}_[0-9]{5} [0-9]+ ([0-9]+) .\},ABSDF2517\{.+?,ssrA_dr\{ssrA_dr ([0-9]+) [0-9]+ .\}/ || 
     $line =~ /(ssrA_dr)\{ssrA_dr [0-9]+ ([0-9]+) .\},.+?,ABSDF2517\{[^\}]+\},$G\{[a-z][a-z][0-9]{5}_[0-9]{5} ([0-9]+) [0-9]+ .\}/) {
    
    my $first = $1;
    my $p1 = $2 + 1;
    my $p2 = $3 - 1;
    my $len = $p2 - $p1 + 1;

    next if $len > 90000;
    #open fna, "/mnt/data/pangenomes/ab/fna/$ab.fna";
    open fna, "/home/ajperez/databases/ab/fna/$ab.fna";
    my $flag; my $seq;
    while (<fna>) {
      chomp;

      if (/^>$id$/) {
        $flag = 1;
      } elsif (/^>/ && $flag) {
        last;
      } elsif ($flag) {
        $seq .= $_;
      }
    }
    my $sseq = substr $seq, $p1, $len;
    if ($first eq "ssrA_dr") {
      my $sseq = reverse $sseq;
      $sseq =~ tr/ACGT/TGCA/;
    }

    next unless $sseq;
    print ">P1virus\_$ab\n$sseq\n";
  }
} 
close in;

exit;

