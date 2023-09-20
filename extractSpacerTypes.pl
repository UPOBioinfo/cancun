#!/usr/bin/perl
use strict;

# I-Fa: grep -wf ../types/ifa.ab ../genes/aba_genes_ssrA_dr2.tsv | cut -f3 | grep "[^,]+,CRISPR,[^,]+" -Eo | grep -v thiE | sus
my $VIR = "unknown5757";
my $BET = "BQ1898_06609";
my $TRA = "unknown1879";
my $IFB  = "thiE|csy4";
my $IFA1 = "cas1__2|ssrA";
my $IFA2 = "cas3__2|unknown975";

my %G;
open in, "../types/crispr.ab";
while (<in>) {
  chomp;
  
  $G{$_} = 1;
}
close in;

open NS, ">nspacers.tsv";
open REP, ">rep_ab.fasta";
open SP, ">spacers_ab.fasta";
print NS "Array\tFrequency\n";

my %nrep; my %nsp;
open in, "../genes/aba_genes.tsv";
while (<in>) {
  chomp;

  next unless /CRISPR/;
  my ($ab, $id) = split/\t/;
  next unless $G{$ab};
  $id =~ s/\.\d+$//;
 
  my $type;
  my $x = 0;
  my @ns; my @type; my @rep; my @sp;
  while (/((\w+)\{[a-z][a-z][0-9]{5}_[0-9]{5} \d+ \d+ [+-]\},)?CRISPR\{CRISPR (\d+) (\d+) .\},?(\w+)?/g) {
    my ($p1, $p2, $g1, $g2) = ($3, $4, $2, $5);

    #print "$ab $id - $p1 $p2 - $g1 $g2\n"; 

    # Classify
    if ($g1 eq $VIR || $g2 eq $VIR) {
      if ($type eq "ifa2") {
        $type[$x-1] = "vir2";
      } else {
        $type = "vir1";
      }
    } elsif ($g1 eq $BET || $g2 eq $BET) {
      if ($type eq "ifa2") {
        $type[$x-1] = "bet2";
      }
      $type = "bet1";
    } elsif ($g1 eq $TRA || $g2 eq $TRA) {
      if ($type eq "ifa1") {
        $type[$x-1] = "tra1";
        $type = "tra2";
      } else {
        $type = "tra2";
      }
    } elsif ($IFB =~ /$g1/ || $IFB =~ /$g2/) {
      $type =  "ifb";
    } elsif ($IFA1 =~ /$g1/ || $IFA1 =~ /$g2/) {
      if ($type eq "tra1") {
        $type[$x-1] = "tra2";
        $type = "tra1";
      } elsif ($type eq "tra2") {
        $type = "tra1"; 
      } else {
        $type =  "ifa1";
      }
    } elsif ($IFA2 =~ /$g1/ || $IFA2 =~ /$g2/) {
      if ($type eq "vir1") {
        $type = "vir2";
      } elsif ($type eq "bet1") {
        $type = "bet2";
      } else {
        $type = "ifa2";
      }
    } 

    next unless $type;

    my $flag = 0;      
    #open crispr, "/mnt/data/pangenomes/ab/ccfinder/$ab/$id.gff";
    open crispr, "/home/ajperez/databases/ab/ccfinder/$ab/$id.gff";
    my (@lines) = <crispr>;
    close crispr;

    foreach my $l (@lines) {
      chomp $l;

      next if $l =~ /^##/ || $l =~ /^$/;
      
      my ($att, $n1, $n2, $line) = (split/\t/, $l)[2,3,4,8];
      if ($att eq "CRISPR") {
        next unless $n1 == $p1 && $n2 == $p2;
        $flag = 1;
        $l =~ m/DR=([ATCG]+);.+Number_of_spacers=(\d+)/;
        my ($rep, $ns) = ($1, $2);
        #print "# $ns\n";
        #print "$rep\n";
        $ns[$x] = $ns; $type[$x] = $type; $rep[$x] = $rep;
      } elsif ($att eq "CRISPRspacer" && $flag == 1) {
        $l =~ m/sequence=([ATCGN]+);/; # tomamos tb N, pero en el _uniq los descartaremos
        my $spacer = $1;
        #print "$spacer\n";
        push @{$sp[$x]}, $spacer;
      } elsif ($att eq "RightFLANK") {
        $flag = 0;
      }
    }    
    
    $x++;
    #print "*** $x\n";
  }
  
  # Print
  for (my $i = 0; $i <= $#ns; $i++) {
    #print "$type[$i] ** $ns[$i]\n";
    print NS "$type[$i]\t$ns[$i]\n";

    my $seq_id = "$type[$i]_$ab";
    $nrep{$seq_id}++;
    print REP ">$seq_id\_$nrep{$seq_id}\n";
    print REP "$rep[$i]\n";

    foreach my $seq_sp (@{$sp[$i]}) {
      $nsp{$seq_id}++;
      print SP ">$seq_id\_$nsp{$seq_id}\n";
      print SP "$seq_sp\n";
    }
  }
  
}
close in;
close NS;
close REP;
close SP;

exit;
