#!/usr/bin/perl
use strict;

# I-Fa: grep -wf ../types/ifa.ab ../genes/aba_genes_ssrA_dr2.tsv | cut -f3 | grep "[^,]+,CRISPR,[^,]+" -Eo | grep -v thiE | sus
my %C = ("AO964_23250" => "C1", "CP54" => "C2", # "cas2" => "C" 
         "ygbF" => "E1", "cas3" => "E2",
         "cas6f" => "F1", "cas1__3" => "F2",
         "yceJ_3" => "V1"); # aparece en I-C, I-E e I-F

my %G;
open in, "../../types/crispr.ab";
while (<in>) {
  chomp;
  
  $G{$_} = 1;
}
close in;

open NS, ">nspacers.tsv";
open REP, ">rep_ab.fasta";
open RP, ">rp_ab.fasta";
open SP, ">spacers_ab.fasta";
print NS "Array\tFrequency\n";

my %nrep; my %nrp; my %nsp;
open in, "../../genes/aba_genes_ssrA_dr.tsv";
while (<in>) {
  chomp;

  next unless /CRISPR/;
  my ($ab, $id) = split/\t/;
  next unless $G{$ab};
  $id =~ s/\.\d+$//;
 
  my $type;
  my $x = 0;
  my @ns; my @type; my @rep; my @sp; my @rp;
  while (/((\w+)\{[a-z][a-z][0-9]{5}_[0-9]{5} \d+ \d+ [+-]\},)?CRISPR\{CRISPR (\d+) (\d+) .\},?(\w+)?/g) {
    my ($p1, $p2, $g1, $g2) = ($3, $4, $2, $5);

    #print "$ab $id - $p1 $p2 - $g1 $g2\n"; 

    # Classify
    if ($C{$g1}) {
      $type = $C{$g1};
    } elsif ($C{$g2}) {
      $type = $C{$g2};
    }

    next unless $type;

    my $strand;
    my $flag = 0;
    open crispr, "/mnt/data/pangenomes/backup/pa/pa/ccfinder/$ab/GFF/$id.gff";
    my (@lines) = <crispr>;
    close crispr;

    foreach my $l (@lines) {
      chomp $l;

      next if $l =~ /^##/ || $l =~ /^$/;
      
      my ($att, $n1, $n2, $line) = (split/\t/, $l)[2,3,4,8];
      if ($att eq "CRISPR") {
        next unless $n1 == $p1 && $n2 == $p2;
        $flag = 1;
        if ($l =~ m/DR=([ATCGN]+);.+Number_of_spacers=(\d+);.+potential_direction=(\w+);/) {
          my ($rep, $ns, $dir) = ($1, $2, $3);
          next if $rep =~ /N/;
          #print "# $ns\n";
          #print "$rep\n";

          # Direction
          #if ($rep =~ /^(GTT|TTC|AGT|CAT)/) {
          #  $strand = 0; # forward
          #} elsif ($rep =~ /^(TTT|ATT)/) {
          #  $strand = 1; # reverse
          #} else {
          #  $strand = 2; # unknown
          #}

          print "$type\t$rep\t$ns\t$dir\n"; # to control repeat sequences
          $ns[$x] = $ns; $type[$x] = $type; $rep[$x] = $rep;
        }
      } elsif ($att eq "CRISPRspacer" && $flag == 1) {
        if ($l =~ m/sequence=([ATCGN]+);/) { # tomamos tb N, pero en el _uniq los descartaremos (tb 2 lineas más abajo)
          my $spacer = $1;
          next if $spacer =~ /N/;
          #print "$spacer\n";
          if ($strand == 1) {
            $spacer = reverse $spacer;
            $spacer =~ tr/ATCG/TAGC/;
          } elsif ($strand == 2) {
            next;
          }
          push @{$sp[$x]}, $spacer;
        }
      } elsif ($att eq "CRISPRdr" && $flag == 1) {
        if ($l =~ m/sequence=([ATCGN]+);/) { # tomamos tb N, pero en el _uniq los descartaremos (tb 2 lineas más abajo)
          my $repeat = $1;
          next if $repeat =~ /N/;
          if ($strand == 1) {
            $repeat = reverse $repeat;
            $repeat =~ tr/ATCG/TAGC/;
          } elsif ($strand == 2) {
            next;
          }
          push @{$rp[$x]}, $repeat;
        }
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
    
    foreach my $seq_rp (@{$rp[$i]}) {
      $nrp{$seq_id}++;
      print RP ">$seq_id\_$nrp{$seq_id}\n";
      print RP "$seq_rp\n";
      #my $len_rp = length $seq_rp; # to measure the frequency of different lengths
      #print "$type[$i]\t$len_rp\n";
    }

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

