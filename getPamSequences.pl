#!/usr/bin/perl
#use strict;

my $BLAST = $ARGV[0] || "spacers_vs_p1virus_95_100.blast";
my $VIRUS = $ARGV[1] || "pivirus"; # if cds is for cds(ffn)
my $R = $ARGV[2]; # 12 if binding ifa1 tra1 vir1 bet1
my $W = 10; # Pam length
my $id;
my %seq;

# Output files
my $ext;
if ($VIRUS eq "cds") {
  $W = 2;
  $VIRUS = "p1virus_cds";
}
$ext = "_12" if $R eq "12";
my @T = ("ifa1", "ifa2", "ifb", "vir1", "vir2", "bet1", "bet2", "tra1", "tra2");
@GT = ("ifa1", "ifa2", "ifb") if $R eq "12";
foreach my $t (@T) {
  my $F = $t;
  open $F, ">pam_$VIRUS\_$F$ext.seq";
}
open FILE, ">pam_$VIRUS$ext.seq";
open FILE2, ">pam_$VIRUS$ext\_protosp.seq";

# Gather fasta
my $FVIRUS = "../$VIRUS.fasta";
$FVIRUS = "./sat_long_variants.fasta" if $VIRUS eq "p1virus_cds";
open in, $FVIRUS;
while (<in>) {
  chomp;
  
  if (/^>(.+)/) {
    $id = $1;
  } else {
    $seq{$id} .= $_;
  }
}
close in;

# Gather types if 12 activated
my %t;
if ($R eq "12") {
  my $TYPE_FILE = "./spacers_ab_uniq_types2.tsv";
  open in, $TYPE_FILE;
  while (<in>) {
    chomp;

    my ($id, $t) = split/\t/;
    push @{$t{$id}}, $t; # can be various types
  }
  close in;
}

# Gather clustered_proteins if cds activated
my %ab;
if ($VIRUS eq "p1virus_cds") {
  open in, "../../roary/ab2/clustered_proteins";
  while (<in>) {
    chomp;

    my ($c1, @c) = split/\t/;
    my ($gn, $id1) = split/: /, $c1;
    push @c, $id1;
  
    foreach my $c (@c) {
      $ab{$c} = $id1;
    }
  }
  close in;
}

# Check Blast and extract PAM sequence
open in, $BLAST;
while (<in>) {
  chomp;

  my $pam;
  my $protosp;
  my ($sp, $vir, $p1, $p2) = (split/\t/)[0,1,4,5];
  my $len_protosp = abs($p2 - $p1) + 1;
  $vir = $ab{$vir} if $VIRUS eq "p1virus_cds";
  if ($p1 < $p2) {
    $protosp = substr $seq{$vir}, $p1, $len_protosp if $VIRUS eq "p1virus_cds";
    my $rr = 11;
    $rr = 3 if $VIRUS eq "p1virus_cds";
    $p1 = $p1 - $rr;
    $pam = substr $seq{$vir}, $p1, $W;
  } else {
    if ($VIRUS eq "p1virus_cds") {
      $protosp = substr $seq{$vir}, $p2, $len_protosp;
      $protosp = reverse $protosp;
      $protosp =~ tr/ATCG/TAGC/;
    }
    $pam = substr $seq{$vir}, $p1, $W;
    $pam = reverse $pam;
    $pam =~ tr/ATCG/TAGC/;
  }
  
  next unless length($pam) == $W; # remove if we want protospacers
  #next unless length($protosp) == $len_protosp; # activate if we want protospacers

  my $type;
  if ($R eq "12" && $VIRUS eq "p1virus_cds") {
    foreach $type (@{$t{$sp}}) {
      print $type "$ab{$vir}\t$pam\n";

      next if $type eq "ifb"; # jump if ifb
      print FILE "$ab{$vir}\t$pam\n";
      print FILE2 "$ab{$vir}\t$protosp\n";
    }
  } elsif ($R eq "12") {
    foreach $type (@{$t{$sp}}) {
      print $type "$pam\n";
      print FILE "$vir\t$pam\n";
    }
  } else {
    ($type) = split/_/, $sp;
    print $type "$pam\n";
    print FILE "$pam\n";
  }
  #print "$sp\t$vir\t$pam\n";
}
close in;

foreach my $t (@GT) {
  close $t;
}
close FILE;

exit;

