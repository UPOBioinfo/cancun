#!/usr/bin/perl
use strict;

my $R = $ARGV[0]; # if '12', if1 will be the same than vir1, tra1 and bet1

my %ab;
my %nvar;
my %r;
my %t;
my %n;
my %s;
my %sp;
my %spab;
my %ifab;
my %nsptype;
my %nab = ("ifa1" => 253, "ifa2" => 190, "ifb" => 906);
my $TYPE_FILE = "./spacers_ab_uniq_types.tsv";
if ($R eq "12") {
  %nab = ("ifa1" => 314, "ifa2" => 225, "ifb" => 906);
  $TYPE_FILE = "./spacers_ab_uniq_types2.tsv";
}
#my %nab = ("ifa1" => 212, "ifa2" => 176, "ifb" => 906); # sin ifa+ifb

# Avoid strains ifa+ifb
open in, "./ifa_ifb.ab";
while (<in>) {
  chomp;
  $ifab{$_} = 1;
}

# Gather virus variants
foreach my $f ("ab_phigaro.asoc", "p1virus.asoc", "helper_assemblies.asoc") {
#foreach my $f ("p1virus.asoc") {
  open in, $f;
  while (<in>) {
    chomp;

    my ($r, $v) = (split/\t/)[0,2];
    my (@v) = split/,/, $v;
    my $nvar = scalar(@v);
  
    next unless $nvar >= 100 || $r eq "Helper";
    $nvar{$r} = $nvar;
    foreach my $v (@v) {
      $r{$v} = $r;
      #my $ab;
      #if ($v =~ /P1virus/) {
      #  ($ab) = (split/_/, $v)[1];
      #} else {
      #  ($ab) = (split/_/, $v)[0];
      #}
      #push @{$ab{$r}}, $ab;
    }
  }
  close in;
}

# Put numbers to 0
foreach my $c (keys %nvar) {
  foreach my $t ("ifa1", "ifa2", "ifb", "vir1", "vir2", "bet1", "bet2", "tra1", "tra2") {
    @{$n{$c}{$t}} = ();
    @{$s{$c}{$t}} = ();
  }
}

# Gather types
open in, $TYPE_FILE;
while (<in>) {
  chomp;

  my ($id, $t) = split/\t/;
  push @{$t{$id}}, $t; # can be various types
  $nsptype{$t}++;
}
close in;

# Gather abs by spacer
open in, "./spacers_ab_uniq.fasta";
while (<in>) {
  chomp;

  next unless /^>/;
  my ($sp, $a1, @sps) = split/ /;
  $sp =~ s/>//;
  foreach my $e (@sps) {
    my ($type, $ab) = (split/_/, $e)[0,1];

    #next if $ifab{$ab}; # avoid ifa+ifb strains
    push @{$spab{$sp}{$type}}, $ab; # assign by spacer and type
  }
}
close in;

# Read Blast
foreach my $f ("spacers_vs_phigaros_95_100.blast", "spacers_vs_p1virus_95_100.blast", "spacers_vs_helper2_95_100.blast") {
#foreach my $f ("spacers_vs_p1virus_95_100.blast") {
  open in, $f;
  while (<in>) {
    chomp;

    my ($sp, $var) = split/\t/;
    next unless $r{$var};
    foreach my $t (@{$t{$sp}}) { 
      push @{$n{$r{$var}}{$t}}, $var; # absolute number of variants
      push @{$sp{$r{$var}}{$t}}, $sp;
      my %hash = map { $_, 1 } (@{$s{$r{$var}}{$t}}, @{$spab{$sp}{$t}}); @{$s{$r{$var}}{$t}} = keys %hash; # variants
      #print "# $t - @{$s{$r{$c2}}{$t}}}\n";
    }
  }
  close in;
}

# Output
print "Phage_cluster\tnvariants\tfa_ifa1\tfr_ifa1\tfasp_ifa1\tfrsp_ifa1\tfas_ifa1\tfrs_ifa1\tfa_ifa2\tfr_ifa2\tfasp_ifa2\tfrsp_ifa2\t";
print "fas_ifa2\tfrs_ifa2\tfa_ifb\tfr_ifb\tfasp_ifb\tfrsp_ifb\tfas_ifb\tfrs_ifb\n"; #\tfa_vir1\tfr_vir1\tfa_vir2\tfr_vir2\tfa_bet1\tfr_bet1\tfa_bet2\tfr_bet2\n";
foreach my $c (keys %nvar) {
  my $c1 = $c;
  $c1 =~ s/Cluster_/Phage_/;
  $c1 =~ s/Helper/Phage-plasmid/;
  $c1 =~ s/P1virus_cluster/P1virus/;
  print "$c1\t$nvar{$c}";
  #foreach my $t ("ifa1", "ifa2", "ifb", "vir1", "vir2", "bet1", "bet2") {
  foreach my $t ("ifa1", "ifa2", "ifb") {
    # Variants
    my %hash = map { $_, 1 } @{$n{$c}{$t}}; @{$n{$c}{$t}} = keys %hash;
    my $fa = scalar(@{$n{$c}{$t}});
    my $fr = $fa / $nvar{$c};
    print "\t$fa\t";
    printf "%.2f", $fr;
    
    # Spacers
    my %hash = map { $_, 1 } @{$sp{$c}{$t}}; @{$sp{$c}{$t}} = keys %hash;
    my $fasp = scalar(@{$sp{$c}{$t}});
    my $frsp = $fasp / $nsptype{$t};
    print "\t$fasp\t";
    printf "%.2f", $frsp;
    
    # Strains with spacers
    my %hash = map { $_, 1 } @{$s{$c}{$t}}; @{$s{$c}{$t}} = keys %hash;
    #my %hash = map { $_, 1 } @{$ab{$c}}; @{$ab{$c}} = keys %hash;
    my $fas = scalar(@{$s{$c}{$t}});
    my $frs = $fas / $nab{$t};
    print "\t$fas\t";
    printf "%.2f", $frs;
  }
  print "\n";
}

exit;

