#!/usr/bin/perl
use strict;

my $F = $ARGV[0] || "spacers_ab_uniq.fasta";
my $B = $ARGV[1]; # 12 if bind ifa1, vir1, tra1 and bet1
my %r; # to avoid redundancy in final results
my %ty;

print "Type\tFrequency\n";

open in, $F;
while (<in>) {
  chomp;
 
  next unless /^>/;
  my ($id, $t) = split/ /;
  $id =~ s/>//;
  my ($n, $t) = split/:/, $t;
  my (@t) = split/,/, $t;

  my @out;
  for my $t (@t) {
    my ($type, $nt) = split/\(/, $t;
    $nt =~ s/\)//;
    my $f = $nt/$n;
    next unless $f > 0.01; # minimal frequency to classify to this type

    # bind to ifa1 and ifa2
    if ($B eq "12") {
      $type = "ifa1" if $type =~ /^bet1|vir1|tra1$/;
      $type = "ifa2" if $type =~ /^bet2|vir2|tra2$/;
    }

    # output
    next if $r{$id}{$type};
    push @out, $type;
    $r{$id}{$type} = 1;
  }

  # Output
  my $ty = join ",", @out;
  $ty{$ty}++;
}
close in;

# Final output
foreach my $ty (sort keys %ty) {
  print "$ty\t$ty{$ty}\n";
}

exit;
