#!/usr/bin/perl -w
use strict;
use List::Util qw(sum);
sub commify($);

die "Run as a filter on the filtered table\n"
  if @ARGV > 0;

# locusId1 => locusId2 => list of #reads for this strain (barcode pair)
my %gp = ();
my $nGenicPairs = 0;
my $nCentralPairs = 0;
while(my $line = <STDIN>) {
  chomp $line;
  my ($bc1, $sc1, $strand1, $pos1, $locus1, $f1, $bc2, $sc2, $strand2, $pos2, $locus2, $f2, $n) = split /\t/, $line;
  die "Wrong number of fields in input\n" unless defined $n;
  next if $bc1 eq "barcode1" && $n eq "n";
  die "Invalid barcode $bc1" unless $bc1 =~ m/^[ACGT]+$/;
  die "Invalid barcode $bc2" unless $bc2 =~ m/^[ACGT]+$/;
  die "Invalid count $n" unless $n =~ m/^\d+$/;
  next if $locus1 eq "" || $locus2 eq "";
  die "Invalid f1 $f1" unless $f1 >= 0 || $f1 <= 1;
  die "Invalid f2 $f2" unless $f2 >= 0 || $f2 <= 1;
  $nGenicPairs++;
  if ($f1 >= 0.1 && $f1 <= 0.9 && $f2 >= 0.1 && $f2 <= 0.9) {
    push @{ $gp{$locus1}{$locus2} }, $n;
    $nCentralPairs++;
  }
}
print STDERR sprintf("Read %s genic pairs and %s central pairs\n",
                     commify($nGenicPairs), commify($nCentralPairs));


my %reads1 = ();
my %strains1 = ();
my %reads2 = ();
my %strains2 = ();
while (my ($locus1, $hash) = each %gp) {
  while (my ($locus2, $list) = each %$hash) {
    # Should self hits be ignored? Decide not to
    $strains1{$locus1} += scalar(@$list);
    $reads1{$locus1} += sum(@$list);
    $strains2{$locus2} += scalar(@$list);
    $reads2{$locus2} += sum(@$list);
  }
}

print join("\t", qw{locusId1 locusId2 nStrains nReads revStrains revReads totStrains1 totReads1 totStrains2 totReads2})."\n";
foreach my $locus1 (sort keys %reads1) {
  next unless exists $reads2{$locus1};
  foreach my $locus2 (sort keys %reads2) {
    next unless exists $reads1{$locus2};
    my $list = $gp{$locus1}{$locus2} || [];
    my $rev = $gp{$locus2}{$locus1} || [];
    die $list unless defined $list;
    die $rev unless defined $rev;
    print join("\t", $locus1, $locus2,
               # sums need || 0 because sum of an empty list is NA, not 0
               scalar(@$list), sum(@$list) || 0,
               scalar(@$rev), sum(@$rev) || 0,
               $strains1{$locus1}, $reads1{$locus1},
               $strains2{$locus2}, $reads2{$locus2})."\n";
  }
}

sub commify($) {
  local $_  = shift;
  1 while s/^(-?\d+)(\d{3})/$1,$2/;
  return $_;
}
