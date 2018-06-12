#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$Bin};
use lib "$Bin/../lib";
use FEBA_Utils; # for ReadFastaEntry

my $usage = <<END
Usage: fastare.pl [ -both ] -query pattern < fasta > out.table
The pattern must be in IUPAC form, such as ACGTNNNNACGT.
Multiple values of query can be given.
By default only searches on the + strand (because patterns are often palindromic anyway).
END
;

my %codes = ("R" => "AG",
             "Y" => "CT",
             "S" => "GC",
             "W" => "AT",
             "K" => "GT",
             "M" => "AC",
             "B" => "CGT",
             "D" => "AGT",
             "H" => "ACT",
             "V" => "ACG",
             "N" => "ACGT");

my ($debug, $both);
my @query;
die $usage
  unless GetOptions('debug' => \$debug, 'both' => \$both, 'query=s{1,}' => \@query)
  && @ARGV == 0
  && @query > 0;
my @regexp;
foreach my $query (@query) {
  $query =~ m/^[ACGTURYSWKMBDHVN]+$/ || die "Invalid query $query -- must use IUPAC codes only\n";
  my $regexp = $query;
  $regexp =~ s/U/T/g;
  while (my ($code, $chars) = each %codes) {
    $regexp =~ s/$code/[$chars]/g;
  }
  print STDERR "Processed query $query to regexp $regexp\n" if $debug;
  push @regexp, $regexp;
}

my $state = {};
while (my ($header, $sequence) = ReadFastaEntry(\*STDIN, $state)) {
  $sequence =~ m/^[ACGTRYSWKMBDHVNX]+$/i ||
    die "Not a nucleotide sequence for $header\n";
  my $rc = reverseComplement($sequence) if defined $both;
  foreach my $i (0..(scalar(@query)-1)) {
    while ($sequence =~ /$regexp[$i]/gi) {
      # @- holds the zero-based offset of the match
      print join("\t", $header, $query[$i], "+", 1 + $-[0])."\n";
    }
    if (defined $both) {
      while ($rc =~ /$regexp[$i]/gi) {
        my $pos = $+[0];
        print join("\t", $header, $query[$i], "-", length($sequence) + 1 - $pos)."\n";
      }
    }
  }
}
