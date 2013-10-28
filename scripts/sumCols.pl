#!/usr/bin/perl -w
use strict;

{
    die "Run as a filter" unless @ARGV == 0;
    my $header = <STDIN>;
    chomp $header;
    my @colnames = split /\t/, $header;
    my @n = map { 0 } @colnames;
    my $nLines = 0;
    while (my $line = <STDIN>) {
	chomp $line;
	my @F = split /\t/, $line;
	foreach my $i (0 .. (scalar(@F)-1)) {
	    $n[$i] += $F[$i] if $F[$i] =~ m/^-?[0-9]+[0-9.]*$/;
	}
	$nLines++;
    }
    print STDERR "Read $nLines data lines\n";
    if ($nLines > 0) {
	print join("\t", @colnames)."\n";
	print join("\t", @n)."\n";
    }
}
