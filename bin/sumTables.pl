#!/usr/bin/perl -w
use strict;
my $usage = <<END
sumTables.pl file1 file2 ... fileN > summedfile\n
  Sums the entries in the files, with the 1st column as the shared key
  Input files must be tab-delimited
END
;

{
    die $usage if @ARGV == 0;
    my @files = @ARGV;

    my %counts = (); # index column to list of counts
    my @colNames = ();
    my $nValues = undef;

    foreach my $file (@files) {
	open(FILE, "<", $file) || die "Cannot read $file";
	my $headerLine = <FILE>;
	$headerLine =~ s/[\r\n]+$//;
	my @colNamesThis = split /\t/, $headerLine;
	if (!defined $nValues) {
	    # this is the first file
	    die "Not enough columns in $file" unless @colNamesThis > 1;
	    @colNames = @colNamesThis;
	    $nValues = scalar(@colNames) - 1;
	} else {
	    # verify that column names match
	    die "Wrong number of columns in $file" unless scalar(@colNames) == scalar(@colNamesThis);
	    foreach my $i (0..$nValues) {
		die "Column name $i in $file does not match"
		    unless $colNamesThis[$i] eq $colNames[$i];
	    }
	}
	my $nRows = 0;
	while (my $line = <FILE>) {
	    $line =~ s/[\r\n]+$//;
	    my @F = split /\t/, $line, -1;
	    die "Wrong number of columns" unless scalar(@F) == scalar(@colNames);
	    $nRows++;
	    my $key = shift @F;
	    my $values = $counts{$key};
	    if (!defined $values) {
		$counts{$key} = \@F;
	    } else {
                # $values is a reference so $counts{$key} is being updated
		for (my $i = 0; $i < $nValues; $i++) {
		    $values->[$i] += $F[$i];
		}
	    }
	}
	close(FILE) || die "Error reading $file";
	print STDERR "Read $nRows rows from $file\n";
    }

    print join("\t", @colNames)."\n";
    while (my ($key,$values) = each %counts) {
	print join("\t", $key, @$values)."\n";
    }
}
