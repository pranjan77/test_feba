#!/usr/bin/perl -w
use strict;
# Given a bunch of files with counts for barcodes, pick out the
# ones from the pool and make a small table of their counts (out.poolcount).
# Also makes a table of total counts (out.colsum)
# and a file of ignored lines (out.codes.ignored).
# Expects either that all columns are in all files,
# or that each file has just one column

{
    die "Usage: combineBarSeq.pl out pool_file codesfiles...\n"
	unless @ARGV >= 3;
    my $out = shift @ARGV;
    my $poolfile = shift @ARGV;
    my @codesFiles = @ARGV;

    my %pool = (); #  rcbarcode to barcode,scaffold,strand,pos
    open(POOL, "<", $poolfile) || die "Cannot read $poolfile";
    while(<POOL>) {
	chomp;
	my ($barcode,$rcbarcode,undef,undef,$scaffold,$strand,$pos) = split /\t/, $_;
	next if $barcode eq "barcode";
	die "Invalid barcode $barcode" unless $barcode =~ m/^[ACGT]+$/;
	die "Invalid rcbarcode $rcbarcode" unless $rcbarcode =~ m/^[ACGT]+$/;
	die "Invalid position $pos" unless $pos eq "" || $pos =~ m/^\d+$/;
	die "Invalid strain $strand" unless $strand eq "+" || $strand eq "-" || $strand eq "";
	die "Duplicate rcbarcode $rcbarcode" if exists $pool{$rcbarcode};
	$pool{$rcbarcode} = [$barcode,$scaffold,$strand,$pos];
    }
    close(POOL) || die "Error reading $poolfile";
    die "No entries in pool file $poolfile" unless scalar(keys %pool) > 0;
    print STDERR "Read " . scalar(keys %pool) . " pool entries from $poolfile\n";

    open(IGNORE, ">", "$out.codes.ignored") || die "Cannot write to $out.codes.ignored";

    my %counts = (); # rcbarcode to vector of counts
    my @indexes = (); # vector of names of samples
    my @colSums = (); # samples to total number of counts
    my @colSumsUsed = (); # samples to total number of counts for used barcodes
    my $nSamples = 0;
    my $oneperfile = 0;
    my $thisIndex = undef; # if in one-per-file mode, which index is this file?

    my $nUsed = 0;
    my $nIgnore = 0;
    foreach my $file (@codesFiles) {
	open(IN, "<", $file) || die "Cannot read $file";
	my $header = <IN>;
	chomp $header;
	my @cols = split /\t/, $header;
	my $first = shift @cols;
	die "Not a barcode counts file -- $file" unless $first eq "barcode";
	die "No index columns in $file" unless @cols > 0;
	if ($nSamples == 0) { # first file
	    @indexes = @cols;
	    $nSamples = scalar(@indexes);
	    @colSums = (0) x $nSamples;
	    @colSumsUsed = (0) x $nSamples;
	    if (scalar(@cols) == 1) {
		$oneperfile = 1;
		$thisIndex = 0;
	    }
	} elsif ($oneperfile) {
	    die "More than one data column in $file" unless scalar(@cols) == 1;
	    my %oldcol = map { $indexes[$_] => $_ } (0..(scalar(@indexes)-1));
	    my $thisCol = $cols[0];
	    if (exists $oldcol{ $thisCol }) {
		$thisIndex = $oldcol{ $thisCol };
	    } else {
		$nSamples++;
		push @indexes, $thisCol;
		$thisIndex = scalar(@indexes)-1;
		push @colSums, 0;
		push @colSumsUsed, 0;
	    }
	} else {
	    die "Wrong number of columns in $file" unless scalar(@cols) == scalar(@indexes);
	    foreach my $i (0..(scalar(@cols)-1)) {
		die "Index mismatch in $file vs. $codesFiles[0] -- $cols[$i] vs. $indexes[$i]"
		    unless $cols[$i] eq $indexes[$i];
	    }
	}
	my $nThisFile = 0;

	while(<IN>) {
	    $nThisFile++;
	    chomp;
	    my @F = split /\t/, $_;
	    my $barcode = shift @F; # actually rcbarcode
	    # note am allowing N in barcode but not in pool
	    die "Invalid barcode: $barcode" unless $barcode =~ m/^[ACGTN]+$/; 
	    die "Wrong number of columns in $file" unless scalar(@F) == ($oneperfile ? 1 : $nSamples);
	    if (exists $pool{$barcode}) {
		$nUsed++;
		if ($oneperfile) {
		    $colSumsUsed[$thisIndex] += $F[0];
		    $counts{$barcode}[$thisIndex] += $F[0];
		} else {
		    for (my $i = 0; $i < $nSamples; $i++) {
			$colSumsUsed[$i] += $F[$i];
		    }
		    if (exists $counts{$barcode}) {
			my $row = $counts{$barcode};
			for (my $i = 0; $i < $nSamples; $i++) {
			    $row->[$i] += $F[$i];
			}
		    } else {
			$counts{$barcode} = \@F;
		    }
		}
	    } else {
		print IGNORE join("\t",$barcode,@F)."\n";
		$nIgnore++;
	    }
	    if ($oneperfile) {
		$colSums[$thisIndex] += $F[0];
	    } else {
		for (my $i = 0; $i < $nSamples; $i++) {
		    $colSums[$i] += $F[$i];
		}
	    }
	}
	close(IN) || die "Error reading from $file";
	print STDERR "Warning: no entries in $file\n" if $nThisFile == 0;
    }
    print STDERR "Read $nUsed lines for codes in pool and $nIgnore ignored lines across " . scalar(@codesFiles) . " files\n";
    close(IGNORE) || die "Error writing to $out.codes.ignored";

    print STDERR "Fraction Used:";
    for(my $i = 0; $i < $nSamples; $i++) {
	print STDERR sprintf(" %s %.3f", $indexes[$i], ($colSumsUsed[$i])/(1.0+$colSums[$i]));
    }
    print STDERR "\n";

    open(COUNT, ">", "$out.poolcount") || die "Cannot write to $out.poolcount";
    print COUNT join("\t", "barcode", "rcbarcode", "scaffold", "strand", "pos", @indexes)."\n";
    while (my ($rcbarcode, $strainrow) = each %pool) {
	my ($barcode,$scaffold,$strand,$pos) = @$strainrow;
	my $counts = $counts{$rcbarcode};
	my @out = ($barcode,$rcbarcode,$scaffold,$strand,$pos);
	if (defined $counts) {
	    for (my $i = 0; $i < $nSamples; $i++) {
		$counts->[$i] += 0; # ensure vector is full length, not undef;
	    }
	    push @out, @$counts;
	} else {
	    push @out, (0) x $nSamples;
	}
	print COUNT join("\t", @out)."\n";
    }
    close(COUNT) || die "Error writing to $out.poolcount";
    print STDERR "Wrote $out.poolcount\n";

    open(SUM, ">", "$out.colsum") || die "Cannot write to $out.colsum";
    print SUM join("\t", @indexes)."\n";
    print SUM join("\t", @colSums)."\n";
    close(SUM) || die "Error writing to $out.colsum";
    print STDERR "Wrote $out.colsum\n";
}
