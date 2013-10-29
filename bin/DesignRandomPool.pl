#!/usr/bin/perl -w
# Given the output of MapTnSeq.pl, possibly from multiple runs, choose the reliable and unambiguous tags
# Makes a tab-delimited table
#
# Input is assumed to contain
# read, barcode, scaffold, position, strand, unique flag, qBeg, qEnd, bit score, %identity
use strict;
use Getopt::Long;
sub reverseComplement($);
sub Variants($); # return all 1-nt variants for a sequence

my $minN = 5;
my $minFrac = 0.75;
my $minRatio = 8.0;
my $maxQBeg = 3;

my $usage = <<END
Usage: DesignRandomPool.pl [ -minN $minN ] [ -minFrac $minFrac ]
	     [ -minRatio $minRatio ] [ -maxQBeg $maxQBeg ]
	     MapTnSeq_files > pool_file
    The input must be tab-delimited with fields
    read,barcode,scaffold,pos,strand,uniq,qBeg,qEnd,score,identity
    where qBeg and qEnd are positions in the read that match the genome after trimming the transposon.

minN is the minimum number of "good" reads for a barcode supporting
its mapping.  (Good means a unique hit to the genome and qBeg=1.)
minFrac is the minimum fraction of input reads for the barcode that
agree with the preferred mapping.  minRatio is the minimum ratio of
reads for preferred mapping over the2nd-most-frequent mapping.
END
    ;

{
    GetOptions('minN=i' => \$minN,
	       'minFrac=f' => \$minFrac, 'minRatio=f' => \$minRatio,
	       'maxQBeg=i' => \$maxQBeg )
	|| die $usage;

    my %barHits = (); # barcode => list of (scaffold,position,strand,uniqueness,qBeg,qEnd)

    my $nMapped = 0;
    my $nHitsRepeat = 0;
    my $nUsable = 0;
    my $nSkipQBeg = 0;

    my %pastEnd = (); # number of reads for that barcode mapped past the end

    while(<>) {
	chomp;
	my ($read,$barcode,$scaffold,$pos,$strand,$uniq,$qBeg,$qEnd,$score,$identity) = split /\t/, $_;
	if ($scaffold eq "pastEnd") {
	    $pastEnd{$barcode}++;
	    $nMapped++;
	} elsif ($maxQBeg >= 1 && $qBeg <= $maxQBeg)  {
	    push @{ $barHits{$barcode} }, [$scaffold,$pos,$strand,$uniq,$qBeg,$qEnd,$score,$identity];
	    $nMapped++;
	} else {
	    $nSkipQBeg++;
	}
    }
    print STDERR "Read $nMapped mapped reads for " . scalar(keys %barHits) . " distinct barcodes\n";
    print STDERR "(Skipped $nSkipQBeg reads with qBeg > $maxQBeg)\n" if $nSkipQBeg > 0;

    print join("\t", "barcode", "rcbarcode", "nTot",
	       "n","scaffold","strand","pos",
	       "n2","scaffold2","strand2","pos2","nPastEnd")."\n";

    my %nInCategory = (); # classification of barcodes
    my ($SCAFFOLD,$POS,$STRAND,$UNIQ,$QBEG,$QEND,$SCORE,$IDENTITY) = (0,1,2,3,4,5,6,7);
    my $nMulti = 0; # barcodes seen >once
    my %barcodeAt = (); # barcode to list of nTot, nMax, at, nNext, nextAt
    my $nReadsForUsable = 0;
    while(my ($barcode, $hits) = each %barHits) {
	my $nPastEnd = $pastEnd{$barcode} || 0;
	my $nTot = scalar(@$hits) + $nPastEnd;
	next unless $nTot >= $minN;
	$nMulti++;
	

	# Is there a location that accounts for at least 90% of the reads, including the past-end reads?
	my @ats = map { join("\t", $_->[$SCAFFOLD], $_->[$STRAND], $_->[$POS]) } @$hits;
	my %nAts = ();
	foreach my $at (@ats) { $nAts{$at}++; }
	my @nAt = sort {$a<=>$b} values(%nAts);
	my $nMax = $nAt[-1];

	if ($nPastEnd >= $nTot/2 || $nPastEnd >= $nMax) {
	    $nInCategory{"PastEnd"}++;
	    my $n2 = $nMax || 0; # note we do not report secondary location (doubt it matters)
	    print join("\t", $barcode, reverseComplement($barcode),
		       $nTot, $nPastEnd, "pastEnd","","",$n2,"","","",$nPastEnd)."\n";
	    next;
	}
	unless($nMax >= $minN && $nMax/$nTot >= $minFrac) {
	    $nInCategory{"NoCons"}++;
	    next;
	}

	my $maxAt = (grep { $nAts{$_} == $nMax } @ats)[0];
	my @iModal = grep { $ats[$_] eq $maxAt } (0..(scalar(@ats)-1));

	# checking unique & qbeg=1 -- but the latter part may be redundant with maxQBeg
	my $nGood = scalar( grep { $_->[$UNIQ] eq "1" && $_->[$QBEG] eq "1" } @$hits[@iModal] );
	unless ($nGood >= $minN && $nGood/$nTot >= $minFrac) {
	    $nInCategory{"FewGood"}++;
	    next;
	}

	my $nNext = @nAt > 1 ? $nAt[-2] : 0;
	my $nextAt = join("\t","","",""); # second-most-common at
	if ($nNext >= 2) {
	    my @nextAt = grep {$nAts{$_} == $nNext} @ats;
	    $nextAt = $nextAt[0];
	}
	unless($nMax >= $minRatio * $nNext) {
	    $nInCategory{"LoRatio"}++;
	    next;
	}
	$barcodeAt{$barcode} = [ $nTot, $nMax, $maxAt, $nNext, $nextAt ];
	$nInCategory{"Usable"}++;
	$nReadsForUsable += $nTot;
    }

    print STDERR "$nMulti barcodes seen $minN or more times, map $nInCategory{Usable} (minFrac $minFrac minRatio $minRatio)\n";
    foreach my $category (sort keys %nInCategory) {
	print STDERR sprintf("%s\t%d\t%.4f\n", $category, $nInCategory{$category}, $nInCategory{$category} / $nMulti);
    }

    my $nOut = 0;
    my $nMasked = 0;
    my $nMaskedReads = 0;
    while (my ($barcode,$row) = each %barcodeAt) {
	my ($nTot,$nMax,$maxAt,$nNext,$nextAt) = @$row;
	my @variants = Variants($barcode);
	my $mask = 0;
	foreach my $variant (@variants) {
	    if (exists $barcodeAt{$variant}
		&& $barcodeAt{$variant}[0] > $nTot
		&& $barcodeAt{$variant}[2] eq $maxAt) {
		$nMasked++;
		$nMaskedReads += $nMax;
		$mask = 1;
		next;
		last;
	    }
	}
	next if $mask;
	my @atSplit = split /\t/, $maxAt;
	my @nextAtSplit = split /\t/, $nextAt, -1; # keep empty entries
	my $nPastEnd = $pastEnd{$barcode} || 0;
	print join("\t", $barcode, reverseComplement($barcode),
		   $nTot, $nMax, @atSplit, $nNext, @nextAtSplit,$nPastEnd)."\n";
	$nOut++;
    }
    print STDERR "Masked $nMasked off-by-1 barcodes ($nMaskedReads reads) leaving $nOut barcodes\n";
    print STDERR sprintf("Reads for those barcodes: %d of %d (%.1f%%)\n",
			 $nReadsForUsable, $nMapped, 100*$nReadsForUsable/($nMapped + 1e-6));
}

sub Variants($) {
    my ($baseseq) = @_;
    my @out = ();
    $baseseq = uc($baseseq);
    my $len = length($baseseq);
    foreach my $i (0..($len-1)) {
	my $pre = substr($baseseq,0,$i);
	my $char = substr($baseseq,$i,1);
	my $post = substr($baseseq,$i+1);
	next unless $char eq "A" || $char eq "C" || $char eq "G" || $char eq "T";
	foreach my $newchar (qw{A C G T}) {
	    push @out, $pre . $newchar . $post unless $newchar eq $char;
	}
    }
    return(@out);
}

# no longer used
sub uniq(@) {
    my %seen = ();
    foreach my $value (@_) { $seen{$value} = 1; }
    return keys(%seen);
}

sub reverseComplement($) 
{
    my $seq = shift;
    chomp $seq;
	my $origSeq=$seq;

    die "Invalid sequence \"$origSeq\" in reverseComplement" if ( ($seq =~ 
tr/RYKMSWBDHVNATCGXrykmswbdhvnatcg-/YRMKWSVHDBNTAGCXyrmkwsvhdbntagc-/) != 
length($seq) );
    $seq = reverse $seq;
    return $seq;
}
