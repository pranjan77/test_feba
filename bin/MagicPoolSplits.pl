#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $nTop = 5;
my $ratio = 50;
my $minN = 3;
my $magic = "/usr2/people/mprice/data/FEBA/magic_pool";

my $usage = <<END
MagicPoolSplits.pl -mapped mapping_file [ -magic magic_pool ]
                   [ -out prefix ] [-nTop $nTop] [-ratio $ratio] [-minN $minN]

    Given the output of MapTnSeq.pl (mapping_file), and a magic pool,
    which defined which transposons had a different (fixed) bracode,
    report which magic pool transposons worked and demultiplex the
    results into a set of non-magic-pool reads (prefix.nomagic) and a
    pool design-like file for the top nTop magic pool items
    (prefix.promotername).

    For each barcode (magic or not), filters out positions that are
    seen ratio x times more often with another barcode -- these are
    presumed to be artefacts of chimeric PCR.

    The default magic pool is in $magic
    The default output prefix is the mapping file name
END
    ;

sub Variants($); # return all 1-nt variants for a sequence
sub IsChimeric($$$$$); # should ignore barcode?
{
    my $input = undef;
    my $out = undef;

    die $usage unless GetOptions('mapped=s' => \$input,
				 'magic=s' => \$magic,
				 'out=s' => \$out,
				 'nTop=i' => \$nTop,
	                         'ratio=f' => \$ratio)
	&& @ARGV == 0;
    die "No mapping file specified:\n$usage" unless defined $input;
    $out = $input if !defined $out;
    die "No magic pool file $magic" unless -e $magic;

    my %magic_barcodes = (); # barcode sequence => name
    open(MAGIC, "<", $magic) || die "Cannot read $magic";
    while(<MAGIC>) {
	chomp;
	my ($promoter,$well,$plasmid,$barcode) = split /\t/, $_, -1;
	next if $barcode eq "barcode";
	die "Invalid barcode in magic pool: $barcode" unless $barcode =~ m/^[ACGT]+$/;
	$magic_barcodes{$barcode} = [$well,$promoter];
    }
    close(MAGIC) || die "Error reading $magic";

    # magic barcode => list of non past-end count, past-end count
    my %magicCounts = map { $_ => [0,0] } (keys %magic_barcodes);
    my $nMagic = 0;
    my $nOff1 = 0;
    my $nOff = 0;

    my %variants = ();
    foreach my $barcode (keys %magic_barcodes) {
	foreach my $variant (Variants($barcode)) {
	    $variants{$variant} = 1;
	}
    }

    my %pos = (); # scaffold => strand => position => barcode => n

    my $nonMagicFile = "$out.nomagic";
    open(IN, "<", $input) || die "Cannot read $input";
    while(<IN>) {
	chomp;
	my @F = split /\t/, $_, -1;
	my (undef,$barcode,$scaffold,$pos,$strand,$uniq,$qbeg,$bits,$identity) = @F;
	die $_ if !defined $identity;
	my $ignore = !$uniq || $qbeg > 3 || $scaffold eq "pastEnd";
	if (exists $magic_barcodes{$barcode}) {
	    $nMagic++;
	    $magicCounts{$barcode}[$scaffold eq "pastEnd" ? 1 : 0]++;
	} elsif (exists $variants{$barcode}) {
	    $nOff1++;
	    $ignore = 1;
	} else {
	    # do not check off by 1 yet
	    $nOff++;
	}
	$pos{$scaffold}{$strand}{$pos}{$barcode}++ unless $ignore;
    }
    close(IN) || die "Error reading $input";

    print STDERR join("\t","Total reads","Magic:",$nMagic,"Off1:",$nOff1,"Other:",$nOff)."\n";
    my @barcodesSorted = sort {$magicCounts{$b}[0] <=> $magicCounts{$a}[0]} (keys %magicCounts);
    for (my $i = 0; $i < scalar(@barcodesSorted) && $i < $nTop; $i++) {
	my $barcode = $barcodesSorted[$i];
	print STDERR join("\t", @{$magic_barcodes{$barcode}},
			  "genomic", $magicCounts{$barcode}[0], "pastEnd", $magicCounts{$barcode}[1],
			  $barcode)."\n";
    }

    open(IN, "<", $input) || die "Cannot read $input";
    open(OUT, ">", $nonMagicFile) || die "Cannot write to $nonMagicFile";
    my $nNonMagicReadsUsed = 0;
    my $nNonMagicReadsChimeric = 0;
    while(<IN>) {
	chomp;
	my @F = split /\t/, $_, -1;
	my (undef,$barcode,$scaffold,$pos,$strand,$uniq,$qbeg,$bits,$identity) = @F;
	die $_ if !defined $identity;
	next unless $uniq && $qbeg <= 3 && $scaffold ne "pastEnd";
	next if exists $magic_barcodes{$barcode} || exists $variants{$barcode};
	if (IsChimeric($barcode,$scaffold,$strand,$pos,\%pos)) {
	    $nNonMagicReadsChimeric++;
	} else {
	    $nNonMagicReadsUsed++;
	    print OUT join("\t",@F)."\n";
	}
    }
    close(IN) || die "Error reading $input";
    close(OUT) || die "Error writing $nonMagicFile";
    print STDERR "Non-magic reads used $nNonMagicReadsUsed (chimeric $nNonMagicReadsChimeric) written to $nonMagicFile\n";

    for (my $i = 0; $i < scalar(@barcodesSorted) && $i < $nTop; $i++) {
	my $barcode = $barcodesSorted[$i];
	my $file = "$out." . $magic_barcodes{$barcode}[0];
	open(OUT, ">", $file) || die "Cannot write to $file";
	print OUT join("\t",qw{barcode rcbarcode nTot n scaffold strand pos n2 scaffold2 strand2 pos2 nPastEnd})."\n";
	my $nRow = 0;
	my $nChimeric = 0;
	my $nRare = 0;
	my $nUsed = 0;
	while (my ($scaffold, $hash1) = each %pos) {
	    while (my ($strand, $hash2) = each %$hash1) {
		while (my ($pos, $hash3) = each %$hash2) {
		    next unless exists $hash3->{$barcode};
		    my $n = $hash3->{$barcode};
		    if (IsChimeric($barcode,$scaffold,$strand,$pos,\%pos)) {
			$nChimeric += $hash3->{$barcode};
		    } elsif ($n < $minN) {
			$nRare += $n;
		    } else {
			$nUsed += $hash3->{$barcode};
			print OUT join("\t","b".$nRow, "r".$nRow, $n, $n, $scaffold, $strand, $pos, 0, 0, 0, 0, 0)."\n";
		    }
		}
	    }
	}
	close(OUT) || die "Error writing $file";
	print STDERR "Magic reads for $barcode used $nUsed (chimeric $nChimeric, rare $nRare) written to $file\n";
    }
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

sub IsChimeric($$$$$) {
    my ($barcode,$scaffold,$strand,$pos,$poshash) = @_;
    my $counthash = $poshash->{$scaffold}{$strand}{$pos};
    die "Invalid input $scaffold $strand $pos" unless defined $counthash;
    die "Unexpected barcode $barcode" unless exists $counthash->{$barcode};
    my $maxcount = 0;
    foreach my $val (values %$counthash) {
	$maxcount = $val if $val > $maxcount;
    }
    return ($maxcount >= $counthash->{$barcode} * $ratio) ? 1 : 0;
}
