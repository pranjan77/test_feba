#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $nTop = 5;
my $magic = "/usr2/people/mprice/data/FEBA/magic_pool";

my $usage = <<END
MagicPoolSplits.pl -mapped mapping_file [ -magic magic_pool ]
                   [ -out prefix ] [-nTop $nTop]

    Given the output of MapTnSeq.pl (mapping_file), and a magic pool,
    report which magic pool transposons worked and demultiplex the
    results into a set of non-magic-pool reads (prefix.nomagic) and a
    pool design-like file for the top nTop magic pool items
    (prefix.promotername)...

    The default magic pool is in $magic
    The default output prefix is the mapping file name
END
    ;

sub Variants($); # return all 1-nt variants for a sequence

{
    my $input = undef;
    my $out = undef;

    die $usage unless GetOptions('mapped=s' => \$input,
				 'magic=s' => \$magic,
				 'out=s' => \$out,
				 'nTop=i' => \$nTop)
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
    my $nonMagicFile = "$out.nomagic";
    open(OUT, ">", $nonMagicFile) || die "Cannot write to $nonMagicFile";
    open(IN, "<", $input) || die "Cannot read $input";
    while(<IN>) {
	chomp;
	my @F = split /\t/, $_, -1;
	my (undef,$barcode,$scaffold,$pos,$strand,$uniq,$qbeg,$bits,$identity) = @F;
	die $_ if !defined $identity;
	if (exists $magic_barcodes{$barcode}) {
	    $nMagic++;
	    $magicCounts{$barcode}[$scaffold eq "pastEnd" ? 1 : 0]++;
	} elsif (exists $variants{$barcode}) {
	    $nOff1++;
	} else {
	    # do not check off by 1 yet
	    $nOff++;
	    print OUT join("\t",@F)."\n";
	}
    }
    close(IN) || die "Error reading $input";
    close(OUT) || die "Error writing $nonMagicFile";
    print STDERR join("\t","Total reads","Magic:",$nMagic,"Off1:",$nOff1,"Other:",$nOff)."\n";
    my @barcodesSorted = sort {$magicCounts{$b}[0] <=> $magicCounts{$a}[0]} (keys %magicCounts);
    my %barcodesTop = (); # top barcode => scaffold => strand => pos => n
    for (my $i = 0; $i < scalar(@barcodesSorted) && $i < $nTop; $i++) {
	my $barcode = $barcodesSorted[$i];
	print STDERR join("\t", @{$magic_barcodes{$barcode}},
			  "genomic", $magicCounts{$barcode}[0], "pastEnd", $magicCounts{$barcode}[1],
			  $barcode)."\n";
	$barcodesTop{$barcode} = {};
    }

    print STDERR "Wrote $nonMagicFile\n";
    open(IN, "<", $input) || die "Cannot read $input";
    while(<IN>) {
	chomp;
	my @F = split /\t/, $_, -1;
	my (undef,$barcode,$scaffold,$pos,$strand,$uniq,$qbeg,$bits,$identity) = @F;
	die $_ if !defined $identity;
	next unless $uniq && $qbeg <= 5;
	$barcodesTop{$barcode}{$scaffold}{$strand}{$pos}++ if exists $barcodesTop{$barcode};
    }
    close(IN) || die "Error reading $input";

    while (my ($barcode, $hash) = each %barcodesTop) {
	my $file = "$out." . $magic_barcodes{$barcode}[0];
	open(OUT, ">", $file) || die "Cannot write to $file";
	print OUT join("\t",qw{barcode rcbarcode nTot n scaffold strand pos n2 scaffold2 strand2 pos2 nPastEnd})."\n";
	my $nRow = 0;
	while (my ($scaffold, $hash2) = each %$hash) {
	    while (my ($strand, $hash3) = each %$hash2) {
		while (my ($pos, $n) = each %$hash3) {
		    $nRow++;
		    print OUT join("\t","b".$nRow, "r".$nRow, $n, $n, $scaffold, $strand, $pos, 0, 0, 0, 0, 0)."\n";
		}
	    }
	}
	close(OUT) || die "Error writing $file";
	print STDERR "Wrote $file\n";
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
