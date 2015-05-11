#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FindGene; # for LocationToGene()
use FEBA_Utils; # for ReadTable(), ReadColumnNames(), ReadFasta()

my $minReads = 2; # minimum #reads at each location to consider reliable
my $minRegion = 100; # minimum size of region with no insertions to consider

my $usage = <<END
Usage: Essentiality.pl [ -minReads $minReads ]
	-out prefix -genome fasta_file -genes tab_file mapped1.tab ... mappedN.tab

The mapped.tab files are from MapTnSeq.pl. The genes table must include the fields
    locusId, scaffoldId, begin, end, strand

By default, there must be at least minReads insertions at a location
or it is ignored.
END
    ;

{
    my $out = undef;
    my $genomeFile = undef;
    my $genesFile = undef;
    die $usage unless GetOptions('minReads=i' => \$minReads,
				 'out=s' => \$out,
				 'genome=s' => \$genomeFile,
				 'genes=s' => \$genesFile)
	&& defined $out && defined $genomeFile && defined $genesFile
	&& @ARGV > 0;
    my @mapfiles = @ARGV;
    my $seqs = &ReadFasta($genomeFile); # scaffoldId => nucleotide sequence
    my $totlen = 0;
    foreach my $seq (values %$seqs) {
	$totlen += length($seq);
    }
    print STDERR "Genome length: $totlen scaffolds: " . scalar(keys %$seqs) . "\n";

    my @genes = &ReadTable($genesFile, qw{locusId scaffoldId begin end strand});
    print STDERR "Read " . scalar(@genes) . " genes\n";
    my %genesSorted = (); # scaffold to list of genes sorted by begin
    foreach my $gene (@genes) {
	push @{ $genesSorted{$gene->{scaffoldId}} }, $gene;
    }
    foreach my $scaffold (keys %genesSorted) {
	my @sorted = sort { $a->{begin} <=> $b->{begin} } @{ $genesSorted{$scaffold} };
	$genesSorted{$scaffold} = \@sorted;
    }
    CheckGeneLocations(\%genesSorted); # writes to STDERR

    my %locusToGene = ();
    # modify begin and end so that no genes overlap.
    while (my ($scaffold, $genes) = each %genesSorted) {
	my $ngenesScaffold = scalar(@$genes);
	foreach my $i (0..($ngenesScaffold-1)) {
	    my $gene = $genes->[$i];
	    my $prevGene = $i==0 ? undef : $genes->[$i-1];
	    my $nextGene = $i==$ngenesScaffold-1 ? undef : $genes->[$i+1];
	    if ($gene->{begin} > $gene->{end}) {
		print STDERR "Warning: skipping wraparound ORF " . $gene->{locusId} . "\n";
		next; # these genes will be skipped
	    }
	    # begin2 and end2 are the non-overlapping regions
	    die "Duplicate locusId " . $gene->{locusId} if exists $locusToGene{$gene->{locusId}};
	    $locusToGene{$gene->{locusId}} = $gene;
	    $gene->{begin2} = $gene->{begin};
	    $gene->{end2} = $gene->{end};
	    $gene->{begin2} = $prevGene->{end}+1
		if defined $prevGene && $prevGene->{end} >= $gene->{begin};
	    $gene->{end2} = $nextGene->{begin}-1
		if defined $nextGene && $nextGene->{begin} < $nextGene->{end} && $nextGene->{begin} <= $gene->{end};
	}
    }

    my $blat = "$Bin/blat";
    if (! -x $blat) {
        $blat = which("blat");
	die "blat is not in $Bin or on path" if ! defined $blat;
    }

    my $nReads = 0;
    my $nSkip = 0; # high qBeg
    my $nUsed = 0; # (genomic and OK qBeg only)
    my %counts = (); # scaffoldId => strand => pos => number of reads
    foreach my $mapfile (@mapfiles) {
	open(MAP, "<", $mapfile) || die "Cannot read $mapfile";
	while(<MAP>) {
	    $nReads++;
	    chomp;
	    my ($readname, $barcode, $scaffoldId, $pos, $strand, $uniq, $readBeg, $readEnd, $score, $identity) = split /\t/, $_, -1;
	    die "Not enough columns in\n$_" unless defined $identity;
	    next if $scaffoldId eq "pastEnd";
	    die "Unrecognized scaffold $scaffoldId" unless exists $seqs->{$scaffoldId};
	    die "Invalid position $pos" unless $pos =~ m/^\d+$/ && $pos >= 1 && $pos <= length($seqs->{$scaffoldId});
	    die "Invalid strand $strand" unless $strand eq "+" || $strand eq "-";
	    if ($readBeg > 3) { # should be rare
		$nSkip++;
		next;
	    }
	    # accept possibly non-unique reads ss the bias is to avoid false positives
	    $counts{$scaffoldId}{$strand}{$pos}++;
	    $nUsed++;
	}
	close(MAP) || die "Error reading $mapfile";
    }

    # Count the "good" reads and write tehm out along with what gene they are in (if any)
    # This does not address the non-unique genome issue
    open(POS, ">", "$out.pos") || die "Cannot write to $out.pos";
    print POS join("\t", qw{scaffoldId strand pos nReads locusId f})."\n";
    my $nUsed2 = 0; # above minReads only
    my $nStrandedLoc = 0;
    while (my ($scaffoldId, $scaffoldHash) = each %counts) {
	while (my ($strand, $strandHash) = each %$scaffoldHash) {
	    while (my ($pos, $count) = each %$strandHash) {
		if ($count >= $minReads) {
		    $nUsed2 += $count;
		    $nStrandedLoc++;
		    my ($locusId, $f) = &LocationToGene($scaffoldId, $pos, \%genesSorted);
		    print POS join("\t", $scaffoldId, $strand, $pos, $count, $locusId, $f)."\n";
		    if ($locusId ne "") {
			my $gene = $locusToGene{$locusId};
			if ($pos >= $gene->{begin2} && $pos <= $gene->{end2}) {
			    $gene->{nPos}++;
			    $gene->{nReads} += $count;
			    if ($f >= 0.1 && $f <= 0.9) {
				$gene->{nPosCentral}++;
				$gene->{nReadsCentral} += $count;
			    }
			}
		    }
		}
	    }
	}
    }
    close(POS) || die "Error writing to $out.pos";
    print STDERR "Wrote $out.pos\n";

    # And write out reads per gene
    my $tmpFna = "$out.$$.fna";
    my $tmpPsl = "$out.$$.psl";
    open(FNA, ">", $tmpFna) || die "Cannot write to $tmpFna";
    while (my ($scaffoldId, $genes) = each %genesSorted) {
	my $scseq = $seqs->{$scaffoldId};
	foreach my $gene (@$genes) {
	    next unless defined $gene->{begin2} && $gene->{begin2} < $gene->{end2};
	    # the non-overlapping part only; ignore strandedness
	    my $ntseq = substr($scseq, $gene->{begin2}-1, $gene->{end2} - $gene->{begin2} + 1);
	    print FNA ">".$gene->{locusId}."\n".$ntseq."\n";
	}
    }
    close(FNA) || die "Error writing to $tmpFna";
    # ntoe -fastMap means no introns high %identity
    my $blatcmd = "$blat -t=dna -out=blast8 -fastMap -minScore=50 $genomeFile $tmpFna $tmpPsl";
    system($blatcmd) == 0 || die "Error running $blatcmd\n$!";

    open(PSL, "<", $tmpPsl) || die "Error reading $tmpPsl";
    while(<PSL>) {
	chomp;
	my @F = split /\t/, $_, -1;
	my ($query, $subject, $identity, $len, $mm, $gaps, $qBeg, $qEnd, $sBeg, $sEnd, $evalue, $score) = @F;
	my $gene = $locusToGene{$query};
	die "Invalid locusId $query in BLAT results" unless defined $gene;
	# ignore self hits
	next if $subject eq $gene->{scaffoldId} && $sBeg >= $gene->{begin2} && $sEnd <= $gene->{end2};
	$gene->{dupScore} = $score unless exists $gene->{dupScore} && $gene->{dupScore} > $score;
    }
    close(PSL) || die "Error reading $tmpPsl";
    unlink($tmpFna);
    unlink($tmpPsl);

    my @geneHeader = ReadColumnNames($genesFile);
    push @geneHeader, qw{ntLenNoOverlap dupScore nPos nReads nPosCentral nReadsCentral};
    open(GENES, ">", "$out.genes") || die "Cannot write to $out.genes";
    print GENES join("\t", @geneHeader)."\n";
    while (my ($scaffold, $genes) = each %genesSorted) {
	foreach my $gene (@$genes) {
	    next unless $gene->{begin} < $gene->{end} && exists $gene->{begin2} && $gene->{begin2} < $gene->{end2};
	    $gene->{ntLenNoOverlap} = $gene->{end2} - $gene->{begin2} + 1;
	    foreach my $field (qw{dupScore nPos nReads nPosCentral nReadsCentral}) {
		$gene->{$field} = 0 if !exists $gene->{$field};
	    }
	    my @fields = map $gene->{$_}, @geneHeader;
	    print GENES join("\t", @fields)."\n";
	}
    }
    close(GENES) || die "Error writing to $out.genes";
    print STDERR "Wrote $out.genes\n";

    # Identify regions with no reads, and test if they are unique
    my %regions = (); # scaffoldId => begin => list of end, dupScore, GC
    open(FNA, ">", $tmpFna) || die "Cannot write to $tmpFna";
    foreach my $scaffoldId (keys %counts) {
	my %insAt = ();
	foreach my $strand (qw{+ -}) {
	    my $hash = $counts{$scaffoldId}{$strand};
	    foreach my $pos (keys %$hash) { $insAt{$pos} = 1; }
	}
	my @insAt = sort { $a <=> $b } keys %insAt;
	foreach my $i (1..(scalar(@insAt)-1)) {
	    my $beg = $insAt[$i-1] + 1;
	    my $end = $insAt[$i] - 1;
	    if ($end - $beg + 1 >= 100) {
		my $seq = substr($seqs->{$scaffoldId}, $beg-1, $end-$beg+1);
		my @hits;
		@hits = $seq =~ m/[GC]/g;
		my $nGC = scalar(@hits);
		@hits = $seq =~ m/[AT]/g;
		my $nAT = scalar(@hits);
		my $GC = $nGC/($nGC + $nAT);
		$regions{$scaffoldId}{$beg} = [ $end, 0, $GC ];
		print FNA ">$scaffoldId:$beg:$end\n$seq\n";
	    }
	}
    }
    close(FNA) || die "Error writing to $tmpFna";
   
    system($blatcmd) == 0 || die "Error running $blatcmd\n$!";
    open(PSL, "<", $tmpPsl) || die "Cannot read $tmpPsl";
    while(<PSL>) {
	chomp;
	my @F = split /\t/, $_, -1;
	my ($query, $subject, $identity, $len, $mm, $gaps, $qBeg, $qEnd, $sBeg, $sEnd, $evalue, $score) = @F;
	my ($scaffoldId,$begin,$end) = split /:/, $query;
	die "Invalid query $query" unless defined $end && $regions{$scaffoldId}{$begin}[0] == $end;
	next if $subject eq $scaffoldId && $sBeg == $begin && $sEnd == $end; # ignore hit to itself
	$regions{$scaffoldId}{$begin}[1] = $score unless $regions{$scaffoldId}{$begin}[1] > $score;
    }
    close(PSL) || die "Error reading $tmpPsl";
    unlink($tmpFna);
    unlink($tmpPsl);

    open(REGIONS, ">", "$out.regions") || die "Cannot write to $out.regions";
    print REGIONS join("\t", qw{scaffoldId begin end GC dupScore})."\n";
    while (my ($scaffoldId, $hash) = each %regions) {
	foreach my $begin (sort {$a <=> $b} keys %$hash) {
	    my ($end,$score,$GC) = @{ $hash->{$begin} };
	    print REGIONS join("\t", $scaffoldId, $begin, $end, $GC, $score)."\n";
	}
    }
    close(REGIONS) || die "Error writing to $out.regions\n";
    print STDERR "Wrote $out.regions\n";

    print STDERR sprintf("Usable genomic reads %d (skipped %d), or average %.3f per nucleotide\n",
			 $nUsed, $nSkip, $nUsed/$totlen);
    print STDERR sprintf("With at least %d reads: %d reads (%.3f per nt) at %d stranded locations (1 per %.3f nt per strand)\n",
			 $minReads, $nUsed2, $nUsed2/$totlen, $nStrandedLoc, $totlen/(2*$nStrandedLoc));
}
