#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);

my $usage = <<END
SetupOrg.pl [ -prefix ABC_ ] -gbk file.gbk [ -name NickName | -out outdir ]
	(where file.gbk is from IMG, RAST, or RefSeq)
SetupOrg.pl -name NickName -gff gff_file -fna genome.fna [ -name NickName | -out outdir ]
	(for assemblies from the JGI)

Creates the output directory (if it does not exist already) and within
it, it creates the genome sequence file (genome.fna), the genes files
(genes.tab, genes.GC), and the protein sequence file (aaseq).

The default output directory is g/NickName

If using IMG, please make sure that you get a genbank file that
includes all scaffolds. From the genome page, use the "Generate
genbank file" option, not the download options.

Use the -aaseq option to write the aaseq file only

Use the -prefix argument to create locusIds like ABC_ not GFFnnnn
(this is mostly relevant if reading genbank files from RAST).

Dependencies: Many of the scripts that this calls require BioPerl.
Side effects: deletes aaseq2, pfam.tab, and tigrfam.tab from the output directory (if they exist)
END
    ;

{
    my ($name, $gbkFile, $aaseq, $fnaFile, $gffFile, $prefix, $outdir);
    $prefix = "GFF";

    die $usage unless GetOptions('aaseq' => \$aaseq,
				 'name=s' => \$name,
				 'gbkFile=s' => \$gbkFile,
				 'gffFile=s' => \$gffFile,
				 'fnaFile=s' => \$fnaFile,
				 'prefix=s' => \$prefix,
				 'out=s' => \$outdir )
	&& @ARGV==0;
    die $usage unless (defined $gbkFile) xor (defined $gffFile);
    die $usage unless (defined $gffFile && defined $fnaFile) || (!defined $gffFile && !defined $fnaFile);
    die $usage unless (defined $name) xor (defined $outdir);
    die "Invalid prefix $prefix" unless $prefix =~ m/^[A-Za-z][A-Za-z0-9_]*$/;

    if (!defined $outdir) {
	die "No g directory" unless -d "g";
	$outdir = "g/$name";
    }
    if (! -d $outdir) {
	mkdir($outdir) || die "Cannot mkdir $outdir";
	system("chmod","g+w",$outdir);
    }
    if (defined $gbkFile) {
	die "No such file: $gbkFile" unless -e $gbkFile;
	print STDERR "Creating $outdir from $gbkFile\n";
    }
    print STDERR "Creating aaseq file only\n" if defined $aaseq;

    # these files should be created from aaseq, so delete them before writing to aaseq
    unlink("$outdir/aaseq2");
    unlink("$outdir/pfam.tab");
    unlink("$outdir/tigrfam.tab");

    if (defined $gbkFile && !defined $aaseq) {
        die "Please install genbank2gff.pl in $Bin and make it executable\n"
            . "This script is available from Ian Holmes' gfftools repository\n"
            . "See https://github.com/ihh/gfftools\n"
            unless -x "$Bin/genbank2gff.pl";
	system("$Bin/gbkToSeq.pl $gbkFile > $outdir/genome.fna") == 0
            || system("$Bin/gbkToSeq2.pl $gbkFile > $outdir/genome.fna") == 0
            || die "Both gbkToSeq.pl and gbkToSeq2.pl failed";
	system("$Bin/genbank2gff.pl < $gbkFile > $outdir/genes.gff") == 0 || die "genbank2gff.pl failed";
	system("$Bin/gffToGenes.pl -prefix $prefix < $outdir/genes.gff > $outdir/genes.tab") == 0 || die "gffToGenes.pl failed";
    }
    if (defined $gbkFile) {
	my $code = system("$Bin/gbkToFaa.pl $outdir/genes.tab $gbkFile > $outdir/aaseq");
	if ($code != 0) {
	    print STDERR "Warning: gbkToFaa.pl failed (common with RAST genbank files) -- translating the nt sequences instead\n";
	    system("$Bin/genesTabTranslation.pl $outdir/genes.tab $outdir/genome.fna > $outdir/aaseq") == 0
		|| die "genesTabTranslation.pl failed";
	}
    }

    if (defined $gffFile && !defined $aaseq) {
	system("cp", $fnaFile, "$outdir/genome.fna") == 0 || die $!;
	system("$Bin/gffToGenes.pl -prefix $prefix < $gffFile > $outdir/genes.tab") == 0 || die "gffToGenes.pl failed";
    }
    if (defined $gffFile) {
	system("$Bin/genesTabTranslation.pl $outdir/genes.tab $outdir/genome.fna > $outdir/aaseq") == 0
	    || die "genesTabTranslation.pl failed";
    }

    if (!defined $aaseq) {
	system("$Bin/RegionGC.pl $outdir/genome.fna $outdir/genes.tab > $outdir/genes.GC") == 0 || die "RegionGC.pl failed";
    }
}

