#!/usr/bin/perl -w
use strict;

{
    my $prefix = glob("~mnprice/FEBA_2013");
    my $scratch = glob("~/scratch");
    my $linesPerPiece = 50*1000*1000;
    die "Usage: RunTnSeq.pl organism library model fastq.gz\n"
        . "   i.e., RunTnSeq.pl Phaeo Phaeo_ML1 pKMW7 fastafile\n"
        unless @ARGV==4;
    my ($organism, $libname, $model, $fastq) = @ARGV;
    die "No such file: $fastq" unless -e $fastq;
    die "No such directory: $prefix/g/$organism" unless -d "$prefix/g/$organism";
    my $fnafile = "$prefix/g/$organism/genome.fna";
    die "No genome file: $fnafile" unless -e $fnafile;
    my $modelfile = "$prefix/model_$model";
    die "No model file: $modelfile" unless -e $modelfile;

    my @parts = glob("$scratch/${libname}_TnSeq.part*[0-9]");
    my $mapGlob = "$scratch/${libname}_TnSeq.part*[0-9].mapped";
    my @mapped = glob($mapGlob);

    print STDERR "See " . scalar(@parts) . " parts files and " . scalar(@mapped) . " mapped files\n";

    if (scalar(@parts ) != scalar(@mapped) && @mapped > 0) {
        die "Warning: mapped files do not match part files for\n$scratch/${libname}_TnSeq.part*[0-9]*\nDelete codes files (or both) and rerun\n";
    }
    if (@parts > 0 && @parts == @mapped && -e $mapped[0]) {
        print STDERR "Reading mapping files $mapGlob\n";
        my $cmd = "(~mnprice/scripts/DesignRandomPool.pl -minN 10 $mapGlob > $prefix/g/$organism/pool.n10) >& $prefix/g/$organism/pool.n10.log";
        system($cmd) == 0 || die "Failed: $cmd";
        exit(0);
    }

    # else, need to make the codes files (slow)
    system("rm $scratch/${libname}_TnSeq.part*") if scalar(@parts) > 0 || scalar(@mapped) > 0;

    print STDERR "gunzipping and splitting $fastq\n";
    system("gunzip -c $fastq | split -l $linesPerPiece -d - $scratch/${libname}_TnSeq.part") == 0 || die "Split failed";

    my $cmdsfile = "$prefix/cmds/${libname}_TnSeq.codecmds";
    open(CMDS, ">", $cmdsfile) || die "Cannot write to $cmdsfile";
    foreach my $i (glob "$scratch/${libname}_TnSeq.part*[0-9]") {
        print CMDS "~mnprice/scripts/MapTnSeq.pl -genome $fnafile -blat ~mnprice/scripts/blat -first $i -model $modelfile > $i.mapped"."\n";
    }
    close(CMDS) || die "Error writing to $cmdsfile";

    system("~mnprice/qsub.pl -d $prefix/cmds -N ${organism}_TnSeq_map < $cmdsfile") == 0 || die "qsub.pl of jobs ${organism}_TnSeq_map failed";
    print STDERR "Submitted jobs for ${organism}_TnSeq_map -- rerun when they complete\n";
}
