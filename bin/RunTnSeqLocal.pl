#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use File::Which;

my $usage = <<END
Usage: RunTnSeqLocal.pl [ -limit nMaxReadsPerFile ] [ -blat blatcmd ]
    nickname library_name modelsuffix fastq_directory
Must be run in the parent directory of g/nickname
The file model_modelsuffix must exist in the feba primers/ directory.
Writes to g/nickname/library_name
-limit and -blat are passed on to MapTnSeq.pl
END
    ;

{
    my $limit = undef;
    my $blat = -e "$Bin/blat" ? "$Bin/blat" : "blat";
    die $usage unless GetOptions('limit=i' => \$limit,
				 'blat=s' => \$blat)
	&& @ARGV==4;
    my ($nickname, $library, $modelsuffix, $fastqdir) = @ARGV;

    die "blat is not on the path!" if ($blat eq "blat" && !defined which("blat"));

    die "Not a directory: $fastqdir" unless -d $fastqdir;
    die "No such file: g/$nickname" unless -d "g/$nickname";
    my $modelfile = "$Bin/../primers/model_".$modelsuffix;
    die "No such file: $modelfile" unless -e $modelfile;

    my $cmdsfile = "$fastqdir/RunTnSeqLocal.cmds";
    system("rm ${cmdsfile}*") if -e $cmdsfile;

    open(CMDS, ">", $cmdsfile) || die "Error writing to $cmdsfile";
    my @inputs = glob("$fastqdir/*.fastq.gz");
    my @mapped = ();
    foreach my $file (@inputs) {
	my $base = $file;
	$base =~ s/.fastq.gz$//;
	push @mapped, "$base.mapped";
	my $cmd = "$Bin/MapTnSeq.pl -genome g/$nickname/genome.fna -model $modelfile -first $file";
	$cmd .= " -limit $limit" if defined $limit;
	$cmd .= " -blat $blat" if defined $blat;
	$cmd .= " > $base.mapped.tmp && mv $base.mapped.tmp $base.mapped";
	unlink("$base.mapped") if -e "$base.mapped";
	print CMDS "$cmd\n";
    }
    close(CMDS) || die "Error writing to $cmdsfile";
    print STDERR "Wrote " . scalar(@inputs) . " mapping jobs to $cmdsfile";
    system("$Bin/submitter.pl",$cmdsfile) == 0
	|| die "Error running $Bin/submitter.pl on $cmdsfile\n";

    system("(grep -h Reads $cmdsfile-*.log; grep -h Proportion $cmdsfile-*.log) >& g/$nickname/$library.reads");

    my $cmd = "$Bin/DesignRandomPool.pl -minN 10 -genes g/$nickname/genes.tab -pool g/$nickname/$library.pool "
	. join(" ",@mapped) . " >& g/$nickname/$library.pool.stats";
    print STDERR "Running: $cmd\n";
    system($cmd) == 0
	|| die "DesignRandomPool.pl failed";
    print STDERR "Wrote g/$nickname/$library.pool and .pool.stats\n";
}
