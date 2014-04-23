#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use lib "$Bin/../lib";
use FEBA_utils;

{
    my $debug = 0;
    my $num_lines = 50*1000*1000;
    my $gp_project = "gentech-rqc";
    my $help = 0;
    GetOptions('debug'       => \$debug,
               'project=s'   => \$gp_project,
               'help|h'      => \$help,
               'lines|l=i'   => \$num_lines);
    die "Usage: RunTnSeq.pl [--project project] ref.fa model_file outdir outbase fastq1.gz fastq2.gz ... fastqN.gz\n"
        . "   i.e., RunTnSeq.pl Phaeo Phaeo_ML1 pKMW7 fastafile reads1.fastq.gz reads2.fastq.gz reads3.fastq.gz\n"
        if $help or @ARGV < 5;
    my $fna_file = shift;
    my $model_file = shift;
    my $outdir = shift; 
    my $outbase = shift;

    my $fastq_files = join(" ",@ARGV);
    
    for my $fastq (@ARGV){
        die "No such file: $fastq" unless -e $fastq;
    }
    
    mkdir $outdir unless -d $outdir;
    mkdir "$outdir/cmds" unless -d "$outdir/cmds";
    die "No genome file: $fna_file" unless -e $fna_file;
    
    die "No model file: $model_file" unless -e $model_file;

    print STDERR "gunzipping and splitting $fastq_files\n";
    system("gunzip -c $fastq_files | split -l $num_lines -d - $outdir/${outbase}_TnSeq.part") == 0 || die "Split failed";
    
    my $part_glob = "$outdir/${outbase}_TnSeq.part[0-9][0-9]";
    my @fq_parts = glob $part_glob; 
    my $map_glob = "$outdir/${outbase}_TnSeq.part[0-9][0-9].mapped";
    my @mapped = glob($map_glob);

    print STDERR "See " . scalar(@fq_parts) . " parts files and " . scalar(@mapped) . " mapped files\n";

    if (scalar(@fq_parts ) != scalar(@mapped) && @mapped > 0) {
        die "Warning: mapped files do not match part files for\n$outdir/${outbase}_TnSeq.part*[0-9]*\nDelete codes files (or both) and rerun\n";
    }
    if (scalar(@mapped) == 0){ 
        my $map_tn_seq_cmd = "$Bin/MapTnSeq.pl -genome $fna_file ".
                             "-blat blat ".
                             "-first $outdir/${outbase}_TnSeq.partSUB ".
                             "-model $model_file ".
                             "> $outdir/${outbase}_TnSeq.partSUB.mapped";
        my $qsub_script = "$outdir/cmds/${outbase}_TnSeq.MapTnSeq.sh";
        open (my $qsub_script_fh,">",$qsub_script) or die "Can't open $qsub_script\n";
        my @modules = ( "blat" );
        my @cmds = ( $map_tn_seq_cmd );
        FEBA_utils::create_job_array_script($qsub_script_fh,1,"5G","12:00:00", $gp_project, \@cmds, \@modules );
        close $qsub_script_fh;
        
        my $job_id = FEBA_utils::submit_job_array($qsub_script,"runTnSeq",1,scalar(@fq_parts),1,$debug);
        
        exit if $debug;
        print STDERR "MapTnSeq job running under job id $job_id\n";
        FEBA_utils::wait($job_id);
        @mapped = glob($map_glob);
    } else {
        print STDERR "Found ".scalar(@mapped)." mapped files\n";
        if (scalar(@fq_parts) != scalar(@mapped)){
            print STDERR "Found different number of parts files and mapped files. Stopping to prevent overwriting wanted files\n";
            exit(-1);
        }
    }

    # if mapping has already been done, and it looks like everything ran, run DesignRandomPool
    if (@fq_parts > 0 && scalar(@fq_parts) == scalar(@mapped) && -e $mapped[0]) {
        print STDERR "Reading mapping files:\n";
        print STDERR join("\n",@mapped)."\n";
        my $qsub_script = "$outdir/cmds/${outbase}_TnSeq.DesignRandomPool.sh";
        open(my $qsub_script_fh,">",$qsub_script) or die "Can't open $qsub_script\n";
        my $design_random_pool_cmd = "$Bin/DesignRandomPool.pl -minN 10 ".join(" ",@mapped)." > $outdir/$outbase.pool.n10 2> $outdir/$outbase.pool.n10.log";
        my @cmds = ( $design_random_pool_cmd );
        FEBA_utils::create_job_script($qsub_script_fh,1,"5G","12:00:00",\@cmds);
        my $job_id = FEBA_utils::submit_job($qsub_script,"designRandomPool",$debug);
        print STDERR "DesignRandomPool job running under job id $job_id\n";
        exit if $debug;
        FEBA_utils::wait($job_id);
    }

}
