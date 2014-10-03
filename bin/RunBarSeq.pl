#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use lib "$Bin/../lib";
use FEBA_utils;

my $usage = <<END
Usage: RunBarSeq.pl [-nosplit] [-debug] [-minQuality 0 ] [ -restart ] [ -pieceLines 20000000 ]
          organism setname barcodes [ fastq.gz or directory ]
   e.g.:
   RunBarSeq.pl psRCH2 psRCH2_ML7_set1 ~mnprice/FEBA_2013/BarSeqPrimersH48 /global/dna/dm_archive/sdm/illumina/00/73/41/7341.4.68147.fastq.gz
   RunBarSeq.pl -minQuality 10 psRCH2 psRCH2_ML7_set1_seqUC BarSeqPrimersH48 fastq/HiSeq_barcodes/psRCH2_ML7_set1

RunBarSeq.pl has two phases -- the first phase splits the fastq file
(if necessary) and launches jobs to the cluster. Wait for these jobs
to finish, then invoke RunBarSeq.pl again to assemble the results for
the second phase. Use -restart to force it to go to the first phase.
To save memory, the second phase ignores barcodes that do not
match the pool file, which is expected to be in
    prefix/g/organism_name/pool or pool.n10

Output files are in FEBA_2013/g/organism/ --
setname.colsum -- total parsed reads per index
setname.poolcount -- counts for strains in the pool
    based on joining to the pool or pool.n10 file.
setname.codes.ignored -- all the counts for the barcodes
	that did not match the pool.

END
    ;

sub maybeRun($); # run command unless $debug is defined
 
my $debug = undef;
{
    my $prefix = glob("~mnprice/FEBA_2013");
    my $scriptdir = "~mnprice/scripts";
    my $scratch = glob("~/scratch");
    my $linesPerPiece = 20*1000*1000;
    my $restart = undef;
    my $nosplit = undef;
    my $minQuality = 0;

    die $usage unless GetOptions('debug' => \$debug,
                                 'nosplit' => \$nosplit,
                                 'restart' => \$restart,
                                 'minQuality=i' => \$minQuality,
                                 'pieceLines=i' => \$linesPerPiece)
        && @ARGV == 4;
    my ($organism, $setname, $barcodes, $fastq) = @ARGV;
    die "No such file: $fastq" unless -e $fastq;
    die "No such file: $barcodes" unless -e $barcodes;
    die "No such directory: $prefix/g/$organism" unless -d "$prefix/g/$organism";

    my @parts;
    my @codes;
    my $partGlob;
    my $codeGlob;
    if (-d $fastq) {
        $partGlob = "$fastq/*fastq.gz";
        @parts = glob($partGlob);
        $codeGlob = "$fastq/*.codes";
        @codes = map { s/.fastq.gz$/.codes/; $_; } @parts;
        @codes = grep { -e $_ } @codes; # only the files that were made
    } else {
        $partGlob ="$scratch/${setname}_BarSeq.part*[0-9]";
        @parts = glob($partGlob);
        $codeGlob = "$scratch/${setname}_BarSeq.part*[0-9].codes";
        @codes = glob($codeGlob);
    }
    print STDERR "See " . scalar(@parts) . " parts and " . scalar(@codes) . " codes files\n";

    my $poolfile = "$prefix/g/$organism/pool";
    $poolfile = "$poolfile.n10" if !-e $poolfile;

    if (scalar(@parts) != scalar(@codes) && @codes > 0 && !defined $restart) {
        die "Warning: codes files do not match part files\n$partGlob\nDelete codes files (or both) and rerun\n";
    }
    if (!defined $restart && @parts > 0 && @parts == @codes && -e $codes[0]) {
        # phase 2 -- combining
        if (!defined $debug) {
            # Parse all the log files
            my @logs = map { s/[.]codes$/.log/; $_; } @codes;
            my $nReads = 0;
            my $nMulti = 0;
            my $nFiles = 0;
            my $nUsable = 0;
            foreach my $log (@logs) {
                if (-e $log) {
                    $nFiles++;
                    open(LOG, "<", $log) || die "Cannot read $log";
                    my $nLines = 0; # numbers of lines parsed with counts of reads
                    while (<LOG>) {
                        if (m/^Reads\s+(\d+)\s+Multiplexed\s+(\d+)\s+Usable\S*\s+(\d+)\b/) {
                            $nReads += $1;
                            $nMulti += $2;
                            $nUsable += $3;
                            $nLines++;
                        }
                    }
                    close(LOG) || die "Error reading $log";
                    print STDERR "Warning: no Reads/Multiplexed/Usable line in $log\n"
                        if $nLines == 0;
                } else {
                    print STDERR "Warning: no log file $log\n";
                }
            }
            print STDERR sprintf("Total reads %d Multi %d (%.1f%%) Usable %d (%.1f%%) for %s from %d files\n",
                                 $nReads,
                                 $nMulti, 100.0*$nMulti/($nReads+0.1),
                                 $nUsable, 100.0*$nUsable/($nReads+0.1),
                                 $setname,
                                 $nFiles);

        }
        if (! -e $poolfile) {
            print STDERR "No pool for $organism -- not combining codes files $codeGlob\n";
        } else {
            print STDERR "Combining codes files $codeGlob\n";
            my $cmd = "$scriptdir/combineBarSeq.pl $prefix/g/$organism/${setname} $poolfile $codeGlob";
            maybeRun($cmd);
        }            
        exit(0);
    }

    # else
    # phase 1, splitting (if necessary) and submitting the jobs to make codes files
    if (! -d $fastq) {
        if (defined $nosplit && scalar(@parts) > 0) {
            print STDERR "Skipping the split step, just redoing the analysis\n";
            maybeRun("rm $codeGlob")
                if scalar(@codes) > 0 && ! defined $debug;
        } else {
            # do the splitting
            system("rm $scratch/${setname}_BarSeq.part*") if scalar(@parts) > 0 || scalar(@codes) > 0 && ! defined $debug;
            my $cmd = "gunzip -c $fastq | split -l $linesPerPiece -d - $scratch/${setname}_BarSeq.part";
            maybeRun($cmd);
            system("touch $scratch/${setname}_BarSeq.part00") if defined $debug;
        }
    }
    my $cmdsfile = "$prefix/cmds/${setname}_BarSeq.codecmds";
    open(CMDS, ">", $cmdsfile) || die "Cannot write to $cmdsfile";
    foreach my $i (glob $partGlob) {
        my $corecmd = "$scriptdir/MultiCodes.pl -primers $barcodes -minQuality $minQuality";
        if ($i =~ m/[.]gz$/) {
            my $out = $i;
            $out =~ s/.fastq.gz$//;
            print CMDS "zcat $i | $corecmd -out $out >& $out.log"."\n";
        } else {
            print CMDS "$corecmd -out $i < $i >& $i.log"."\n";
        }
    }
    close(CMDS) || die "Error writing to $cmdsfile";

    maybeRun("~mnprice/qsub.pl -d $prefix/cmds -N $setname. < $cmdsfile");
    print STDERR "Submitted jobs for $setname -- rerun when they complete\n" if !defined $debug;
    print STDERR "Cannot find $poolfile -- poolcounts will not be created in second pass\n" if ! -e $poolfile;
}

sub maybeRun($) {
    my ($cmd) = @_;
    if (defined $debug) {
        print STDERR "Would run: $cmd\n";
    } else {
        system($cmd) == 0 || die "script failed: $cmd";
    }
}
