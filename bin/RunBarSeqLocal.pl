#!/usr/bin/perl -w
# RunBarSeqLocal.pl -- process a lane of barseq into a table of counts for each strain in each sample
# issuing jobs locally using submitter.pl
use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use lib "$Bin/../lib";

my $usage = <<END
Usage: RunBarSeqLocal.pl  [-minQuality 0 ] [ -pieceLines 20000000 ]
	  [-nosplit] [-debug] [ -limit 1000 ]
          [ -indexes BarSeqPrimersH48 ]
          organism_directory setname fastq.gz_file_or_directory_with_fastq.gz_files
   e.g.:
   feba/bin/RunBarSeqLocal.pl g/Keio Keio_ML9_set1 fastq/Keio_ML9_set1
   feba/bin/RunBarSeqLocal.pl -indexes feba/primers/BarSeqPrimersH48 g/MR1 MR1_ML3_set1 fastq/MR1_ML3_set1
   feba/bin/RunBarSeqLocal.pl -indexes feba/primers/BarSeqPrimersH48 g/psRCH2 psRCH2_ML7_set1 fastq/7341.4.68147.fastq.gz 

RunBarSeq.pl has two phases -- the first phase counts the barcodes in
the fastq file(s).  If the BarSeq was run with the newer style of
primers that have names like IT001, then do not use the -indexes
argument.  In this case it will only accept a directory as input, and
it will assume that each fastq or fastq.gz file in that directory is
already demultiplexed, with a code in the filename such as _IT094_ or
_10_TAGCTT_ to indicate which sample it is from. There can be more
than one file per sample.  For older experiments, where the primers
have names like H01 or M01, use the -indexes argument to describe
which primers were used.  RunBarSeq.pl will use this information to
demultiplex the reads.  If -indexes is used, then RunBarSeq.pl can
accept as input either a directory that contains fastq or fastq.gz
files or a single large_fastq.gz file. Given a single fastq.gz file,
it will split it and process each piece in parallel.

The second phase aggregates the counts of these barcodes. To save
memory, it ignores barcodes that do not match the pool file, which is
expected to be in

    organism_directory/pool or pool.n10

Output files are in organism_directory/ --
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
    my $linesPerPiece = 20*1000*1000;
    my $nosplit = undef;
    my $minQuality = 0;
    my $limitReads = undef;
    my $barcodes = undef;

    die $usage unless GetOptions('debug' => \$debug,
                                 'nosplit' => \$nosplit,
                                 'minQuality=i' => \$minQuality,
                                 'pieceLines=i' => \$linesPerPiece,
				 'limit=i' => \$limitReads,
				 'indexes=s' => \$barcodes)
        && @ARGV == 3;
    my ($gdir, $setname, $fastq) = @ARGV;
    die "No such file: $fastq" unless -e $fastq;
    die "No such directory: $gdir" unless -d $gdir;
    die "No such file: $barcodes" unless !defined $barcodes || -e $barcodes;

    my $poolfile = "$gdir/pool";
    $poolfile = "$poolfile.n10" if !-e $poolfile;
    die "Cannot find $poolfile (or withut .n10)" unless -e $poolfile;

    # design the parts
    my $prefix = $fastq;
    my @parts;
    my @codes;
    my $codeGlob;
    if (-d $fastq) {
	@parts = glob("$fastq/*fastq.gz");
	@parts = grep { !m/^[.]/ }  @parts;
	if (scalar(@parts) == 0) {
	    die "No *.gz files in $fastq" if defined $barcodes;
	    # Indexed runs sometimes have one directory per sample, with fastq.gz file(s) within each directory
	    @parts = glob("$fastq/*/*.fastq.gz");
	    @parts = grep !m"[.]/", @parts; # remove hidden files
	    die "Cannot find *.gz files in $fastq or its subdirectories" unless scalar(@parts) > 0;
	    $codeGlob = "$fastq/*/*.codes";
	} else {
	    $codeGlob = "$fastq/*.codes";
	}
	@codes = map { my $i = $_; $i =~ s/.fastq.gz$/.codes/; $i; } @parts;
	@codes = grep { -e $_ } @codes; # only the files that were previously made
    } else {
	die "Cannot demultiplex the file $fastq without the -indexes option" unless defined $barcodes;
	$prefix = $fastq; $prefix =~ s/fastq[.]gz$//;
	$prefix = "$prefix.parts";
	if (! -d $prefix) {
	    mkdir($prefix) || die "Cannot mkdir $prefix; $!";
	}
        @parts = glob("$prefix/${setname}_BarSeq.part*[0-9]");
        $codeGlob = "$prefix/${setname}_BarSeq.part*[0-9].codes";
        @codes = glob($codeGlob); # only the preexisting files
    }
    print STDERR "See " . scalar(@parts) . " parts and " . scalar(@codes) . " codes files; codes files will be overwritten\n"
	if scalar(@codes) > 0;

    # splitting if necessary
    if (! -d $fastq) {
        if (defined $nosplit && scalar(@parts) > 0) {
            print STDERR "Skipping the split step, just redoing the analysis with " . scalar(@parts) . " existing pieces\n";
            maybeRun("rm $codeGlob")
                if scalar(@codes) > 0 && ! defined $debug;
        } else {
            # do the splitting
            maybeRun("rm $prefix/${setname}_BarSeq.part*") if (scalar(@parts) > 0 || scalar(@codes) > 0);
            my $cmd = "gunzip -c $fastq | split -l $linesPerPiece -d - $prefix/${setname}_BarSeq.part";
            maybeRun($cmd);
            system("touch $prefix/${setname}_BarSeq.part00") if defined $debug;
	    @parts = glob("$prefix/${setname}_BarSeq.part*[0-9]"); # now @parts is actually there
        }
    }

    # build the list of commands, and update @codes to be what we'll make
    @codes = (); # will update with expected results
    my $cmdsfile = "$prefix/${setname}_BarSeq.codecmds";
    maybeRun("rm $cmdsfile*") if -e $cmdsfile;
    maybeRun("rm -f $codeGlob");
    open(CMDS, ">", $cmdsfile) || die "Cannot write to $cmdsfile";
    foreach my $i (@parts) {
	print STDERR "Considering part $i\n" if $debug;

        my $corecmd = "$Bin/MultiCodes.pl -minQuality $minQuality";
	$corecmd .= " -limit $limitReads" if defined $limitReads;
	if (defined $barcodes) {
	    $corecmd .= " -primers $barcodes";
	} else {
	    my @path = split "/", $i;
	    my $name = pop @path;
	    my @pieces = split /[._]/, $name;
	    my @indexes = grep m/^IT\d+$/, @pieces;
	    my $index = undef;
	    # sometimes ITO (capital O) not IT0 (numeral zero)
	    my @indexes2 = grep m/^ITO\d+$/, @pieces;

	    if (@indexes == 1) {
		$index = $indexes[0];
	    } elsif (@indexes2 == 1){
		# e.g. HiSeq_barcodes/FEBA_BS_81/Sample_FEBA_BS_81_ITO88/FEBA_BS_81_ITO88_GGCCTG_L001_R1_001.fastq.gz
		$index = $indexes2[0];
		$index =~ s/ITO/IT0/;
	    } elsif ($name =~ m/_(\d+)_[ACGT][ACGT][ACGT][ACGT][ACGT][ACGT]_/
		     || $name =~ m/_Index(\d+)_[ACGT][ACGT][ACGT][ACGT][ACGT][ACGT]_/i
                     || $name =~ m/_Index(\d+)_S\d+_/) {
                # e.g. FEBA_BS_60_10_TAGCTT_L001_R1_001.fastq.gz
                # e.g. FEBA_BS_125_Index10_S10_L001_R1_001.fastq.gz
		$index = sprintf("IT%03d", $1);
            } elsif ($name =~ m/_(\d+)_S\d+_L\d+_/) {
                # e.g. FEBA_BS_117_24_S120_L002_R1_001.fastq.gz is IT024
                $index = sprintf("IT%03d", $1);
	    } elsif ($name =~ m/undetermined/i) {
		print STDERR "Skipping $name\n";
		next;
	    } else {
		die "Cannot identify the index ITnnn from file $i";
	    }
	    $corecmd .= " -index $index";
	}
        if ($i =~ m/[.]gz$/) {
            my $out = $i;
            $out =~ s/.fastq.gz$//;
	    push @codes, "$out.codes";
            print CMDS "zcat $i | $corecmd -out $out >& $out.log"."\n";
        } else {
	    push @codes, "$i.codes";
            print CMDS "$corecmd -out $i < $i >& $i.log"."\n";
        }
    }
    close(CMDS) || die "Error writing to $cmdsfile";

    # run them all
    maybeRun("$Bin/submitter.pl $cmdsfile");

    # combine the results
    if (!defined $debug) {
	# Parse all the log files
	my @logs = map { my $out = $_; $out =~ s/[.]codes$/.log/; $out; } @codes;
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
    print STDERR "Combining codes files $codeGlob\n";
    my $cmd = "$Bin/combineBarSeq.pl $gdir/${setname} $poolfile $codeGlob";
    maybeRun($cmd);
}

sub maybeRun($) {
    my ($cmd) = @_;
    if (defined $debug) {
        print STDERR "Would run: $cmd\n";
    } else {
        system($cmd) == 0 || die "script failed: $cmd";
    }
}
