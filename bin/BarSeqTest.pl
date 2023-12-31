#!/usr/bin/perl -w
# Process a small set of barseq test samples, including MultiCodes.pl, combineBarSeq.pl, and BarSeqR.pl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);

my $usage =<<END
Usage: BarSeqTest.pl -org organism [ -n25 | -bs3 ] -index S1:S2:S3 -desc Time0:Lactate:Glucose
	    -fastqdir directory

Assumes that g/organism includes all the typical files, including a pool.n10
or pool file, and writes to g/organism/barseqtest by default.

The fastq directory should include file(s) named _index_*.fastq.gz or
index_*.fastq.gz for each index. BarSeqTest looks in subdirectories,
not just in the directory itself.

Optional arguments:
  -pool g/organism/pool.n10
  -outdir g/organism/barseqtest
  -hascodes -- use the *.codes files that already exist in the fastq directory
  -test -- do not actually run any commands
END
    ;

sub maybeRun($); # run command unless $test is defined
my $test = undef; # check for files, but do no work

{
    my ($org, $indexSpec, $descSpec, $fastqdir, $pool, $outdir);
    my ($n25, $bs3, $hascodes);

    GetOptions('org=s' => \$org,
	       'index=s' => \$indexSpec,
	       'desc=s' => \$descSpec,
	       'fastqdir=s' => \$fastqdir,
	       'pool=s' => \$pool,
	       'test' => \$test,
               'n25' => \$n25,
               'bs3' => \$bs3,
               'hascodes' => \$hascodes,
	       'outdir=s' => \$outdir) || die $usage;
    @ARGV == 0 || die $usage;
    die $usage unless defined $org && defined $indexSpec && defined $descSpec && defined $fastqdir;
    die "No such directory: g/$org" unless -d "g/$org";
    die "No such directory: $fastqdir" unless -d $fastqdir;
    die "Cannot specify both -n25 and -bs3" if defined $n25 && defined $bs3;
    
    if (!defined $pool) {
	$pool = "g/$org/pool.n10";
	unless (-e $pool) {
	    $pool = "g/$org/pool";
	    die "No such file: g/$org/pool.n10 or g/$org/pool" unless -e $pool;
	}
    } else {
	die "No such file: $pool" unless -e $pool;
    }

    # parse the index and description fields
    my @indexes = split /:/, $indexSpec, -1;
    die "Invalid -index $indexSpec -- must have at least two, colon-separated " if @indexes < 2;
    my @desc = split /:/, $descSpec, -1;
    die "Invalid -desc $descSpec -- must match -index" if scalar(@desc) != scalar(@indexes);
    my @time0s = grep { $_ eq "Time0" } @desc;
    die "Must have both Time0 and non Time0 samples" if scalar(@time0s) < 1 || scalar(@time0s) == scalar(@desc);

    my %index;
    my %numToIndex;
    foreach my $index (@indexes) {
      $index{$index} = 1;
      if ($index =~ m/^[A-Z]+0*(\d+)$/) {
        $numToIndex{$1} = $index;
      }
    }
    die "Non-unique indexes\n" unless scalar(keys %index) == scalar(@indexes);
    if (scalar(@indexes) != scalar(keys %numToIndex)) {
      print STDERR "Warning: cannot identify unique sample numbers for each index; not matching by number\n";
      %numToIndex = ();
    }

    # find the fastq files
    my %fastqFiles = (); # index => list of matching files within $fastqdir
    # -H means follow symbolic link if it is the argument
    my @files = `find -H $fastqdir -name '*.fastq.gz'`;
    foreach my $file (@files) {
      chomp $file;
      my $name = $file;
      $name =~ s/[.]fastq[.]gz$//;
      my ($sampleNum, $index);
      if ($name =~ m/_(\d+)_S\d+_L\d+_/
          || $name =~ m/_Index(\d+)_S\d+/
          || $name =~ m/_Index(\d+)_[ACGT][ACGT][ACGT][ACGT][ACGT][ACGT]_/i) {
        # e.g. FEBA_BS_117_24_S120_L002_R1_001.fastq.gz is #24
        # e.g. FEBA_BS_125_Index10_S10_L001_R1_001.fastq.gz is #10
        # e.g. FEBA_BS_60_10_TAGCTT_L001_R1_001.fastq.gz is #10
        $sampleNum = $1;
      } elsif ($name =~ m/(IT\d+)_/ || $name =~ m/_(IT\d+)_/
              || $name =~ m!_(IT\d+)/!) {
        # e.g. FEBA_BS_195_IT001/FEBA_BS_195_S97_L002_R1_001.fastq.gz is IT001
        $index = $1;
      } elsif ($name =~ m/^[A-Z]*\d*_S(\d+)_*L*\d*_R\d+/) {
        # e.g. A10_S106_L003_R1_001.fastq.gz is #106
        # e.g. A10_S10_R1_001.fastq.gz is #10
        $sampleNum = $1;
      } else {
        die "Cannot identify sample number of index for file $file\n";
      }
      if (defined $sampleNum) {
        $sampleNum =~ s/^0*// if defined $sampleNum;
        $index = $numToIndex{$sampleNum} if exists $numToIndex{$sampleNum};
      }
      if (defined $index && exists $index{$index}) {
        push @{ $fastqFiles{$index} }, $file;
        print STDERR "fastq for index $index: $file\n";
      }
    }

    $outdir = "g/$org/barseqtest" if !defined $outdir;
    mkdir($outdir) if ! -d $outdir && !defined $test;
    my @codesFiles = ();
    foreach my $index (@indexes) {
      die "No fastq files for $index in $fastqdir\n" if !exists $fastqFiles{$index};
      my @filenames = @{ $fastqFiles{$index} };
      if (defined $hascodes) {
        my @codesThis = @filenames;
        s/[.]fastq[.]gz$/.codes/ for @codesThis;
        foreach my $file (@codesThis) {
          die "No such file: $file\n" unless -e $file;
          push @codesFiles, $file;
        }
      } else {
        my $pathSpec = join(" ", @filenames);
        push @codesFiles, "$outdir/$index.codes";
        my $extraopt = "";
        $extraopt = "-n25" if defined $n25;
        $extraopt = "-bs3" if defined $bs3;
        &maybeRun("zcat $pathSpec | $Bin/MultiCodes.pl $extraopt -minQuality 0 -index $index -out $outdir/$index");
      }
    }
    &maybeRun("$Bin/combineBarSeq.pl -all $outdir/test $pool " . join(" ", @codesFiles));
    unless (defined $test) {
	# Write the metadata file
	open(EXPS, ">", "$outdir/exps_table") || die "Cannot write to $outdir/exps_table";
	print EXPS join("\t", qw{SetName Index Person Date_pool_expt_started Description})."\n";
	foreach my $i (0..(scalar(@indexes)-1)) {
	    print EXPS join("\t", "test", $indexes[$i], "someone", "somedate", $desc[$i])."\n";
	}
	close(EXPS) || die "Error writing to $outdir/exps_table";
    }
    &maybeRun("$Bin/BarSeqR.pl -org $org -exps $outdir/exps_table  -pool $pool -indir $outdir -outdir $outdir -genes g/$org/genes.GC test");
}

sub maybeRun($) {
    my ($cmd) = @_;
    if (defined $test) {
        print STDERR "Would run: $cmd\n";
    } else {
	print STDERR "Running: $cmd\n";
        system($cmd) == 0 || die "script failed: $cmd";
    }
}
