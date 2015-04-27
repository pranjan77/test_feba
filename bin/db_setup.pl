#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils; # for ReadTable()
use DBI;

my $usage = <<END
Usage: db_setup.pl [ -db db_file_name ] -orginfo orginfo -indir htmldir nickname1 ... nicknameN
    Sets up a sqlite database suitable for the cgi scripts, reading from
    the html directories indir/nickname1 ... indir/nicknameN
    (created by BarSeqR.pl)
    Does not set up the fasta database for the CGI.

    The orginfo file should include all of the columns for the Orginfo
    table, but may contain organisms that are not included in the
    command line -- these will not be removed before loading.

    Creates intermediate files db.tablename.org in the working directory.
    If no db is specified, it just creates these and writes to stdout the
    commands to be used to import them into a database (that already has
    the tables set up by using lib/db_setup_tables.sql; otherwise,
    it uses them to load the database and deletes the files.
END
;

# global variables for maintaining the work list
my %workRows = (); # table name to number of rows
my $workFile;
my $workTable;
my @workFiles = ();
my @workCommands = ();

sub StartWorkFile($$); # table, file
sub StartWork($$); # table, organism
sub EndWork(); # 
sub WorkPutRow(@); # list of fields to be written as a line
sub WorkPutHash($@); # hash and list of fields to use

{
    my ($dbfile,$indir,$orgfile);
    (GetOptions('db=s' => \$dbfile,
		'orginfo=s' => \$orgfile,
		'indir=s' => \$indir)
     && defined $indir && defined $orgfile)
	|| die $usage;
    my @orgs = @ARGV;
    die "No such directory: $indir" unless -d $indir;
    die "No such file: $orgfile" unless -e $orgfile;
    die "No organism nicknames specified\n" if scalar(@orgs) < 1;

    foreach my $org (@orgs) {
	die "No such directory: $indir/$org" unless -d "$indir/$org";
	foreach my $file (qw{.FEBA.success genes expsUsed fit_quality.tab fit_logratios_good.tab fit_t.tab}) {
	    die "Missing file: $indir/$org/$file" unless -e "$indir/$org/$file";
	}
    }

    print STDERR "Reading " . scalar(@orgs) . " organisms from $indir\n";
    if (defined $dbfile) {
	my $sqlfile = "$Bin/../lib/db_setup_tables.sql";
	die "No such file: $sqlfile" unless -e $sqlfile;
	print STDERR "Creating seqlite3 database $dbfile\n";
	unlink($dbfile) if -e $dbfile; # start with an empty database
	system("sqlite3 $dbfile < $sqlfile") == 0 || die $!;
    }
    else {
	print STDERR "Writing test files\n";
    }

    # Load orginfo
    my @orgfields = qw{name division genus species strain taxonomyId};
    my @orginfo = &ReadTable($orgfile, @orgfields);
    my %orgSeen = map { $_->{name} => 1 } @orginfo;
    my %orgUse = map { $_ => 1 } @orgs;
    foreach my $org (@orgs) {
	die "Organism $org is not described in $orgfile"
	    unless exists $orgSeen{$org};
    }

    # Create db.Organism
    StartWorkFile("Organism", "db.Organism");
    foreach my $row (@orginfo) {
	$row->{orgId} = $row->{name};
	WorkPutHash($row, qw{orgId division genus species strain taxonomyId})
	    if exists $orgUse{$row->{orgId}};
    }
    EndWork();

    # Create db.Gene.*
    foreach my $org (@orgs) {
	my @genes = &ReadTable("$indir/$org/genes",
			       qw{locusId sysName scaffoldId begin end desc name GC nTA});
	die "No genes for $org" unless @genes > 0;
	StartWork("Gene", $org);
	foreach my $row (@genes) {
	    $row->{orgId} = $org;
	    $row->{gene} = $row->{name};
	    $row->{type} = 1 if !exists $row->{type};
	    WorkPutHash($row, qw{orgId locusId sysName scaffoldId begin end type strand gene desc GC nTA});
	}
	EndWork();
    }

    # Create db.Experiment.*
    foreach my $org (@orgs) {
	my @q = &ReadTable("$indir/$org/fit_quality.tab",
			   qw{name short t0set num nMapped nPastEnd nGenic nUsed gMed gMedt0 gMean
				cor12 mad12 mad12c mad12c_t0 opcor adjcor gccor maxFit u});
	# exps has additional fields; not included in schema yet; they can be lined up by name
	my @exps = &ReadTable("$indir/$org/expsUsed",
			      qw{name SetName Date_pool_expt_started Person Mutant.Library Description
                                 gDNA.plate gDNA.well Index Media Growth.Method Group
                                 Condition_1 Units_1 Concentration_1
                                 Condition_2 Units_2 Concentration_2});
	my %exps = map {$_->{name} => $_} @exps;
	StartWork("Experiment",$org);
	foreach my $row (@q) {
	    $row->{"orgId"} = $org;
	    $row->{"expName"} = $row->{"name"};
	    $row->{"expDesc"} = $row->{"short"};
	    $row->{"timeZeroSet"} = $row->{"t0set"};

	    my $id = $row->{name};
	    my $exp = $exps{$id} || die "No matching metadata for experiment $org $id";
	    $row->{expDescLong} = $exp->{Description};
	    $row->{mutantLibrary} = $exp->{"Mutant.Library"};
	    $row->{expGroup} = $exp->{Group};
	    $row->{dateStarted} = $exp->{Date_pool_expt_started};
	    $row->{setName} = $exp->{SetName};
	    $row->{seqindex} = $exp->{Index};
	    foreach my $field (qw{Person Media Condition_1 Units_1 Concentration_1
                                 Condition_2 Units_2 Concentration_2}) {
		$row->{lc($field)} = $exp->{$field} eq "NA" ? "" : $exp->{$field};
	    }
	    WorkPutHash($row, qw{orgId expName expDesc timeZeroSet num nMapped nPastEnd nGenic
                                 nUsed gMed gMedt0 gMean cor12 mad12 mad12c mad12c_t0
                                 opcor adjcor gccor maxFit
				 expGroup expDescLong mutantLibrary person dateStarted setName seqindex media
				 condition_1 units_1 concentration_1
				 condition_2 units_2 concentration_2})
		if $row->{u} eq "TRUE";
	}
	EndWork();
    }

    # Create db.GeneFitness.*
    foreach my $org (@orgs) {
	my $fit_file = "$indir/$org/fit_logratios_good.tab";
	my @fitNames = &ReadColumnNames($fit_file);
	my @colNamesFull = grep m/^set/, @fitNames;
	my @colNames = @colNamesFull;
	map { s/ .*$//; } @colNames; # skip annotations after the name
	my @fit = &ReadTable($fit_file, "locusId");
	my @t = &ReadTable("$indir/$org/fit_t.tab",grep {$_ ne "comb"} @fitNames);
	StartWork("GeneFitness",$org);
	foreach my $iRow (0..$#fit) {
	    my $fit_row = $fit[$iRow];
	    my $t_row = $t[$iRow];
	    my $locusId = $fit_row->{locusId};
	    die "Mismatched locusIds in $fit_file and fit_t.tab, row $iRow"
		unless $locusId eq $t_row->{locusId};
	    foreach my $iCol (0..$#colNamesFull) {
		my $colFull = $colNamesFull[$iCol];
		die "Missing from fit: $iCol $colFull" unless defined $fit_row->{$colFull};
		die "Missing from t: $iCol $colFull" unless defined $t_row->{$colFull};
		die "Missing from colNames: $iCol $colFull" unless defined $colNames[$iCol];
		WorkPutRow($org, $locusId, $colNames[$iCol],
			   $fit_row->{$colFull}, $t_row->{$colFull});
	    }
	}
	EndWork();
    }

    # Create db.Cofit.*
    foreach my $org (@orgs) {
	my $cofitFile = "$indir/$org/cofit";
	if (! -e $cofitFile) {
	    print STDERR "No cofit files for $org (usually means too few experiments)\n";
	    next;
	}
	my @cofit = &ReadTable($cofitFile, qw{locusId hitId rank cofit});
	StartWork("Cofit",$org);
	foreach my $row (@cofit) {
	    $row->{orgId} = $org;
	    WorkPutHash($row, qw{orgId locusId hitId rank cofit});
	}
	EndWork();
    }

    if (defined $dbfile) {
	# Run the commands
	open(SQLITE, "|-", "sqlite3", "$dbfile") || die "Cannot run sqlite3 on $dbfile";
	print SQLITE ".bail on\n";
	print SQLITE ".mode tabs\n";
	print STDERR "Loading tables\n";
	foreach my $workCommand (@workCommands) {
	    print SQLITE "$workCommand\n";
	}
	close(SQLITE) || die "Error running sqlite3 commands\n";

	# Check #rows in each table
	my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", "", "", { RaiseError => 1 }) || die $DBI::errstr;
	while (my ($table, $nRowsExpect) = each %workRows) {
	    my ($nRowsActual) = $dbh->selectrow_array("SELECT COUNT(*) FROM $table")
		|| die "counting rows in $table failed";
	    die "Failed to load $table: expect $nRowsExpect rows but see $nRowsActual rows instead\n"
		unless $nRowsActual == $nRowsExpect;
	}
	$dbh->disconnect();
	print STDERR "Successfully created $dbfile\n";
	print STDERR "Cleaning up\n";

	# Delete the files
	foreach my $file (@workFiles) {
	    unlink($file);
	}
    }
}

sub StartWorkFile($$) {
    my ($table,$file) = @_;
    die "Already working on $workTable $workFile when calling StartWorkFile with $table $file"
	if defined $workTable || defined $workFile;
    $workTable = $table;
    $workFile = $file;
    push @workFiles, $workFile;
    push @workCommands, ".import $workFile $workTable";
    open(WORK, ">", $workFile) || die "Cannot write to $workFile";
    $workRows{$workTable} = 0 if !exists $workRows{$workTable};
}

sub StartWork($$) {
    my ($table,$org) = @_;
    StartWorkFile($table,"db.$table.$org");
}

sub EndWork() {
    die "Illegal call to EndWork()" unless defined $workFile && defined $workTable;
    close(WORK) || die "Error writing to $workFile";
    print STDERR "Wrote $workFile\n";
    $workFile = undef;
    $workTable = undef;
}

sub WorkPutRow(@) {
    die "Illegal call to WorkPutRow" unless defined $workFile && defined $workTable;
    print WORK join("\t", @_)."\n" || die "Error writing to $workFile";
    $workRows{$workTable}++;
}

sub WorkPutHash($@) {
    my ($row,@fields) = @_;
    foreach my $field (@fields) {
	die "No such field $field" unless exists $row->{$field};
    }
    WorkPutRow(map $row->{$_}, @fields);
}

