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
    Does not set up the fasta database.

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

{
    my ($dbfile,$indir,$orgfile);
    my @orgs;
    (GetOptions('db=s' => \$dbfile,
		'orginfo=s' => \$orgfile
		'indir=s' => \$indir)
     && defined $indir && defined $orgfile)
	|| die $usage;
    my @orgs = @ARGV;
    die "No such directory: $indir" unless -d $indir;
    die "No such file: $orgfile" unless -e $orgfile;
    die "No organism nicknames specified" if scalar(@orgs) < 1;
    foreach my $org (@orgs) {
	die "No such directory: $indir/$org" unless -d "$indir/$org";
	foreach my $file (qw{.FEBA.success genes expsUsed fit_quality.tab fit_logratios_good.tab fit_t.tab}) {
	    die "Missing file: $indir/$org/$file" unless -e "$indir/$org/$file";
	}
    }
    print STDERR "Reading " . scalar(@orgs) . " organisms from $indir\n";
    if (defined $db) {
	my $sqlfile = "$Bin/../lib/db_setup_tables.sql";
	die "No such file: $sqlfile" unless -e $sqlfile;
	print STDERR "Creating seqlite3 database $dbfile\n";
	unlink($db) if -e $db;
	system("sqlite3 $db < $sqlfile");
    }
    else {
	print STDERR "Writing test files\n";
    }

    my @commands = ();
    my @files = ();
    my $file;

    # Organisms table
    my @orgfields = qw{name division genus species strain taxonomyId};
    my @orginfo = &ReadTable($orgfile, @orgfields);
    my %orgSeen = map { $_->{name} => 1 } @orginfo;
    foreach my $org (@orgs) {
	die "Organism $org is not described in $orgfile"
	    unless exists $orgSeen{$org};
    }

    $file = "db.orginfo";
    open(FILE, ">", $file) || die "Error writing to $file";
    foreach my $row (@orginfo) {
	print FILE join("\t", @$row->{ @orgfields })."\n";
    }
    close(FILE) || die "Error writing to $file";
    print STDERR "Wrote $file\n";
    push @files, $file;
    push @commands, ".import db.orginfo Organism";
}
