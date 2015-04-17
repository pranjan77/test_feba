#!/usr/bin/perl -w
# Add cofit information to a sqlite database.
use strict;
use Getopt::Long;

my $usage = <<END
db_add_cofit.pl [ -db sqlite_database | -out cofit.tab ] -html html_parent_directory nickname1 ... nicknameN
	For each organism nickname, loads the cofit table from parent_directory/nickname/cofit
	into the sqlite_database, or else just saves the import table into the file.
END
    ;

{
    my ($db, $htmldir, $outfile) = undef;
    die $usage
	unless GetOptions('db=s' => \$db, 'out=s' => \$outfile, 'html=s' => \$htmldir)
	&& @ARGV > 0
	&& (defined $db xor defined $outfile) && defined $htmldir;
    my @orgs = @ARGV;
    die "No such file: $db" if defined $db && ! -e $db;
    $outfile = "$db.cofit.$$" if !defined $outfile;
    die "No such directory: $htmldir" unless -d $htmldir;

    foreach my $org (@orgs) {
	die "No such directory: $htmldir/$org" unless -d "$htmldir/$org";
	die "Not a successful output directory, no $htmldir/$org/.FEBA.success"
	    unless -e "$htmldir/$org/.FEBA.success";
    }

    open(OUT, ">", $outfile) || die "Error writing to $outfile";
    # Column order: nickname locusId hitId rank cofit
    my $nOrgSuccess = 0;
    my $nRows = 0;
    foreach my $org (@orgs) {
	my $cofit_file = "$htmldir/$org/cofit";
	if (! -e $cofit_file) {
	    print STDERR "Skipping $org, no cofit file $cofit_file (too few experiments?)\n";
	    next;
	}
	# else
	open(COFIT, "<", $cofit_file) || die "Cannot read $cofit_file";
	my $header = <COFIT>;
	chomp $header;
	my @colnames = split /\t/, $header;
	my %colnames = map { $colnames[$_] => $_ } (0..((@colnames)-1));
	foreach my $field (qw{locusId hitId rank cofit}) {
	    die "Invalid header in $cofit_file\n$header\nNo field named $field"
		unless exists $colnames{$field};
	}
	while(<COFIT>) {
	    chomp;
	    my @F = split /\t/, $_, -1;
	    foreach my $field (qw{locusId hitId rank cofit}) {
		die "Invalid line\n$_\nno value for $field"
		    unless defined $F[$colnames{$field}] && $F[$colnames{$field}] ne "";
	    }
	    print OUT join("\t", $org,
			   $F[$colnames{locusId}],
			   $F[$colnames{hitId}], 
			   $F[$colnames{rank}], 
			   $F[$colnames{cofit}])."\n";
	    $nRows++;
	}
	close(COFIT) || die "Error reading $cofit_file";
	$nOrgSuccess++;
    }
    print STDERR "Read cofit for $nOrgSuccess organisms of " . scalar(@orgs)."\n";
    print STDERR "Wrote $nRows rows to $outfile\n";
    close(OUT) || die "Error writing to $outfile";

    if (defined $db) {
	my $db_cmds =<<END
DROP TABLE IF EXISTS Cofit;
CREATE TABLE Cofit(
	nickname TEXT NOT NULL,
        locusId TEXT NOT NULL,
	hitId TEXT NOT NULL,
        rank INT NOT NULL,
        cofit REAL NOT NULL,
        PRIMARY KEY (nickname,locusId,hitId)
);
END
    ;
	$db_cmds .= qq{.separator "\t"\n};
	$db_cmds .=".import $outfile Cofit\n";

        open(SQLITE, "|-", "sqlite3", "$db") || die "Cannot run sqlite3 on $db";
	print SQLITE $db_cmds;
	close(SQLITE) || die "Error running sqlite3";
        my $db_rows = `echo "SELECT COUNT(*) FROM Cofit;" | sqlite3 $db`;
	chomp $db_rows;
        if ($db_rows == $nRows) {
	    print STDERR "Success: wrote $db_rows rows to Cofit table in $db\n";
	} else {
	    die "Error: Cofit table in $db has $db_rows rows\n";
	}
    }
}





