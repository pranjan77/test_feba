#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use FEBA_Utils; # for ReadTable()
use DBI;

my $usage = <<END
Usage: db_setup.pl -reanno reannotation_file -db db_file_name

    reannotation_file should be tab-delimited and should include the
    fields orgId, locusId, new_annotation, and comment.
END
    ;

{
    my ($reannofile,$db);
    die $usage unless GetOptions('reanno=s' => \$reannofile,
                                 'db=s' => \$db)
        && @ARGV == 0
        && defined $reannofile
        && defined $db;
    my @annos = ReadTable($reannofile, ['orgId', 'locusId', 'new_annotation', 'comment']);
    print STDERR "Read " . scalar(@annos) . " rows from $reannofile\n";
    my $dbh = DBI->connect("dbi:SQLite:dbname=$db", "", "", { RaiseError => 1 }) || die $DBI::errstr;
    my $clear_statement = $dbh->prepare(qq{ DELETE FROM Reannotation; });
    $clear_statement->execute() || die "Cannot clear the Reannotation table in $db";
    my $insert_statement = $dbh->prepare(qq{ INSERT INTO Reannotation VALUES (?,?,?,?); });
    my $orgs = $dbh->selectcol_arrayref("SELECT orgId FROM Organism");
    my %orgs = map { $_ => 1 } @$orgs;
    my $nAdd = 0;
    my %skipOrgs = ();
    foreach my $row (@annos) {
      if (! exists $orgs{ $row->{orgId} }) {
        $skipOrgs{ $row->{orgId} } = 1;
      } else {
        my ($n) = $dbh->selectrow_array("SELECT COUNT(*) FROM Gene WHERE orgId = ? AND locusId = ?",
                                        {}, $row->{orgId}, $row->{locusId});
        if ($n == 0) {
            print STDERR "Warning: reannotation of non-existent gene $row->{locusId} in $row->{orgId}\n";
        } else {
            $insert_statement->execute($row->{orgId}, $row->{locusId}, $row->{new_annotation}, $row->{comment})
                || die "Failed to insert into Reannotation";
            $nAdd++;
        }
      }
    }
    $dbh->disconnect();
    print STDERR "Cleared Reannotation table and added $nAdd rows\n";
    print STDERR "Ignored reannotations for unknown organisms: " . join(" ", sort keys %skipOrgs) . "\n"
      if keys(%skipOrgs) > 0;
}
