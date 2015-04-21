#!/usr/bin/perl -w
#######################################################
## myFitShow.cgi
##
## Copyright (c) 2015 All Right Reserved by UC Berkeley
##
## Authors:
## Wenjun Shao (wjshao@berkeley.edu) and Morgan Price
#######################################################

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use List::Util 'min';

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();

# read the input information

my $species = $cgi->param('species') || "";
my $gene = $cgi->param('gene') || die "no gene name found\n";
my $showAll = $cgi->param('showAll') || 0;

$gene =~ s/ *$//;
$gene =~ s/^ *//;

print $cgi->header;
print $cgi->start_html(
    -title =>"Fitness",
    -style => {-code => $style},
    -author=>'wjshaoATberkeley.edu',
    -meta=>{'copyright'=>'copyright 2015 UC Berkeley'},
);

# check user input

Utils::fail($cgi, "$gene is invalid. Please enter correct gene name!") unless ($gene =~ m/^[A-Za-z0-9_]*$/);

print $cgi->h2("Fitness Profile");

# connect to database

my $dbh = Utils::get_dbh();
my ($stmt, $sth, $rv, $nickname, $locusId, $sys, $desc, $genus, $g);

# When no species has been specified, get the species name and gene locusId, which may not be unique

if ($species eq "All 17 genomes" || $species eq "") {
    $stmt = qq(SELECT nickname, locusId, sysName, desc, gene FROM Gene WHERE locusId = ? OR sysName = ? OR upper(gene) = upper(?););
    $sth = $dbh->prepare( $stmt );
    $rv = $sth->execute($gene,$gene,$gene) or die $DBI::errstr;
    print $DBI::errstr if($rv < 0);

    my $cnt = 0;
    my @hits = ();

    while(my @row = $sth->fetchrow_array()) {
        ($nickname,$locusId,$sys,$desc,$g) = @row;

        my $stmt = qq(SELECT genus, species FROM Organism WHERE name = ?;);
        my $sth2 = $dbh->prepare( $stmt );
        $rv = $sth2->execute($nickname) or die $DBI::errstr;
        print $DBI::errstr if($rv < 0);

        ($genus,$species) = $sth2->fetchrow_array();

        $stmt = qq(SELECT expName, fit, t FROM GeneFitness WHERE species = ? AND locusId = ?;);
        my $sth3 = $dbh->prepare( $stmt );
        $rv = $sth3->execute($nickname, $locusId) or die $DBI::errstr;
        print $DBI::errstr if($rv < 0);

        my $fitness = "no data";
        if (my @row3 = $sth3->fetchrow_array()) {
            my $dest = "myFitShow.cgi?species=$nickname&gene=$locusId";
            $fitness = qq(<a href=$dest>check data</a>);
        }

        my @hit = ($locusId,$sys,$g,$desc,$species,$fitness);
        push @hits, @hit;
        $cnt++;
    }

    if ($cnt > 1) {

        print $cgi->p("Gene(s) found for $gene:");

        my @td = ();
        while ( my @elems = splice @hits, 0, 6 ) {
            push @td, $cgi->td( \@elems );
        }
        print $cgi->table(
            { -border=>1, cellpadding=>3 },
            $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
                $cgi->th( [ 'geneId','sysName','geneName','description','species','fitness' ] ) ),
                $cgi->Tr( \@td )
        );

        $dbh->disconnect();
        Utils::endHtml($cgi);

    } elsif ($cnt == 0) {
        print $cgi->p("No gene found for $gene!");
        $dbh->disconnect();
        Utils::endHtml($cgi);
    }

} else {

    # full name was received, get the nickname for the species

    my @t = split(" ",$species);

    if ($#t >= 1) {
        $stmt = qq(SELECT name FROM Organism WHERE species = ?;);
        $sth = $dbh->prepare( $stmt );
        $rv = $sth->execute($species) or die $DBI::errstr;
        print $DBI::errstr if($rv < 0);

        ($nickname) = $sth->fetchrow_array();
    } else {
        $stmt = qq(SELECT species, name FROM Organism WHERE name = ?;);
        $sth = $dbh->prepare( $stmt );
        $rv = $sth->execute($species) or die $DBI::errstr;
        print $DBI::errstr if($rv < 0);

        ($species, $nickname) = $sth->fetchrow_array();
    }

    # get other information for the gene

    $stmt = qq(SELECT locusId, sysName, gene, desc FROM Gene WHERE (nickname = ? OR nickname = ?) AND (locusId = ? OR sysName = ? OR upper(gene) = upper(?)););
    $sth = $dbh->prepare( $stmt );
    $rv = $sth->execute($nickname, $species, $gene, $gene, $gene) or die $DBI::errstr;
    print $DBI::errstr if($rv < 0);

    ($locusId, $sys, $g, $desc) = $sth->fetchrow_array();

}

# get the fitness data from the database

$stmt = qq(SELECT expName, fit, t FROM GeneFitness WHERE species = ? AND locusId = ?;);
$sth = $dbh->prepare( $stmt );
$rv = $sth->execute($nickname, $locusId) or die $DBI::errstr;
print $DBI::errstr if($rv < 0);

my @outRows = ();
my $fitCnt = 0;
while(my @row = $sth->fetchrow_array()) {
    push @outRows, {exp => $row[0], fit => $row[1], t => $row[2]};
    $fitCnt += 1;
}

if ($fitCnt > 0) {

    print $cgi->p("Fitness data for gene $gene: $sys $locusId $desc in $species");

    # sort by fitness value

    my $limit = $#outRows > 19 ? 19 : $#outRows;
    $limit = $showAll ? $#outRows : $limit;
    @outRows = sort { abs($b->{ 'fit' }) <=> abs($a->{ 'fit' }) } @outRows;
    @outRows = @outRows[0..$limit];
    @outRows = sort { $a->{ 'fit' } <=> $b->{ 'fit' } } @outRows;

    my @absfit = map { abs($_->{ 'fit' }) } @outRows;
    my $sig = sprintf("%.1f",min(@absfit));

    my $len = $limit + 1;

    if (!$showAll) {
        print $cgi->h5("Note: only the top $len experiments with the most significant phenotypes are shown (|fitness| > $sig).");
    } else {
        print $cgi->h5("All experiments are shown here.");
    }

    # get meta info for the experiments

    my @out;
    foreach my $i (0..$#outRows) {
        $stmt = qq(SELECT expDesc FROM QualityExperiment WHERE nickname = ? AND expName = ?;);
        $sth = $dbh->prepare( $stmt );
        $rv = $sth->execute($nickname, $outRows[$i]->{ 'exp' }) or die $DBI::errstr;
        print $DBI::errstr if($rv < 0);

        my $desc;
        while(my @row = $sth->fetchrow_array()) {
            ($desc) = @row;
        }

        push @out, {exp => $desc, fit => sprintf("%.1f",$outRows[$i]->{ 'fit' }), t => sprintf("%.1f",$outRows[$i]->{ 't' })};
    }
    my @td = ();
    foreach my $elem ( @out ) {

        # color code fitness value in range (-3 .. 3)

        my $fit = $elem->{ 'fit' };
        my $perc = ($fit + 3)/6;
        $perc = $perc > 1 ? 1 : $perc;
        $perc = $perc < 0 ? 0 : $perc;

        my ($color) = Utils::get_color($perc);

        # color code only the fitness cell 

        my $td = $cgi->td( [ $elem->{ 'exp' } ] ) . $cgi->td({-bgcolor => "$color" }, [ $elem->{ 'fit' } ] ) . $cgi->td( [ $elem->{ 't' } ] );
        push @td, $td;
    }
    print $cgi->table(
        { -border=>1, cellpadding=>3 },
        $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
            $cgi->th( [ 'experiment','fitness','t score' ] ) ),
            $cgi->Tr( \@td )
    );

    if ($showAll == 0) {
        my $showAllDest = qq(myFitShow.cgi?species=$nickname&gene=$locusId&showAll=1);
        print $cgi->h4(qq(<a href=$showAllDest>Show all fitness data</a>));
    } else {
        my $showFewDest = qq(myFitShow.cgi?species=$nickname&gene=$locusId&showAll=0);
        print $cgi->h4(qq(<a href=$showFewDest>Show fewer fitness data</a>));
    }
    print $cgi->h4($cgi->a({href => "cofit.cgi?species=$nickname&gene=$locusId"}, "Top cofit genes"));

    my $dest = qq(mySeqSearch.cgi?species=$nickname&gene=$locusId);
    print $cgi->h4(qq(<a href=$dest>Check homologs</a>));

} else {

    print $cgi->p("Sorry, no fitness data for the gene $gene.");
    print $cgi->h5(q(Note: if gene locusId or sysName were used, please choose "All 17 genomes" or matched species.));

}

$dbh->disconnect();
Utils::endHtml($cgi);

exit 0;

# END

#----------------------------------------
