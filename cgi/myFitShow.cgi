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
my ($stmt, $sth, $rv, $nickname, $locusId, $sys, $desc, $genus, $strain, $g);

my $orginfo = Utils::orginfo($dbh);

# When no species has been specified, get the species name and gene locusId, which may not be unique

if ($species eq "All" || $species eq "") {
    my $hits = $dbh->selectall_arrayref(qq{SELECT nickname, locusId, sysName, desc, gene FROM Gene
                                           WHERE locusId = ? OR sysName = ? OR upper(gene) = upper(?)},
	                                undef, $gene, $gene, $gene)
	|| die;

    # and add another column for whether there is fitnes data
    foreach my $row (@$hits) {
    	my ($orgId, $locusId, $sysName, $desc, $geneName) = @$row;
	push @$row, Utils::gene_has_fitness($dbh,$orgId,$locusId);
    }

    if (@$hits == 0) {
        print $cgi->p("No gene found for $gene");
        $dbh->disconnect();
        Utils::endHtml($cgi);
    } elsif (@$hits > 1) {
        print $cgi->p("Genes found for $gene:");
	my @trows = ();
	push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
			      $cgi->th( [ 'geneId','sysName','geneName','description','genome','fitness' ] ) );
	foreach my $row (@$hits) {
	    my ($orgId, $locusId, $sysName, $desc, $geneName, $hasData) = @$row;
	    my $dest = "myFitShow.cgi?species=$orgId&gene=$locusId";
	    my $fitString = $hasData ? qq(<a href="$dest">check data</a>) : "no data";
	    my @trow = map $cgi->td($_), ($locusId, $sysName, $geneName, $desc,
					$orginfo->{$orgId}->{genome}, $fitString );
	    push @trows, $cgi->Tr(@trow);
	}

        print $cgi->table( { -border=>1, cellpadding=>3 }, @trows);

        $dbh->disconnect();
        Utils::endHtml($cgi);
    }
    # else fall through to 1-gene case
    ($nickname, $locusId, $sys, $desc, $g, undef) = @{ $hits->[0] };
    $species = $orginfo->{$nickname}->{genome}; # description of organism
} else {
    $nickname = $species;
    $species = $orginfo->{$nickname}->{genome}; # description of organism

    # get other information for the gene
    $stmt = qq(SELECT locusId, sysName, gene, desc FROM Gene WHERE nickname = ? AND (locusId = ? OR sysName = ? OR upper(gene) = upper(?)););
    $sth = $dbh->prepare( $stmt );
    $rv = $sth->execute($nickname, $gene, $gene, $gene) or die $DBI::errstr;
    print $DBI::errstr if($rv < 0);

    ($locusId, $sys, $g, $desc) = $sth->fetchrow_array();
    if (!defined $locusId) {
	print $cgi->p("No gene found for $gene in $species");
        $dbh->disconnect();
        Utils::endHtml($cgi);
    }
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

    print $cgi->p("Fitness data for gene $gene: $sys $locusId $desc in $genus $species $strain");

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
}

$dbh->disconnect();
Utils::endHtml($cgi);

exit 0;

# END

#----------------------------------------
