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
use List::Util 'max';

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();

# read the input information

my $species = $cgi->param('species') || "";
my $gene = $cgi->param('gene') || die "no gene name found\n";

$gene =~ s/ *$//;
$gene =~ s/^ *//;

print $cgi->header;
print $cgi->start_html(
    -title =>"Fitness",
    -style => {-code => $style},
    -author=>'wjshaoATberkeley.edu',
    -meta=>{'copyright'=>'copyright 2015 UC Berkeley'},
#    -BGCOLOR=>'#fffacd'
);

# check user input

Utils::fail("$gene is invalid. Please enter correct gene name!") unless ($gene =~ m/^[A-Za-z0-9_]*$/);

print $cgi->h2("Fitness Profile");

# connect to database

my $dbh = Utils::get_dbh();
my ($stmt, $sth, $rv, $nickname, $locusId, $sys, $desc, $genus, $g);

# When no species has been specified, get the species name and gene locusId, which may not be unique

if ($species eq "All 17 genomes" || $species eq "") {
    $stmt = qq(SELECT nickname, locusId, sysName, desc, gene FROM Gene WHERE locusId = "$gene" OR sysName = "$gene" OR upper(gene) = upper("$gene"););
    $sth = Utils::execute_db($dbh, $stmt);

    my $cnt = 0;
    my @hits = ();

    while(my @row = $sth->fetchrow_array()) {
        ($nickname,$locusId,$sys,$desc,$g) = @row;

        my $stmt = qq(SELECT genus, species FROM Organism WHERE name = "$nickname";);
        my $sth2 = Utils::execute_db($dbh, $stmt);
        ($genus,$species) = $sth2->fetchrow_array();

        $stmt = qq(SELECT expName, fit, t FROM GeneFitness WHERE species = "$nickname" AND locusId = "$locusId";);
        my $sth3 = Utils::execute_db($dbh, $stmt);
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
        $stmt = qq(SELECT name FROM Organism WHERE species = "$species";);
        $sth = Utils::execute_db($dbh, $stmt);
        ($nickname) = $sth->fetchrow_array();
    } else {
        $stmt = qq(SELECT species, name FROM Organism WHERE name = "$species";);
        $sth = Utils::execute_db($dbh, $stmt);
        ($species, $nickname) = $sth->fetchrow_array();
    }

    # get other information for the gene

    $stmt = qq(SELECT locusId, sysName, gene, desc FROM Gene WHERE (nickname = "$nickname" OR nickname = "$species") AND (locusId = "$gene" OR sysName = "$gene" OR upper(gene) = upper("$gene")););
    $sth = Utils::execute_db($dbh, $stmt);
    ($locusId, $sys, $g, $desc) = $sth->fetchrow_array();

}

# get the fitness data from the database

$stmt = qq(SELECT expName, fit, t FROM GeneFitness WHERE species = "$nickname" AND locusId = "$locusId";);
$sth = Utils::execute_db($dbh, $stmt);

my (@exp,@fit,@t);
my $hasFit = 0;
while(my @row = $sth->fetchrow_array()) {
    push @exp, $row[0];
    push @fit, $row[1];
    push @t, $row[2];
    $hasFit = 1;
}

if ($hasFit > 0) {

    print $cgi->p("Fitness data for gene $gene: $sys $locusId $desc in $species");

    # sort by fitness value

    my @sorted_fit_indexes = sort { abs($fit[$b]) <=> abs($fit[$a]) } 0..$#fit;
    my $limit = $#fit > 19 ? 19 : $#fit;
    my @top_indexes = @sorted_fit_indexes[0..$limit];
    @exp = @exp[@top_indexes];
    @fit = @fit[@top_indexes];
    @t = @t[@top_indexes];
    @sorted_fit_indexes = sort { $fit[$a] <=> $fit[$b] } 0..$#fit;
    @exp = @exp[@sorted_fit_indexes];
    @fit = @fit[@sorted_fit_indexes];
    @t = @t[@sorted_fit_indexes];

    my $len = $limit + 1;
    my @absfit = map { abs($_) } @fit;
    my $sig = sprintf("%.1f",max(@absfit));

    print $cgi->h5("Note: only the top $len experiments with the most significant phenotypes are shown (|fitness| > $sig).");

    # get meta info for the experiments

    my @out;
    foreach my $i (0..$#t) {
        $stmt = qq(SELECT expDesc FROM QualityExperiment WHERE nickname = "$nickname" AND expName = "$exp[$i]";);
        $sth = Utils::execute_db($dbh, $stmt);
        my $desc;
        while(my @row = $sth->fetchrow_array()) {
            ($desc) = @row;
        }

        push @out, $desc;
        push @out, sprintf("%.1f",$fit[$i]);
        push @out, sprintf("%.1f",$t[$i]);
    }
    my @td = ();
    while ( my @elems = splice @out, 0, 3 ) {
        my ($exp, $fit, $t) = @elems;

        # color code fitness value in range (-3 .. 3)

        my $perc;
        if ($fit > 3) {
            $perc = 1;
        } elsif ($fit < -3) {
            $perc = 0;
        } else {
            $perc = ($fit + 3)/6;
        }
        my ($color) = Utils::get_color($perc);

        # color code only the fitness cell 

        my @elem1 = ($exp);
        my @elem2 = ($fit);
        my @elem3 = ($t);
        my $td = $cgi->td( \@elem1 ) . $cgi->td({-bgcolor => "$color" }, \@elem2 ) . $cgi->td( \@elem3 );
        push @td, $td;
    }
    print $cgi->table(
        { -border=>1, cellpadding=>3 },
        $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
            $cgi->th( [ 'experiment','fitness','t score' ] ) ),
            $cgi->Tr( \@td )
    );

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

