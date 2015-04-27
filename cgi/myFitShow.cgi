#!/usr/bin/perl -w
#######################################################
## myFitShow.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Wenjun Shao (wjshao@berkeley.edu) and Morgan Price
#######################################################
#
# Required CGI parameters:
# gene -- a locusId, sysName, or gene name to match on
#	(may show multiple hits)
# Optional CGI parameters:
# orgId -- which organism to search in
# optionshowAll -- 1 if showing all fitness values instead of just the most extreme ones

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

my $orgSpec = $cgi->param('orgId') || "";
my $geneSpec = $cgi->param('gene') || die "no gene name found\n";
my $showAll = $cgi->param('showAll') || 0;

$geneSpec =~ s/ *$//;
$geneSpec =~ s/^ *//;

print $cgi->header;
print $cgi->start_html(
    -title =>"Fitness",
    -style => {-code => $style},
    -author=>'wjshaoATberkeley.edu',
    -meta=>{'copyright'=>'copyright 2015 UC Berkeley'},
);

# check user input

Utils::fail($cgi, "$geneSpec is invalid. Please enter correct gene name!") unless ($geneSpec =~ m/^[A-Za-z0-9_]*$/);

# connect to database

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);

my $query = qq{SELECT orgId, locusId, sysName, desc, gene FROM Gene
		WHERE locusId = ? OR sysName = ? OR upper(gene) = upper(?)};
my $hits;
if ($orgSpec) {
    $query .= " AND orgId = ?";
    $hits = $dbh->selectall_arrayref($query, { Slice => {} }, $geneSpec, $geneSpec, $geneSpec, $orgSpec);
} else {
    $hits = $dbh->selectall_arrayref($query, { Slice => {} }, $geneSpec, $geneSpec, $geneSpec);
}

# and add another column for whether there is fitness data
foreach my $gene (@$hits) {
    $gene->{has_fitness} = Utils::gene_has_fitness($dbh,$gene->{orgId},$gene->{locusId});
}

if (@$hits == 0) {
    print $cgi->p("No gene found for $geneSpec");
} elsif (@$hits > 1) {
    print $cgi->p("Genes found for $geneSpec:");
    my @trows = ();
    push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
			  $cgi->th( [ 'geneId','sysName','geneName','description','genome','fitness' ] ) );
    foreach my $gene (@$hits) {
	my $dest = "myFitShow.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}";
	my $fitString = $gene->{has_fitness} ? qq(<a href="$dest">check data</a>) : "no data";
	my @trow = map $cgi->td($_), ($gene->{locusId}, $gene->{sysName}, $gene->{gene}, $gene->{desc},
				      $orginfo->{$gene->{orgId}}->{genome}, $fitString );
	push @trows, $cgi->Tr(@trow);
    }
    
    print $cgi->table( { -border=>1, cellpadding=>3 }, @trows);
    
} else {
    # just 1 hit
    my $gene = $hits->[0];
    my $orgId = $gene->{orgId};
    my $locusId = $gene->{locusId};

    if ($hits->[0]{has_fitness} == 0) {
	print $cgi->p("Sorry, no fitness data for the gene $geneSpec in " . $orginfo->{$gene->{orgId}}{genome});
    } else {
	# show fitness data for gene
	my @fit = @{ $dbh->selectall_arrayref("SELECT expName,fit,t from GeneFitness where orgId=? AND locusId=?",
					      { Slice => {} },
					      $orgId, $locusId) };
	my $nTotValues = scalar(@fit);
	die "Unreachable" if $nTotValues == 0;
	my $limitRows = $showAll ? $nTotValues : 20;
	my $minAbsFit;
	if ($nTotValues > $limitRows) {
	    # subset the rows
	    @fit = sort { abs($b->{fit}) <=> abs($a->{fit}) } @fit;
	    @fit = @fit[0..($limitRows-1)];
	    $minAbsFit = abs($fit[$#fit]{fit});
	}
	@fit = sort { $a->{fit} <=> $b->{fit} } @fit;

	# and get metadata about experiments
	my $expinfo = Utils::expinfo($dbh,$orgId);

	print $cgi->h2("Fitness Profile");
	print $cgi->p("Fitness data for gene $geneSpec: $gene->{sysName} $locusId $gene->{desc} in $orginfo->{$orgId}{genome}");
	if (@fit < $nTotValues) {
	    $minAbsFit = sprintf("%.1f", $minAbsFit);
	    print $cgi->h5("Note: only the top $limitRows experiments with the most significant phenotypes are shown (|fitness| &ge; $minAbsFit)");
	}

	my @out = (); # specifiers for HTML rows for each fitness value
	foreach my $fitrow (@fit) {
	    my $expName = $fitrow->{expName};
	    push @out, join(" ",
			    $cgi->td( $expinfo->{$expName}{expDesc} ),
			    $cgi->td( { -bgcolor => Utils::fitcolor($fitrow->{fit}) },
				      sprintf("%.1f", $fitrow->{fit}) ),
			    $cgi->td( sprintf("%.1f", $fitrow->{t}) ));
	}
	print $cgi->table(
	    { -border=>1, cellpadding=>3 },
	    $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
		     $cgi->th( [ 'experiment','fitness','t score' ] ) ),
            $cgi->Tr( \@out ) );
	# links
	if ($showAll == 0) {
	    my $showAllDest = qq(myFitShow.cgi?orgId=$orgId&gene=$locusId&showAll=1);
	    print $cgi->h4(qq(<a href=$showAllDest>Show all fitness data</a>));
	} else {
	    my $showFewDest = qq(myFitShow.cgi?orgId=$orgId&gene=$locusId&showAll=0);
	    print $cgi->h4(qq(<a href=$showFewDest>Show fewer fitness data</a>));
	}
	print $cgi->h4($cgi->a({href => "cofit.cgi?orgId=$orgId&locusId=$locusId"}, "Top cofit genes"));
    } # end if gene has fitness data
    # check homologs where or not gene has data
    my $dest = qq(mySeqSearch.cgi?orgId=$orgId&locusId=$locusId);
    print $cgi->h4(qq(<a href=$dest>Check homologs</a>));
}
$dbh->disconnect();
Utils::endHtml($cgi);
