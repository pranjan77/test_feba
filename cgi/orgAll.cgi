#!/usr/bin/perl -w

#######################################################
## orgGroup.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo and Morgan Price
#######################################################
#
# This page shows a table of organisms, along with how many conditions
# of fitness data are available.
#
# CGI parameters -- none

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;


my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);

# main table
# gather data and slice it into an array of hashes
# my $cond = $dbh->selectall_arrayref(qq{SELECT orgId, expGroup, COUNT(DISTINCT condition_1) AS nCond, COUNT(*) as nExp FROM Experiment GROUP BY orgId,expGroup ORDER BY orgId; },
    # { Slice => {} });
# my $cond = $dbh->selectall_arrayref(qq{SELECT *, COUNT(DISTINCT expDesc) AS nCond, COUNT(*) as nExp FROM Experiment WHERE expGroup NOT IN ('carbon source', 'nitrogen source', 'stress') GROUP BY orgId ORDER BY orgId; },
#     { Slice => {} });

my $all = $dbh->selectall_arrayref(qq{SELECT e.orgId, o.division
	FROM "Experiment" as e 
	JOIN "Organism" as o on e.orgId = o.orgId
	GROUP BY e.orgId ORDER BY o.genus, o.species, o.strain; },
    { Slice => {} });

my $count = $dbh->selectall_hashref(qq{SELECT orgId, COUNT(DISTINCT expDesc) AS nCond
	FROM "Experiment"
	WHERE expGroup NOT IN ('carbon source', 'nitrogen source', 'stress') 
	GROUP BY orgId  },
    'orgId');

# my $total = $dbh->selectall_hashref(qq{SELECT *, COUNT(*) as nExp FROM Experiment GROUP BY orgId ORDER BY genus, species, strain; },
#     'orgId');

my $total = $dbh->selectall_hashref(qq{SELECT e.*, COUNT(*) as nExp FROM 
	"Experiment" as e 
	JOIN "Organism" as o on e.orgId = o.orgId
	GROUP BY e.orgId ORDER BY o.genus, o.species, o.strain; },
    'orgId');
# my $cond = $dbh->selectall_arrayref(qq{SELECT orgId, expGroup, COUNT(DISTINCT condition_1) AS nCond, COUNT(*) as nExp FROM Experiment GROUP BY orgId ORDER BY orgId, expGroup; },
#     { Slice => {} });

my $nNitrogen = $dbh->selectall_hashref(qq{SELECT orgId, COUNT(DISTINCT Condition_1) AS nCon FROM Experiment WHERE expGroup='nitrogen source' AND NOT Condition_1='' GROUP BY orgId; },
    'orgId');
my $nCarbon = $dbh->selectall_hashref(qq{SELECT orgId, COUNT(DISTINCT Condition_1) AS nCon FROM Experiment WHERE expGroup='carbon source' AND NOT Condition_1='' GROUP BY orgId; },
    'orgId');
my $nStress = $dbh->selectall_hashref(qq{SELECT orgId, COUNT(DISTINCT Condition_1) AS nCon FROM Experiment WHERE expGroup='stress' AND NOT Condition_1='' GROUP BY orgId; },
    'orgId');



# write the title
my $title = scalar(keys %$orginfo) . " Organisms";
my $start = Utils::start_page("$title");

print
    header,
    $start, '<div id="ntcontent">',
    h2("$title");

#exit if no results
# Utils::fail($cgi, "No experiments for this organism.") if @$cond == 0;

#create table
# my @headings = qw{Organism Experiments Carbon Nitrogen Stresses Other};
my @trows = ( Tr({-valign => "top", -align => 'center'}, th([ 'Organism', 'Division', 'Experiments', 'Carbon Sources', 'Nitrogen Sources', 'Stress Compounds', 'Other Conditions'])));
#Tr({ -valign => 'top', -align => 'center' }, map { th($_) } \@headings) );
foreach my $row (@$all) {
	my $orgId = $row->{orgId};
	# my $other = $row->{nCond} - $nCarbon->{$row->{orgId}}->{nCon} - $nNitrogen->{$row->{orgId}}->{nCon} - $nStress->{$row->{orgId}}->{nCon};
    push @trows, Tr({ -valign => 'top', -align => 'right' },
    	# display result row by row
	     td ({-align=>'left'}, [a({href=>"org.cgi?orgId=$row->{orgId}"},$orginfo-> {$row->{orgId}}{genome}), $row->{division}]), #organism
		 	# $row->{expGroup} || "unrecorded", #group
		 td([#$row->{nCond}, #conditions
		 	a( { href => "org.cgi?orgId=$row->{orgId}"},
		    $total->{$row->{orgId}}->{nExp}),#$row->{nExp} ), #experiments
		 	a({href=>"exps.cgi?orgId=$orgId&expGroup=carbon%20source"}, $nCarbon->{$row->{orgId}}->{nCon}),
		 	a({href=>"exps.cgi?orgId=$orgId&expGroup=nitrogen%20source"}, $nNitrogen->{$row->{orgId}}->{nCon}), 
		 	a({href=>"exps.cgi?orgId=$orgId&expGroup=stress"}, $nStress->{$row->{orgId}}->{nCon}),
		 	$count->{$row->{orgId}}->{nCond},
		 ]));
}

print table({cellspacing => 0, cellpadding => 3}, @trows);
print "<br><br>";


# display number of genes that have data out of total genes
# gather number of genes and data
# my $numGenes = $dbh->selectrow_array(qq{SELECT COUNT (DISTINCT locusId) FROM Gene;}, undef);
# my $numData = $dbh->selectrow_array(qq{SELECT COUNT (DISTINCT locusId) FROM GeneFitness;}, undef);
# print $cgi->p("Fitness data for $numData genes of $numGenes genes.");


#display taxonomic information and link
# if ((defined $orginfo->{$orgId}{taxonomyId}) && ($orginfo->{$orgId}{taxonomyId} ne "")) {
	# print $cgi->p($cgi->a({href => "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=$orginfo->{$orgId}{taxonomyId}"}, "NCBI Taxonomy"));
# }

#file of experimental data - generated by createExpData.cgi
# my $data = `./createExpData.cgi orgId=$orgId`;
# print $cgi->p($cgi->a({href => "download.cgi?orgId=$orgId"}, "Download experimental data"), " - Note: May take a minute or so to load once clicked.");
# print $cgi->p($cgi->a({href => "createExpData.cgi?orgId=$orgId"}, "Download experimental data"), " - Note: May take a few seconds to load once clicked.");


$dbh->disconnect();
Utils::endHtml($cgi);
