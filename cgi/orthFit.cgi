#!/usr/bin/perl -w
#######################################################
## orthFit.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Given Group, Condition_1, and an anchor gene, list fitness values
# of the gene and its orthologs for relevant conditions.
# This is useful even if there are no orthologs because it
# shows other replicates or other concentrations.
#
# Required CGI parameters:
# orgId, locusId -- these specify the anchor gene
# expGroup and condition1 -- these specify the condition.
#	Note that condition1 may be empty but it must still be present
#	(and is used to choose which to show).

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();

my $orgId = $cgi->param('orgId') || die "no orgId parameter";
my $locusId = $cgi->param('locusId') || die "no locusId parameter";
my $expGroup = $cgi->param('expGroup') || die "no expGroup parameter";
my $condition1 = $cgi->param('condition1');
die "no condition1 parameter" unless defined $condition1;

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Invalid orgId" unless exists $orginfo->{$orgId};
my $gene = $dbh->selectrow_hashref("SELECT * from Gene where orgId=? AND locusId=?",
				   {}, $orgId, $locusId);
die "No such gene" unless defined $gene->{locusId};

my $orths = Utils::get_orths($dbh, $orgId, $locusId);
my @orths = sort { $b->{ratio} <=> $a->{ratio} } values(%$orths);
my @genes = ($gene);
push @genes, @orths;
# ensure that condition1 is valid
my ($cnt) = $dbh->selectrow_array(qq{SELECT count(*) from Experiment WHERE orgId = ? AND expGroup = ? AND condition_1 = ?},
				  {}, $orgId, $expGroup, $condition1);
die "No experiments for specified expGroup and condition1\n" unless $cnt > 0;

my $showId = $gene->{sysName} || $gene->{locusId};
my $title = "$showId and its orthologs: fitness in $expGroup $condition1";
print header,
    start_html(
	-title =>$title,
	-style => {-code => $style},
	-author=>'Morgan Price',
	-meta=>{'copyright'=>'copyright 2015 UC Berkeley'}),
    h2($title),
    div({-style => "float: right; vertical-align: top;"},
	a({href => "help.cgi#ortholog"}, "Help")),
    div({-style => "clear: right"});

my @headings = ( a({title=>"score ratio for ortholog: blast_score_of_alignment / self_score"},'ratio'),
		 qw{organism gene name description experiment fitness t} );
my @trows = ( $cgi->Tr({ -valign=> 'top', -align => 'center'},
		       map { th($_) } @headings) );
my $nSkipOrth = 0;
foreach my $o (@genes) {
    my $data = $dbh->selectall_arrayref(qq{SELECT * from Experiment JOIN GeneFitness USING  (orgId,expName)
                                           WHERE orgId=? AND locusId=? AND expGroup=? AND condition_1=? ORDER BY fit},
					{ Slice => {} }, $o->{orgId}, $o->{locusId}, $expGroup, $condition1);
    $nSkipOrth++ if @$data == 0 && $o->{orgId} ne $orgId;
    my $first = 1;
    foreach my $row (@$data) {
	my $ratio = $o->{orgId} eq $orgId ? "&mdash;" : sprintf("%.2f",$o->{ratio});
	my $orgShort = "";
	if ($first) {
	    my $d = $orginfo->{$o->{orgId}};
	    my $short = $d->{genome}; $short =~ s/^(.)[^ ]+/\1./;
	    $orgShort = a({title => $d->{genome}}, $short);
	}
	push @trows, $cgi->Tr({ -valign => 'top', -align => 'left'},
			      td($first ?  $ratio : ""),
			      td($orgShort),
			      td($first ? a({href => "myFitShow.cgi?orgId=$o->{orgId}&gene=$o->{locusId}", style => "color: rgb(0,0,0)"},
					    $o->{sysName} || $o->{locusId}) : ""),
			      td($first ? $o->{gene} : ""),
			      td($first ? $o->{desc} : ""),
			      td(a({href => "exp.cgi?orgId=$o->{orgId}&expName=$row->{expName}", style => "color: rgb(0,0,0)"},
				   $row->{expDesc})),
			      td({ -bgcolor => Utils::fitcolor($row->{fit}) }, sprintf("%.1f", $row->{fit})),
			      td(sprintf("%.1f", $row->{t})) );
	$first = 0;
    }
}
print
    table({ cellpadding => 3, cellspacing =>0}, @trows),
    $nSkipOrth ? p("$nSkipOrth orthologs are not shown because they lack fitness data for this condition (or they lack data entirely)") : "",
    p(a({href => "mySeqSearch.cgi?orgId=$orgId&locusId=$locusId"}, "Show all homologs"));


$dbh->disconnect();
Utils::endHtml($cgi);
