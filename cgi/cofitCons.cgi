#!/usr/bin/perl -w
#######################################################
## cofit.cgi -- examine cofitness of orthologs of two genes
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required parameters: orgId, locusId, hitId

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

sub cofitStrings($$$$); # dbh, orgId, locus1, locus2 => cofitness, rank
my $cgi=CGI->new;
my $style = Utils::get_style();
print $cgi->header;

my $orgId = $cgi->param('orgId') || die "no orgId";
my $locusId = $cgi->param('locusId') || die "no locusId";
my $hitId = $cgi->param('hitId') || die "no hitId";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Invalid orgId" unless exists $orginfo->{$orgId};

my $gene1 = $dbh->selectrow_hashref("SELECT * from Gene WHERE orgId = ? AND locusId = ?;",
				    {}, $orgId, $locusId);
die "Invalid locusId" unless defined $gene1->{locusId};
my $gene2 = $dbh->selectrow_hashref("SELECT * from Gene WHERE orgId = ? AND locusId = ?;",
				    {}, $orgId, $hitId);
die "Invalid hitId" unless defined $gene2->{locusId};

my $show1 = $gene1->{sysName} || $gene1->{locusId};
my $show2 = $gene2->{sysName} || $gene2->{locusId};

my $orth1 = Utils::get_orths($dbh,$orgId,$locusId);
my $orth2 = Utils::get_orths($dbh,$orgId,$hitId);
my @orth1 = sort { $b->{ratio} <=> $a->{ratio} } values(%$orth1);
my @orth1b = grep { exists $orth2->{$_->{orgId}} } @orth1;

my @trows = ( Tr({-valign => 'top', -align => 'center'},
		 th(["organism",
		     a({title=>"score ratio for ortholog 1: blast_score_of_alignment / self_score"},'ratio1'),
		     "gene1", "name1", "description1",
		     a({title=>"score ratio for ortholog 2: blast_score_of_alignment / self_score"},'ratio2'),
		     "gene2", "name2", "description2", "cofit", "rank" ])),
	     Tr({-valign => 'top', -align => 'left'},
		 td([ $cgi->a({href => "org.cgi?orgId=". $orginfo->{$orgId}{orgId}}, $orginfo->{$orgId}{genome}),
		      "1.0",
		      a({ href => "myFitShow.cgi?orgId=$orgId&gene=$locusId"}, $show1),
		      $gene1->{gene},
		      $gene1->{desc},
		      "1.0",
		      a({ href => "myFitShow.cgi?orgId=$orgId&gene=$hitId"}, $show2),
		      $gene2->{gene},
		      $gene2->{desc},
		      &cofitStrings($dbh,$orgId,$locusId,$hitId) ])) );

my $nBoth = scalar(@orth1b);
my $n1only = scalar(@orth1) - $nBoth;
my $n2only = scalar(keys %$orth2) - $nBoth;

foreach my $o1 (@orth1b) {
    my $orthOrg = $o1->{orgId};
    my $o2 = $orth2->{$orthOrg};
    push @trows, Tr({-valign => 'top', -align => 'left'},
		    td([ $cgi->a({href => "org.cgi?orgId=". $orginfo->{$orthOrg}{orgId}}, $orginfo->{$orthOrg}{genome}),
			 sprintf("%.2f", $o1->{ratio}),
			 a({ href => "myFitShow.cgi?orgId=$orthOrg&gene=$o1->{locusId}"}, $o1->{sysName} || $o1->{locusId}),
			 $o1->{gene},
			 $o1->{desc},
			 sprintf("%.2f", $o2->{ratio}),
			 a({ href => "myFitShow.cgi?orgId=$orthOrg&gene=$o2->{locusId}"}, $o2->{sysName} || $o2->{locusId}),
			 $o2->{gene},
			 $o2->{desc},
			 &cofitStrings($dbh, $orthOrg, $o1->{locusId}, $o2->{locusId}) ]));
}

my $title = "Conservation of cofitness between $show1 and $show2 in $orginfo->{$orgId}{genome}";
print
    start_html(
	-title => $title,
	-style => {-code => $style},
	-author=>'morgannprice@yahoo.com',
	-meta=>{'copyright'=>'copyright 2015 UC Berkeley'} ),
    h2($title),
    div({-style => "float: right; vertical-align: top;"},
	a({href => "help.cgi#ortholog"}, "Help")),
    h3("$nBoth genomes with putative orthologs of both genes"),
    table({-cellspacing => 0, cellpadding => 3}, @trows);
print p("Not shown: $n1only genomes with orthologs for $show1 only; $n2only genomes with orthologs for $show2 only");

# then print table of the two genes
# then print graceful error message if either lacks orthologs
# then print graceful error message if no organisms with orthologs for both (@orth1b empty)
# then for each in 1 with orth1b, look for cofitness...

Utils::endHtml($cgi);

sub cofitStrings($$$$) {
    my ($dbh,$orgId,$locus1,$locus2) = @_;
    my ($cofit,$rank) = $dbh->selectrow_array("SELECT cofit,rank FROM Cofit where orgId=? AND locusId=? AND hitId=?",
					      {}, $orgId, $locus1, $locus2);
    my $URL = "genesFit.cgi?orgId=$orgId&locusId=$locus1&locusId=$locus2";
    return (a({href => $URL}, sprintf("%.2f", $cofit)), $rank) if defined $cofit;
    # else !defined cofit
    # get min cofit and max rank
    my $hasData1 = Utils::gene_has_fitness($dbh,$orgId,$locus1);
    my $hasData2 = Utils::gene_has_fitness($dbh,$orgId,$locus2);
    if ($hasData1 && $hasData2 ) {
	($cofit,$rank) = $dbh->selectrow_array("SELECT min(cofit), max(rank) FROM Cofit where orgId=? AND locusId=? GROUP BY orgId,locusId;",
					       {}, $orgId, $locus1);
	return (a({href => $URL, title => sprintf("under %.2f",$cofit)}, "low"), "&gt; $rank") if defined $cofit;
	# else
	return (a({title => "no cofitness for this organism"}, "&mdash;"), a({title => "no cofitness for this organism"}, "&mdash;"));
    }
    # else
    return (a({title => "no data for either"}, "&mdash;"), a({title => "no data for either"}, "&mdash;"))
	if (!$hasData1 && !$hasData2);
    return (a({title => "no data for gene 1"}, "&mdash;"), a({title => "no data for gene 1"}, "&mdash;"))
	if !$hasData1 && $hasData2;
    return (a({title => "no data for gene 2"}, "&mdash;"), a({title => "no data for gene 2"}, "&mdash;"))
	if $hasData1 && !$hasData2;
    die "Unreachable";
}

