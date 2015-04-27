#!/usr/bin/perl -w
#######################################################
## cofit.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required parameters: orgId and locusId, for organism and which gene

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();
print $cgi->header;

my $orgId = $cgi->param('orgId') || die "no species identifier\n";
my $locusId = $cgi->param('locusId') || die "no gene identifier\n";

Utils::fail($cgi, "$orgId is invalid. Please enter correct species name!") unless ($orgId =~ m/^[A-Za-z0-9_]*$/);
Utils::fail($cgi, "$locusId is invalid. Please enter correct locusId!") unless ($locusId =~ m/^[A-Za-z0-9_]*$/);

my $dbh = Utils::get_dbh();
my ($sysName,$geneName,$desc) = $dbh->selectrow_array("SELECT sysName,gene,desc FROM Gene WHERE orgId=? AND locusId=?", undef,
					    $orgId, $locusId);
Utils::fail($cgi, "No locus $locusId in species $orgId") unless defined $sysName;
my ($genus,$species,$strain) = $dbh->selectrow_array("SELECT genus,species,strain FROM Organism WHERE orgId=?", undef, $orgId);
Utils::fail($cgi, "No species information for $orgId") unless defined $species;

print $cgi->start_html(
    -title =>"Cofitness for $sysName ($genus $species $strain)",
    -style => {-code => $style},
    -author=>'morgannprice@yahoo.com',
    -meta=>{'copyright'=>'copyright 2015 UC Berkeley'},
);
print $cgi->h2("Cofitness");
print $cgi->h4("Top cofit genes for $geneName $locusId $sysName: $desc from $genus $species $strain");
#print $cgi->p(qq{<A HREF="myFitShow.cgi?orgId=$orgId&gene=$locusId">$sysName: $desc</A>});

my $cofitResults = $dbh->selectall_arrayref(qq{
	SELECT hitId, rank, cofit, gene AS hitName, sysName AS hitSysName, desc AS hitDesc
		FROM Cofit JOIN Gene ON Cofit.hitId=Gene.locusId AND Cofit.orgId=Gene.orgId
		WHERE Cofit.orgId=? AND Cofit.locusId=?
		ORDER BY rank LIMIT 20}, undef, $orgId, $locusId) || die;
if (@$cofitResults == 0) {
    print $cgi->p(qq{Cofitness results are not available for this organism, sorry.});
} else {
    my @headRow = map { $cgi->td($cgi->b($_)) } qw{Rank Cofitness Hit Name Description};

    my @trows = ( $cgi->Tr(@headRow) );
    my @colors = ('#FFFFDD', '#FFFFFF');
    my $iRow = 0;

    foreach  my $row (@$cofitResults) {
	my ($hitId,$rank,$cofit,$hitName,$hitSysName,$hitDesc) = @$row;
	my $showId = $hitSysName || $hitId;
	my $rowcolor = $colors[ $iRow++ % scalar(@colors) ];
	$cofit = sprintf("%.2f",$cofit);
	push @trows, $cgi->Tr({bgcolor => $rowcolor, align => 'left', valign => 'top' },
			      $cgi->td($rank),
			      $cgi->td($cofit),
			      $cgi->td($cgi->a( {href => "myFitShow.cgi?orgId=$orgId&gene=$hitId" },
						$showId )),
			      $cgi->td($hitName),
			      $cgi->td($hitDesc));
    }
    print $cgi->table( {cellpadding => 3, cellspacing => 0 }, @trows );
}
print $cgi->h4($cgi->a({href => "myFitShow.cgi?orgId=$orgId&gene=$locusId"}, "Fitness data"));
print $cgi->h4($cgi->a({href => "mySeqSearch.cgi?orgId=$orgId&locusId=$locusId"}, "Check homologs"));

Utils::endHtml($cgi);
