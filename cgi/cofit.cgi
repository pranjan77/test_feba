#!/usr/bin/perl -w
#######################################################
## cofit.cgi
##
## Copyright (c) 2015 All Right Reserved by UC Berkeley
##
## Authors:
## Morgan Price
#######################################################

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();
print $cgi->header;

my $species = $cgi->param('species') || die "no species identifier\n";
my $locusId = $cgi->param('gene') || die "no gene identifier\n";

Utils::fail($cgi, "$species is invalid. Please enter correct species name!") unless ($species =~ m/^[A-Za-z0-9_]*$/);
Utils::fail($cgi, "$locusId is invalid. Please enter correct locusId!") unless ($locusId =~ m/^[A-Za-z0-9_]*$/);

my $dbh = Utils::get_dbh();
my ($sysName,$desc) = $dbh->selectrow_array("SELECT sysName,desc FROM Gene WHERE nickname=? AND locusId=?", undef,
					    $species, $locusId);
Utils::fail($cgi, "No locus $locusId in species $species") unless defined $sysName;
my ($speciesName) = $dbh->selectrow_array("SELECT species FROM Organism WHERE name=?", undef, $species);
Utils::fail($cgi, "No species information for $species") unless defined $speciesName;

print $cgi->start_html(
    -title =>"Cofitness for $sysName ($speciesName)",
    -style => {-code => $style},
    -author=>'morgannprice@yahoo.com',
    -meta=>{'copyright'=>'copyright 2015 UC Berkeley'},
);
print $cgi->h2("Cofitness");
print $cgi->h4("Top cofit genes for $sysName: $desc from $speciesName");
#print $cgi->p(qq{<A HREF="myFitShow.cgi?species=$species&gene=$locusId">$sysName: $desc</A>});

my $cofitResults = $dbh->selectall_arrayref(qq{
	SELECT hitId, rank, cofit, sysName AS hitSysName, desc AS hitDesc
		FROM Cofit JOIN Gene ON Cofit.hitId=Gene.locusId AND Cofit.nickname=Gene.nickname
		WHERE Cofit.nickname=? AND Cofit.locusId=?
		ORDER BY rank LIMIT 20}, undef, $species, $locusId) || die;
if (@$cofitResults == 0) {
    print $cgi->p(qq{Cofitness results are not available for this organism, sorry.});
} else {
    my @headRow = map { $cgi->td($cgi->b($_)) } qw{Rank Cofitness Hit Description};

    my @trows = ( $cgi->Tr(@headRow) );
    my @colors = ('#FFFFDD', '#FFFFFF');
    my $iRow = 0;

    foreach  my $row (@$cofitResults) {
	my ($hitId,$rank,$cofit,$hitSysName,$hitDesc) = @$row;
	my $rowcolor = $colors[ $iRow++ % scalar(@colors) ];
	$cofit = sprintf("%.2f",$cofit);
	push @trows, $cgi->Tr({bgcolor => $rowcolor, align => 'left', valign => 'top' },
			      $cgi->td($rank),
			      $cgi->td($cofit),
			      $cgi->td($cgi->a( {href => "myFitShow.cgi?species=$species&gene=$hitSysName" },
						$hitSysName )),
			      $cgi->td($hitDesc));
    }
    print $cgi->table( {cellpadding => 3, cellspacing => 0 }, @trows );
}
print $cgi->h4($cgi->a({href => "myFitShow.cgi?species=$species&gene=$locusId"}, "Fitness data"));
print $cgi->h4($cgi->a({href => "mySeqSearch.cgi?species=$species&gene=$locusId"}, "Check homologs"));

Utils::endHtml($cgi);
