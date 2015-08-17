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
# Optional CGI parameters: 
# help -- 1 if on help/tutorial mode

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;
# my $style = Utils::get_style();
print $cgi->header;

my $orgId = $cgi->param('orgId') || die "no species identifier\n";
my $locusId = $cgi->param('locusId') || die "no gene identifier\n";
my $help = $cgi->param('help') || "";

Utils::fail($cgi, "$orgId is invalid. Please enter correct species name!") unless ($orgId =~ m/^[A-Za-z0-9_]*$/);
Utils::fail($cgi, "$locusId is invalid. Please enter correct locusId!") unless ($locusId =~ m/^[A-Za-z0-9_]*$/);

my $dbh = Utils::get_dbh();
my ($sysName,$geneName,$desc,$type) = $dbh->selectrow_array("SELECT sysName,gene,desc,type FROM Gene WHERE orgId=? AND locusId=?", undef,
					    $orgId, $locusId);
Utils::fail($cgi, "No locus $locusId in species $orgId") unless defined $sysName;
my ($genus,$species,$strain) = $dbh->selectrow_array("SELECT genus,species,strain FROM Organism WHERE orgId=?", undef, $orgId);
Utils::fail($cgi, "No species information for $orgId") unless defined $species;

my $showId = $sysName || $locusId;
my $showName = $geneName || $showId;
my $start = Utils::start_page("Cofitness for $showName ($genus $species $strain");
my $tabs = Utils::tabsGene($dbh,$cgi,$orgId,$locusId,0,$type,"cofit");

print
	$start, $tabs,
#     start_html( -title =>"Cofitness for $sysName ($genus $species $strain)",
# 		-style => {-code => $style},
# 		-author=>'morgannprice@yahoo.com',
# 		-meta=>{'copyright'=>'copyright 2015 UC Berkeley'}),
    h2("Top cofit genes for $showId from " . $cgi->a({href => "org.cgi?orgId=$orgId"}, "$genus $species $strain")),
 #    div({-style => "float: right; vertical-align: top;"},
	# a({href => "help.cgi#cofitness"}, "Help")),

	# qq[<div style="position: relative;"><div class="floatbox">],
	# start_form(-name => 'input', -method => 'GET', -action => 'genesFit.cgi'),
	# "<BR>Compare selected genes to $showId $geneName: " . submit(-name=>"heatmap"),
	# end_form,
	# qq[</div>];

    h3({class=>"short"},"$showId $geneName : $desc");
#print $cgi->p(qq{<A HREF="myFitShow.cgi?orgId=$orgId&gene=$locusId">$sysName: $desc</A>});

if ($help == 1) {
		print qq[<BR><BR><div class="helpbox">
		<b><u>About this page:</u></b><BR><ul>
		<li>Show the genes with the most similar fitness patterns as this gene. </li>
		<li>To get to this page, search for any gene and click on the "Cofit" tab.</li> 
		<li>To view a scatterplot of the fitness patterns, click on a cofitness value.</li>
		<li>To make a heatmap, checkmark genes of interest and click the "Heatmap" button at the bottom.</li>
                <li>For more about cofitness, see the <A HREF="help.cgi#cofitness">help page</A>.
		</ul></div>];
	}


my $cofitResults = $dbh->selectall_arrayref(qq{
	SELECT hitId, rank, cofit, gene AS hitName, sysName AS hitSysName, desc AS hitDesc
		FROM Cofit JOIN Gene ON Cofit.hitId=Gene.locusId AND Cofit.orgId=Gene.orgId
		WHERE Cofit.orgId=? AND Cofit.locusId=?
		ORDER BY rank LIMIT 20}, undef, $orgId, $locusId) || die;
if (@$cofitResults == 0) {
    print "<br>", $cgi->p(qq{Cofitness results are not available for this gene, sorry.});
} else {
    my @headRow = map { $cgi->th($cgi->b($_)) } qw{Rank Hit Name Description}, a({title => "Maximum cofitness of orthologs"}, "Conserved?"), qw{Cofitness &nbsp;} ;

    my @trows = ( $cgi->Tr(@headRow) );
    my @colors = ('#FFFFDD', '#FFFFFF');
    my $iRow = 0;

    foreach  my $row (@$cofitResults) {
	my ($hitId,$rank,$cofit,$hitName,$hitSysName,$hitDesc) = @$row;
	my $showId = $hitSysName || $hitId;
	my $rowcolor = $colors[ $iRow++ % scalar(@colors) ];
	$cofit = sprintf("%.2f",$cofit);
	my ( $cofitCons ) = $dbh->selectrow_array(qq{ SELECT max(cofit) FROM Ortholog o1, Ortholog o2, Cofit c
                                 	WHERE o1.orgId1 = ? AND o1.locusId1 = ?
	                                 AND o1.orgId2 = c.orgId AND o1.locusId2 = c.locusId
	                                 AND o2.orgId1 = o1.orgId1 AND o2.orgId2=o1.orgId2
	                                 AND o2.locusId1 = ?
                                         AND o2.locusId2 = c.hitId; },
				    {}, $orgId, $locusId, $hitId);
	my $url2 = "compareGenes.cgi?orgId=$orgId&locus1=$locusId&locus2=$hitId";
	$url2 .= "&help=1" if $help == 1;
	push @trows, $cgi->Tr({bgcolor => $rowcolor, align => 'left', valign => 'top' },
			      $cgi->td($rank),
			      $cgi->td($cgi->a( {href => "myFitShow.cgi?orgId=$orgId&gene=$hitId" },
						$showId )),
			      $cgi->td($hitName),
			      $cgi->td($hitDesc),
                  $cgi->td( $cgi->a({href => "cofitCons.cgi?orgId=$orgId&locusId=$locusId&hitId=$hitId"},
						defined $cofitCons ? sprintf("%.2f", $cofitCons) : "no") ),
                  $cgi->td($cgi->a({title=>"Compare genes via scatterplot", href => "$url2"}, $cofit)),
                  $cgi->td(checkbox('locusId',0,$hitId,'')),
	                      );

    }
    my $name = $geneName || $showId;
    print
	start_form(-name => 'input', -method => 'GET', -action => 'genesFit.cgi'),
	hidden('orgId', $orgId),
	hidden('locusId', $locusId),
	'<br>',
	p("Select other genes you would like to make a comparative heatmap with $name with and click the button below."),
	# qq[<br><div style="position: relative;"><div class="floatbox floatbox2">],
	# submit(-class=>"heatmap", -name=>"heatmap with $name"),
	# qq[</div></div>],
	table( {cellpadding => 3, cellspacing => 0 }, @trows ),
	p(submit(-class=>"heatmap", -name=>"Heatmap of $name with selected genes")),
	# "<BR>Compare selected genes to $showId $geneName: " . 
	end_form,
	"<BR><BR>";
}
# print $cgi->p($cgi->a({href => "myFitShow.cgi?orgId=$orgId&gene=$locusId"}, "Fitness data"));

 #    print
	# p(a({href => "getSeq.cgi?orgId=$orgId&locusId=$locusId"}, "Show sequence"),
	#   "or",
	#   a({href => "mySeqSearch.cgi?orgId=$orgId&locusId=$locusId"}, "Check homologs"))
	# if $type == 1;

Utils::endHtml($cgi);
