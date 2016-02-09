#!/usr/bin/perl -w
#######################################################
## exp.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price, Victoria Lo
#######################################################
#
# Required parameters: orgId and expName, for organism and which experiment
# Optional parameter show:
# specific -- show specific phenotypes
# important or detrimental -- show top 200 genes either way
# quality -- show quality metrics
# (by default, shows all of these optoins)
# help -- 1 if on help/tutorial mode

use strict;
use CGI qw(:standard Vars -nosticky);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
sub CompoundToHTML($);

use lib "../lib";
use Utils;
use URI::Escape;

my $cgi=CGI->new;
print $cgi->header;

my $orgId = $cgi->param('orgId') || die "no species identifier\n";
die "Invalid species name!" unless $orgId =~ m/^[A-Za-z0-9_]*$/;
my $expName = $cgi->param('expName') || die "no experiment identifier\n";
die "Invalid experiment name!" unless $expName =~ m/^[A-Za-z0-9_]*$/;
my $show = $cgi->param('show') || "";
my $help = $cgi->param('help') || "";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "No such orgId: $orgId" unless exists $orginfo->{$orgId};

my $exp = $dbh->selectrow_hashref("SELECT * FROM Experiment WHERE orgId = ? AND expName = ?", {}, $orgId, $expName);
die "No such experiment: $expName" unless defined $exp->{expName};

# Fetch specific phenotypes
my $spec = $dbh->selectall_arrayref(qq{SELECT * from SpecificPhenotype JOIN GeneFitness USING (orgId,locusId,expName)
                                       WHERE orgId = ? AND expName = ?},
				    { Slice => {} },
				    $orgId, $expName);
my $start = Utils::start_page("Experiment $expName for $orginfo->{$orgId}{genome}");
my $tabs;
if ($show eq "") {
	$tabs = Utils::tabsExp($dbh,$cgi,$orgId,$expName,$exp->{expGroup},$exp->{condition_1},"overview");
} elsif ($show eq "specific") {
	$tabs = Utils::tabsExp($dbh,$cgi,$orgId,$expName,$exp->{expGroup},$exp->{condition_1},"specific");
} elsif ($show eq "important"){
	$tabs = Utils::tabsExp($dbh,$cgi,$orgId,$expName,$exp->{expGroup},$exp->{condition_1},"-gene");
} elsif ($show eq "detrimental") {
	$tabs = Utils::tabsExp($dbh,$cgi,$orgId,$expName,$exp->{expGroup},$exp->{condition_1},"+gene");
} elsif ($show eq "quality") {
	$tabs = Utils::tabsExp($dbh,$cgi,$orgId,$expName,$exp->{expGroup},$exp->{condition_1},"metrics");
}

print $start, $tabs,
    h2("Experiment $expName for ". $cgi->a({href => "org.cgi?orgId=$orgId"}, "$orginfo->{$orgId}{genome}")),
	# corner compare box
	qq[<div style="position: relative;"><div class="floatbox">],
    start_form(-name => 'input', -method => 'GET', -action => 'compareExps.cgi'),
    hidden('orgId', $orgId),
    hidden('expName2', $expName),
    "Compare to: ",
    textfield(-name => 'query1', -value => '', -size => 20, -maxlength => 100),
    # submit('Go'),
    "<button type='submit'>Go</button>",
    end_form,
    qq[</P></div></div>],
    h3($exp->{expDescLong});

# print qq[<div style="position: relative;"><div class="floatbox">],
#     start_form(-name => 'input', -method => 'GET', -action => 'compareExps.cgi'),
#     hidden('orgId', $orgId),
#     hidden('expName2', $expName),
#     "Compare to another experiment: ",
#     textfield(-name => 'query1', -value => '', -size => 20, -maxlength => 100),
#     end_form,
#     qq[</P></div></div>],

if ($help == 1) {
        print qq[<div class="helpbox">
        <b><u>About this page:</u></b><BR><ul>
        <li>View the genes that had strong and specific phenotypes in this experiment. </li>
        <li>To get to this page, search for any experiment and click on the "Specific" tab.</li> 
        <li>To compare to another experiment via scatterplot, add another experiment using the box above. (Try "cisplatin".)</li>
        <li>To make a comparative heatmap, check the genes of interest and click the "Heatmap" link at the bottom.</li>
        <li>For more about how we define a specific phenotype, see the <A HREF="help.cgi#specific">help page</A>.
        </ul></div>];
    }


my @fit = (); # sorted list of fitness values to show
my $header = undef;
if ($show eq "") {
    my ($html1,$html2);
    $html1 = CompoundToHTML( $exp->{condition_1} );
    $html2 = CompoundToHTML( $exp->{condition_2} );
    my $cond1 = $exp->{condition_1} ? join(" ", $html1,
                                           $exp->{concentration_1}, $exp->{units_1}) : "";
    my $cond2 = $exp->{condition_2} ? join(" ", $html2,
                                           $exp->{concentration_2}, $exp->{units_2}) : "";
    my $media = $exp->{media};
    if ($cond2) {
        $media = join(" + ", $media, $cond1, $cond2);
    } elsif ($cond1) {
        $media = join(" + ", $media, $cond1);
    }
    $media .= " pH=$exp->{pH}" if $exp->{pH} ne "";
    my @culture = ("Culturing: ". $exp->{mutantLibrary});
    push @culture, $exp->{vessel} if $exp->{vessel} ne "";
    push @culture, $exp->{aerobic} if $exp->{aerobic} ne "";
    push @culture, "at $exp->{temperature} (C)" if $exp->{temperature} ne "";
    push @culture, "shaken=$exp->{shaking}" if $exp->{shaking} ne "";
    push @culture, "($exp->{liquid})" if $exp->{liquid} ne "" && lc($exp->{liquid}) ne "liquid";
    
    my @pieces = ("Media: $media", join(", ",@culture), "By: $exp->{person} on $exp->{dateStarted}");
    
    my $mediaComponents = $dbh->selectall_arrayref(qq{SELECT * from MediaComponents LEFT JOIN Compounds USING (compound)
                                                          WHERE media = ?},
                                                   { Slice => {} },
                                                   $exp->{media});
    if (@$mediaComponents > 0) {
        my @compStrings = ();
        foreach my $row (@$mediaComponents) {
            my $compString = $row->{CAS} ?
                a({-href => "http://commonchemistry.org/ChemicalDetail.aspx?ref=$row->{CAS}"}, $row->{compound})
                :  $row->{compound};
            $compString = "$row->{concentration} $row->{units} $compString" if $row->{concentration} && $row->{units};
            push @compStrings, $compString;
        }
        push @pieces, "Media components: " . small(join(", ", @compStrings));
    }
    if ($exp->{growthPlate} ne "" && $exp->{growthWells} ne "") {
        push @pieces, "Growth plate: $exp->{growthPlate} $exp->{growthWells}";
    }
    print join("<BR>\n", @pieces)."\n";
} elsif ($show eq "specific") {
    $header = "Genes with " . a({href => "help.cgi#specific"}, "specific") . " phenotypes:";
    @fit = @{ $dbh->selectall_arrayref(qq{SELECT * FROM SpecificPhenotype JOIN GeneFitness USING (orgId,expName,locusId)
                                          JOIN Gene USING (orgId,locusId)
					  WHERE orgId=? AND expName=? ORDER BY fit},
				       { Slice => {} },
				       $orgId, $expName) };
} elsif ($show eq "important") {
    $header = "200 most important genes:";
    @fit = @{ $dbh->selectall_arrayref(qq{SELECT * from GeneFitness JOIN GENE USING (orgId,locusId)
                                          WHERE orgId=? AND expName=?
                                          ORDER BY fit LIMIT 200},
				       { Slice => {} },
				       $orgId, $expName) };
} elsif ($show eq "detrimental") {
    $header = "200 most detrimental genes:";
    @fit = @{ $dbh->selectall_arrayref(qq{SELECT * from GeneFitness JOIN GENE USING (orgId,locusId)
                                          WHERE orgId=? AND expName=?
                                          ORDER BY fit DESC LIMIT 200},
				       { Slice => {} },
				       $orgId, $expName) };
} elsif ($show eq "quality") {
    $header = "Quality Metrics:";
}

print $cgi->h3($header) if defined $header;

if (@fit > 0) { # show the table
    my @trows = ();
    foreach my $row (@fit) {
	my $geneId = $row->{sysName} || $row->{locusId};
    my $strainUrl = "strainTable.cgi?orgId=$orgId&locusId=$row->{locusId}&expName=$expName";
    $strainUrl .= "&help=1" if $help == 1;
    my $orthUrl = "orthFit.cgi?orgId=$orgId&locusId=$row->{locusId}"
                          . "&expGroup=" . uri_escape($exp->{expGroup})
                                              . "&condition1=" . uri_escape($exp->{condition_1});
    $orthUrl .= "&help=1" if $help == 1;
	push @trows, $cgi->Tr({align => 'left', valign => 'top'},
			      td(checkbox('locusId',0,$row->{locusId},'')),
			      td(a({href => "myFitShow.cgi?orgId=$orgId&gene=$row->{locusId}",},
					       # style => "color:rgb(0,0,0)"},
						$geneId)),
			      td($row->{gene}),
			      td({ -bgcolor => Utils::fitcolor($row->{fit}) },
                                 a({ -href => $strainUrl,
                                     -title => "per-strain data",
                                     -style => "color:rgb(0,0,0)" },
                                     sprintf("%.1f", $row->{fit})) ),
			      td( sprintf("%.1f", $row->{t}) ),
			      td($row->{desc}),
			      td(a({ 
				     title => "Compare to data from similar experiments or orthologs",
				     href => $orthUrl},
				 "compare")) );
    }
    print
	start_form(-name => 'input', -method => 'GET', -action => 'genesFit.cgi'),
	hidden('orgId', $orgId),
	$cgi->table(
	    { cellspacing => 0, cellpadding => 3},
	    $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
		     $cgi->th(['&nbsp;', 'gene','name','fitness','t score','description', '&nbsp;']),
		     @trows) ),
	"<BR><BR>",
	submit(-name => 'Heatmap of selected genes'),
	end_form;
} elsif ($show eq "quality") {
    my @trows = ();
    push @trows, $cgi->Tr({align => 'left', valign => 'top' },
			  $cgi->td([ "Time0", $exp->{timeZeroSet}, "which Time0s the sample was compared to"]));
    my %expl = ("cor12" => "rank correlation(fit1, fit2), where fit1 is fitness for the first half (10-50%) and fit2 is fitness for the second half (50-90%) of each gene",
		"maxFit" => "The maximum fitness value",
		"opcor" => "rank correlation(upstream gene, downstream gene) over pairs that are adjacent and likely to be in the same operon",
		"adjcor" => "like opcor but for adjacent genes that are not on the same strand",
		"gccor" => "linear correlation of gene fitness and gene GC content",
		"mad12" => "median absolute difference of fit1, fit2",
		"mad12c" => "median absolute difference of log count for 1st and 2nd half of genes in this sample",
		"mad12c_t0" => "like mad12c but for the Time0s",
		"gMed" => "median reads per gene in this sample",
		"gMedt0" => "median reads per gene in the Time0 sample",
		"gMean" => "mean reads per gene in this sample",
		"nMapped" => "#reads for this sample that corresponded to a known strain (in millions)",
		"nPastEnd" => "#reads that corresponded to a strain that has an insertion within the suicide vector instead of within the genome.",
		"nGenic" => "#reads that lie within central 10-90% of a gene",
		"nUsed" => "#reads used to estimate gene fitness (genic and enough coverage for strain and for gene)" );
    foreach my $key (qw{cor12 maxFit opcor adjcor gccor mad12 mad12c mad12c_t0}) {
	my $val = $exp->{$key};
	push @trows, $cgi->Tr({align => 'left', valign => 'top' },
			      $cgi->td([$key, sprintf("%.2f", $exp->{$key}), $expl{$key}]));
    }
    foreach my $key (qw{gMed gMedt0 gMean}) {
	push @trows, $cgi->Tr({align => 'left', valign => 'top' },
			      $cgi->td([$key, sprintf("%.0f", $exp->{$key}), $expl{$key}]));
    }
    foreach my $key (qw{nMapped nPastEnd nGenic nUsed}) {
	push @trows, $cgi->Tr({align => 'left', valign => 'top' },
			      $cgi->td([$key, sprintf("%.3f M", $exp->{$key}/1e6), $expl{$key}]));
    }
    print $cgi->table({cellpadding => 3, cellspacing => 0}, @trows);
}

if ($show ne "specific") {
    if (@$spec > 0) {
	print p(a({href => "exp.cgi?orgId=$orgId&expName=$expName&show=specific"},
			      "Specific phenotypes for " . scalar(@$spec). " genes in this experiment"));
        print p(a({href => "spec.cgi?orgId=$orgId&expGroup=".uri_escape($exp->{expGroup})."#".$exp->{condition_1} },
                  "Specific phenotypes for $exp->{expGroup} $exp->{condition_1} in $orginfo->{$orgId}{genome}"))
            if $exp->{expGroup} && $exp->{condition_1};
    }  else {
	print $cgi->p("No genes had specific phenotypes in this experiment.");
    }
}


print
    p(a({href => "orthCond.cgi?expGroup=" . uri_escape($exp->{expGroup})
         . "&condition1=" . uri_escape($exp->{condition_1})},
	"Specific phenotypes for $exp->{expGroup} $exp->{condition_1} across organisms"));
    
$dbh->disconnect();
Utils::endHtml($cgi);

sub CompoundToHTML($) {
    my ($compound) = @_;
    return $compound if ! $compound;
    my ($cas) = $dbh->selectrow_array("SELECT CAS FROM Compounds WHERE compound=?", {}, $compound);
    return $cas ? a({-href => "http://commonchemistry.org/ChemicalDetail.aspx?ref=$cas"}, $compound)
        : $compound;
}
