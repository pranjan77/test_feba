#!/usr/bin/perl -w
#######################################################
## exp.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required parameters: orgId and expName, for organism and which experiment
# Optional parameter show:
# specific -- show specific phenotypes
# important or detrimental -- show top 200 genes either way
# quality -- show quality metrics
# (by default, shows all of these optoins)

use strict;
use CGI qw(:standard Vars -nosticky);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();
print $cgi->header;

my $orgId = $cgi->param('orgId') || die "no species identifier\n";
die "Invalid species name!" unless $orgId =~ m/^[A-Za-z0-9_]*$/;
my $expName = $cgi->param('expName') || die "no experiment identifier\n";
die "Invalid experiment name!" unless $expName =~ m/^[A-Za-z0-9_]*$/;
my $show = $cgi->param('show') || "";

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

print $cgi->start_html(
    -title =>"Experiment $expName for $orginfo->{$orgId}{genome}",
    -style => {-code => $style},
    -author=>'morgannprice@yahoo.com',
    -meta=>{'copyright'=>'copyright 2015 UC Berkeley'},
);

print $cgi->h2("Experiment $expName for $orginfo->{$orgId}{genome}");
print $cgi->h3($exp->{expDescLong});
my $cond1 = $exp->{condition_1} ? join(" ", $exp->{condition_1}, $exp->{concentration_1}, $exp->{units_1}) : "";
my $cond2 = $exp->{condition_2} ? join(" ", $exp->{condition_2}, $exp->{concentration_2}, $exp->{units_2}) : "";
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

print $cgi->p(join("<BR>", "Media: $media", join(", ",@culture), "By: $exp->{person} on $exp->{dateStarted}"));

my @fit = (); # sorted list of fitness values to show
my $header = undef;
if ($show eq "specific") {
    $header = "Genes with specific phenotypes:";
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
	push @trows, $cgi->Tr({align => 'left', valign => 'top'},
			      $cgi->td(checkbox('locusId',0,$row->{locusId},'')),
			      # $cgi->td(qq{<INPUT TYPE="CHECKBOX" name="locusId" value="$row->{locusId}" unchecked}),
			      $cgi->td($cgi->a({href => "myFitShow.cgi?orgId=$orgId&gene=$row->{locusId}",
					       style => "color:rgb(0,0,0)"},
						$geneId)),
			      $cgi->td($row->{gene}),
			      $cgi->td({ -bgcolor => Utils::fitcolor($row->{fit}) },
				       sprintf("%.1f", $row->{fit})),
			      $cgi->td( sprintf("%.1f", $row->{t}) ),
			      $cgi->td($row->{desc}));
    }
    print
	start_form(-name => 'input', -method => 'GET', -action => 'genesFit.cgi'),
	hidden('orgId', $orgId),
	$cgi->table(
	    { cellspacing => 0, cellpadding => 3},
	    $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
		     $cgi->th(['&nbsp;', 'gene','name','fitness','t score','description']),
		     @trows) ),
	submit(-name => 'Top fitness of selected genes'),
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
	print $cgi->p($cgi->a({href => "exp.cgi?orgId=$orgId&expName=$expName&show=specific"},
			      "See specific phenotypes for " . scalar(@$spec). " genes"));
    }  else {
	print $cgi->p("No genes had specific phenotypes in this experiment.");
    }
}
print $cgi->p($cgi->a({href => "exp.cgi?orgId=$orgId&expName=$expName&show=important"}, "Show important genes"))
    if $show ne "important";
print $cgi->p($cgi->a({href => "exp.cgi?orgId=$orgId&expName=$expName&show=detrimental"}, "Show detrimental genes"))
    if $show ne "detrimental";
print $cgi->p($cgi->a({href => "exp.cgi?orgId=$orgId&expName=$expName&show=quality"}, "Show quality metrics"))
    if $show ne "quality";

$dbh->disconnect();
Utils::endHtml($cgi);
