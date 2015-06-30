#!/usr/bin/perl
#######################################################
## Utils.pm -- for cgi code
##
## Copyright (c) 2015 All Right Reserved by UC Berkeley
##
## Authors:
## Wenjun Shao (wjshao@berkeley.edu) and Morgan Price
#######################################################

package Utils;
require Exporter;

use strict;
use warnings;

use DBI;
use Time::HiRes;
use Carp;

sub start_page($);
sub get_style();
sub tabsGene($$$$$$$);
sub tabsExp($$$$$$$);
sub fail($$);
sub formatFASTA($$);
sub crc64($);
sub get_color($);
sub get_dbh(); # open the sqlite db and return the resulting database handle
sub blast_db();
sub tmp_dir();
sub orginfo($);
sub get_orths($$$);

#--------------------------------------------------------

sub start_page($) {
    my ($title) = @_;
    my $header = <<"EOT";
    <!DOCTYPE html>
    <html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US">
    <head>
    <meta http-equiv="X-UA-Compatible" content="IE=edge"> 
    <link rev="made" href="mailto:morgannprice%40yahoo.com" />
    <title>$title</title>
    <meta name="copyright" content="copyright 2015 UC Berkeley" />
    <link rel="shortcut icon" href="../images/favicon.ico" type="image/x-icon">
    <link rel="icon" href="../images/favicon.ico" type="image/x-icon">
    <link rel="stylesheet" href="../images/feba.css">
    <link href='http://fonts.googleapis.com/css?family=Montserrat:700' rel='stylesheet' type='text/css'>
    <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
        </head>
    <body>
    <div id="page">
    <div id="nav"> <div class="box">
        <li class="header">Fitness Browser</li>
        <li><a href="myFrontPage.cgi">Home</a></li>
        <li><a href="geneSearch.cgi">Find Gene</a></li>
        <li><a href="blastSearch.cgi">BLAST</a></li>
        <li><a href="orgAll.cgi">Organisms</a></li>
        <li><a href="expSearch.cgi">Experiments</a></li>
        <li><a href="help.cgi">Help</a></li>
        </div></div>
      
    <div id="main">
EOT
    return $header;
}

sub get_style() {
    my $style = <<"EOT";
    body {
        font-family: verdana, arial, sans-serif;
    }
    H2 {
        color: red;
        /*border-bottom: 1pt solid;*/
    }
    H4 {
        color: blue;
    }
    table {
        border: black thin solid;
    }
    th,td {
        border: lightgrey thin solid;
    }

.axis path,
.axis line {
  fill: none;
  stroke: #000;
  shape-rendering: crispEdges;
}

.dot {
  stroke: none;
}

.line {
  position: relative;
  margin: 0 auto;
  height: 7px;
  z-index: 0;
  margin-top: 5px;
  margin-bottom: 0 auto;
  padding-left: 10px;
  padding-right: 10px;
}

.line2 {
  position: absolute;
  top: 0px;
  height: 7px;
  z-index: 2;
}

EOT
    return $style;
}

sub fail($$) {
    my ($cgi, $notice) = @_;
    print $cgi->h3(qq(Sorry: $notice));
    print $cgi->h4(qq(<a href="myFrontPage.cgi">Go back to front page</a>));
    print $cgi->end_html;
    exit;
}

sub endHtml($) {
    my ($cgi) = @_;
    # print $cgi->p(qq(<a href="myFrontPage.cgi">Go back to front page</a>));
    print $cgi->end_html;
    exit 0;
}

#--------------------------------------------------------

sub formatFASTA($$)
{
    my ($title, $sequence) = @_;
    return undef unless ($sequence);
    my $s = qq/>$title\n/;
    my @whole = $sequence =~ /.{1,60}/g;
    my $line = join "\n", @whole;
    return $s."$line\n";
}

#--------------------------------------------------------

# color code fitness value in range (-3 .. 0 .. 3) to blue => white => yellow
sub fitcolor($) {
    my ($fit) = @_;
    return " #BFBFBF" if !defined $fit;
    my $perc = ($fit + 3)/6;
    $perc = $perc > 1 ? 1 : $perc;
    $perc = $perc < 0 ? 0 : $perc;
    return fractioncolor($perc);
}

# from a fraction 0:1 to a color
sub fractioncolor($) {
#RED (255, 0, 0)
#YELLOW (255, 255, 0)
#BLUE (0, 0, 255)
#WHITE (255, 255, 255)

    my ($perc) = shift @_;
    my ($red, $blue, $green);

# I need "BLUE-WHITE-YELLOW" scheme
    if($perc >= 0.5) {
        $red = 255;
        $green = $red;
        $blue = int(255 * 2 * (1 - $perc) + 0.5);
    } else {
        $red = int(($perc/0.5)*255 + 0.5);
        $green = $red;
        $blue = 255;
    }

    my $string=sprintf ("#%2.2X%2.2X%2.2X",$red,$green,$blue);
    return ($string);
}

sub get_dbh() {
    my $database = "../cgi_data/feba.db";
    return DBI->connect("dbi:SQLite:dbname=$database","","",{ RaiseError => 1 }) || die $DBI::errostr;
}

sub blast_db() {
    return "../cgi_data/aaseqs";
}

sub tmp_dir() {
    return "../tmp";
}

####
sub gene_has_fitness($$$) {
    my ($dbh,$orgId,$locusId) = @_;
    die unless defined $dbh && $orgId && defined $locusId;
    my @result = $dbh->selectrow_array("SELECT expName FROM GeneFitness WHERE orgId=? AND locusId=? LIMIT 1",
				       undef, $orgId, $locusId);
    return scalar(@result) > 0 ? 1 : 0;
}

sub gene_has_cofitness($$$) {
    my ($dbh,$orgId,$locusId) = @_;
    die unless defined $dbh && $orgId && defined $locusId;
    my @result = $dbh->selectrow_array("SELECT cofit FROM Cofit WHERE orgId=? AND locusId=? LIMIT 1",
                       undef, $orgId, $locusId);
    return scalar(@result) > 0 ? 1 : 0;
}

# returns a hash of orgId => attribute => value, including "genome" for the display name
sub orginfo($) {
    my ($dbh) = @_;
    my $orginfo = $dbh->selectall_hashref("SELECT * FROM Organism", "orgId");
    while (my ($orgId,$row) = each %$orginfo) {
	$row->{genome} = join(" ", $row->{genus}, $row->{species}, $row->{strain});
    }
    return $orginfo;
}

# returns a hash of expName => attribute => value
sub expinfo($$) {
    my ($dbh,$orgId) = @_;
    return $dbh->selectall_hashref("SELECT * FROM Experiment WHERE orgId = ?", "expName", {}, $orgId);
}

# Like cmp but for experiments, for ordering by group and condition/concentration, with
# expDesc as a fallback.
# Each argument should be a hash corresponding to a row in the Experiment table
sub CompareExperiments($$) {
    my ($expA,$expB) = @_;
    die unless defined $expA->{expGroup} && defined $expB->{expGroup};
    return $expA->{expGroup} cmp $expB->{expGroup}
           || lc($expA->{condition_1}) cmp lc($expB->{condition_1})
           || $expA->{concentration_1} cmp $expB->{concentration_1}
           || lc($expA->{expDesc}) cmp lc($expB->{expDesc});
}

# Should check that orgId is valid (if it is not empty) before calling.
# Allows partial-word matches in Condition_1 or Condition_2, or full word matches to expDesc or expDescLong, or exact match to Group.
# Note is not case sensitive
sub matching_exps($$$) {
    my ($dbh,$orgId,$expSpec) = @_;
    die if !defined $dbh || !defined $orgId || !defined $expSpec;
    # make the query safe to include in SQL
    $expSpec =~ s/ +$//;
    $expSpec =~ s/^ +$//;
    $expSpec =~ s/[\'\"\n\r]//g;

    my $orgClause = $orgId eq "" ? "" : qq{ AND orgId = "$orgId"};
    my $sql = qq{SELECT * from Organism JOIN Experiment USING (orgId)
             WHERE (expName = "$expSpec"
	            OR Condition_1 LIKE "$expSpec%" OR Condition_1 LIKE "% $expSpec%"  OR Condition_1 LIKE "%-$expSpec%"
		    OR Condition_2 LIKE "$expSpec%" OR Condition_2 LIKE "% $expSpec%"  OR Condition_2 LIKE "%-$expSpec%"
	     	    OR expGroup = '$expSpec'
	            OR expDesc = '$expSpec' OR expDesc LIKE "$expSpec %" OR expDesc LIKE "% $expSpec" OR expDesc LIKE "% $expSpec %"
	            OR expDescLong = '$expSpec' OR expDescLong LIKE "$expSpec %" OR expDescLong LIKE "% $expSpec" OR expDescLong LIKE "% $expSpec %")
	     $orgClause
	     ORDER BY genus, species, strain, expGroup, condition_1, concentration_1, expDesc;};
         # die $sql;
    return $dbh->selectall_arrayref($sql, { Slice => {} });
}

# Returns a reference to a hash from ortholog's orgId => gene
# Each gene is a hash that includes the fields in Gene and the ratio
sub get_orths($$$) {
    my ($dbh,$orgId,$locusId) = @_;
    die unless defined $orgId && defined $locusId;
    return $dbh->selectall_hashref(qq{ SELECT Gene.*, Ortholog.ratio FROM GENE JOIN ORTHOLOG
					ON Ortholog.orgId2=Gene.orgId AND Ortholog.locusId2=Gene.locusId
					AND orgId1=? AND locusId1=? },
	                           "orgId", {}, $orgId, $locusId);
}

# returns a simple summary of the gene's fitness pattern, like "no data", "weak", "strong", or "cofit",
# and a longer piece of text that can be used as a tooltip
sub gene_fit_string($$$) {
    my ($dbh,$orgId,$locusId) = @_;
    die unless defined $orgId && defined $locusId;
    my ($minFit,$maxFit,$minT,$maxT) = $dbh->selectrow_array(qq{ SELECT min(fit), max(fit), min(t), max(t)
                                                                 FROM GeneFitness WHERE orgId = ? AND locusId = ? ; },
							     {}, $orgId, $locusId);
    return ("no data", "no fitness data for this gene ") if !defined $minFit;
    my $tip = sprintf("fit = %.1f to +%.1f (t = %.1f to +%.1f)", $minFit, $maxFit, $minT, $maxT);
    my ($maxCofit) = $dbh->selectrow_array(qq{ SELECT cofit FROM Cofit WHERE orgId = ? AND locusId = ? AND rank = 1 LIMIT 1; },
					   {}, $orgId, $locusId);
    $tip .= sprintf(", max(cofit) = %.2f", $maxCofit) if defined $maxCofit;
    return ("strong",$tip) if ($minT < -5 && $minFit < -2) || ($maxT > 5 && $maxFit > 2);
    return ("cofit",$tip) if defined $maxCofit && $maxCofit > 0.75;
    return ("sig.",$tip) if ($minT < -4 && $minFit < -1) || ($maxT > 4 && $maxFit > 1);
    return ("weak", $tip);
}



sub tabsGene($$$$$$$) {
    my($dbh,$cgi,$orgId,$locusId,$showAll,$type,$curr) = @_;
    my $hasFit = gene_has_fitness($dbh, $orgId, $locusId);
    my $hasCofit = gene_has_cofitness($dbh, $orgId, $locusId);
    # print "hasFit ".$hasFit." hasCofit ".$hasCofit;

    my $code = qq[<div id="tabs">];

    my $gClass = "gene";
    $gClass = "selected" if $curr eq "gene";
    $code .= qq[<li class="$gClass"> <a href="geneOverview.cgi?orgId=$orgId&gene=$locusId">Gene</a></li>];

    my $fClass = "gene";
    $fClass = "selected" if $curr eq "fit";
    $fClass = "disable" if $hasFit == 0;

    $code .= qq[<li class="$fClass"> <a href="myFitShow.cgi?orgId=$orgId&gene=$locusId&showAll=0">Fitness</a></li>];
    # magic switching tabs for fitness!
    # if ($showAll == 0) {
    #     my $showAllDest = qq(myFitShow.cgi?orgId=$orgId&gene=$locusId&showAll=1);
    #     $code .= qq[<li class="$fClass"> <a href="$showAllDest">All Fitness Data</a></li>];
    #     # print $cgi->p(qq(<a href=$showAllDest>All fitness data</a>));
    # } else {
    #     my $showFewDest = qq(myFitShow.cgi?orgId=$orgId&gene=$locusId&showAll=0);
    #     $code .= qq[<li class="$fClass"> <a href="$showFewDest">Strongest Phenotypes</a></li>];
    #     # print $cgi->p(qq(<a href=$showFewDest>Strongest phenotypes</a>));
    # }


    my $nClass = "gene";
    $nClass = "selected" if $curr eq "nearby";
    $code .= qq[<li class="$nClass"> <a href="genesFit.cgi?orgId=$orgId&locusId=$locusId&around=2">Nearby</a></li>];
    my ($maxCofit) = $dbh->selectrow_array(qq{ SELECT cofit FROM Cofit WHERE orgId = ? AND locusId = ? AND rank = 1 LIMIT 1; },
                           {}, $orgId, $locusId);

    my $max = sprintf("(max cofit %.2f)", $maxCofit) if defined $maxCofit;
    my $cClass = "gene";
    $cClass = "selected" if $curr eq "cofit";
    $cClass = "disable" if $hasCofit == 0;
    $code .= qq[<li class="$cClass"><a href="cofit.cgi?orgId=$orgId&locusId=$locusId" title="$max">Cofit</a></li>];

    if ($type == 1) {
        my $hClass, my $pClass = "gene";
        if ($curr eq "protein") {
            $pClass = "selected";
        } elsif ($curr eq "homo") {
            $hClass = "selected";
        }
        $code .= qq[<li class="$pClass"><a href="domains.cgi?orgId=$orgId&locusId=$locusId">Protein</a></li>];
        $code .= qq[<li class="$hClass"><a href="mySeqSearch.cgi?orgId=$orgId&locusId=$locusId">Homologs</a></li>];

    }
    $code .= qq[</div><div id="tabcontent">];

    return $code;
    

#     my @links = ();
#     if ($gene->{locusId} =~ m/^\d+$/) {
#         push @links, $cgi->a({href => "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=$gene->{locusId}"},
#                  "MicrobesOnline");
#     }
#     if ($orgId eq "Keio" && $gene->{sysName} =~ m/^b\d+$/) {
#         push @links, $cgi->a({href => "http://ecocyc.org/ECOLI/search-query?type=GENE&gname=$gene->{sysName}"}, "EcoCyc");
#     }
#     print $cgi->p("Links: " . join(", ", @links)) if (@links > 0);
# } #  end if just 1 hit
    # $code .= q['']

}

sub tabsExp($$$$$$$) {
    my($dbh,$cgi,$orgId,$expName,$expGroup,$cond,$curr) = @_;
    # print "show: $show, curr: $curr";

    my $code = qq[<div id="tabs">];

    my $oClass = "exp";
    $oClass = "selected" if $curr eq "overview";
    $code .= qq[<li class="$oClass"> <a href="exp.cgi?orgId=$orgId&expName=$expName">Overview</a></li>];

    my $sClass = "exp";
    $sClass = "selected" if $curr eq "specific";
    $code .= qq[<li class="$sClass"> <a href="exp.cgi?orgId=$orgId&expName=$expName&show=specific" title="Specific Phenotypes">Specific</a></li>];

    my $nClass = "exp";
    $nClass = "selected" if $curr eq "-gene"; #and $show eq "important");
    $code .= qq[<li class="$nClass"> <a href="exp.cgi?orgId=$orgId&expName=$expName&show=important" title="Important Genes">- Genes</a></li>];

    my $pClass = "exp";
    $pClass = "selected" if $curr eq "+gene";
    $code .= qq[<li class="$pClass"> <a href="exp.cgi?orgId=$orgId&expName=$expName&show=detrimental" title="Detrimental Genes">+ Genes</a></li>];

    my $mClass = "exp";
    $mClass = "selected" if $curr eq "metrics";
    $code .= qq[<li class="$mClass"><a href="exp.cgi?orgId=$orgId&expName=$expName&show=quality" title="Quality Metrics">Metrics</a></li>];

    # my $cClass = "exp";
    # $cClass = "selected" if $curr eq "compare";
    # $code .= qq[<li class="$cClass"><a href="expComp.cgi?orgId=$orgId&expName=$expName">Compare</a>];

    # if ($type == 1) {
    #     my $hClass, my $pClass = "gene";
    #     if ($curr eq "protein") {
    #         $pClass = "selected";
    #     } elsif ($curr eq "homo") {
    #         $hClass = "selected";
    #     }
    #     $code .= qq[<li class="$pClass"><a href="domains.cgi?orgId=$orgId&locusId=$locusId">Protein</a></li>];
    #     $code .= qq[<li class="$hClass"><a href="mySeqSearch.cgi?orgId=$orgId&locusId=$locusId">Homologs</a></li>];
    # }

    $code .= qq[</div><div id="tabcontent">];

    return $code;

}



#END 

return 1;

