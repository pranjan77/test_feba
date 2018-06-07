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
die "Invalid species name!" unless $orgId =~ m/^[A-Za-z0-9_-]*$/;
my $expName = $cgi->param('expName') || die "no experiment identifier\n";
$expName =~ s/^[ \t]+//;
$expName =~ s/^[ \t\r\n]+$//;
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

if ($help) {
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
    $media .= " ($exp->{mediaStrength}x)" if $exp->{mediaStrength} != 1;
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
    if ($exp->{pubId}) {
      my $pub = $dbh->selectrow_hashref("SELECT * from Publication WHERE pubId = ?",
                                        {}, $exp->{pubId});
      if (defined $pub) { # should always exist
        push @pieces, "Reference: " . a({ -href => $pub->{URL}, -title => $pub->{title} },
                                     $pub->{pubId});
      }
    }

    my $mediaComponents = $dbh->selectall_arrayref(qq{SELECT * from MediaComponents LEFT JOIN Compounds USING (compound)
                                                          WHERE media = ? },
                                                   { Slice => {} },
                                                   $exp->{media});
    if (@$mediaComponents > 0) {
        my %compStrings = (); # mix => list of components; note mix = "" for most
        foreach my $row (@$mediaComponents) {
            my $compString = $row->{CAS} ?
                a({-href => "http://commonchemistry.org/ChemicalDetail.aspx?ref=$row->{CAS}"}, $row->{compound})
                :  $row->{compound};
            if ($row->{concentration} && $row->{units}) {
              my $conc = $row->{concentration} * $exp->{mediaStrength};
              $compString = join(" ", $conc, $row->{units}, $compString);
              $row->{mix} = "" if !defined $row->{mix};
              push @{ $compStrings{$row->{mix}} }, $compString;
            }
        }
        my $comp = "Media components: " . join(", ", @{ $compStrings{""} });
        foreach my $mix (sort keys %compStrings) {
            next if $mix eq "";
            $comp .= ", $mix " . small("(" . join(", ", @{ $compStrings{$mix} }) . ")");
        }
        $comp .= " " . small("(final concentrations)") if $exp->{mediaStrength} != 1;
        push @pieces, $comp;
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

my %cons = (); # locusId => nInOG if it has a specific phenotype in this condition
if ($exp->{condition_1} ne "") {
    my $specOG = $dbh->selectall_arrayref("SELECT locusId,nInOG from SpecOG WHERE orgId = ? AND expGroup = ? AND condition = ?",
                                          {}, $orgId, $exp->{expGroup}, $exp->{condition_1});
    foreach my $row (@$specOG) {
        my ($locusId,$nInOG) = @$row;
        $cons{$locusId} = $nInOG;
    }
}

print $cgi->h3($header) if defined $header;

if (@fit > 0) { # show the table
  my @trows = ();
  foreach my $row (@fit) {
    my $geneId = $row->{sysName} || $row->{locusId};
    my $strainUrl = "strainTable.cgi?orgId=$orgId&locusId=$row->{locusId}&expName=$expName";
    $strainUrl .= "&help=1" if $help;
    my $orthUrl = "orthFit.cgi?orgId=$orgId&locusId=$row->{locusId}"
      . "&expGroup=" . uri_escape($exp->{expGroup})
        . "&condition1=" . uri_escape($exp->{condition_1});
    $orthUrl .= "&help=1" if $help;
    my $alt_desc = Utils::alt_descriptions($dbh,$orgId,$row->{locusId});
    push @trows,
      $cgi->Tr( {align => 'left', valign => 'top'},
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
                td( a({ -title => $alt_desc || "no other information",
                        -href => "domains.cgi?orgId=$orgId&locusId=$row->{locusId}" },
                      $row->{desc}) ),
                td(a({ 
                      title => "Compare to data from similar experiments or orthologs",
                      href => $orthUrl},
                     exists $cons{$row->{locusId}} && $cons{$row->{locusId}} > 1?
                     "<i>conserved</i>" : "compare")) );
  }
  print
    start_form(-name => 'input', -method => 'GET', -action => 'genesFit.cgi'),
      hidden('orgId', $orgId),
	$cgi->table( { cellspacing => 0, cellpadding => 3},
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

    print h3("Specific Phenotypes");
    if (@$spec > 0) {
	print p(a({href => "exp.cgi?orgId=$orgId&expName=$expName&show=specific"},
                  "For " . scalar(@$spec). " genes in this experiment"));
        print p(a({href => "spec.cgi?orgId=$orgId&expGroup=".uri_escape($exp->{expGroup})."#".$exp->{condition_1} },
                  "For $exp->{expGroup} $exp->{condition_1} in $orginfo->{$orgId}{genome}"))
            if $exp->{expGroup} && $exp->{condition_1};
    }  else {
	print $cgi->p("None in this experiment");
        if ($exp->{expGroup}) {
            print
                p(a({href => "spec.cgi?orgId=$orgId&expGroup=" . uri_escape($exp->{expGroup})
                         . ($exp->{condition_1} eq "" ? "" : "#" . uri_escape($exp->{condition_1}))},
                    "For $orginfo->{$orgId}{genome} in $exp->{expGroup} experiments"));
        }
    }
}


if ($exp->{condition_1} ne "") {
    print
        p(a({href => "orthCond.cgi?expGroup=" . uri_escape($exp->{expGroup})
                 . "&condition1=" . uri_escape($exp->{condition_1})},
            ($show eq "specific" ? "Specific phenotypes for" : "For")
            . " $exp->{expGroup} $exp->{condition_1} across organisms"));
}

if ($show ne "specific") {
    print
        h3("Metabolic Maps"),
        p("Color code by fitness: see",
          a({href => "keggmap.cgi?mapId=01100&orgId=$orgId&expName=$expName"}, "overview map"),
          "or",
          a({href => "keggmaplist.cgi?orgId=$orgId&expName=$expName"}, "list of maps")."." );

    if (@$spec > 0) {
        # try to highlight useful maps that contain genes with specific phenotypes
        my $locusIn = join(",", map "'".$_->{locusId}."'", @$spec);
        my %specEc = (); # ec => 1 for those locusIds
        my $ecTIGR = $dbh->selectcol_arrayref(
            "SELECT ec FROM GeneDomain WHERE orgId = ? AND locusId IN ( $locusIn );",
            {}, $orgId);
        foreach my $ec (@$ecTIGR) { $specEc{$ec} = 1; }
        my $ecKEGG = $dbh->selectcol_arrayref(
            qq{ SELECT ecnum FROM BestHitKEGG
                JOIN KEGGMember USING (keggOrg,keggId)
                JOIN KgroupEC USING (kgroup)
                WHERE orgId = ? AND locusId IN ( $locusIn ); },
            {}, $orgId);
        foreach my $ec (@$ecKEGG) { $specEc{$ec} = 1; }
        my $ecSEED = $dbh->selectcol_arrayref(
            "SELECT num FROM SEEDClass WHERE orgId = ? AND locusId IN ( $locusIn ) AND type = 1;",
            {}, $orgId);
        foreach my $ec (@$ecSEED) { $specEc{$ec} = 1; }
        if (keys(%specEc) > 0) {
            my $ecIn = join(",", map "'".$_."'", keys %specEc);
            my $maps = $dbh->selectall_arrayref(
                qq{ SELECT mapId, title, COUNT(DISTINCT objectId) nEc
                    FROM KEGGConf JOIN KEGGMap USING (mapId)
                    WHERE objectId IN ( $ecIn ) AND type=1
                    GROUP BY mapId
                    ORDER BY nEc DESC });
            if (scalar(@$maps) > 0) {
                my @mapShow = ();
                foreach my $map (@$maps) {
                    my ($mapId,$title,$nEc) = @$map;
                    push @mapShow, li(a({href => "keggmap.cgi?orgId=$orgId&expName=$expName&mapId=$mapId"},
                                        $title));
                }
                print p("Maps containing gene(s) with specific phenotypes:",
                        ul(@mapShow));
            }
        }
    }
}
$dbh->disconnect();
Utils::endHtml($cgi);

sub CompoundToHTML($) {
    my ($compound) = @_;
    return $compound if ! $compound;
    $compound =~ s/^supernat[ea]nt; /supernatant from /i;
    my ($cas) = $dbh->selectrow_array("SELECT CAS FROM Compounds WHERE compound=?", {}, $compound);
    return $cas ? a({-href => "http://commonchemistry.org/ChemicalDetail.aspx?ref=$cas"}, $compound)
        : $compound;
}
