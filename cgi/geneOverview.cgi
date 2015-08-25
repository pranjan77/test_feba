#!/usr/bin/perl -w
#######################################################
## myFitShow.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo, Wenjun Shao (wjshao@berkeley.edu) and 
## Morgan Price
#######################################################
#
# Required CGI parameters:
# gene -- a locusId, sysName, or gene name to match on
#	(may show multiple hits)
# Optional CGI parameters:
# orgId -- which organism to search in
# showAll -- 1 if showing all fitness values instead of just the most extreme ones

use strict;

use CGI qw(:standard Vars -nosticky);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;
print $cgi->header;
my $style = Utils::get_style();

my $orgSpec = $cgi->param('orgId') || "";
my $geneSpec = $cgi->param('gene');
my $showAll = $cgi->param('showAll') ? 1 : 0;
my $start = Utils::start_page("Overview for $geneSpec ($orgSpec)");

$geneSpec =~ s/ *$//;
$geneSpec =~ s/^ *//;

if (!defined $geneSpec || $geneSpec eq "") {
    Utils::fail($cgi, "you must enter the gene name or locus tag");
}

# check user input
Utils::fail($cgi, "$geneSpec is invalid. Please enter correct gene name!") unless ($geneSpec =~ m/^[A-Za-z0-9_-]*$/);

# connect to database
my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);

my $query = qq{SELECT * FROM Gene
		WHERE ( locusId = ? OR sysName = ? OR upper(gene) = upper(?) )};
my $hits;
if ($orgSpec) {
    $query .= " AND orgId = ?";
    $hits = $dbh->selectall_arrayref($query, { Slice => {} }, $geneSpec, $geneSpec, $geneSpec, $orgSpec);
} else {
    $hits = $dbh->selectall_arrayref($query, { Slice => {} }, $geneSpec, $geneSpec, $geneSpec);
}

# and add another column for whether there is fitness data
foreach my $gene (@$hits) {
    $gene->{has_fitness} = Utils::gene_has_fitness($dbh,$gene->{orgId},$gene->{locusId});
}

if (@$hits == 0) {
	print $start,
    print $cgi->h3("No gene found for $geneSpec",
		   (exists $orginfo->{$orgSpec}{genome} ? " in " . $orginfo->{$orgSpec}{genome} : ""));
} elsif (@$hits > 1) {
	print $start,

	h3("Genes found for $geneSpec:");
    my @trows = ();
    push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
			  $cgi->th( [ 'geneId','sysName','geneName','description','genome','fitness' ] ) );
    foreach my $gene (@$hits) {
	my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $gene->{orgId}, $gene->{locusId});
	my @trow = map $cgi->td($_), ($gene->{locusId}, $gene->{sysName}, $gene->{gene}, $gene->{desc},
				      $cgi->a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId}}, "$orginfo->{$gene->{orgId}}->{genome}"),
				      a( {href => "myFitShow.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}", title => $fittitle, },
					 $fitstring));
	push @trows, $cgi->Tr(@trow);
    }
    
    print $cgi->table( { cellspacing=>0, cellpadding=>3 }, @trows);
    # PRINT TABLE

} else {
    # just 1 hit
    my $gene = $hits->[0];
    my $orgId = $gene->{orgId};
    my $locusId = $gene->{locusId};
    my $scaffold = $gene->{scaffoldId};
    my $begin = $gene->{begin};
    my $end = $gene->{end};
    my $strand = $gene->{strand};
    my $type = $gene->{type};
    my $typeName = "";
    $typeName = "Protein-coding gene" if $type == 1;
    $typeName = "23S (large subunit) ribosomal RNA" if $type == 2;
	$typeName = "16S (small subunit) ribosomal RNA" if $type == 3;
	$typeName = "5S ribosomal RNA" if $type == 4;
	$typeName = "Transfer RNA" if $type == 5;
	$typeName = "Other non-coding RNA" if $type == 6;
	$typeName = "Pseudogene derived from a protein-coding gene" if $type == 7;
	$typeName = "Pseudogene derived from an RNA gene" if $type == 8;
	$typeName = "CRISPR" if $type == 9;
	$typeName = "CRISPR spacer" if $type == 10;
	$typeName = "Antisense RNA" if $type == 11;
	$typeName = "Unclassified feature (possibly a pseudogene)" if $type == 99;

	#nearby 5 genes
	my @locusIds = $geneSpec;
	my $idShow = $gene->{sysName} || $gene->{locusId};
	my %spacingDesc = (); # locusId => spacing description
	my $type;
    die "Cannot specify nearby with multiple locusIds or with addgene" unless @locusIds == 1;
    my $centralId = $locusIds[0];
    my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId = ? AND locusId = ?",
				       {}, $orgId, $centralId);
    $type = $gene->{type};
    my $scgenes = $dbh->selectall_arrayref("SELECT * from Gene where orgId = ? AND scaffoldId = ? ORDER BY begin",
					   { Slice => {} }, $orgId, $gene->{scaffoldId});
    die "Cannot find genes for $gene->{scaffoldId}" unless scalar(@$scgenes) > 0;
    my ($iCentral) = grep { $scgenes->[$_]{locusId} eq $centralId } (0..(scalar(@$scgenes)-1));
    die if !defined $iCentral;
    my $i1 = $iCentral - 5;
    $i1 = 0 if $i1 < 0;
    my $i2 = $iCentral + 5;
    $i2 = scalar(@$scgenes)-1 if $i2 > scalar(@$scgenes);
    @locusIds = map { $scgenes->[$_]{locusId} } ($i1..$i2);

 # 	my @genes = ();
	# foreach my $locusId (@locusIds) {
 #    	my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId = ? AND locusId = ?",
	# 				       {}, $orgSpec, $locusId);
 #    	# print $gene->{desc};
	#     die "No such locus $locusId in org $orgId" if !defined $gene->{locusId};
	# }


    my @genes = ();
	foreach my $locusId (@locusIds) {
	    my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId = ? AND locusId = ?",
					       {}, $orgSpec, $locusId);
	    die "No such locus $locusId in org $orgId" if !defined $gene->{locusId};
	    # expName => "fit" => fitness value
	    $gene->{fit} = $dbh->selectall_hashref(qq{SELECT expName,fit,t FROM GeneFitness
	                                               WHERE orgId = ? AND locusId = ?},
						   "expName", {}, $orgId, $locusId);
	    foreach my $expName (keys %{ $gene->{fit} }) {
		# die "No such experiment: $expName" unless exists $expinfo->{$expName};
	    }
	    $gene->{nExps} = scalar(keys %{ $gene->{fit} });
	    push @genes, $gene;
	}



    if ($hits->[0]{has_fitness} == 0) {
		my $start = Utils::start_page("Fitness data for $idShow in $orginfo->{$orgId}{genome}");
		my $tabs = Utils::tabsGene($dbh,$cgi,$orgId,$locusId,0,$gene->{type},"gene");
		
		print 
			$start, $tabs, 
			h2("Gene $idShow in " . $cgi->a({href => "org.cgi?orgId=$orgId"}, "$orginfo->{$orgId}{genome}")),
			$cgi->h3("$idShow $gene->{gene}: $gene->{desc} in " . $orginfo->{$gene->{orgId}}{genome}),
			# $cgi->p("Sorry, no fitness data for $idShow.");
    } else {
		# show fitness data for gene
			my @fit = @{ $dbh->selectall_arrayref("SELECT expName,fit,t from GeneFitness where orgId=? AND locusId=?",
						      { Slice => {} },
						      $orgId, $locusId) };
			my $nTotValues = scalar(@fit);
			die "Unreachable" if $nTotValues == 0;
			my $limitRows = $showAll ? $nTotValues : 20;
			my $minAbsFit;
		if ($nTotValues > $limitRows) {
		    # subset the rows
		    @fit = sort { abs($b->{fit}) <=> abs($a->{fit}) } @fit;
		    @fit = @fit[0..($limitRows-1)];
		    $minAbsFit = abs($fit[$#fit]{fit});
		}

		# and get metadata about experiments
		my $expinfo = Utils::expinfo($dbh,$orgId);

		if ($showAll) {
		    @fit = sort { Utils::CompareExperiments($expinfo->{$a->{expName}}, $expinfo->{$b->{expName}}) } @fit;
		} else {
		    @fit = sort { $a->{fit} <=> $b->{fit} } @fit;
		}

		my $start = Utils::start_page("Overview for $idShow ($orginfo->{$orgId}{genome})");
		my $title = "Gene Info for $idShow in $orginfo->{$orgId}{genome}";
		my $tabs = Utils::tabsGene($dbh,$cgi,$orgSpec,$geneSpec,$showAll,$gene->{type},"gene");


		print
			$start, $tabs,
		    h2("Gene info for $idShow in " . $cgi->a({href => "org.cgi?orgId=$orgId"}, "$orginfo->{$orgId}{genome}")),
		    h3("$idShow $gene->{gene}: $gene->{desc}");
	}

   	print "Type $type: $typeName<BR>
   	Located on scaffold $scaffold, $strand strand, nucleotides $begin - $end";


	my @links = ();
	if ($gene->{locusId} =~ m/^\d+$/) {
	push @links, $cgi->a({href => "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=$gene->{locusId}"},
			     "MicrobesOnline");
	}
	if ($orgId eq "Keio" && $gene->{sysName} =~ m/^b\d+$/) {
	push @links, $cgi->a({href => "http://ecocyc.org/ECOLI/search-query?type=GENE&gname=$gene->{sysName}"}, "EcoCyc");
	}
	print $cgi->p("Links: " . join(", ", @links)) if (@links > 0);

	print $cgi->h3({style=>'text-align:center'},"Nearby Genes");

	my $lastTot = ($genes[-1]->{end} - $genes[0]->{begin});
	my $width = 800; #'100%'; #input the max-width according to the .arrow class in css
	my $factor = $width/$lastTot;
	my $xScale = 10 * $factor;
	my $scale = 510 * $factor;

	my $svg = qq[
	<center><svg class="arrows" viewBox="0 -30 $width 100">
	<defs>
	<marker id='right' orient="auto"
	    markerWidth='2' markerHeight='4'
	    refX='0.1' refY='2'>
	    <!-- triangle pointing right (+x) -->
	    <path d='M0,0 V4 L2,2 Z' fill="black"/>
	</marker>
	<marker id='right2' orient="auto"
	    markerWidth='2' markerHeight='4'
	    refX='0.1' refY='2'>
	    <!-- triangle pointing right (+x) -->
	    <path d='M0,0 V4 L2,2 Z' fill="red"/>
	</marker>
	<marker id='left' orient="auto"
	    markerWidth='2' markerHeight='4'
	    refX='0.5' refY='2'>
	    <!-- triangle pointing left (-x) -->
	    <path d='M0,2 L2,0 L2,4 Z' fill="black"/>
	</marker>
	<marker id='left2' orient="auto"
	    markerWidth='2' markerHeight='4'
	    refX='0.5' refY='2'>
	    <!-- triangle pointing left (-x) -->
	    <path d='M0,2 L2,0 L2,4 Z' fill="red"/>
	</marker>
	<marker id="scaleEnd" orient="auto" 
		markerWidth='2' markerHeight='4'
		refX='0.5' refY='2'>
		<line x1="0" y1="5" x2="0" y2="-5" style="stroke: black;"/>
	</marker>
	</defs>

	<line 
	    id='scale'
	    marker-start='url(#scaleEnd)'
	    marker-end='url(#scaleEnd)'
	    x1="$xScale" y1="55" x2="$scale" y2="55" style="stroke:black;stroke-width:2" />
	<text x="3" y="53" font-family="Verdana" font-size="10" fill="black">500 nt</text>

	<!--<line 
	    id='arrow-line'
	    marker-end='url(#head)'
	    x1="0" y1="0" x2="100" y2="0" style="stroke:black;stroke-width:2" />
	</svg>-->
	];

	# my @headings = qw{Locus Name Description Strand Distance(nt) Phenotype};
	# my @trows = ( Tr({ -valign => 'top', -align => 'center' }, map { th($_) } \@headings) );
	my @trows = ( Tr({-valign => "top", -align => 'center'}, th([ 'Locus', 'Name', 'Description', 'Strand', a({title=>'Distance between the end of previous gene and beginning of current gene', color=>"#011B47"},'Distance (nt)'), 'Phenotype', ])));

	my $diff = "";
	my $prevrow;
	my $total = 0;
	my $newDist = 0;
	my $newDistAdj = 0;
	my $totalAdj = 0;
	my $num = "15";
	my $pos = 0;
	foreach my $row (@genes) {
		$diff = $row->{begin} - $prevrow->{end} if defined $prevrow;		# $prevrow->{end} = 0 if !defined $prevrow;
		$newDist = $total + $diff;
		$total = $newDist + ($row->{end} - $row->{begin});
		$newDistAdj = $newDist * $factor;
		$totalAdj = $total * $factor;
		my $textX = $newDist + ($total - $newDist)/2;
		my $textXAdj = $textX * $factor - 11-3;
		my $textY = $pos - 2-3;
		my $textYAdj = $textY + 15;
		my $posAdj = $pos + 15;
		my $color = "black";
		my $head = "marker-end='url(#right)'";
		$head = "marker-start='url(#left)'" if $row->{strand} eq "-";
		my $bgcolor = undef;
		my $text = "#00A8C6";
		if ($row->{locusId} eq $geneSpec) {
			$color = "red";
			$text = "red";
			$head = "marker-end='url(#right2)'";
			$head = "marker-start='url(#left2)'" if $row->{strand} eq "-";
			$bgcolor = "#FFFFFF";
		}
		my $label = $row->{gene} || $row->{sysName}; #|| $row->{locusId};
		my $label2 =  $row->{gene} || $row->{sysName} || $row->{locusId};
		my $label3 = $label2;
		$label3 =~ s/^.*_/_/ if $row->{sysName} || $row->{locusId};
		$label2 = $row->{sysName}. ": " . $label2 if $row->{sysName};

		$svg .= qq[
		<g class="click" onclick="window.location.href='singleFit.cgi?orgId=$row->{orgId}&locusId=$row->{locusId}'"><title>$label2 - $row->{desc}, $row->{begin} - $row->{end}</title><line id='arrow-line' $head x1="$newDistAdj" y1="$posAdj" x2="$totalAdj" y2="$posAdj" style="stroke:$color;stroke-width:2" />
        <text x="$textXAdj" y="$textYAdj" font-family="Verdana" font-size="13" fill="$text" onmouseover="this.style.fill='#CC0024'" onmouseout="this.style.fill='#00A8C6'">$label3</text></g>];

		my ($phen, $tip) = Utils::gene_fit_string($dbh,$orgSpec,$row->{locusId});
		 # if $row->{locusId} eq $geneSpec; #f4f3e4
	    push @trows, Tr({ -valign => 'top', -align => 'left', -bgcolor=>"$bgcolor"},
	    	# display result row by row
		    td([ a({href => "geneOverview.cgi?orgId=$orgId&gene=$row->{locusId}"},$row->{sysName}||$row->{locusId}), #locus
			 	a({href => "geneOverview.cgi?orgId=$orgId&gene=$row->{locusId}"},$row->{gene} || $row->{sysName}), 
			 	$row->{desc}, 
			 	$row->{strand},
			 	a({title=>"From $prevrow->{end} to $row->{begin}"},$diff), # $row->{begin},
			 	a({href => "myFitShow.cgi?orgId=$orgId&gene=$row->{locusId}", title=>$tip},$phen),
			 	# a( { href => "exps.cgi?orgId=$orgId&expGroup=$row->{expGroup}"},
			  #   $row->{nExp} ), #experiments
			 ]));
	    $prevrow = $row;
	    if ($pos == 0) {
		    $num = $num * -1;
		    $pos = $num;

	    } else {
	    	$pos = 0;
	    }
	}
	$svg .= "</svg></center>";
	print $svg;
	print table({cellspacing => 0, cellpadding => 3}, @trows);
	   

	# links
# 	if ($showAll == 0) {
# 	    my $showAllDest = qq(myFitShow.cgi?orgId=$orgId&gene=$locusId&showAll=1);
# 	    print $cgi->p(qq(<a href=$showAllDest>All fitness data</a>));
# 	} else {
# 	    my $showFewDest = qq(myFitShow.cgi?orgId=$orgId&gene=$locusId&showAll=0);
# 	    print $cgi->p(qq(<a href=$showFewDest>Strongest phenotypes</a>));
# 	}
# 	print $cgi->p($cgi->a( { href => "genesFit.cgi?orgId=$orgId&locusId=$locusId&around=2" }, "Fitness for nearby genes"));
# 	my ($maxCofit) = $dbh->selectrow_array(qq{ SELECT cofit FROM Cofit WHERE orgId = ? AND locusId = ? AND rank = 1 LIMIT 1; },
# 					       {}, $orgId, $locusId);
# 	print $cgi->p($cgi->a({href => "cofit.cgi?orgId=$orgId&locusId=$locusId"}, "Top cofit genes"),
# 		      sprintf("(max cofit %.2f)", $maxCofit)) if defined $maxCofit;
#     } # end else unique hit has data

#     if ($gene->{type} == 1) {
#     	print
# 			p(a({href => "getSeq.cgi?orgId=$orgId&locusId=$locusId"}, "Show sequence"),
# 	  		"or",
# 	  		a({href => "mySeqSearch.cgi?orgId=$orgId&locusId=$locusId"}, "Check homologs")
# 	  	);
# 		print p(a({href => "domains.cgi?orgId=$orgId&locusId=$locusId"}, "See Domains"));
#     }
    

} #  end if just 1 hit

print "<br><br>";

$dbh->disconnect();
Utils::endHtml($cgi);
