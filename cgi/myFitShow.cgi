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
# my $style = Utils::get_style();

my $orgSpec = $cgi->param('orgId') || "";
my $geneSpec = $cgi->param('gene');
my $showAll = $cgi->param('showAll') ? 1 : 0;
my $start = Utils::start_page("Fitness for $geneSpec ($orgSpec)");

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

# my $query = qq{SELECT orgId, locusId, sysName, desc, gene, type FROM Gene
		# WHERE ( locusId = ? OR sysName = ? OR upper(gene) = upper(?) )};
my $hits = Utils::matching_exact($dbh, $orgSpec, $geneSpec);


#sort by exact name/ID match, descriptor match, domain match - hisg
# if ($orgSpec) {
#     $query .= " AND orgId = ?";
#     $hits = $dbh->selectall_arrayref($query, { Slice => {} }, $geneSpec, $geneSpec, $geneSpec, $orgSpec);
# } else {
#     $hits = $dbh->selectall_arrayref($query, { Slice => {} }, $geneSpec, $geneSpec, $geneSpec);
# }


# and add another column for whether there is fitness data + duplicates + > 100
my %used;
my $count = 0; 
foreach my $gene (@$hits) {
	next if (exists $used{$gene->{locusId}}) or ($count > 100);
    $gene->{has_fitness} = Utils::gene_has_fitness($dbh,$gene->{orgId},$gene->{locusId});
    $used{$gene->{orgId}}->{$gene->{locusId}} = 1;
	$count += 1;
}


if (@$hits == 1) {
    # just 1 hit
    my $gene = $hits->[0];
    my $orgId = $gene->{orgId};
    my $locusId = $gene->{locusId};

    print redirect(-url=>"singleFit.cgi?orgId=$orgId&locusId=$locusId&showAll=0");

}

my $descs = Utils::matching_descs($dbh, $orgSpec, $geneSpec);
@$descs = grep { !exists $used{ $_->{orgId} }{ $_->{locusId} } } @$descs;
foreach my $gene (@$descs) {
	next if (exists $used{$gene->{locusId}}) or ($count > 100);
    $gene->{has_fitness} = Utils::gene_has_fitness($dbh,$gene->{orgId},$gene->{locusId});
    $used{$gene->{orgId}}->{$gene->{locusId}} = 1;
	$count += 1;
}


my $domains = Utils::matching_domains($dbh, $orgSpec, $geneSpec);
@$domains = grep { !exists $used{ $_->{orgId} }{ $_->{locusId} } } @$domains;



if (@$domains == 1) {
    # just 1 hit
    my $gene = $domains->[0];
    my $orgId = $gene->{orgId};
    my $locusId = $gene->{locusId};

    print redirect(-url=>"singleFit.cgi?orgId=$orgId&locusId=$locusId&showAll=0");
}

if (@$hits == 0 and @$descs == 0 and @$domains == 0) {
	print $cgi->header;
 #    print $cgi->start_html(
	# -title =>"Gene Search",
	# -style => {-code => $style},
	# -author=>'wjshaoATberkeley.edu',
	# -meta=>{'copyright'=>'copyright 2015 UC Berkeley'} );
	print $start,'<div id="ntcontent">';
    print $cgi->h3("No gene found for $geneSpec",
		   (exists $orginfo->{$orgSpec}{genome} ? " in " . $orginfo->{$orgSpec}{genome} : ""));

} elsif (@$hits > 1 or @$descs > 1 or @$domains > 1) {
	print $cgi->header;
 #    print $cgi->start_html(
	# -title =>"Gene Search",
	# -style => {-code => $style},
	# -author=>'wjshaoATberkeley.edu',
	# -meta=>{'copyright'=>'copyright 2015 UC Berkeley'} ),
	# div({-style => "float: right; vertical-align: top;"},
	#     a({href => "help.cgi"}, "Help")),
	

	print $start, qq[<div id="ntcontent">];

	my @trows = ();
    if (@$hits > 1) {
    	print h3(b("Match by gene/locus for $geneSpec:"));
    	push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
			  $cgi->th( [ 'Gene ID','Gene Name','Description','Genome','Fitness' ] ) );

	    foreach my $gene (@$hits) {
		my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $gene->{orgId}, $gene->{locusId});
		my @trow = map $cgi->td($_), (
			a( {href => "geneOverview.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}"}, $gene->{sysName}||$gene->{locusId}), 
			# $gene->{sysName}, 
			$gene->{gene}, 
			$gene->{desc},
					      $cgi->a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId}}, "$orginfo->{$gene->{orgId}}->{genome}"),
					      a( {href => "myFitShow.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}", title => $fittitle, },
						 $fitstring));
		push @trows, $cgi->Tr(@trow);
	    }
	}
	print $cgi->table( { cellspacing=>0, cellpadding=>3 }, @trows);
    # print qq[<br><br>];

    my @trows = ();
    if (@$descs > 1) {
    	print h3(b("Match by description for $geneSpec:"));
    	push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
			  $cgi->th( [ 'Gene ID','Gene Name','Description','Genome','Fitness' ] ) );

	    foreach my $gene (@$descs) {
		my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $gene->{orgId}, $gene->{locusId});
		my @trow = map $cgi->td($_), (
			a( {href => "geneOverview.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}"}, $gene->{sysName}||$gene->{locusId}), 
			# $gene->{sysName}, 
			$gene->{gene}, 
			$gene->{desc},
					      $cgi->a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId}}, "$orginfo->{$gene->{orgId}}->{genome}"),
					      a( {href => "myFitShow.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}", title => $fittitle, },
						 $fitstring));
		push @trows, $cgi->Tr(@trow);
	    }
	}
	print $cgi->table( { cellspacing=>0, cellpadding=>3 }, @trows);
    # print qq[<br><br>];

	my @trows = ();
    $count = 0;    
	if (@$domains > 1) {
		print h3(b("Match by domain for $geneSpec:"));
	    push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
			  $cgi->th( [ 'Gene ID','Gene Name','Description','Genome', 'Domain ID', 'Domain Name', 'Fitness' ] ) );

		foreach my $gene (@$domains) {
			next if (exists $used{$gene->{locusId}}) or ($count > 100);
		my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $gene->{orgId}, $gene->{locusId});
		my @trow = map $cgi->td($_), (
			a( {href => "geneOverview.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}"}, $gene->{sysName}||$gene->{locusId}), 
			# $gene->{sysName}, 
			$gene->{gene}, 
			$gene->{desc},
	      	$cgi->a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId}}, "$orginfo->{$gene->{orgId}}->{genome}"),
			$gene->{domainId},
			$gene->{domainName},
	      	a( {href => "myFitShow.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}", title => $fittitle, },
						 $fitstring));
		push @trows, $cgi->Tr(@trow);
		$used{$gene->{locusId}} = 1;
		$count += 1;
		}
	}
    
    print $cgi->table( { cellspacing=>0, cellpadding=>3 }, @trows);
    print qq[<br><br>];
    # PRINT TABLE

} else {
    # just 1 hit
    my $gene = $hits->[0];
    my $orgId = $gene->{orgId};
    my $locusId = $gene->{locusId};

    print redirect(-url=>"singleFit.cgi?orgId=$orgId&locusId=$locusId&showAll=0");
	
#     if ($hits->[0]{has_fitness} == 0) {
# 		my $idShow = $gene->{sysName} || $gene->{locusId};
# 		my $start = Utils::start_page("Fitness data for $idShow in $orginfo->{$orgId}{genome}");
# 		my $tabs = Utils::tabsGene($dbh,$cgi,$orgId,$locusId,0,$gene->{type},"fitness");
		
# 		print 
# 			# $cgi->div({-id=>"ntcontent"},
# 			$start, $tabs, 
# 			$cgi->h3("$idShow $gene->{gene}: $gene->{desc} in " . $orginfo->{$gene->{orgId}}{genome}), 
# 			$cgi->p("Sorry, no fitness data for $idShow");
#     } else {
# 	# show fitness data for gene
# 		my @fit = @{ $dbh->selectall_arrayref("SELECT expName,fit,t from GeneFitness where orgId=? AND locusId=?",
# 					      { Slice => {} },
# 					      $orgId, $locusId) };
# 		my $nTotValues = scalar(@fit);
# 		die "Unreachable" if $nTotValues == 0;
# 		my $limitRows = $showAll ? $nTotValues : 20;
# 		my $minAbsFit;
# 	if ($nTotValues > $limitRows) {
# 	    # subset the rows
# 	    @fit = sort { abs($b->{fit}) <=> abs($a->{fit}) } @fit;
# 	    @fit = @fit[0..($limitRows-1)];
# 	    $minAbsFit = abs($fit[$#fit]{fit});
# 	}

# 	# and get metadata about experiments
# 	my $expinfo = Utils::expinfo($dbh,$orgId);

# 	if ($showAll) {
# 	    @fit = sort { Utils::CompareExperiments($expinfo->{$a->{expName}}, $expinfo->{$b->{expName}}) } @fit;
# 	} else {
# 	    @fit = sort { $a->{fit} <=> $b->{fit} } @fit;
# 	}

# 	my $idShow = $gene->{sysName} || $gene->{locusId};
# 	my $title = "Fitness data for $idShow in $orginfo->{$orgId}{genome}";
# 	#my $tabs = Utils::tabsGene($dbh,$cgi,$orgSpec,$geneSpec,$showAll,$gene->{type},"fit");

# 	print
# 	   #  start_html( -title => $title, -style => {-code => $style}, -author=>'wjshaoATberkeley.edu',
# 			 # -meta=>{'copyright'=>'copyright 2015 UC Berkeley'} ),
# 		$start, #$tabs,
# 		'<div id="ntcontent">',
# 	    h2("Fitness data for $idShow in " . $cgi->a({href => "org.cgi?orgId=$orgId"}, "$orginfo->{$orgId}{genome}")),
# 	    qq[<div style="position: relative;"><div class="floatbox">],
# 	    start_form(-name => 'input', -method => 'GET', -action => 'genesFit.cgi'),
# 	    # "<P><center><b>Compare</b></center> 
# 	    "<br>Add gene: ",
# 	    hidden( -name => 'orgId', -value => $orgId, -override => 1 ),
# 	    hidden( -name => 'showAll', -value => $showAll, -override => 1  ),
# 	    hidden( -name => 'locusId', -value => $locusId, -override => 1 ),
# 	    textfield( -name => 'addgene', -default => "", -override => 1, -size => 20, -maxLength => 100 ),
# 	    end_form,
# 	    qq[</P></div></div>],
# 	 #    div({-style => "float: right; vertical-align: top;"},
# 		# a({href => "help.cgi#fitness"}, "Help")),
# 	    h3("$idShow $gene->{gene}: $gene->{desc}");

# 	    # Option to add a gene (links to genesFit.cgi)
# 	# print
# 	#     qq[<div style="position: relative;"><div class="floatbox">],
# 	#     start_form(-name => 'input', -method => 'GET', -action => 'genesFit.cgi'),
# 	#     # "<P><center><b>Compare</b></center> 
# 	#     "<br>Add gene: ",
# 	#     hidden( -name => 'orgId', -value => $orgId, -override => 1 ),
# 	#     hidden( -name => 'showAll', -value => $showAll, -override => 1  ),
# 	#     hidden( -name => 'locusId', -value => $locusId, -override => 1 ),
# 	#     textfield( -name => 'addgene', -default => "", -override => 1, -size => 20, -maxLength => 100 ),
# 	#     end_form,
# 	#     qq[</P></div></div>];



# 	if ($showAll) {
# 	    print $cgi->p("All " . scalar(@fit) . " fitness values, sorted by group and condition");
# 	} else {
# 	    if (defined $minAbsFit) {
# 		$minAbsFit = sprintf("%.1f", $minAbsFit);
# 		print $cgi->p("Top $limitRows experiments with the strongest phenotypes (|fitness| &ge; $minAbsFit)");
# 	    } else {
# 		print $cgi->p("All " . scalar(@fit) . " fitness values, sorted by value");
# 	    }
# 	}

# 	my $option = "Or see ";
# 	if ($showAll == 0) {
#         my $showAllDest = qq(myFitShow.cgi?orgId=$orgId&gene=$locusId&showAll=1);
#         $option .= qq[<a href="$showAllDest">all fitness data</a>];
#         # print $cgi->p(qq(<a href=$showAllDest>All fitness data</a>));
#     } else {
#         my $showFewDest = qq(myFitShow.cgi?orgId=$orgId&gene=$locusId&showAll=0);
#         $option .= qq[<a href="$showFewDest">strongest phenotypes</a>];
#     }
#     print $option . '.';
	    
# 	my @out = (); # specifiers for HTML rows for each fitness value
# 	my $lastGroup = ""; # only enter the first time
# 	foreach my $fitrow (@fit) {
# 	    my $expName = $fitrow->{expName};
# 	    my $exp = $expinfo->{$expName};
# 	    my $group = $exp->{expGroup};
# 	    push @out, join(" ",
# 			    td($group eq $lastGroup ? "" : $group),
# 			    td(a({ 
# 			    	# style => "color:rgb(0,0,0)",
# 					       title => "$expName: $exp->{expDescLong}",
# 					       href => "exp.cgi?orgId=$orgId&expName=$expName" },
# 					       $exp->{expDesc})),
# 			    td( { -bgcolor => Utils::fitcolor($fitrow->{fit}) },
# 				      sprintf("%.1f", $fitrow->{fit}) ),
# 			    td( sprintf("%.1f", $fitrow->{t}) ),
# 			    td(a({ 
# 			    	# style => "color:rgb(0,0,0)",
# 				   title => "Compare to data from similar experiments or orthologs",
# 				   href => "orthFit.cgi?orgId=$orgId&locusId=$locusId&expGroup=$exp->{expGroup}&condition1=$exp->{condition_1}" }),
# 				 "compare") );
# 	    $lastGroup = $group if $showAll;
# 	}
# 	my $relsize = $showAll ? "70%" : "100%";
# 	print $cgi->table(
# 	    { cellspacing => 0, cellpadding => 3, },
# 	    $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
# 		     $cgi->th( [ 'group', 'condition','fitness','t score', '&nbsp;' ] ) ),
#             $cgi->Tr({-align=>'left',-valign=>'top',-style=>"font-size: $relsize"}, \@out ) );

	

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
    

 #    my @links = ();
 #    if ($gene->{locusId} =~ m/^\d+$/) {
	# push @links, $cgi->a({href => "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=$gene->{locusId}"},
	# 		     "MicrobesOnline");
 #    }
 #    if ($orgId eq "Keio" && $gene->{sysName} =~ m/^b\d+$/) {
	# push @links, $cgi->a({href => "http://ecocyc.org/ECOLI/search-query?type=GENE&gname=$gene->{sysName}"}, "EcoCyc");
 #    }
 #    print $cgi->p("Links: " . join(", ", @links)) if (@links > 0);
    print "<br><br>";
} #  end if just 1 hit



$dbh->disconnect();
Utils::endHtml($cgi);
