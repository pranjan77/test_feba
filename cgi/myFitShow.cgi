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
use IO::Handle;

use lib "../lib";
use Utils;

my $cgi=CGI->new;

my $orgSpec = $cgi->param('orgId') || "";
my $geneSpec = $cgi->param('gene');
my $help = $cgi->param('help') || "";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my $orgTitle = $cgi->param('orgId') ? "in $orginfo->{$orgSpec}{genome}" : "";
my $start = Utils::start_page("Genes matching $geneSpec $orgTitle");

$geneSpec =~ s/[ \t]*$//;
$geneSpec =~ s/^[ \t]*//;

# if no gene or locus specified
if (!defined $geneSpec || $geneSpec eq "") {
	print $cgi->header;
	print $start,'<div id="ntcontent">';
    Utils::fail($cgi, "you must enter the gene name or locus tag");
}

# match for exact gene or locus ID
my $hits = Utils::matching_exact($dbh, $orgSpec, $geneSpec);

my %used;
my $count = 0; 
# index results to avoid duplicates
foreach my $gene (@$hits) {
	next if (exists $used{$gene->{locusId}});
    $gene->{has_fitness} = Utils::gene_has_fitness($dbh,$gene->{orgId},$gene->{locusId});
    $used{$gene->{orgId}}->{$gene->{locusId}} = 1;
}


if (scalar(@$hits) == 1) {
    # just 1 hit, redirect
    my $gene = $hits->[0];
    my $orgId = $gene->{orgId};
    my $locusId = $gene->{locusId};
    my $url = "singleFit.cgi?orgId=$orgId&locusId=$locusId&showAll=0";
    $url = "singleFit.cgi?orgId=$orgId&locusId=$locusId&showAll=0&help=1" if $help == 1;

    print redirect(-url=>"$url");
    exit(0);

}

#start page
print $cgi->header, $start,'<div id="ntcontent">';

# make table for exact gene/locus matches
if (scalar(@$hits) > 0) {
	my @trows = ();
	print h3(b("Match by gene/locus for $geneSpec:"));
	push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
		  $cgi->th( [ 'Gene ID','Gene Name','Description','Genome','Fitness' ] ) );

    foreach my $gene (@$hits) {
    	if ($count < 100) {
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
			$count += 1;
		}
	}
	print $cgi->table( { cellspacing=>0, cellpadding=>3 }, @trows);
        print "\n"; # to allow flushing
}

autoflush STDOUT 1; # so preliminary results appear

# make table for description matches if total matches < 100 so far, filtering for repeats
my $descs;
$descs = Utils::matching_descs($dbh, $orgSpec, $geneSpec) if $count < 100;
@$descs = grep { !exists $used{ $_->{orgId} }{ $_->{locusId} } } @$descs;
foreach my $gene (@$descs) {
	if ($count < 100) {
	    $gene->{has_fitness} = Utils::gene_has_fitness($dbh,$gene->{orgId},$gene->{locusId});
	    $used{$gene->{orgId}}->{$gene->{locusId}} = 1;
	}
}

my @trows = ();
    if (@$descs >= 1) {
    	print h3(b("Match by description for $geneSpec:"));
    	push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
			  $cgi->th( [ 'Gene ID','Gene Name','Description','Genome','Fitness' ] ) );

	    foreach my $gene (@$descs) {
	    	if ($count < 100) {
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
				$count += 1;
	    	} 
	}
	print $cgi->table( { cellspacing=>0, cellpadding=>3 }, @trows);
        print "\n"; # to allow flush
}

# make table for domain matches if total matches < 100 so far, filtering for repeats
my $domains;
$domains = Utils::matching_domains($dbh, $orgSpec, $geneSpec) if $count < 100;

@$domains = grep { !exists $used{ $_->{orgId} }{ $_->{locusId} } } @$domains;

@trows = ();
    # $count = 0;    
	if (@$domains >= 1) {
		print h3(b("Match by domain for $geneSpec:"));
	    push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
			  $cgi->th( [ 'Gene ID','Gene Name','Description','Genome', 'Domain ID', 'Domain Name', 'Fitness' ] ) );

		foreach my $gene (@$domains) {
			if ($count < 100) {
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
}


if ($count >= 100) {
	print "<BR> Only the first 100 results are shown. <br>";
}
# if (@$domains == 1) {
#     # just 1 hit
#     my $gene = $domains->[0];
#     my $orgId = $gene->{orgId};
#     my $locusId = $gene->{locusId};

#     print redirect(-url=>"singleFit.cgi?orgId=$orgId&locusId=$locusId&showAll=0");
# }


# if no hits
if (@$hits == 0 and @$descs == 0 and @$domains == 0) {
    print $cgi->h3("No gene found for $geneSpec",
		   (exists $orginfo->{$orgSpec}{genome} ? " in " . $orginfo->{$orgSpec}{genome} : ""));

} 

print qq[<br><br>];

$dbh->disconnect();
Utils::endHtml($cgi);
