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
sub commify;

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

    my $beginC = &commify($begin);
    my $endC = &commify($end);
    print "Type $type: $typeName<BR>
   	Located on scaffold $scaffold, $strand strand, nucleotides $beginC to $endC";


    my @links = ();
    if ($gene->{locusId} =~ m/^\d+$/) {
	push @links, $cgi->a({href => "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=$gene->{locusId}"},
			     "MicrobesOnline");
    }
    if ($orgId eq "Keio" && $gene->{sysName} =~ m/^b\d+$/) {
	push @links, $cgi->a({href => "http://ecocyc.org/ECOLI/search-query?type=GENE&gname=$gene->{sysName}"}, "EcoCyc");
    }

    # Linking to NCBI is tricky -- for example, Cupriavidus basilensis 4G11 is in NCBI but the locus tags are not indexed yet.
    # Similarly PS417_05815 from P. fluorescens WCS417 -- hits P. simaie recA instead if search Protein, wierd.
    # So assume that it will work if the genome has an NCBI taxonomyId. This is a hack.
    push @links, $cgi->a({href => "http://www.ncbi.nlm.nih.gov/gene/?term=$gene->{sysName}#see-all"}, "NCBI")
        if $gene->{type} == 1 && $gene->{sysName} ne "" && $orginfo->{$orgId}{taxonomyId} ne "";

    # Use LocusXRef to try to link to IMG
    my $imgBase = "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=GeneDetail&page=geneDetail&gene_oid=";
    my ($imgId) = $dbh->selectrow_array(qq{SELECT xrefId FROM LocusXref
                                           WHERE orgId=? AND locusId=? AND xrefDb='IMG' LIMIT 1},
                                        {}, $orgId, $locusId);
    push @links, $cgi->a({href => $imgBase . $imgId}, "IMG") if defined $imgId;

    print $cgi->p("Links: " . join(", ", @links)) if (@links > 0);

    print $cgi->h3({style=>'text-align:center'},"Nearby Genes");

    print Utils::geneArrows(\@genes, $locusId, undef, undef);

    my @trows = ( Tr({-valign => "top", -align => 'center'},
                     th([ 'Locus', 'Name', 'Description', 'Strand',
                          a({title=>'Distance between the end of previous gene and beginning of current gene', color=>"#011B47"},'Distance (nt)'),
                          'Phenotype', ])));

    my $diff = "";
    my $prevrow;
    foreach my $row (@genes) {
        $diff = $row->{begin} - $prevrow->{end} if defined $prevrow;		# $prevrow->{end} = 0 if !defined $prevrow;
        my $label = $row->{gene} || $row->{sysName}; #|| $row->{locusId};
        my $label2 =  $row->{gene} || $row->{sysName} || $row->{locusId};
        my $label3 = $label2;
        $label3 =~ s/^.*_/_/ if $row->{sysName} || $row->{locusId};
        $label2 = $row->{sysName}. ": " . $label2 if $row->{sysName};
        my $bgcolor = undef;
        $bgcolor = "#FFFFFF" if $row->{locusId} eq $locusId;
        my ($phen, $tip) = Utils::gene_fit_string($dbh,$orgSpec,$row->{locusId});
        push @trows, Tr({ -valign => 'top', -align => 'left', -bgcolor=>$bgcolor},
                        td([ a({href => "geneOverview.cgi?orgId=$orgId&gene=$row->{locusId}"},$row->{sysName}||$row->{locusId}), #locus
                             a({href => "geneOverview.cgi?orgId=$orgId&gene=$row->{locusId}"},$row->{gene} || $row->{sysName}), 
                             $row->{desc}, 
                             $row->{strand},
                             a({title=>"From $prevrow->{end} to $row->{begin}"},$diff), # $row->{begin},
                             a({href => "myFitShow.cgi?orgId=$orgId&gene=$row->{locusId}", title=>$tip},$phen),
                           ]));
        $prevrow = $row;
    }
    print table({cellspacing => 0, cellpadding => 3}, @trows);
    
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

sub commify($) {
    local $_  = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}
