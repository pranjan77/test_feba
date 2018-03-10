#!/usr/bin/perl -w

#######################################################
## org.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo and Morgan Price
#######################################################
#
# This page shows on overview  for a single 
# organism, with an overview of experiment groups nad
# links to specific phenotypes and downloads.
#
# Required CGI parameters:
# orgId -- which organism

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use IO::Handle; # for autoflush

use lib "../lib";
use Utils;
use URI::Escape;

my $cgi=CGI->new;

my $orgId = $cgi->param('orgId') || die "No orgId parameter";
my $expGroup = $cgi->param('expGroup');
my $condition_1 = $cgi->param('condition_1'); 

my $style = Utils::get_style();
my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
Utils::fail($cgi, "Unknown organism: $orgId") unless $orgId eq "" || exists $orginfo->{$orgId};

# main table
# gather data and slice it into an array of hashes
my $cond = $dbh->selectall_arrayref(qq{SELECT orgId, expGroup, COUNT(DISTINCT condition_1) AS nCond, COUNT(*) as nExp FROM Experiment WHERE orgId=? GROUP BY expGroup ORDER BY expGroup; },
    { Slice => {} },
    $orgId);


# gather number of genes and data
my $numGenes = $dbh->selectrow_array(qq{SELECT COUNT (DISTINCT locusId) FROM Gene WHERE orgId = ?;}, undef, $orgId);
my $numData = $dbh->selectrow_array(qq{SELECT COUNT (DISTINCT locusId) FROM GeneFitness WHERE orgId = ?;}, undef, $orgId);


# write the title
my $title = scalar(@$cond) > 0 ? "$orginfo->{$orgId}{division}: $orginfo->{$orgId}{genome}" : "No experiments for this organism.";
my $start = Utils::start_page("$title");

print
    header,
    $start, '<div id="ntcontent">',
    h2($title);

#exit if no results
Utils::fail($cgi, "No experiments for this organism.") if @$cond == 0;

#create table
my @headings = qw{Group Conditions Experiments};
my @trows = ( Tr({ -valign => 'top', -align => 'center' }, map { th($_) } \@headings) );
foreach my $row (@$cond) {
    push @trows, Tr({ -valign => 'top', -align => 'left' },
                    td([
                         $row->{expGroup} || "unrecorded", #group
                         $row->{nCond}, #conditions
                         a( { href => "exps.cgi?orgId=$orgId&expGroup=" . uri_escape($row->{expGroup}) },
                            $row->{nExp} ), #experiments
                       ]));
}

#display taxonomic information and link
if ((defined $orginfo->{$orgId}{taxonomyId}) && ($orginfo->{$orgId}{taxonomyId} ne "")) {
	print $cgi->p("NCBI taxonomy id:", $cgi->a({href => "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=$orginfo->{$orgId}{taxonomyId}"}, $orginfo->{$orgId}{taxonomyId}));
}

# not sure this plus use of newlines below is actually working to make top of page load faster
autoflush STDOUT 1; # show preliminary results

print
    p("Fitness data for $numData genes of $numGenes genes."),
    h3("Fitness Experiments"),
    table({cellspacing => 0, style => "margin-left: 50px;", cellpadding => 3}, @trows),
    "\n";

# Overview of specific phenotypes
my $specN = $dbh->selectall_hashref(
    qq{ SELECT expGroup, count(DISTINCT condition_1) AS nCond, count(DISTINCT locusId) AS nGene
        FROM SpecificPhenotype JOIN Experiment USING (orgId,expName)
        WHERE orgId = ?
        GROUP BY expGroup },
    "expGroup", { Slice => {} }, $orgId);
if (scalar(keys %$specN) > 0) {
    my ($nSpecGene) = $dbh->selectrow_array(
        "SELECT COUNT(DISTINCT locusId) FROM SpecificPhenotype WHERE orgId=?",
        {}, $orgId);
    my ($nConsSpec) = $dbh->selectrow_array(
        "SELECT COUNT(DISTINCT locusId) FROM SpecOG WHERE orgId=? AND nInOG > 1",
        {}, $orgId);
    print
        h3("Specific Phenotypes"),
        p("$nSpecGene genes with specific phenotypes,  and $nConsSpec with conserved-specific phenotypes.");
    
    my @trows = ();
    push @trows, $cgi->Tr({-valign => "top"}, $cgi->th([ 'Group', '# Conditions', '# Genes' ]));
    foreach my $expGroup (sort keys %$specN) {
        my $row = $specN->{$expGroup};
        push @trows, $cgi->Tr({-valign => "top"}, $cgi->td([ 
                                  a({ href => "spec.cgi?orgId=$orgId&expGroup=".uri_escape($expGroup) }, $expGroup),
                                  $row->{nCond},
                                  $row->{nGene} ]));
    }
    print table({cellspacing => 0, style => "margin-left: 50px;", cellpadding => 3}, @trows);
}
print "\n";

my $reanno = $dbh->selectall_arrayref("SELECT * from Reannotation WHERE orgId=?",
                                      {}, $orgId);
if (scalar(@$reanno) > 0) {
    print
        h3("Updated Annotations"),
        a({-href => "orgReanno.cgi?orgId=$orgId"},
          "Updated annotations for " . scalar(@$reanno) . " genes");
}

print
  h3("Metabolic Maps from KEGG"),
  p("See",
    a( {href => "keggmap.cgi?mapId=01100&orgId=$orgId"}, "overview map"),
    "or",
    a( {href => "keggmaplist.cgi?orgId=$orgId"}, "list of maps")."."),
  h3("SEED Subsystems"),
  p(a({href => "seedsubsystemsOrg.cgi?orgId=$orgId"}, "See list")),
  h3("MetaCyc Pathways"),
  p(a{href => "pathwaysOrg.cgi?orgId=$orgId"}, "See list"),
  "\n";
      
# No COUNT DISTINCT in sqlite3 so use GROUP BY as workaround
my $loci = $dbh->selectcol_arrayref("SELECT locusId FROM ConservedCofit WHERE orgId = ? GROUP BY locusId",
                                          {}, $orgId);
my $nConsCofit = scalar(@$loci);
#file of experimental data - generated by createExpData.cgi
# my $data = `./createExpData.cgi orgId=$orgId`;
print h3("Downloads"),
    ul(
        li(a({href => "createFitData.cgi?orgId=$orgId"}, "Fitness values"), "(tab-delimited, slow)"), 
        li(a({href => "createCofitData.cgi?orgId=$orgId"}, "Top cofitness for each gene"), "(tab-delimited)"),
        ul(li("includes $nConsCofit genes with conserved cofitness")),
        li(a({href => "spec.cgi?orgId=$orgId&download=1"}, "Specific phenotypes"), "(tab-delimited)"),
        li(a({href => "createExpData.cgi?orgId=$orgId"}, "Experiment meta-data"), "(tab-delimited)"),
        li(a({href => "orgSeqs.cgi?orgId=$orgId"}, "Protein sequences"), "(fasta)"),
        li(a({href => "orgGenes.cgi?orgId=$orgId"}, "List of genes"), "(tab-delimited)")
    );

$dbh->disconnect();
Utils::endHtml($cgi);
