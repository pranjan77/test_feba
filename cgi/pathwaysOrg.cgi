#!/usr/bin/perl -w
#######################################################
## pathwaysOrg.cgi -- list MetaCyc pathways for an organism
##
## Copyright (c) 2018 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId -- which organism

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;

my $orgId = $cgi->param('orgId') || die "No orgId parameter";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Unknown organism" unless $orgId eq "" || exists $orginfo->{$orgId};

print
  header,
  Utils::start_page("MetaCyc Pathways for $orginfo->{$orgId}{genome}"),
  h2("Metacyc Pathways for",
     a({ -href => "org.cgi?orgId=$orgId" }, $orginfo->{$orgId}{genome}));


# For each reaction that is in a pathway, find potential genes
# by BestHitMetacyc, EC numbers, or SEED annotations (via links to KEGG reactions)

my $rxn2path = $dbh->selectall_arrayref("SELECT rxnId, pathwayId FROM MetacycPathwayReaction");
my %rxnInPath = map { $_->[0] => 1 } @$rxn2path;
my $rxn2locus = $dbh->selectall_arrayref("SELECT rxnId, locusId FROM BestHitMetacyc WHERE orgId = ?",
                                  {}, $orgId);

my %rxnFound = (); # rxnId => 1 if found
foreach my $row (@$rxn2locus) {
  my ($rxnId,$locusId) = @$row;
  $rxnFound{$rxnId} = 1;
}

my $rxnSponteaneous = $dbh->selectcol_arrayref("SELECT rxnId FROM MetacycReaction WHERE isSpontaneous = 1");
# rxnId => 1 if spontaenous
my %rxnSponteaneous = map { $_ => 1 } @$rxnSponteaneous;

# EC mappings
my $rxn2ec = $dbh->selectall_arrayref("SELECT rxnId, ecnum FROM MetacycReactionEC");
my %ec = ();
foreach my $row (@$rxn2ec) {
  my ($rxnId,$ecnum) = @$row;
  $ec{$ecnum} = 1 if exists $rxnInPath{$rxnId};
}
my @ecs = keys %ec;
my $ecGenes = Utils::EcToGenesAll($dbh, $orgId); # ec => locusId => 1

foreach my $row (@$rxn2ec) {
  my ($rxnId,$ec) = @$row;
  $rxnFound{$rxnId} = 1 if exists $ecGenes->{$ec};
}

# SEED mappings
my $rxnFoundSEED = $dbh->selectcol_arrayref(qq{ SELECT DISTINCT rxnId FROM SEEDAnnotation
                                                 JOIN SEEDAnnotationToRoles USING (seed_desc)
                                                 JOIN SEEDRoleReaction USING (seedrole)
                                                 JOIN SEEDReaction USING (seedrxnId)
                                                 JOIN MetacycReaction USING (keggrxnId)
                                                 WHERE orgId = ? AND keggrxnId <> "" }, {}, $orgId);
foreach my $rxnId (@$rxnFoundSEED) {
  $rxnFound{$rxnId} = 1;
}

print p("Found candidate genes for ", scalar(keys %rxnFound), "reactions");

my %pathRxns = ();
foreach my $row (@$rxn2path) {
  my ($rxnId,$pathId) = @$row;
  push @{ $pathRxns{$pathId} }, $rxnId;
}

my %pathFound = ();
my %pathSize = ();
my %pathScore = ();
foreach my $pathId (keys %pathRxns) {
  my $rxns = $pathRxns{$pathId};
  $pathSize{$pathId} = scalar(@$rxns);
  my $nFound = 0;
  foreach my $rxnId (@$rxns) {
    $nFound++ if exists $rxnFound{$rxnId} || exists $rxnSponteaneous{$rxnId};
  }
  $pathFound{$pathId} = $nFound;
  $pathScore{$pathId} = Utils::MetacycPathwayScore($pathSize{$pathId}, $nFound);
}

my $rxnname = $dbh->selectall_arrayref("SELECT pathwayId, pathwayName FROM MetacycPathway");
my %pathname = map { $_->[0] => $_->[1] } @$rxnname;

my @path = sort { $pathScore{$b} <=> $pathScore{$a}
                    || $pathSize{$a} <=> $pathSize{$b}
                    || $pathname{$a} cmp $pathname{$b} }
  keys %pathScore;

my @trows = ();
push @trows, Tr(th("Pathway"), th("Steps Found"));

foreach my $pathId (@path) {
  next if $pathFound{$pathId} == 0;
  push @trows, Tr( td( a( { -href => "pathway.cgi?orgId=$orgId&pathwayId=$pathId" },
                          $pathname{$pathId} ) ),
                   td( $pathFound{$pathId} . " / " . $pathSize{$pathId} ) );
}
print table( { cellspacing => 0, cellpadding => 3 }, @trows);

$dbh->disconnect();
Utils::endHtml($cgi);
