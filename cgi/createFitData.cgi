#!/usr/bin/perl -w

#######################################################
## createExpData.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo and Morgan Price
#######################################################
#
# This page generates output to STDOUT and presents it as
# a download, which includes all of the experimental FITNESS data
# for a single organism. It is run from org.cgi. 
# 
#
# Required CGI parameters:
# orgId -- which organism to search for
#
# Optional CGI parameters:
# expName -- which experiments to include. If not specified, include all experiments.

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi = CGI->new;

my $orgId = $cgi->param('orgId');
$orgId = "" if !defined $orgId;
my @expNames = $cgi->param('expName');

my $dbh = Utils::get_dbh();

# gather all of the data needed

my $gene = $dbh->selectall_arrayref("SELECT orgId, locusId, sysName, gene, desc FROM Gene WHERE orgId = ?", {Slice => {}}, $orgId);
die "Illegal orgId" if scalar(@$gene) == 0;

my $exp = $dbh->selectall_arrayref("SELECT expName, expDesc FROM Experiment WHERE orgId = ?", {Slice => {}}, $orgId);
my %exp = map { $_->{expName} => $_ } @$exp;
foreach my $expName (@expNames) {
  die "Invalid expName $expName in org $orgId\n"
    unless exists $exp{$expName};
}

my $fitness = $dbh->selectall_arrayref("SELECT locusId, expName, fit FROM GeneFitness WHERE orgId = ?", {}, $orgId);

my %fit=();
foreach my $frow(@$fitness) {
	my ($locusId, $expName, $fit) = @$frow;
	$fit{$locusId}{$expName} = $fit;
}

my $outfile = @expNames == 0 ? "fit_organism_${orgId}.txt" : "fit_organism_${orgId}_selected.txt";
print ("Content-Type:application/x-download\n");
print "Content-Disposition: attachment; filename=$outfile\n\n";

# open my $logFile, '>', "organism_$orgId.txt" or die "error trying to (over)write: $!";

@expNames = map { $_->{expName} } @$exp
  if @expNames == 0;
my @expsUse = map $exp{$_}, @expNames;

# print the header row
print join("\t",
           qw{orgId locusId sysName geneName desc},
           map { $_->{expName} . " " . $_->{expDesc} } @expsUse) . "\n";

# print the data row by row
foreach my $row (@$gene) {
  my $locusId = $row->{locusId};
  next if !exists $fit{$locusId};
  my @out = ($row->{orgId}, $locusId, $row->{sysName}, $row->{gene}, $row->{desc});
  push @out, map { sprintf("%.3f", $fit{$locusId}{$_}) } @expNames;
  print join("\t", @out) . "\n";
}
