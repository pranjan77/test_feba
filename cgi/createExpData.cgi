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
# a download, which includes all of the experimental data
# for a single organism. It is run from org.cgi. 
# 
#
# Required CGI parameters:
# orgId -- which organism to search for

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi = CGI->new;

my $orgId = $cgi->param('orgId');
$orgId = "" if !defined $orgId;

my $dbh = Utils::get_dbh();

# gather all of the data needed
print STDERR "getting data...\n";

my $gene = $dbh->selectall_arrayref("SELECT orgId, locusId, sysName, gene FROM Gene WHERE orgId = ?", {Slice => {}}, $orgId);
die "Illegal orgId" if scalar(@$gene) == 0;

my $exp = $dbh->selectall_arrayref("SELECT * FROM Experiment WHERE orgId = ?", {Slice => {}}, $orgId);



print ("Content-Type:application/x-download\n");
print "Content-Disposition: attachment; filename=exp_organism_$orgId.txt\n\n";

# open my $logFile, '>', "organism_$orgId.txt" or die "error trying to (over)write: $!";

# print the header row
print STDERR "writing headers...\n";
print "orgId" . "\t" . "expName" . "\t" . "expDesc" . "\t" . "timeZeroSet" . "\t" . "num" . "\t" . "nMapped" . "\t" . "nPastEnd" . "\t" . "nGenic" . "\t" . "nUsed" . "\t" . "gMed" . "\t" . "gMedt0" . "\t" . "gMean" . "\t" . "cor12" . "\t" . "mad12" . "\t" . "mad12c" . "\t" . "mad12c_t0" . "\t" . "opcor" . "\t" . "adjcor" . "\t" . "gccor" . "\t" . "maxFit" . "\t" . "expGroup" . "\t" . "expDescLong" . "\t" . "mutantLibrary" . "\t" . "person" . "\t" . "dateStarted" . "\t" . "steName" . "\t" . "seqindex" . "\t" . "media" . "\t" . "temperature" . "\t" . "pH" . "\t" . "vessel" . "\t" . "aerobic" . "\t" . "liquid" . "\t" . "shaking" . "\t" . "condition_1" . "\t" . "units_1" . "\t" . "concentration_1" . "\t" . "condition_2" . "\t" . "units_2" . "\t" . "concentration_2" . "\t" . "growthPlate" . "\t" . "growthWells" . "\n";
print "\n";

# print the data row by row
print STDERR "writing data...\n";
foreach my $row(@$exp) {
	print $row->{orgId} . "\t" . $row->{expName} . "\t" . $row->{expDesc} . "\t" . $row->{timeZeroSet} . "\t" . $row->{num} . "\t" . $row->{nMapped} . "\t" . $row->{nPastEnd} . "\t" . $row->{nGenic} . "\t" . $row->{nUsed} . "\t" . $row->{gMed} . "\t" . $row->{gMedt0} . "\t" . $row->{gMean} . "\t" . $row->{cor12} . "\t" . $row->{mad12} . "\t" . $row->{mad12c} . "\t" . $row->{mad12c_t0} . "\t" . $row->{opcor} . "\t" . $row->{adjcor} . "\t" . $row->{gccor} . "\t" . $row->{maxFit} . "\t" . $row->{expGroup} . "\t" . $row->{expDescLong} . "\t" . $row->{mutantLibrary} . "\t" . $row->{person} . "\t" . $row->{dateStarted} . "\t" . $row->{steName} . "\t" . $row->{seqindex} . "\t" . $row->{media} . "\t" . $row->{temperature} . "\t" . $row->{pH} . "\t" . $row->{vessel} . "\t" . $row->{aerobic} . "\t" . $row->{liquid} . "\t" . $row->{shaking} . "\t" . $row->{condition_1} . "\t" . $row->{units_1} . "\t" . $row->{concentration_1} . "\t" . $row->{condition_2} . "\t" . $row->{units_2} . "\t" . $row->{concentration_2} . "\t" . $row->{growthPlate} . "\t" . $row->{growthWells} . "\n";
}

print STDERR "done\n";
