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

my $gene = $dbh->selectall_arrayref("SELECT orgId, locusId, sysName, gene, desc FROM Gene WHERE orgId = ?", {Slice => {}}, $orgId);
die "Illegal orgId" if scalar(@$gene) == 0;

# my $exp = $dbh->selectall_arrayref("SELECT expName, expDesc FROM Experiment WHERE orgId = ?", {Slice => {}}, $orgId);
my $cofitness = $dbh->selectall_arrayref("SELECT * FROM Cofit WHERE orgId = ?", {Slice => {}}, $orgId);


my %genes = ();
foreach my $frow(@$gene) {
	my $locus = $frow->{locusId};
	$genes{$locus} = $frow;
}

print ("Content-Type:application/x-download\n");
print "Content-Disposition: attachment; filename=cofit_organism_$orgId.txt\n\n";

# open my $logFile, '>', "organism_$orgId.txt" or die "error trying to (over)write: $!";

# print the header row
print STDERR "writing headers...\n";
print 'orgID' ."\t". 'geneName'. "\t". 'locusID' ."\t" . 'sysName' . "\t" . "desc" . "\t" . "hitId\thitSysName\thitName\trank\tcofit";
print "\n";

# print the data row by row
print STDERR "writing data...\n";
foreach my $row (@$cofitness) {
	my $locus = $row->{locusId};
	print $row->{orgId} ."\t" . $genes{$locus}{gene} . "\t". $locus ."\t". $genes{$locus}{desc} . "\t" . $row->{hitId} . "\t" . $genes{$row->{hitId}}{sysName} . $genes{$row->{hitId}}{gene} . "\t" . $row->{rank} . "\t" . $row->{cofit};
	print "\n";
};

print STDERR "done\n";