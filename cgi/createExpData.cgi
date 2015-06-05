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
# for a single organism. It is run from orgGroup.cgi. 
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

# my $style = Utils::get_style();
# print $cgi->header;

my $orgId = $cgi->param('orgId');
$orgId = "" if !defined $orgId;

my $dbh = Utils::get_dbh();

# gather all of the data needed
print STDERR "getting data...\n";

my $gene = $dbh->selectall_arrayref("SELECT * FROM Gene WHERE orgId = ?", {Slice => {}}, $orgId);
die "Illegal orgId" if scalar(@$gene) == 0;

my $exp = $dbh->selectall_arrayref("SELECT expName, expDesc FROM Experiment WHERE orgId = ?", {Slice => {}}, $orgId);
my $fitness = $dbh->selectall_arrayref("SELECT * FROM GeneFitness WHERE orgId = ?", {Slice => {}}, $orgId);
my %fit=();
foreach my $frow(@$fitness) {
	my $expName = $frow->{expName};
	my $locusId = $frow->{locusId};
	$fit{$locusId}{$expName} = $frow->{fit};
}

print ("Content-Type:application/x-download\n");
print "Content-Disposition: attachment; filename=organism_$orgId.txt\n\n";

# open my $logFile, '>', "organism_$orgId.txt" or die "error trying to (over)write: $!";

# print the header row
print STDERR "writing headers...\n";
print 'orgID' ."\t". 'geneName'. "\t". 'locusID' ."\t" . 'sysName';
	foreach my $titlerow (@$exp) {
		print "\t" . $titlerow->{expName} ." ". $titlerow->{expDesc};
	};
print "\n";

# print the data row by row
print STDERR "writing data...\n";
foreach my $row (@$gene) {
	my $locus = $row->{locusId};
	next if !exists $fit{$locus};
	print $row->{orgId} ."\t" . $row->{gene} . "\t". $row->{locusId} ."\t". $row->{sysName};
	foreach my $subrow(@$exp) {
		my $expr = $subrow->{expName};
		my $rounded = "";
		if (exists $fit{$locus}{$expr}) { 
			$rounded = sprintf("%.3f", $fit{$locus}{$expr});
		};
		print "\t" . $rounded;
	};
	print "\n";
	# };
};

print STDERR "done\n";