#!/usr/bin/perl -w

#######################################################
## orgSeqs.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# This page generates an amino acid fasta sequence file
# for an organism.
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
die "No orgId" if !defined $orgId || $orgId eq "";

my $dbh = Utils::get_dbh();

# gather all of the data needed
my $genes = $dbh->selectall_hashref("SELECT * FROM Gene WHERE orgId = ?",
                                    "locusId", # field to index hash by
                                    {Slice => {}}, $orgId);
die "Illegal orgId" if scalar(keys %$genes) == 0;

print "Content-Type:application/x-download\n";
print "Content-Disposition: attachment; filename=organism_$orgId.faa\n\n";

my $aaseqs = Utils::blast_db();
open(FAA, "<", $aaseqs) || die "Cannot read $aaseqs";
my $printing = 0;
while(<FAA>) {
    chomp;
    if (m/^>(.*)$/) {
        my ($orgLocus,$locusId) = split /:/, $1;
        if ($orgLocus eq $orgId) {
            $printing = 1;
            my $gene = $genes->{$locusId};
            die "Unrecognized gene $locusId in $orgId" if !defined $gene;
            print ">$orgId:$locusId $gene->{sysName} $gene->{desc}\n";
        } else {
            $printing = 0;
        }
    } else {
        print $_."\n" if $printing;
    }
}
close(FAA) || die "Error reading $aaseqs";
