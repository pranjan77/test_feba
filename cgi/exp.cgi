#!/usr/bin/perl -w
#######################################################
## exp.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required parameters: orgId and expName, for organism and which experiment

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();
print $cgi->header;

my $orgId = $cgi->param('orgId') || die "no species identifier\n";
die "Invalid species name!" unless $orgId =~ m/^[A-Za-z0-9_]*$/;
my $expName = $cgi->param('expName') || die "no experiment identifier\n";
die "Invalid experiment name!" unless $expName =~ m/^[A-Za-z0-9_]*$/;

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "No such orgId: $orgId" unless exists $orginfo->{$orgId};

my $exp = $dbh->selectrow_hashref("SELECT * FROM Experiment WHERE orgId = ? AND expName = ?", {}, $orgId, $expName);
die "No such experiment: $expName" unless defined $exp->{expName};

print $cgi->start_html(
    -title =>"Experiment $expName for $orginfo->{$orgId}{genome}",
    -style => {-code => $style},
    -author=>'morgannprice@yahoo.com',
    -meta=>{'copyright'=>'copyright 2015 UC Berkeley'},
);

my @trows = ();
foreach my $key (sort keys %$exp) {
    my $val = $exp->{$key};
    push @trows, $cgi->Tr({align => 'left', valign => 'top'}, $cgi->td($key), $cgi->td($val));
}
print $cgi->table({cellpadding => 3, cellspacing => 0}, @trows);
$dbh->disconnect();
Utils::endHtml($cgi);
