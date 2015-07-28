#!/usr/bin/perl -w
#######################################################
## cond.cgi -- list conditions for a given group
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Key parameters: expGroup -- which experiment group to summarize

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;
use URI::Escape;

my $cgi=CGI->new;

my $expGroup = $cgi->param('expGroup') || die "No expGroup parameter";

# my $style = Utils::get_style();
my $dbh = Utils::get_dbh();

my $cond = $dbh->selectall_arrayref(qq{SELECT condition_1, COUNT(DISTINCT orgId) AS nOrg, COUNT(*) AS nExp
                                       FROM Experiment WHERE expGroup = ?
				       GROUP BY condition_1 ORDER BY condition_1 },
				    { Slice => {} },
				    $expGroup);
my $title = scalar(@$cond) > 0 ? "Overview of $expGroup experiments" : "No conditions for this group";
my $start = Utils::start_page("$title");

print
    header, $start, '<div id="ntcontent">',
  #   start_html( -title => $title, -style => { -code => $style }, -author => 'Morgan Price',
		# -meta => { 'copyright' => 'copyright 2015 UC Berkeley' }),
    h2($title);
    # div({-style => "float: right; vertical-align: top;"}, a({href => "help.cgi#specific"}, "Help"));

Utils::fail($cgi, "no conditions in this group") if @$cond == 0;

my @headings = qw{Condition Organisms Experiments &nbsp;};
my @trows = ( Tr({ -valign => 'top', -align => 'center' }, map { th($_) } \@headings) );
foreach my $row (@$cond) {
    push @trows, Tr({ -valign => 'top', -align => 'left' },
		    td([ $row->{condition_1} || "unrecorded",
			 $row->{nOrg},
			 a( { href => "exps.cgi?expGroup=" . uri_escape($expGroup)
                                  . "&condition1=" . uri_escape($row->{condition_1}) },
			    $row->{nExp} ),
			 a( { href => "orthCond.cgi?expGroup=" . uri_escape($expGroup)
                                  . "&condition1=" . uri_escape($row->{condition_1}) },
			    "compare" ) ]));
}

print table({cellspacing => 0, cellpadding => 3}, @trows);
print "<br><br>";

$dbh->disconnect();
Utils::endHtml($cgi);
