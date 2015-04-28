#!/usr/bin/perl -w
#######################################################
## exps.cgi -- search for experiments, or list all in one organism
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Key parameters: orgId and query, for which organism and to look up experiments
# At least one must be meaningful (present and not empty)

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();
print $cgi->header;

my $orgId = $cgi->param('orgId');
$orgId = "" if !defined $orgId;

my $dbh = Utils::get_dbh();

# Make sure both parameters are safe
my $orginfo = Utils::orginfo($dbh);
Utils::fail($cgi, "Unknown organism: $orgId") unless $orgId eq "" || exists $orginfo->{$orgId};
my $expSpec = $cgi->param('query');
$expSpec = "" if !defined $expSpec;

Utils::fail($cgi, "Must specify organism or query by experiment") if $orgId eq "" && $expSpec eq "";
# make the query safe to include in SQL
$expSpec =~ s/ +$//;
$expSpec =~ s/^ +$//;
$expSpec =~ s/[\'\"\n\r]//g;
# allow partial-word matches in Condition_1 or Condition_2, or full word matches to expDesc or expDescLong, or match to Group
# note is not case sensitive
my $orgClause = $orgId eq "" ? "" : qq{ AND orgId = "$orgId"};
my $sql = qq{SELECT * from Organism JOIN Experiment USING (orgId)
             WHERE (expName = "$expSpec"
	            OR Condition_1 LIKE "$expSpec%" OR Condition_1 LIKE "% $expSpec%"  OR Condition_1 LIKE "%-$expSpec%"
		    OR Condition_2 LIKE "$expSpec%" OR Condition_2 LIKE "% $expSpec%"  OR Condition_2 LIKE "%-$expSpec%"
	     	    OR expGroup = "$expSpec"
	            OR expDesc = "$expSpec" OR expDesc LIKE "$expSpec %" OR expDesc LIKE "% $expSpec" OR expDesc LIKE "% $expSpec %"
	            OR expDescLong = "$expSpec" OR expDescLong LIKE "$expSpec %" OR expDescLong LIKE "% $expSpec" OR expDescLong LIKE "% $expSpec %")
	     $orgClause
	     ORDER BY genus, species, strain, expGroup, condition_1, concentration_1, expDesc};
my $exps = $dbh->selectall_arrayref($sql, { Slice => {} });

print $cgi->start_html(
    -title =>"Experiments for $expSpec",
    -style => {-code => $style},
    -author=>'Morgan Price',
    -meta=>{'copyright'=>'copyright 2015 UC Berkeley'},
);

if (@$exps == 0) {
   print $cgi->h3(qq{No experiment found matching "$expSpec"});
} else {
  my $heading = "Experiments";
  $heading .= " in $orginfo->{$orgId}{genome}" if $orgId ne "";
  $heading .= qq{ matching "$expSpec"} if $expSpec ne "";
  print $cgi->h3($heading);
  my @trows = ();
  push @trows, $cgi->Tr({-valign => "top"}, $cgi->th([ 'organism', 'name', 'group', 'condition', 'description' ]));
  foreach my $row (@$exps) {
      push @trows, $cgi->Tr({-valign => "top"},
             $cgi->td([ $orginfo->{$row->{orgId}}{genome},
	              $cgi->a({href => "exp.cgi?orgId=$row->{orgId}&expName=$row->{expName}"}, $row->{expName}),
		      $row->{expGroup}, $row->{condition_1}, $row->{expDesc} ]));
  }
  print $cgi->table({cellspacing => 0, cellpadding => 3}, @trows);
}

$dbh->disconnect();
Utils::endHtml($cgi);
