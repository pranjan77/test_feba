#!/usr/bin/perl
#######################################################
## Utils.pm -- for cgi code
##
## Copyright (c) 2015 All Right Reserved by UC Berkeley
##
## Authors:
## Wenjun Shao (wjshao@berkeley.edu) and Morgan Price
#######################################################

package Utils;
require Exporter;

use strict;
use warnings;

use DBI;
use Time::HiRes;
use Carp;

sub get_style();
sub fail($$);
sub formatFASTA($$);
sub crc64($);
sub get_color($);
sub get_dbh(); # open the sqlite db and return the resulting database handle
sub blast_db();
sub tmp_dir();
sub orginfo($);
sub get_orths($$$);

#--------------------------------------------------------

sub get_style() {
    my $style = <<"EOT";
    body {
        font-family: verdana, arial, sans-serif;
    }
    H2 {
        color: red;
        border-bottom: 1pt solid;
    }
    H4 {
        color: blue;
    }
    table {
        border: black thin solid;
    }
    th,td {
        border: lightgrey thin solid;
    }

.axis path,
.axis line {
  fill: none;
  stroke: #000;
  shape-rendering: crispEdges;
}

.dot {
  stroke: none;
}

EOT
    return $style;
}

sub fail($$) {
    my ($cgi, $notice) = @_;
    print $cgi->h3(qq(Sorry: $notice));
    print $cgi->h4(qq(<a href="myFrontPage.cgi">Go back to front page</a>));
    print $cgi->end_html;
    exit;
}

sub endHtml($) {
    my ($cgi) = @_;
    print $cgi->p(qq(<a href="myFrontPage.cgi">Go back to front page</a>));
    print $cgi->end_html;
    exit 0;
}

#--------------------------------------------------------

sub formatFASTA($$)
{
    my ($title, $sequence) = @_;
    return undef unless ($sequence);
    my $s = qq/>$title\n/;
    my @whole = $sequence =~ /.{1,60}/g;
    my $line = join "\n", @whole;
    return $s."$line\n";
}

#--------------------------------------------------------

# color code fitness value in range (-3 .. 0 .. 3) to blue => white => yellow
sub fitcolor($) {
    my ($fit) = @_;
    return " #BFBFBF" if !defined $fit;
    my $perc = ($fit + 3)/6;
    $perc = $perc > 1 ? 1 : $perc;
    $perc = $perc < 0 ? 0 : $perc;
    return fractioncolor($perc);
}

# from a fraction 0:1 to a color
sub fractioncolor($) {
#RED (255, 0, 0)
#YELLOW (255, 255, 0)
#BLUE (0, 0, 255)
#WHITE (255, 255, 255)

    my ($perc) = shift @_;
    my ($red, $blue, $green);

# I need "BLUE-WHITE-YELLOW" scheme
    if($perc >= 0.5) {
        $red = 255;
        $green = $red;
        $blue = int(255 * 2 * (1 - $perc) + 0.5);
    } else {
        $red = int(($perc/0.5)*255 + 0.5);
        $green = $red;
        $blue = 255;
    }

    my $string=sprintf ("#%2.2X%2.2X%2.2X",$red,$green,$blue);
    return ($string);
}

sub get_dbh() {
    my $database = "../cgi_data/feba.db";
    return DBI->connect("dbi:SQLite:dbname=$database","","",{ RaiseError => 1 }) || die $DBI::errostr;
}

sub blast_db() {
    return "../cgi_data/aaseqs";
}

sub tmp_dir() {
    return "../tmp";
}

sub gene_has_fitness($$$) {
    my ($dbh,$orgId,$locusId) = @_;
    die unless defined $dbh && $orgId && defined $locusId;
    my @result = $dbh->selectrow_array("SELECT expName FROM GeneFitness WHERE orgId=? AND locusId=? LIMIT 1",
				       undef, $orgId, $locusId);
    return scalar(@result) > 0 ? 1 : 0;
}

# returns a hash of orgId => attribute => value, including "genome" for the display name
sub orginfo($) {
    my ($dbh) = @_;
    my $orginfo = $dbh->selectall_hashref("SELECT * FROM Organism", "orgId");
    while (my ($orgId,$row) = each %$orginfo) {
	$row->{genome} = join(" ", $row->{genus}, $row->{species}, $row->{strain});
    }
    return $orginfo;
}

# returns a hash of expName => attribute => value
sub expinfo($$) {
    my ($dbh,$orgId) = @_;
    return $dbh->selectall_hashref("SELECT * FROM Experiment WHERE orgId = ?", "expName", {}, $orgId);
}

# Like cmp but for experiments, for ordering by group and condition/concentration, with
# expDesc as a fallback.
# Each argument should be a hash corresponding to a row in the Experiment table
sub CompareExperiments($$) {
    my ($expA,$expB) = @_;
    die unless defined $expA->{expGroup} && defined $expB->{expGroup};
    return $expA->{expGroup} cmp $expB->{expGroup}
           || lc($expA->{condition_1}) cmp lc($expB->{condition_1})
           || $expA->{concentration_1} cmp $expB->{concentration_1}
           || lc($expA->{expDesc}) cmp lc($expB->{expDesc});
}

# Should check that orgId is valid (if it is not empty) before calling.
# Allows partial-word matches in Condition_1 or Condition_2, or full word matches to expDesc or expDescLong, or exact match to Group.
# Note is not case sensitive
sub matching_exps($$$) {
    my ($dbh,$orgId,$expSpec) = @_;
    die if !defined $dbh || !defined $orgId || !defined $expSpec;
    # make the query safe to include in SQL
    $expSpec =~ s/ +$//;
    $expSpec =~ s/^ +$//;
    $expSpec =~ s/[\'\"\n\r]//g;

    my $orgClause = $orgId eq "" ? "" : qq{ AND orgId = "$orgId"};
    my $sql = qq{SELECT * from Organism JOIN Experiment USING (orgId)
             WHERE (expName = "$expSpec"
	            OR Condition_1 LIKE "$expSpec%" OR Condition_1 LIKE "% $expSpec%"  OR Condition_1 LIKE "%-$expSpec%"
		    OR Condition_2 LIKE "$expSpec%" OR Condition_2 LIKE "% $expSpec%"  OR Condition_2 LIKE "%-$expSpec%"
	     	    OR expGroup = '$expSpec'
	            OR expDesc = '$expSpec' OR expDesc LIKE "$expSpec %" OR expDesc LIKE "% $expSpec" OR expDesc LIKE "% $expSpec %"
	            OR expDescLong = '$expSpec' OR expDescLong LIKE "$expSpec %" OR expDescLong LIKE "% $expSpec" OR expDescLong LIKE "% $expSpec %")
	     $orgClause
	     ORDER BY genus, species, strain, expGroup, condition_1, concentration_1, expDesc};
    return $dbh->selectall_arrayref($sql, { Slice => {} });
}

# Returns a reference to a hash from ortholog's orgId => gene
# Each gene is a hash that includes the fields in Gene and the ratio
sub get_orths($$$) {
    my ($dbh,$orgId,$locusId) = @_;
    die unless defined $orgId && defined $locusId;
    return $dbh->selectall_hashref(qq{ SELECT Gene.*, Ortholog.ratio FROM GENE JOIN ORTHOLOG
					ON Ortholog.orgId2=Gene.orgId AND Ortholog.locusId2=Gene.locusId
					AND orgId1=? AND locusId1=? },
	                           "orgId", {}, $orgId, $locusId);
}

#END 

return 1;

