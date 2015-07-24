#!/usr/bin/perl -w
#######################################################
## strainTable.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId -- which organism
# and either
#	scaffoldId, begin, pos -- the range of interest (limited to 50 kb)
# or
#	locusId -- which gene to show (by default, range is exactly the width of the gene)
# Optional CGI parameters:
# expName -- which experiments to show. (Can be more than one.)
# addexp -- additional experiments (i.e. a setname or a condition)
# zoom -- in or out
# pan -- left or right

use strict;
use CGI qw(-nosticky :standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;
use StrainFitness;
sub commify($);

my $cgi=CGI->new;
my $style = Utils::get_style();
my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my $orgId = $cgi->param('orgId') || die "No orgId found";
my $genome = $orginfo->{$orgId}{genome} || die "Invalid orgId $orgId";

my @expNames = $cgi->param('expName');
my $scaffoldId = $cgi->param('scaffoldId');
my $begin = $cgi->param('begin');
my $end = $cgi->param('end');
my $locusSpec = $cgi->param('locusId');
my $locusSpecShow;

if (defined $locusSpec && $locusSpec ne "") {
    my $sysName;
    ($scaffoldId,$begin,$end,$sysName) = $dbh->selectrow_array(
        "SELECT scaffoldId,begin,end,sysName FROM Gene WHERE orgId = ? AND locusId = ?",
        {}, $orgId, $locusSpec);
    die "Unknown locus $locusSpec in $orgId" if !defined $end;
    $locusSpecShow = $sysName || $locusSpec;
    my $widen = int(1 + 0.2 * ($end-$begin+1));
    $begin -= $widen;
    $end += $widen;
} elsif (defined $scaffoldId && defined $begin && defined $end) {
    die "Invalid scaffold parameter" if $scaffoldId eq "";
    die "Invalid begin parameter" unless $begin =~ m/^-?\d+$/;
    die "Invalid end parameter" unless $end =~ m/^\d+$/;
}
my $zoom = $cgi->param('zoom');
my $initwidth = $end - $begin + 1;
if ($zoom eq "in") {
    $begin += 0.2 * $initwidth;
    $end -= 0.2 * $initwidth;
} elsif ($zoom eq "out") {
    $begin -= 0.4 * $initwidth;
    $end += 0.4 * $initwidth;
}
my $pan = $cgi->param('pan');
if ($pan eq "left") {
    
    $begin -= 0.4 * $initwidth;
    $end -= 0.4 * $initwidth;
} elsif ($pan eq "right") {
    $begin += 0.4 * $initwidth;
    $end += 0.4 * $initwidth;
}
$begin = int($begin);
$end = int($end);

$end = $begin + 1 if $begin eq $end;
unless ($begin < $end) {
    print $cgi->header;
    Utils::fail($cgi, "Invalid begin/end $begin $end")
}
my $maxWidth = 25*1000;
if ($end - $begin >= $maxWidth) {
    print $cgi->header;
    Utils::fail($cgi, "Too wide, the limit is to show a range of " . &commify($maxWidth) . " nucleotides");
}

my $addexp = $cgi->param('addexp');

# additional experiments?
if (defined $addexp && $addexp ne "") {
    my @expsNew = @{ Utils::matching_exps($dbh,$orgId,$addexp) };
    if (@expsNew == 0) {
        print header;
        Utils::fail($cgi, qq{No experiment matching "$addexp"}); # XXX
    }
    # else
    my %oldExps = map { $_ => 1 } @expNames;
    push @expNames, grep {!exists $oldExps{$_} } map { $_->{expName} } @expsNew;
}

my $expinfo = Utils::expinfo($dbh,$orgId);
foreach my $expName (@expNames) {
    die "No such experiment: $expName" unless exists $expinfo->{$expName};
}

my $begComma = &commify($begin);
my $endComma = &commify($end);
print
    header,
    Utils::start_page("Strain Fitness in $genome"),
    q{<div id="ntcontent">},
    h2("Strain Fitness in ",
       a({-href => "org.cgi?orgId=$orgId"}, "$genome"),
       defined $locusSpecShow ? "around " . a({-href => "singleFit.cgi?orgId=$orgId&locusId=$locusSpec"}, $locusSpecShow)
       : " at $scaffoldId : $begComma to $endComma"),
    start_form(-name => 'input', -method => 'GET', -action => 'strainTable.cgi'),
    hidden( -name => 'orgId', -value => $orgId, -override => 1),
    hidden( -name => 'scaffoldId', -value => $scaffoldId, -override => 1),
    hidden( -name => 'begin', -value => $begin, -override => 1),
    hidden( -name => 'end', -value => $end, -override => 1),
    join("\n", map { hidden( -name => 'expName', -value => $_, -override => 1) } @expNames),
    p({-class => "inline"},
      "Add experiment(s): ",
      textfield(-name => 'addexp', -default => "", -override => 1, -size => 20, -maxLength => 100),
      "Zoom", submit('zoom','in'), submit('zoom','out'),
      "Pan:", submit('pan','left'), submit('pan','right')),
    end_form,
    p(small("Only strains with sufficient reads to estimate fitness are shown, but the strain fitness values are still rather noisy. Strains near the edge of a gene are not shown as being associated with that gene (the Gene column will be empty)."));

# should I add zoom in/out and pan left/right buttons??
my $rows = StrainFitness::GetStrainFitness("../cgi_data", $dbh, $orgId, $scaffoldId, $begin, $end);

if (@$rows == 0) {
    print "No fitness data for strains in range " . commify($begin) . " to " . commify($end) . "\n";
}
my @trows = (); # the table
# header row
my @headings = qw{Position Strand Gene};
push @headings, a({-title => "Fractional position within gene"}, "fraction");
foreach my $expName (@expNames) {
    push @headings, a({-href => "exp.cgi?orgId=$orgId&expName=$expName", -title => $expName},
                      $expinfo->{$expName}{expDesc});
}
push @trows, Tr({-align => 'CENTER', -valign=>'TOP'}, th(\@headings));

# leave out gene if not a used strain
foreach my $row (@$rows) {
    $row->{locusId} = "" unless $row->{used} eq "TRUE";
}
my %locusIds = map { $_->{locusId} => 1 } @$rows;
my %genes = (); # locusId => row
foreach my $locusId (keys %locusIds) {
    next if $locusId eq "";
    my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId = ? AND locusId = ?",
                                       {}, $orgId, $locusId);
    die "Unknown locusId $locusId" unless exists $gene->{locusId};
    $genes{$locusId} = $gene;
}

foreach my $row (@$rows) {
    my $locusId = $row->{locusId};
    my $locusShow = "";
    my $gene = undef;
    if ($locusId ne "") {
        $gene = $genes{$locusId};
        $locusShow = $gene->{sysName} || $gene->{locusId};
    }
    my @row = ( a({-title => "barcode $row->{barcode}"}, &commify($row->{pos})),
                $row->{strand},
                $locusId eq "" ? "" : a({-title => $gene->{desc}, -href => "singleFit.cgi?orgId=$orgId&locusId=$locusId"},
                                        $locusShow),
                $locusId eq "" ? "" : sprintf("%.2f",
                                              ($row->{pos} - $gene->{begin}) / ($gene->{end} - $gene->{begin} + 1))
        );
    @row = map { td($_) } @row;
    foreach my $expName (@expNames) {
        my $fit = $row->{ $expName };
        push @row, td( { -bgcolor => Utils::fitcolor($fit) }, sprintf("%.1f",$fit));
    }
    push @trows, Tr({-align => 'CENTER', -valign=>'TOP'}, @row);
}

if (scalar(@expNames) > 0) {
    # add row for removing items
    my @row = (td(""),td(""),td(""),td(""));
    my $baseURL = "strainTable.cgi?orgId=$orgId&scaffoldId=$scaffoldId&begin=$begin&end=$end";
    foreach my $expName (@expNames) {
        my @otherExps = grep { $_ ne $expName } @expNames;
        my @otherExpSpec = map { "expName=$_" } @otherExps;
        push @row, td( a({ -title => "$expName : $expinfo->{$expName}{expDesc}",
                           -href => join("&", $baseURL, @otherExpSpec) },
                         "remove") );
    }
    push @trows, Tr({-align => 'CENTER', -valign=>'TOP'}, @row);
}
    
print small(table({ cellspacing => 0, cellpadding => 3, }, @trows));

$dbh->disconnect();
Utils::endHtml($cgi);

# i.e., 1234567 => 1,234,567
sub commify($) {
    local $_  = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}
