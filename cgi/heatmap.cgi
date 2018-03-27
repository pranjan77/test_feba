#!/usr/bin/perl -w
#######################################################
## heatmap.cgi
##
## Copyright (c) 2018 University of California
##
## Author: Morgan Price
#######################################################

# Required CGI parameters:
# orgId -- which organism
#
# Common CGI parameters:
# r -- 1 or more row specifiers, in order from top to bottom. Each is either a locusId
#	or a label specifier like _l1 or similar label.
# rt.X for various row specifiers (i.e., rtext,_l1) -- text to show instead of
#	the usual gene annotation. (The original annotation will be in the popup text.)
# c -- 1 or more column specifiers (experiments), from left to right
# view -- set to 1 if in view-only mode # not implemented
# 
# CGI parameters used to update the view:
# addrow -- search for matching gene(s) in this organism [multiple space-delimited terms]
# addrowAt -- -1 or 1 for top or bottom
# addcol -- ditto but may add multiple experiments
# addcolAt -- -1 or 1 for left or right

use strict;

use CGI qw(-nosticky :standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use List::Util 'sum';
use HTML::Entities;

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();
my $orgId = $cgi->param('orgId') || die "No orgId found";
my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my $genome = $orginfo->{$orgId}{genome} || die "Invalid orgId $orgId";
my $expinfo = Utils::expinfo($dbh,$orgId); # expName => attribute => value

my @r = $cgi->param('r');
my @c = $cgi->param('c');
my $view = $cgi->param('view') || 0;

my @errors = ();

my $addrow = $cgi->param('addrow') || "";
$addrow =~ s/[ \t]+$//;
$addrow =~ s/^[ \t]+//;
my @addrow = split /[ \t]+/, $addrow;
my $addrowAt = $cgi->param('addrowAt') || 0;

foreach my $add (@addrow) {
  my $r = undef; # the row to add, if set
  if (grep { $_ eq $addrow} @r) {
    push @errors, "Gene $addrow is already included";
  } elsif ($add =~ m/^_l([0-9]+)$/) {
    $r = $add;
  } elsif ($add ne "") {
    if ($add !~ m/^[A-Za-z90-9_-]*$/) {
      push @errors, "Invalid gene to add";
    } else {
      my ($locusId) = $dbh->selectrow_array(qq{ SELECT locusId FROM Gene WHERE orgId = ?
                                                  AND (locusId = ? OR sysName = ? OR gene = ? COLLATE NOCASE) LIMIT 1 },
                                            {}, $orgId, $add, $add, $add);
      if (!defined $locusId) {
        push @errors, qq{Cannot find gene "$addrow"};
      } elsif (sum(map { $_ eq $locusId } @r) > 0) {
        push @errors, qq{Gene "$add" (locus $locusId) is already included};
      } else {
        $r = $locusId;
      }
    }
  }
  if (defined $r) {
    if ($addrowAt eq "-1") {
      unshift @r, $r;
    } else {
      push @r, $r;
    }
  }
}

my $addcol = $cgi->param('addcol') || "";
$addcol =~ s/[ \t]+$//;
$addcol =~ s/^[ \t]+//;
if ($addcol) {
  my $exps = Utils::matching_exps($dbh,$orgId,$addcol);
  if (@$exps == 0) {
    push @errors, qq{No experiments matching "$addcol"};
  } else {
    my %c = map { $_ => 1 } @c;
    my @keep = grep !exists $c{$_->{expName}}, @$exps;
    my @ignore = grep exists $c{$_->{expName}}, @$exps;
    if (@keep == 0) {
      push @errors, scalar(@ignore) . " matching experiments are already shown";
    } else {
      my $addcolAt = $cgi->param('addcolAt') || "";
      if ($addcolAt eq "-1") {
        unshift @c, map $_->{expName}, @keep;
      } else {
        push @c, map $_->{expName}, @keep;
      }
    }
  }
}

my @locusIds = grep !m/^_/, @r;

my %genes = (); # locusId => gene
foreach my $locusId (@locusIds) {
  my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId = ? AND locusId = ?",
                                     {}, $orgId, $locusId);
  die "No such locus $locusId in org $orgId" if !defined $gene->{locusId};
  # expName => "fit" => fitness value
  $gene->{fit} = $dbh->selectall_hashref(qq{SELECT expName,fit,t FROM GeneFitness
                                               WHERE orgId = ? AND locusId = ?},
                                         "expName", {}, $orgId, $locusId);
  foreach my $expName (keys %{ $gene->{fit} }) {
    die "No such experiment: $expName" unless exists $expinfo->{$expName};
  }
  $genes{$locusId} = $gene;
}

my $title1 = @locusIds > 0 ? "Fitness Heatmap for " . scalar(@locusIds) . " genes in "
  : "Build Fitness Heatmap for";

print $cgi->header,
  Utils::start_page($title1 . " " . $genome),
  h2($title1, a({-href => "org.cgi?orgId=$orgId"}, $genome)),
  qq[<div id="ntcontent">];

foreach my $error (@errors) {
  print $cgi->h3($error);
}

my @trows = ();
my @headings = qw{Gene Description};
foreach my $expName (@c) {
  die "Invalid column: $expName" unless exists $expinfo->{$expName};
  my $exp = $expinfo->{$expName};
  push @headings, a( {-href => "exp.cgi?orgId=$orgId&expName=$expName", -title => $expName},
                     $expinfo->{$expName}{expDesc});
}

my $ncol = scalar(@headings);
push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'}, $cgi->th(\@headings));

# remove-column row
if (@c > 0 && ! $view) {
  my @row = ();
  push @row, td(""), td("");
  foreach my $expName (@c) {
    my @args = ("orgId=$orgId", "view=$view");
    push @args, map "r=$_", @r;
    push @args, map "c=$_", grep $_ ne $expName, @c;
    foreach my $r (@r) {
      my $rt = $cgi->param("rt.$r");
      push @args, "rt.$r=$rt" if defined $rt;
    }
    my $URL = "heatmap.cgi?" . join("&", @args);
    push @row, td(a({ -href => $URL, title => "remove $expName" }, "remove"));
  }
  push @trows, $cgi->Tr({-align=>'CENTER', -valign=>'TOP'}, @row);
}

foreach my $rowspec (@r) {
  if ($rowspec =~ m/^_l\d+$/) {
    my $text = $cgi->param("rt.$rowspec");
    $text = " " if !defined $text;
    push @trows, $cgi->Tr(td( { -colspan => $ncol }, $text));
  } else {
    my $locusId = $rowspec;
    die "Invalid row $rowspec" unless exists $genes{$locusId};
    my $gene = $genes{$locusId};
    my @out = (); # fields for each column
    my $name =  Utils::gene_link($dbh, $gene, "name", "myFitShow.cgi");
    $name .= " " . "(" . $gene->{gene} . ")" if $gene->{gene};
    push @out, td($name);
    push @out, td(Utils::gene_link($dbh, $gene, "desc", "domains.cgi"));
    foreach my $expName (@c) {
      my $showId = $gene->{sysName} || $gene->{locusId};
      my $fit = $gene->{fit}{$expName}{fit};
      my $strainUrl = "strainTable.cgi?orgId=$orgId&locusId=$gene->{locusId}&expName=$expName";
      my $t = $gene->{fit}{$expName}{t};
      my $show = "&nbsp;";
      if (defined $fit) {
        $show = a({ -href => $strainUrl,
                    -title => "$showId: t = " . sprintf("%.1f",$t),
                    -style => "color:rgb(0,0,0)" },
                  sprintf("%.1f", $fit));
      }
      push @out, td({ -bgcolor => Utils::fitcolor($fit) }, $show);
    }
    push @trows, Tr(@out);
  }
}
print table( { cellspacing => 0, cellpadding => 3 }, @trows) if @r > 0;

my @hidden = (); # all of the common (static) arguments as hidden fields
# Do not save view as it never needs to be set by a hidden field (a form can only turn it off)
push @hidden, hidden('orgId', $orgId );
foreach my $r (@r) {
  push @hidden, hidden(-name => 'r', -default => $r, -override => 1);
}
foreach my $r (@r) {
  my $arg = "rt.$r";
  my $value = $cgi->param($arg);
  push @hidden, hidden($arg, $value) if defined $value;
}
foreach my $c (@c) {
  push @hidden, hidden(-name => 'c', -default => $c, -override => 1);
}

my $nlabel = 0;
while (grep($_ eq "_l$nlabel", @r)) {
  $nlabel++;
  die if $nlabel > 1000;
}

my $selectRowAt = qq{
<SELECT name="addrowAt">
<OPTION value=-1 >top</OPTION>
<OPTION value=1 SELECTED >bottom</OPTION>
</SELECT>
};

my $selectColAt = qq{
<SELECT name="addcolAt">
<OPTION value=-1 >left</OPTION>
<OPTION value=1 SELECTED >right</OPTION>
</SELECT>
};

my $go = "<BUTTON type='submit'>Go</BUTTON>";

if ($view) {
  # this form need not set view back to 0 as that is the default
  print
    start_form(-name => 'input', -method => 'GET', -action => 'heatmap.cgi'),
      p(join("",@hidden),
        "<BUTTON type='submit'>Edit</BUTTON>"),
          end_form;
} else {
    # I should really use div and text-style: width=12em or whatever to make the text have a consistent width...
  print
    start_form(-name => 'input', -method => 'GET', -action => 'heatmap.cgi'),
      p("Add gene(s):",
        join("", @hidden),
        textfield( -name => 'addrow', -default => "", -override => 1, -size => 10, -maxLength => 1000 ),
        " at $selectRowAt $go"),
          end_form;

  print
    start_form(-name => 'input', -method => 'GET', -action => 'heatmap.cgi'),
      p("Add label:",
        join("", @hidden),
        hidden(-name => 'addrow', -default => "_l$nlabel", -override => 1),
        textfield( -name => "rt._l$nlabel", -default => "", -override => 1, -size => 10, -maxLength => 100 ),
        " at $selectRowAt $go"),
          end_form;

  print
    start_form(-name => 'input', -method => 'GET', -action => 'heatmap.cgi'),
      p("Add experiment:",
        join("", @hidden),
        textfield( -name => 'addcol', -default => "", -override => 1, -size => 15, -maxLength => 100 ),
        " at $selectColAt $go"),
          end_form;

  print
    start_form(-name => 'input', -method => 'GET', -action => 'heatmap.cgi'),
      p(join("", @hidden),
        hidden(-name => 'view', -default => 1, -override => 1),
        "<BUTTON type='submit'>Hide controls</BUTTON>"),
          end_form;

  my $nloci = scalar(@locusIds);
  print p("Or see",
          a({ -href => "genesFit.cgi?orgId=$orgId&" . join("&", map "locusId=$_", @locusIds) },
            "top conditions for $nloci genes"))
    if $nloci > 1;
}

$dbh->disconnect();
Utils::endHtml($cgi);

