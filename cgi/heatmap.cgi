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
# r -- list of more row specifiers, in order from top to bottom. Each is either a locusId
#	or xref or a comment/label specifier like _l1 but varying the number
# rt.X for various row specifiers (i.e., rt._l1) -- text to show instead of
#	the usual gene annotation. (The original annotation will be in the popup text.)
#	This is also how the text for labels is set.
# c -- list of column specifiers (experiments), from left to right
# view -- set to 1 if in view-only mode
# edit -- set to a row (0-indexed) if editing the text for this row
#
# CGI parameters used to update the view:
# addrow -- search for matching gene(s) in this organism [multiple space-delimited terms]
# addrowAt -- -1, 0, or 1 for top, middle, or bottom
# addcol -- ditto but may add multiple experiments
# addcolAt -- -1 or 1 for left or right
# ca -- additional columns that may be added at left or right

use strict;
sub SetupHidden();

# use oldstyle urls to ignore ; as a potential separator of queries.
use CGI qw(-nosticky -oldstyle_urls :standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use List::Util 'sum';
use HTML::Entities;
use URI::Escape;

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();
my $orgId = $cgi->param('orgId') || die "No orgId found";
my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my $genome = $orginfo->{$orgId}{genome} || die "Invalid orgId $orgId";
my $expInfo = Utils::expinfo($dbh,$orgId); # expName => attribute => value
my $view = $cgi->param('view') || 0;

# To be set by SetupHidden()
# Do not save $view as it never needs to be set by a hidden field (a form can only turn it off)
my @hidden = (); # all of the common (static) arguments as hidden fields
my %hiddenrt = (); # row specifier => the hidden field for its text (if any)

# Process the parameters that define the veiw
# First, r and c specify the genes (rows) and experiments (columns)
# addrow / addrowAt adds a gene
# rt.locusId sets the labels for the genes
#   this must be processed before addcol or else the labels might be lost
#   in the experiment selection form
# addcol / addcolAt adds experiment(s)
#	and may need to show a form to select which ones to include
#	ca (short for column add) is set used by that form

my @loadErrors = ();
my @r = $cgi->param('r');
my @c = $cgi->param('c');

my $addrow = $cgi->param('addrow') || "";
$addrow =~ s/[ ,\t\r\n]+$//;
$addrow =~ s/^[ ,\t\r\n]+//;
my @addrow = split /[ ,\t\r\n]+/, $addrow;
my $addrowAt = $cgi->param('addrowAt') || 0;

foreach my $add (@addrow) {
  my $r = undef; # the row to add, if set
  if (grep { $_ eq $addrow} @r) {
    push @loadErrors, "Gene $addrow is already included";
  } elsif ($add =~ m/^_l([0-9]+)$/) {
    $r = $add;
  } elsif ($add ne "") {
    if ($add !~ m/^[A-Za-z90-9_.-]*$/) {
      push @loadErrors, "Invalid gene to add";
    } else {
      my $list = Utils::matching_genes($dbh, $orgId, $add);
      my $locusId = undef;
      $locusId = $list->[0]{locusId} if @$list == 1;
      if (!defined $locusId) {
        push @loadErrors, qq{Cannot find gene "$add"};
      } elsif (sum(map { $_ eq $locusId } @r) > 0) {
        push @loadErrors, qq{Gene "$add" (locus $locusId) is already included};
      } else {
        $r = $locusId;
      }
    }
  }
  if (defined $r) {
    if ($addrowAt eq "-1") {
      unshift @r, $r; # at end (bottom)
    } elsif ($addrowAt eq "1") {
      push @r, $r; # at beginning (top)
    } else { # add at middle
      my $mid = int((scalar(@r)-1)/2); # 0 if empty
      $mid++ if $mid < scalar(@r)-1;
      splice @r, $mid, 0, ($r);
    }
  }
}

my %genes = (); # locusId => gene
my @locusIds = ();
my %labels = (); # locusId or _{labelNumber} to label
# Convert locus specifiers to actual genes, and fetch all the genes and fitness data;
# also save the labels (which are input by locus specifier, which may change)
foreach my $i (0..(scalar(@r)-1)) {
  my $locusId = $r[$i];
  my $rt = param("rt.$locusId");
  unless ($locusId =~ m/^_/) {
    # not a comment line
    my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId = ? AND locusId = ?",
                                     {}, $orgId, $locusId);
    if (!defined $gene->{locusId}) {
      # Try again using xrefs
      # This allows these heatmaps to work if the assembly
      # changes (if xrefs were added for the old ids)
      my $list = $dbh->selectall_arrayref(qq{ SELECT * FROM LocusXref JOIN Gene USING (orgId,locusId)
                                              WHERE orgId = ? AND xrefId = ? },
                                          { Slice => {} }, $orgId, $locusId);
      if (@$list == 1) {
        $gene = $list->[0];
        $locusId = $gene->{locusId};
        $r[$i] = $gene->{locusId};
      } else {
        push @loadErrors, "Locus $locusId not found";
        $locusId = "";
        $r[$i] = ""; # will delete later
      }
    }
    if ($locusId ne "") {
      $genes{$locusId} = $gene;
      push @locusIds, $locusId;
    }
  }
  # Save label under the new locusId
  $labels{$locusId} = $rt if defined $rt && $rt ne "";
}
# Remove unknown loci
@r = grep { $_ ne "" } @r;

# Set up @rtargs for saving the labels and comments
my @rtargs;
foreach my $r (@r) {
  push @rtargs, "rt.$r=$labels{$r}" if defined $labels{$r};
}

my $addcol = $cgi->param('addcol') || "";
$addcol =~ s/[ \t]+$//;
$addcol =~ s/^[ \t]+//;
my $addcolAt = $cgi->param('addcolAt') || "";

if ($addcol) {
  my $exps = Utils::matching_exps($dbh,$orgId,$addcol);
  if (@$exps == 0) {
    push @loadErrors, qq{No experiments matching "$addcol"};
  } else {
    my %c = map { $_ => 1 } @c;
    my @keep = grep !exists $c{$_->{expName}}, @$exps;
    my @ignore = grep exists $c{$_->{expName}}, @$exps;
    if (@keep == 0) {
      push @loadErrors, scalar(@ignore) . " matching experiments are already shown";
    } elsif (@keep >= 5) {
      # Show checkbox form to choose which of these to include
      # I can just reuse the c parameter for each checkbox -- it only gets
      # added to the list if it is specified
      my $title1 = "Select experiments to include in heatmap for";
      print $cgi->header,
        Utils::start_page("$title1 $genome"),
            h2($title1, a({ -href => "org.cgi?orgId=$orgId" }, $genome));
      if (scalar(@keep) == scalar(@$exps)) {
        print p(scalar(@keep), qq{experiments matched "$addcol".},
                "Check the ones you want to add to the heatmap:");
      } else {
        print p("Of", scalar(@$exps), qq{experiments that matched "$addcol",},
                scalar(@keep), "are not in the heatmap.",
                "Check the ones you want to add:");
      }
      SetupHidden();
      print
        start_form(-name => 'select', -method => 'GET', -action => 'heatmap.cgi'),
          join("", @hidden),
            hidden(-name => 'addcolAt');
      my @th = qw{Name Description};
      my @trows = ();
      push @trows, Tr(th(\@th));
      foreach my $exp (@keep) {
        my @td = (checkbox(-name => "ca", -checked => 0, -value => $exp->{expName},
                           -label => $exp->{expName}),
                  a({ -href => "exp.cgi?orgId=$orgId&expName=$exp->{expName}" },
                    $exp->{expDesc}));
        push @trows, Tr(td(\@td));
      }
      print table({cellspacing => 0, cellpadding => 3 }, @trows),
        p(qq{<BUTTON type='submit'>Add</BUTTON>}),
        p(small("You can use % as a wild card when searching for experiments.")),
        end_form;
      $dbh->disconnect();
      Utils::endHtml($cgi); # exits
    } else {
      # small #experiments matched, just add them
      if ($addcolAt eq "-1") {
        unshift @c, map $_->{expName}, @keep;
      } else {
        push @c, map $_->{expName}, @keep;
      }
    }
  }
}

# Add the selected subset of experiments, at left or right
my @ca = $cgi->param('ca');
if (@ca) {
  if ($addcolAt eq "-1") {
    unshift @c, @ca;
  } else {
    push @c, @ca;
  }
}

my $maxC = 50;
my $maxR = 200;
if (@c > $maxC) {
  $#c = $maxC-1;
  push @loadErrors, "This page limits the number of columns to $maxC";
}
if (@r > $maxR) {
  $#r = $maxR-1;
  push @loadErrors, "This page limits the number of rows to $maxR";
}


my $title1 = @locusIds > 0 ? "Heatmap for " . scalar(@locusIds) . " genes x " . scalar(@c) . " experiments in "
  : "Build Heatmap in";

print $cgi->header,
  Utils::start_page($title1 . " " . $genome),
  h2($title1, a({-href => "org.cgi?orgId=$orgId"}, $genome));

foreach my $error (@loadErrors) {
  print $cgi->h3($error);
}

SetupHidden();

# Load the fitness data (if any) for each gene
while (my ($locusId, $gene) = each %genes) {
  # expName => "fit" => fitness value
  $gene->{fit} = $dbh->selectall_hashref(qq{SELECT expName,fit,t FROM GeneFitness
                                            WHERE orgId = ? AND locusId = ?},
                                         "expName", {}, $orgId, $locusId);
  foreach my $expName (keys %{ $gene->{fit} }) {
    die "No such experiment: $expName" unless exists $expInfo->{$expName};
  }
}

my $rmmark = "&#10799;"; # Unicode Character 'VECTOR OR CROSS PRODUCT' (U+2A2F)
my $editmark = span({ -style => "display: inline-block; transform: rotateZ(90deg);" }, "&#9998;");

my $edit = $cgi->param("edit");
if (defined $edit && $edit ne "") {
  # If in edit mode, just show a box to edit this text
  die "Invalid edit argument" unless $edit =~ m/^\d+$/ && $edit >= 0 && $edit < scalar(@r);
  my $r = $r[$edit];
  my $orig = $labels{$r};
  my $label;
  my $isGene = $r !~ m/^_l/;
  if ($isGene) {
    my $gene = $genes{$r};
    die unless defined $gene;
    $orig = $gene->{desc} unless defined $orig && $orig ne "";
    $label = "Revise description for " . Utils::gene_link($dbh, $gene, "name", "myFitShow.cgi");
    $label .= " " . "(" . $gene->{gene} . ")" if $gene->{gene};
    $label .= " to ";
  } else {
    $label = "Revise comment: ";
  }

  print
    p(start_form(-name => 'edit', -method => 'GET', -action => 'heatmap.cgi'),
      join("", grep { $_ ne $hiddenrt{$r} } @hidden),
      $label,
      br(),
      textfield( -name => "rt.$r", -default => $orig, -override => 1, -size => 50, -maxLength => 500 ),
      " <BUTTON type='submit'>Go</BUTTON> ",
      end_form);
  $dbh->disconnect();
  Utils::endHtml($cgi); # exits
}

my @trows = ();
my @headings = ();
push @headings, "&nbsp;" unless $view; # delete/edit/up/down controls
push @headings, qw{Gene Description};
foreach my $expName (@c) {
  print p({ -style => "color: red;" }, "Ignoring invalid experiment $expName (no longer in the Fitness Browser?)")
    unless exists $expInfo->{$expName};
}
@c = grep { exists $expInfo->{$_} } @c;
SetupHidden(); # updated @c

foreach my $expName (@c) {
  die "Invalid column: $expName" unless exists $expInfo->{$expName};
  my $exp = $expInfo->{$expName};
  my $title = "$expName by $exp->{person} on $exp->{dateStarted} - $exp->{media}";
  foreach my $i (1..4) {
    my $cond = $exp->{"condition_${i}"};
    my $units = $exp->{"units_${i}"};
    my $conc = $exp->{"concentration_${i}"};
    $title .= " + $conc $units $cond" if $cond ne "";
  }
  my $show = a( {-href => "exp.cgi?orgId=$orgId&expName=$expName", -title => $title},
                     $exp->{expDesc});
  unless ($view) {
    my @args = ("orgId=$orgId", "view=$view");
    push @args, map "r=$_", @r;
    push @args, @rtargs;
    push @args, map "c=$_", grep $_ ne $expName, @c;
    my $rmURL = "heatmap.cgi?" . join("&", @args);
    $show .= br() . a({ -href => $rmURL, title => "remove $expName" }, $rmmark);
  }
  push @headings, $show;
}
push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'}, $cgi->th(\@headings));
my $ncol = scalar(@headings);
my $ncolLabel = $ncol - ($view ? 0 : 1);

my $maxRow = scalar(@r)-1;
my $actionStyle = "color: rgb(103,146,160); font-size: smaller;";
# this includes everything but the rows
# view need not be specified as there are no outbound URLs in view mode
my $baseURL = "heatmap.cgi?" . join("&", "orgId=$orgId", @rtargs, map { "c=$_" } @c);
foreach my $iRow (0..$maxRow) {
  my $editURL = join("&", $baseURL, "edit=$iRow", map { "r=$_" } @r);
  my $rowspec = $r[$iRow];
  my $rmURL = join("&", $baseURL, map { "r=$_" } grep { $_ ne $rowspec } @r);

  my @td = ();
  unless ($view) {
    my @controls = ();
    push @controls, a({ -href => $rmURL, -title => "remove", -style => $actionStyle}, $rmmark);
    push @controls, a({ -href => $editURL, -title => "edit text" , -style => $actionStyle}, $editmark);
    if ($iRow > 0) {
      my @newr = ();
      push @newr, @r[0..($iRow-2)] if $iRow >= 2;
      push @newr, $r[$iRow], $r[$iRow-1];
      push @newr, @r[($iRow+1)..$maxRow] if $iRow < $maxRow;
      my $URL = join("&", $baseURL, map { "r=$_" } @newr);
      push @controls, a({-href => $URL, -title => "move up", -style => $actionStyle}, "&uarr;");
    }
    if ($iRow < $maxRow) {
      my @newr = ();
      push @newr, @r[0..($iRow-1)] if $iRow > 0;
      push @newr, $r[$iRow+1], $r[$iRow];
      push @newr, @r[($iRow+2)..$maxRow] if $iRow+2 <= $maxRow;
      my $URL = join("&", $baseURL, map {"r=$_" } @newr);
      push @controls, a({-href => $URL, -title => "move down", -style => $actionStyle}, "&darr;");
    }
    push @td, td(@controls);
  }

  if ($rowspec =~ m/^_l\d+$/) {
    push @td, td({-colspan => $ncolLabel}, encode_entities($labels{$rowspec}));
  } else {
    my $locusId = $rowspec;
    die "Invalid row $rowspec" unless exists $genes{$locusId};
    my $gene = $genes{$locusId};
    my @out = (); # fields for each column

    my $name =  Utils::gene_link($dbh, $gene, "name", "myFitShow.cgi");
    $name .= " " . "(" . $gene->{gene} . ")" if $gene->{gene};
    push @td, td($name);
    my $rt = $labels{$locusId};
    if (defined $rt && $rt ne "") {
      push @td, td( a({ -href => "domains.cgi?orgId=$orgId&locusId=$locusId",
                         -title => Utils::gene_link($dbh, $gene, "text") },
                      encode_entities($rt)));
    } else {
      push @td, td(Utils::gene_link($dbh, $gene, "desc", "domains.cgi"));
    }
    foreach my $expName (@c) {
      my $showId = $gene->{sysName} || $gene->{locusId};
      $showId .= " (" . $gene->{gene} . ")" if $gene->{gene};
      my $fitObj = $gene->{fit}{$expName};
      push @td, Utils::fittd(fit => $fitObj->{fit},
                             t => $fitObj->{t},
                             gene => $gene,
                             expName => $expName,
                             expDesc => $expInfo->{$expName}{expDesc});
    }
  }
  push @trows, Tr(@td);
}


# End the divs early so that any horizontal scroll-bar is visible at the bottom of the table instead of hidden farther down
# Also note that most pages have div #page, div #main, and div #ntcontent,
# but this page does not have the ntcontent div
print "</div></div>"; # end page, end main

print table( { cellspacing => 0, cellpadding => 3,
               -style => "table-layout: auto; overflow-x: scroll" },
             @trows) if @r > 0;

my $nlabel = 0;
while (grep($_ eq "_l$nlabel", @r)) {
  $nlabel++;
  die if $nlabel > 1000;
}

my $selectRowAt = qq{
<SELECT name="addrowAt">
<OPTION value=-1 >top</OPTION>
<OPTION value=0 >middle</OPTION>
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

my $tinyform = "";
if ($ENV{SERVER_NAME}) { # not set if testing from command line
  my $page_URL = "http://" . $ENV{SERVER_NAME} . $ENV{REQUEST_URI};
  $tinyform = join("",
                   start_form(-style => "display: inline;",
                              -name => 'tiny', -method => 'POST', target => "blank",
                              -action => 'https://tinyurl.com/create.php'),
                   hidden(-name => 'url', -default => $page_URL, -override => 1),
                   " <BUTTON type='submit'>TinyURL</BUTTON> ",
                   end_form);
}

if ($view) {
  # this form need not set view back to 0 as that is the default
  print
    p(start_form(-style => "display: inline;",
                 -name => 'input', -method => 'GET', -action => 'heatmap.cgi'),
      join("",@hidden),
      " <BUTTON type='submit'>Edit</BUTTON> ",
      end_form,
      $tinyform);
} else {
    # I should really use div and text-style: width=12em or whatever to make the text have a consistent width...
  print
    p(start_form(-style => "display: inline;",
                 -name => 'input', -method => 'GET', -action => 'heatmap.cgi'),
      "Add genes:",
      join("", @hidden),
      textfield( -name => 'addrow', -default => "", -override => 1, -size => 10, -maxLength => 100 ),
      " at $selectRowAt $go",
      end_form),
    p(start_form(-style => "display: inline;",
                 -name => 'input', -method => 'GET', -action => 'heatmap.cgi'),
      "Add comment:",
      join("", @hidden),
      hidden(-name => 'addrow', -default => "_l$nlabel", -override => 1),
      textfield( -name => "rt._l$nlabel", -default => "", -override => 1, -size => 10, -maxLength => 500 ),
      " at $selectRowAt $go",
        end_form),
    p(start_form(-name => 'input', -method => 'GET', -action => 'heatmap.cgi'),
      "Add experiments:",
      join("", @hidden),
      textfield( -name => 'addcol', -default => "", -override => 1, -size => 15, -maxLength => 100 ),
      " at $selectColAt $go",
      end_form),
    p(start_form(-style => "display: inline;", -name => 'input', -method => 'GET', -action => 'heatmap.cgi'),
      join("", @hidden),
      hidden(-name => 'view', -default => 1, -override => 1),
      "<BUTTON type='submit'>Hide controls</BUTTON>",
      end_form,
      $tinyform);

  my $nloci = scalar(@locusIds);
  print p("Or see",
          a({ -href => "genesFit.cgi?orgId=$orgId&" . join("&", map "locusId=$_", @locusIds) },
            "top conditions for $nloci genes"))
    if $nloci > 1;
  print p("Or see",
          a({ -href => "expCor.cgi?orgId=$orgId&" . join("&", map { "e=$_" } @c) },
            "experiment correlations"))
    if scalar(@c) > 1;
  print p("Or download",
          a({ -href => "createFitData.cgi?" . join("&", "orgId=$orgId", map { "expName=$_" } @c) },
            "fitness data"),
          "for these", scalar(@c), "experiments and all genes")
    if @c  > 0;
  print p("Or try the",
          a({ -href => "cmpbrowser.cgi?anchorOrg=$orgId&" . join("&", map {"anchorLoci=$_"} @locusIds) },
            "comparative fitness browser"))
    if $nloci > 0;
}

$dbh->disconnect();
Utils::endHtml($cgi);

sub SetupHidden() {
  # Empty these in case this isn't the first invocation
  @hidden= ();
  %hiddenrt = ();
  push @hidden, hidden('orgId', $orgId );
  foreach my $r (@r) {
    push @hidden, hidden(-name => 'r', -default => $r, -override => 1);
  }
  foreach my $r (@r) {
    if (defined $labels{$r}) {
      my $arg = "rt.$r";
      my $h = hidden($arg, $labels{$r});
      push @hidden, $h;
      $hiddenrt{$r} = $h;
    }
  }
  foreach my $c (@c) {
    push @hidden, hidden(-name => 'c', -default => $c, -override => 1);
  }
}
