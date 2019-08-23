#!/usr/bin/perl -w
use strict;
#######################################################
## cmpbrowser.cgi -- show regions of scaffolds,
## as defined by genes, and color code
## similar or orthologous genes
##
## Copyright (c) 2019 University of California
##
## Author: Morgan Price
#######################################################
#
# Required arguments:
# anchorOrg and anchorLoci
#   These are used to easily add genomes to the comparative browser.
#   They are also used to set orgId, locusIdLeft, and locusIdRight,
#   if those are not set.
#   anchorLoci may have multiple values
#
# Common arguments defining a multi-track view:
# o for orgId, b for locusIdBeg, e for locusIdEnd, and s for strand
#   These may have multiple values (1 per track)
#   They define the loci shown on each track
#   For shorter encoding, strand may be "1" or "+" for + strand; otherwise means - strand
# e.orgId -- experiments for that orgId (multiply valued)
#
# Arguments for modifying the view:
# addOrg -- (select orthologous genes of the anchor(s) from this organism
# addRange -- a locusId or sysName, or left:right
# flipTrack (0-based) -- changes which strand is used for that track
# removeTrack -- a track to remove
# changeTrack and changeLeft  (+1 for another gene, -1 for fewer genes)
# changeTrack and changeRight (+1 for another gene, -1 for fewer genes)
# upTrack -- a track to move up
# downTrack -- a track to move down
# expTrack (0-based) -- a track to add experiments to. All tracks for this organism are modified.
# addExp -- a search term for experiments
#
# Arguments for coloring:
# colorOff -- offset to the color rotation
#
# Argument for modifying anchors:
# resetAnchor -- reset anchor to the first (top) track

use strict;
use CGI qw(-nosticky :standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use lib "../lib";
use Utils;
use List::Util qw{min max sum};
use HTML::Entities;

sub GetGene($$$);
sub AddGeneToTracks($$);

my $maxNGenes = 50; # on one track
my $newTrackDistance = 10*1000; # how far a gene is before we put it on a new track
my $kbWidth = 150; # svg units/kb
my $trackHeight = 50; # units per track
my $barHeight = $trackHeight * 0.7; # height of scale bar region
my $arrowSize = 30; # size of the arrow part of a gene
my $padding = 30; # at left only

{
  my $cgi=CGI->new;
  my $dbh = Utils::get_dbh();
  my $orginfo = Utils::orginfo($dbh);

  # Parse anchor parameters
  my $anchorOrg = param('anchorOrg') || die "Must specify anchorOrg";
  die "Invalid organism $anchorOrg" unless exists $orginfo->{$anchorOrg};
  $CGI::LIST_CONTEXT_WARN = 0; # no warnings for lixt context use
  my @anchorLoci = param('anchorLoci');
  die "Must specify anchorLoci" unless @anchorLoci > 0;
  my @anchorGenes = map GetGene($dbh, $anchorOrg, $_), @anchorLoci;

  my @orgId = param('o');
  my @locusIdBeg = param('b');
  my @locusIdEnd = param('e');
  my @trackStrand = param('s');

  die "o b e and s must be the same lengths"
    unless scalar(@orgId) == scalar(@locusIdBeg)
      && scalar(@orgId) == scalar(@locusIdEnd)
      && scalar(@orgId) == scalar(@trackStrand);

  # Build tracks from anchor parameters or track parameters
  my @tracks = (); # each is a hash including orgId, geneBeg, geneEnd, strand
  # and additional fields defined later:
  # genes -- the list of genes, from left to right
  if (@orgId == 0) {
    # Design a layout from the anchor genes
    foreach my $gene (@anchorGenes) {
      AddGeneToTracks(\@tracks, $gene);
    }
  } else {
    for (my $i = 0; $i < scalar(@orgId); $i++) {
      my $orgId = $orgId[$i];
      die "Invalid orgId $orgId" unless exists $orginfo->{$orgId};
      my $geneBeg = GetGene($dbh, $orgId, $locusIdBeg[$i]);
      my $geneEnd = GetGene($dbh, $orgId, $locusIdEnd[$i]);
      die "Mismatched scaffolds for $geneBeg $geneEnd"
        unless $geneBeg->{scaffoldId} eq $geneEnd->{scaffoldId};
      push @tracks, { orgId => $orgId,
                      geneBeg => $geneBeg,
                      geneEnd => $geneEnd,
                      strand => $trackStrand[$i] eq "+" || $trackStrand[$i] eq "1" ? "+" : "-" };
    }
  }

  # Load experiment parameters
  # orgId => list of experiment objects
  my %orgExps = map { $_->{orgId} => [] } @tracks;
  foreach my $orgId (keys %orgExps) {
    my @expNames = param("e.$orgId");
    my %seen = ();
    foreach my $expName (@expNames) {
      next if exists $seen{$expName};
      $seen{$expName} = 1;
      my $exp = $dbh->selectrow_hashref("SELECT * FROM Experiment WHERE orgId = ? AND expName = ?",
                                        {}, $orgId, $expName);
      die "Unknown experiment $expName for $orgId" unless defined $exp;
      push @{ $orgExps{$orgId} }, $exp;
    }
  }

  # Handle view-modification parameters
  my @warnings = ();
  if (param('addOrg')) {
    # Add tracks by organism & orthology
    my $orgA = param('addOrg');
    die "Cannot add anchor" if $orgA eq $anchorOrg;
    die "Invalid addOrg $orgA" unless exists $orginfo->{$orgA};
    # Find orthologs of all the anchor genes
    my @orthIds = ();
    foreach my $anchorGene (@anchorGenes) {
      my ($orthId) = $dbh->selectrow_array(qq{ SELECT locusId2 FROM Ortholog
                                               WHERE orgId1 = ? AND locusId1 = ? AND orgId2 = ? },
                                           {}, $anchorOrg, $anchorGene->{locusId}, $orgA);
      push @orthIds, $orthId if defined $orthId;
    }
    if (@orthIds == 0) {
      push @warnings, "Sorry, no orthologs of " . scalar(@anchorGenes). " genes from $orginfo->{$anchorOrg}{genome}"
        . " were found in $orginfo->{$orgA}{genome}.";
    } else {
      push @warnings, "Only " . scalar(@orthIds). " of " . scalar(@anchorGenes) . "genes from $orginfo->{$anchorOrg}{genome}"
        . " were found in $orginfo->{$orgA}{genome}."
          unless scalar(@orthIds) == scalar(@anchorGenes);
      my $nOldTracks = scalar(@tracks);
      my $nChange = 0;
      foreach my $orthId (@orthIds) {
        my $orth = GetGene($dbh, $orgA, $orthId);
        $nChange += AddGeneToTracks(\@tracks, $orth);
      }
      if ($nChange == 0) {
        push @warnings, "All " . scalar(@orthIds) . " orthologs in $orginfo->{$orgA}{genome} were already shown.";
      } elsif ($nOldTracks == scalar(@tracks)) {
        push @warnings, "Track(s) for $orginfo->{$orgA}{genome} were expanded to include " . scalar(@orthIds) . " orthologs.";
      }
    }
  } elsif (param('addRange')) {
    # Add tracks by organism & genes
    my $addSpec = param('addRange');
    # Removing leading white space
    $addSpec =~ s/^ +//;
    $addSpec =~ s/ +$//;
    # turn whitespace comma or : into :
    $addSpec =~ s/[, \t:]+/:/g;
    if ($addSpec !~ m/^[0-9a-zA-Z_.:-]+/) {
      push @warnings, "Invalid add. Must be a locus tag or locus1:locus2";
    }
    my @sysNames = ( $addSpec );
    @sysNames = split /:/, $addSpec if $addSpec =~ m/:/;
    foreach my $sysName (@sysNames) {
      my $genes = $dbh->selectall_arrayref("SELECT * from Gene WHERE sysName = ?",
                                           { Slice => {} }, $sysName);
      $genes = $dbh->selectall_arrayref("SELECT * from Gene WHERE locusId = ?",
                                           { Slice => {} }, $sysName)
        if @$genes == 0;
      if (@$genes == 0) {
        push @warnings, "Locus $sysName was not found.";
      } elsif (@$genes > 1) {
        push @warnings, "Locus $sysName was ambiguous.";
      } else {
        my $change = AddGeneToTracks(\@tracks, $genes->[0]);
        push @warnings, "Locus $sysName was already shown." unless $change;
      }
    }
  } elsif (defined param('flipTrack')) {
    my $iTrack = param('flipTrack');
    die "Invalid flip $iTrack" unless $iTrack >= 0 && $iTrack < @tracks;
    my $track = $tracks[$iTrack];
    $track->{strand} = $track->{strand} eq "+" ? "-" : "+";
  } elsif (defined param('removeTrack')) {
    my $iTrack = param('removeTrack');
    die "Invalid remove $iTrack" unless $iTrack >= 0 && $iTrack < @tracks;
    my @tracksNew = ();
    for (my $i = 0; $i < scalar(@tracks); $i++) {
      push @tracksNew, $tracks[$i] unless $i == $iTrack;
    }
    @tracks = @tracksNew;
  } elsif (defined param('changeTrack') &&
           (defined param('changeLeft') || defined param('changeRight'))) {
    my $iTrack = param('changeTrack');
    die "Invalid change $iTrack" unless $iTrack >= 0 && $iTrack < @tracks;
    my $track = $tracks[$iTrack];
    my $orgId = $track->{orgId};

    my ($side, $nAdd);
    if (defined param('changeLeft')) {
      $side = "left";
      $nAdd = param('changeLeft');
    } else {
      $side = "right";
      $nAdd = param('changeRight');
    }
    die "Invalid changeLeft or changeRight $nAdd" unless $nAdd == 1 || $nAdd == -1;
    my $modBeg = $side eq ($track->{strand} eq "+" ? "left" : "right");
    my $geneToChange = $modBeg ? "geneBeg" : "geneEnd";
    if ($nAdd == -1 && $track->{geneBeg}{locusId} eq $track->{geneEnd}{locusId}) {
      # No genes left!
      $track->{geneBeg} = undef;
      @tracks = grep defined $_->{geneBeg}, @tracks;
      push @warnings, "Removed track for $orginfo->{$orgId}{genome} with no genes left.";
    } else {
      my $scGenes = $dbh->selectall_arrayref(qq{ SELECT * FROM Gene
                                                 WHERE orgId = ? AND scaffoldId = ?
                                                 ORDER BY begin },
                                             { Slice => {} }, $orgId, $track->{geneBeg}{scaffoldId} );
      die "No genes on scaffold" unless @$scGenes > 0;
      my $iGene;
      ($iGene) = grep { $scGenes->[$_]{locusId} eq $track->{ $geneToChange }{locusId} } (0..(@$scGenes)-1);
      die unless defined $iGene;
      if ($modBeg) {
        $iGene -= $nAdd;
      } else {
        $iGene += $nAdd;
      }
      if ($iGene < 0 || $iGene > scalar(@$scGenes)) {
        push @warnings, "Sorry, no additional genes are available for this track\n";
      } else {
        $track->{ $modBeg ? "geneBeg" : "geneEnd" } = $scGenes->[$iGene];
      }
    }
  } elsif (defined param('upTrack')) {
    my $iTrack = param('upTrack');
    die "Invalid move track $iTrack" unless $iTrack >= 0 && $iTrack < @tracks;
    if ($iTrack > 0) {
      my @tracksNew = ();
      push @tracksNew, @tracks[0..($iTrack-2)] if $iTrack >= 2;
      push @tracksNew, $tracks[$iTrack], $tracks[$iTrack-1];
      push @tracksNew, @tracks[($iTrack+1)..(scalar(@tracks)-1)] if $iTrack+1 < scalar(@tracks);
      @tracks = @tracksNew;
    }
  } elsif (defined param('downTrack')) {
    my $iTrack = param('downTrack');
    die "Invalid move track $iTrack" unless $iTrack >= 0 && $iTrack < @tracks;
    if ($iTrack < scalar(@tracks)-1) {
      my @tracksNew = ();
      push @tracksNew, @tracks[0..($iTrack-1)] if $iTrack >= 1;
      push @tracksNew, $tracks[$iTrack+1], $tracks[$iTrack];
      push @tracksNew, @tracks[($iTrack+2)..(scalar(@tracks)-1)] if $iTrack+2 < scalar(@tracks);
      @tracks = @tracksNew;
    }
  } elsif (defined param('expTrack') && defined param('addExp')) {
    my $iTrack = param('expTrack');
    die "Invalid add-exp track $iTrack" unless $iTrack >= 0 && $iTrack < @tracks;
    my $expSpec = param('addExp');
    $expSpec =~ s/^ +//;
    $expSpec =~ s/ +$//;
    if ($expSpec ne "") {
      my $orgId = $tracks[$iTrack]{orgId};
      my $exps = Utils::matching_exps($dbh, $orgId, $expSpec);
      if (@$exps == 0) {
        push @warnings, "No experiments matching '" . encode_entities($expSpec) . "' were found"
          . " in $orginfo->{$orgId}{genome}";
      } else {
        my %known = ();
        foreach my $exp (@{ $orgExps{$orgId} }) {
          $known{$exp->{expName}} = 1;
        }
        my @add = grep !exists $known{$_->{expName}}, @$exps;
        if (@add == 0) {
          push @warnings, "All experiments matching '" . encode_entities($expSpec)
            . " for $orginfo->{$orgId}{genome} are already shown";
        } else {
          my $maxAdd = 5;
          if (@add > $maxAdd) {
            push @warnings, "Too many experiments match '" . encode_entities($expSpec) . "', only adding $maxAdd";
            @add = splice(@add, 0, $maxAdd);
          } else {
            push @warnings, "Adding " . scalar(@add) . " experiments for $orginfo->{$orgId}{genome}";
          }
          foreach my $exp (@add) {
            push @{ $orgExps{$orgId} }, $exp;
          }
        }
      }
    }
  }

  # Fetch all the genes for each track
  # Note that genes overlapping the 1st or last gene could be excluded from the view
  foreach my $track (@tracks) {
    die unless $track->{geneBeg}{orgId} eq $track->{orgId};
    die unless $track->{geneEnd}{orgId} eq $track->{orgId};
    die unless $track->{geneBeg}{scaffoldId} eq $track->{geneEnd}{scaffoldId};
    my $scaffoldId = $track->{geneBeg}{scaffoldId};
    my $minEnd = $track->{geneBeg}{end};
    my $maxBeg = $track->{geneEnd}{begin};
    my $genes = $dbh->selectall_arrayref(qq{ SELECT * FROM Gene
                                             WHERE orgId = ? AND scaffoldId = ?
                                              AND end >= $minEnd AND begin <= $maxBeg },
                                         { Slice => {} },
                                         $track->{orgId}, $scaffoldId);
    die "No genes in range for $track->{orgId} $scaffoldId $minEnd $maxBeg"
      unless @$genes > 0;
    die "geneBeg $track->{geneBeg}{locusId} not included in $track->{orgId} $scaffoldId $minEnd $maxBeg"
      unless scalar(grep $_->{locusId} eq $track->{geneBeg}{locusId}, @$genes) == 1;
    die "geneEnd $track->{geneEnd}{locusId} not included in $track->{orgId} $scaffoldId $minEnd $maxBeg"
      unless scalar(grep $_->{locusId} eq $track->{geneEnd}{locusId}, @$genes) == 1;
    $track->{genes} = $genes;
    die "Too many genes for $track->{orgId}" if @$genes > $maxNGenes;
  }
  my $nTracks = scalar(@tracks);
  my $anchorShow = $anchorGenes[0]{sysName} || $anchorGenes[0]{locusId};
  print header,
    Utils::start_page("Comparative Browser for $anchorShow ($nTracks tracks)"),
    q{<div id="ntcontent">};
  print p({-style => "color: red;" }, join(br(), @warnings))
    if @warnings > 0;

  # Compute colors -- first, compute groups
  my %geneGroup = (); # orgId => locusId => iGroup
  my $nGeneGroups = 0;
  my %orth = (); # orgId => locusId => orgId => orthId
  foreach my $track (@tracks) {
    my $orgId = $track->{orgId};
    foreach my $gene (@{ $track->{genes} }) {
      $orth{$orgId}{ $gene->{locusId} } = {};
    }
  }
  while (my ($orgId, $hash) = each %orth) {
    while(my ($locusId, $ohash) = each %$hash) {
      my $rows = $dbh->selectall_arrayref(qq{ SELECT orgId2,locusId2 FROM Ortholog
                                              WHERE orgId1 = ? AND locusId1 = ? },
                                          {}, $orgId, $locusId);
      foreach my $row (@$rows) {
        my ($orgId2,$locusId2) = @$row;
        # Ignore orthologs that are not shown
        $ohash->{$orgId2} = $locusId2 if exists $orth{$orgId2}{$locusId2};
      }
    }
  }
  foreach my $track (@tracks) {
    my $orgId = $track->{orgId};
    foreach my $gene (@{ $track->{genes} }) {
      my $locusId = $gene->{locusId};
      my $ohash = $orth{$orgId}{$locusId};
      if (keys(%$ohash) > 0) {
        while (my ($orgId2, $locusId2) = each %$ohash) {
          if (exists $geneGroup{$orgId2}{$locusId2}) {
            $geneGroup{$orgId}{$locusId} = $geneGroup{$orgId2}{$locusId2};
            last;
          }
        }
        if (!exists $geneGroup{$orgId}{$locusId}) {
          $geneGroup{$orgId}{$locusId} = $nGeneGroups++;
          while (my ($orgId2,$locusId2) = each %$ohash) {
            $geneGroup{$orgId2}{$locusId2} = $geneGroup{$orgId}{$locusId};
          }
        }
      }
    }
  }

  if (param('resetAnchor') && @tracks > 0) {
    $anchorOrg = $tracks[0]{orgId};
    @anchorGenes = @{ $tracks[0]{genes} };
    @anchorLoci = map $_->{locusId}, @anchorGenes;
  }

  my @colors = qw{Red Green Yellow Blue Orange Purple Cyan Magenta Lime Pink Teal Lavender Brown Beige Maroon MediumSpringGreen Olive Coral};
  my $colorOff = param('colorOff') || 0;

  # Save the values to use for a URL or hidden form elements
  my @hidden = (['anchorOrg', $anchorOrg]);
  push @hidden, map ['anchorLoci', $_], @anchorLoci;
  foreach my $track (@tracks) {
    push @hidden, ['o', $track->{orgId}],
      ['b', $track->{geneBeg}{locusId}],
      ['e', $track->{geneEnd}{locusId}],
      ['s', $track->{strand} eq "+" ? 1 : 0];
  }
  while (my ($orgId, $exps) = each %orgExps) {
    foreach my $exp (@$exps) {
      push @hidden, [ "e.$orgId", $exp->{expName} ];
    }
  }
  push @hidden, [ 'colorOff', $colorOff ];
  my $tracksURL = "cmpbrowser.cgi?" . join("&", map $_->[0] . "=" . $_->[1], @hidden);
  my $tracksHidden = join("\n",
                          map qq{<input type="hidden" name="$_->[0]" value="$_->[1]">}, @hidden);

  # Output the tracks
  for (my $iTrack = 0; $iTrack < @tracks; $iTrack++) {
    my $track = $tracks[$iTrack];
    my $orgId = $track->{orgId};
    my $genes = $track->{genes};
    my $invert = $track->{strand} ne "+" && $track->{strand} ne "1";
    my $bAddScale = $iTrack == scalar(@tracks-1);

    # header line shows organism name and link to fitness data
    my @loci = map $_->{locusId}, @$genes;
    @loci = reverse @loci if $invert;
    my $baseCGI = @loci > 1 ? "genesFit.cgi" : "singleFit.cgi";
    my $URL = join("&", "${baseCGI}?orgId=${orgId}", map "locusId=$_", @loci);
    my $g = $orginfo->{$orgId};
    my $linkColorStyle = "color: darkblue;";
    my $removeGlyph = "&#10799;"; # looks like an x
    my $removeLink = a({-href => "$tracksURL&removeTrack=$iTrack",
                        -style => $linkColorStyle,
                       -title => "Remove this track"},
                       $removeGlyph);
    my $arrowStyle = "$linkColorStyle font-weight: bold; font-size: 150%;";
    my $flipLink = a({ -href => "$tracksURL&flipTrack=$iTrack", -style => $arrowStyle,
                       -title => "Flip strand for this track" },
                     "&harr;"); # left-right arrow
    my $linkUp = a({ -href => "$tracksURL&upTrack=$iTrack", -style => $arrowStyle,
                     -title => "Move this track up" },
                   "&uarr;");
    $linkUp = "" if $iTrack == 0;
    my $linkDn = a({ -href => "$tracksURL&downTrack=$iTrack", -style => $arrowStyle,
                     -title => "Move this track down" },
                   "&darr;");
    $linkDn = "" if $iTrack == scalar(@tracks)-1;
    my $nGenes = scalar(@$genes);

    # For each organism, there's a div with the top line containing
    # the organism name, controls for adding experiments, and the "arrow" controls
    # (And then there's the svg, the (optional) heatmap, and the spacer div)
    print
      div({-style => "width: 100%; padding-top: 0.5em; padding-bottom: 1em; position: relative;"},
          start_form(-name => 'input', -method => 'GET', -action => 'cmpbrowser.cgi'),
          $tracksHidden,
          qq{<input type="hidden" name="expTrack" value="$iTrack">},
          a({-href => $URL, -title => "see all fitness data for $nGenes genes",
             style => "$linkColorStyle; font-size: 110%;" },
            i($g->{genus}, $g->{species}), $g->{strain}),
          span({-style => "text-align: right; padding-left: 10%; position: absolute; left: 40%; top: 0px;"},
             "Experiments:",
             textfield(-name => 'addExp', -size => 15, -maxlength => 50, -default => '', -override => 1),
             submit(-style => 'text-align: left; float: none; display:inline; padding: 2px 2px; ',
                    -value => "Add", -name => ""),
             "&nbsp;",
             $linkUp . $linkDn . $flipLink . $removeLink),
          end_form),
      "\n";

    # svg shows the genes in order
    my $xmin = min(map $_->{begin}, @$genes);
    my $xmax = max(map $_->{end}, @$genes);
    my $xdiff = $xmax - $xmin;
    # min. 1 kb for scale bar
    my $xdiffUse = $xdiff;
    $xdiffUse = 1000 if $xdiffUse < 1000 && $bAddScale;
    my $svg_width = $padding + $xdiffUse * $kbWidth / 1000.0;

    # SVG has +y axis going down, not up
    my $top = 0;
    my $genetop = $top;
    my $bottom = $trackHeight;
    my $geneymid = ($genetop+$bottom)/2;
    my $right = $padding + $xdiff * $kbWidth/1000.0;
    my $svgHeight = $trackHeight;
    $svgHeight += $barHeight if $bAddScale;

    # Build the svg (but do not output it)
    my $svg = qq{<svg width="${svg_width}" height="${svgHeight}" style="position: relative; left: 1em;">\n};
    # center line
    $svg .= qq{<line x1="$padding" y1="$geneymid" x2="$right" y2="$geneymid" style="stroke:black; stroke-width:1;"/>\n};

    foreach my $gene (@$genes) {
      my $start = $gene->{strand} eq "+" ? $gene->{begin} : $gene->{end};
      my $stop = $gene->{strand} eq "+" ? $gene->{end} : $gene->{begin};
      if ($invert) {
        $start = $xmax - $start;
        $stop = $xmax - $stop;
      } else {
        $start = $start - $xmin;
        $stop = $stop - $xmin;
      }
      my $xstart = $padding + $start * $kbWidth/1000.0;
      my $xstop = $padding + $stop * $kbWidth/1000.0;
      my @points; # list of x,y pairs
      if (abs($xstop-$xstart) > $arrowSize) {
        my $xmid = $xstart < $xstop ? $xstop - $arrowSize : $xstop + $arrowSize;
        @points = ([$xstart,$bottom], [$xstart,$genetop], [$xmid,$genetop], [$xstop,$geneymid], [$xmid, $bottom]);
      } else {
        @points = ([$xstart,$bottom], [$xstart,$genetop], [$xstop,$geneymid]);
      }
      my $pointstr = join(" ", map { $_->[0].",".$_->[1] } @points);
      $gene->{color} = "white"; # default color
      if (exists $geneGroup{$orgId}{$gene->{locusId}}) {
        my $iColor = ($colorOff + $geneGroup{$orgId}{$gene->{locusId}}) % scalar(@colors);
        $gene->{color} = $colors[$iColor];
      }
      my $poly = qq{<polygon points="$pointstr" style="fill:$gene->{color}; stroke:black; stroke-width:1;" />};

      my $showId = $gene->{gene} || $gene->{sysName} || $gene->{locusId};
      $showId =~ s/^.*_/_/ if defined $showId;
      my $xlabel = ($xstart+$xstop)/2;
      my $label = qq{<text x="$xlabel" y="$geneymid" fill="black" alignment-baseline="middle" text-anchor="middle">$showId</text>};
      $baseCGI = $gene->{type} eq 1 ? "domains.cgi" : "geneOverview.cgi";
      my $URL = encode_entities( "${baseCGI}?orgId=$orgId&gene=$gene->{locusId}" );
      $svg .= qq{<a xlink:href="$URL">};
      $svg .= "<title>" . encode_entities( Utils::gene_link($dbh, $gene, "text") ) . "</title>\n";
      $svg .= $poly. $label . qq{</a>\n};
    }
    if ($iTrack == scalar(@tracks) - 1) { # add scale bar for 1 kb at bottom
      my @bary = ($trackHeight + $barHeight * 0.7, $trackHeight + $barHeight * 0.95);
      my $barAt = ($bary[0] + $bary[1])/2;
      my $barright = $padding + $kbWidth;
      $svg .= qq{<line x1="$padding" y1="$barAt" x2="$barright" y2="$barAt" style="stroke:black; stroke-width:1"/>\n};
      foreach my $x ($padding, $barright) {
        $svg .= qq{<line x1="$x" y1="$bary[0]" x2="$x" y2="$bary[1]" style="stroke:black; stroke-width:1"/>\n};
      }
      my $barcenter = ($padding + $barright)/2;
      my $barLabelY = $barAt - $barHeight/10;
      $svg .= qq{<text x="$barcenter" y="$barLabelY">1 kb</text>\n};
    }
    $svg .= qq{</svg>\n};

    # the track (svg object) has controls to the left and right, so use a div container
    # Use position relative for the div so that position "absolute" for the pieces
    # will be relative to it
    my @arrowsLeft = ( a({ -href => "$tracksURL&changeTrack=${iTrack}&changeLeft=1",
                           -title => "Add the next gene on the left",
                           -style => "position: absolute; top: -0.250em; $arrowStyle" }, "&larr;"),
                       a({ -href => "$tracksURL&changeTrack=${iTrack}&changeLeft=-1",
                           -title => "Remove the left-most gene",
                           -style => "position: absolute; top: 0.75em; $arrowStyle" }, "&rarr;") );
    my @arrowsRight = ( a({ -href => "$tracksURL&changeTrack=${iTrack}&changeRight=1",
                           -title => "Add the next gene on the right",
                           -style => "position: absolute; top: -0.25em; $arrowStyle" }, "&rarr;"),
                       a({ -href => "$tracksURL&changeTrack=${iTrack}&changeRight=-1",
                           -title => "Remove the right-most gene",
                           -style => "position: absolute; top: 0.75em; $arrowStyle" }, "&larr;") );
    print div({-style => "width:100%; position: relative;"},
              join("\n", @arrowsLeft, ""),
              $svg,
              span({-style => "display: inline-block; position: relative; left: 2em; top: 0em; vertical-align: top;"},
                   join("\n", @arrowsRight, "")));

    # Add fitness data heatmap, if any
    my $exps = $orgExps{$orgId};
    if (@$exps > 0) {
      # Fetch data
      foreach my $gene (@$genes) {
        $gene->{fit} = $dbh->selectall_hashref(qq{SELECT expName,fit,t FROM GeneFitness
                                                  WHERE orgId = ? AND locusId = ?},
                                               "expName", {}, $orgId, $gene->{locusId} );
      }

      # Build table -- Group, Condition, 1 column per gene, and the remove control
      my @trows = ();
      my @headings = map th($_), qw{Group Condition};
      my @genesInOrder = @$genes;
      @genesInOrder = reverse @genesInOrder if $invert;
      foreach my $gene (@genesInOrder) {
        # Gene headings have background color as in svg and black foreground color
        my $link = Utils::gene_link($dbh, $gene, "name", "myFitShow.cgi");
        $link =~ s/<A /<A style="color:black;" /i;
        push @headings, th({-style => "background-color: $gene->{color};"}, $link);
      }
      push @headings, "";
      push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'}, @headings);
      foreach my $exp (@$exps) {
        my $expName = $exp->{expName};
        my @td = ( td( $exp->{expGroup} ),
                   td( a({ -href => "exp.cgi?orgId=$orgId&expName=$expName",
                           -title => $expName }, $exp->{expDesc}) )
                 );
        foreach my $gene (@genesInOrder) {
          my $showId = $gene->{sysName} || $gene->{locusId};
          my $fit = $gene->{fit}{$exp->{expName} }{fit};
          my $strainUrl = "strainTable.cgi?orgId=$orgId&locusId=$gene->{locusId}&expName=$expName";
          my $t = $gene->{fit}{$expName}{t};
          my $show = "&nbsp;";
          if (defined $fit) {
            $show = a({ -href => $strainUrl,
                        -title => "$showId: t = " . sprintf("%.1f",$t),
                        -style => "color:rgb(0,0,0)" },
                      sprintf("%.1f", $fit));
          }
          push @td, td({ -bgcolor => Utils::fitcolor($fit) }, $show);
        }
        my @hidden2 = grep { $_->[0] ne "e.$orgId" || $_->[1] ne $expName } @hidden;
        push @td, td(a({-href => "cmpbrowser.cgi?" . join("&", map $_->[0] . "=" . $_->[1], @hidden2),
                        -title => "Remove this experiment",
                        -style => $linkColorStyle }, $removeGlyph));
        push @trows, Tr(@td);
      }
      print br(), table( { cellspacing => 0, cellpadding => 3 }, @trows) . "\n";
    }

    # Add space at the bottom of each track
    print qq{<div style="width:100%; position: relative; height: 2.5em;"></div>},
  } # end loop to show tracks

  #print p(a({-href => $tracksURL }, "self"));
  # Form for adding more tracks
  my @orgOptions = ("");
  my %orgLabels = ("" => "Choose a genome:");
  my @orginfo = sort { $a->{genome} cmp $b->{genome} } values(%$orginfo);
  foreach my $hash (@orginfo) {
    my $orgId = $hash->{orgId};
    if ($orgId ne $anchorOrg) {
      push @orgOptions, $orgId;
      $orgLabels{$orgId} = $hash->{genome};
    }
  }
  print start_form(-name => 'input', -method => 'GET', -action => 'cmpbrowser.cgi'),
    $tracksHidden,
    p("Add genes by locus tags",
      textfield(-name => 'addRange', -size => 15, -maxlength => 50, -default => '', -override => 1),
      "or find orthologs in",
      popup_menu(-style => 'width: 12em',
                 -name => 'addOrg',
                 -values => \@orgOptions, -labels => \%orgLabels,
                 -default => $orgOptions[0], -override => 1),
      submit(-style => 'text-align: left; float: none; display:inline; padding: 2px 2px;',
              -value => "Add genes", -name => "")),
    end_form;
  print p("Orthologs are searched against", scalar(@anchorGenes), "anchor genes from",
           a({-href => "org.cgi?orgId=$anchorOrg"}, $orginfo->{$anchorOrg}{genome}) . ".",
          "Or", a({-href => "$tracksURL&resetAnchor=1" }, "reset anchor genes"),
          "to the top track."
         )
    if @tracks > 0;
  print q{</div>}; # end of content
  $dbh->disconnect();
  Utils::endHtml($cgi);
}

sub GetGene($$$) {
  my ($dbh, $orgId, $locusId) = @_;
  die "Empty locusId in $orgId" if !defined $locusId || $locusId eq "";
  my $gene = $dbh->selectrow_hashref("SELECT * from Gene WHERE orgId = ? AND locusId = ?",
                                     {}, $orgId, $locusId);
  die "Invalid locus $locusId in $orgId" unless defined $gene;
  return $gene;
}


sub AddGeneToTracks($$) {
  my ($tracks, $gene) = @_;
  my $mid = ($gene->{begin} + $gene->{end}) / 2;
  foreach (my $i = 0; $i < scalar(@$tracks); $i++) {
    my $track = $tracks->[$i];
    if ($track->{orgId} eq $gene->{orgId}
        && $track->{geneBeg}{scaffoldId} eq $gene->{scaffoldId}
        && $track->{geneBeg}{scaffoldId} eq $gene->{scaffoldId}
        && (abs($track->{geneBeg}{begin} - $mid) < $newTrackDistance
            || abs($track->{geneEnd}{end} - $mid) < $newTrackDistance)) {
      # Modify the track to include this gene (if needed)
      # In case of overlapping genes -- the logic used to load genes
      # uses minEnd and maxBeg, so maintain those
      if ($gene->{end} < $track->{geneBeg}{end}) {
        $track->{geneBeg} = $gene;
      } elsif ($gene->{begin} > $track->{geneEnd}{begin}) {
        $track->{geneEnd} = $gene;
      } else {
        return 0; # no change needed
      }
      return 1;
    }
  }
  # New track needed
  push @$tracks, { orgId => $gene->{orgId},
                   geneBeg => $gene,
                   geneEnd => $gene,
                   strand => $gene->{strand} };
  return 1;
}

