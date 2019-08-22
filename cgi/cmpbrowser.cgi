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
#
# Arguments for modifying the view:
# addOrg -- (select orthologous genes of the anchor(s) from this organism
# addRange -- a locusId or sysName, or left:right
# TODO implement:
# flipTrack (0-based) -- changes which strand is used for that track
# removeTrack -- a track to remove
# upTrack -- a track to move up
# downTrack -- a track to move down
# changeOrg and changeLeft  (+1 for another gene, -1 for fewer genes)
# changeOrg and changeRight (+1 for another gene, -1 for fewer genes)

use strict;
use CGI qw(-nosticky :standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use lib "../lib";
use Utils;
use List::Util qw{min max};
use HTML::Entities;

sub GetGene($$$);
sub AddGeneToTracks($$);

my $maxNGenes = 50; # on one track
my $newTrackDistance = 10*1000; # how far a gene is before we put it on a new track
my $kbWidth = 150; # svg units/kb
my $trackHeight = 50; # units per track
my $barHeight = $trackHeight * 0.7; # height of scale bar region
my $arrowSize = 30; # size of the arrow part of a gene
my $padding = 20; # at left and right

{
  my $cgi=CGI->new;
  my $dbh = Utils::get_dbh();
  my $orginfo = Utils::orginfo($dbh);

  my $anchorOrg = param('anchorOrg') || die "Must specify anchorOrg";
  die "Invalid organism $anchorOrg" unless exists $orginfo->{$anchorOrg};
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
        unless $geneBeg->{scaffold} == $geneEnd->{scaffold};
      push @tracks, { orgId => $orgId,
                      geneBeg => $geneBeg,
                      geneEnd => $geneEnd,
                      strand => $trackStrand[$i] eq "+" || $trackStrand[$i] eq "1" ? "+" : "-" };
    }
  }

  my @warnings = ();
  # Add new tracks?
  if (param('addOrg')) {
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
    my $addSpec = param('addRange');
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

  # Save the values to use for a URL or hidden form elements
  my @hidden = (['anchorOrg', $anchorOrg]);
  push @hidden, map ['anchorLoci', $_], @anchorLoci;
  foreach my $track (@tracks) {
    push @hidden, ['o', $track->{orgId}],
      ['b', $track->{geneBeg}{locusId}],
      ['e', $track->{geneEnd}{locusId}],
      ['s', $track->{strand} eq "+" ? 1 : 0];
  }
  my $trackURL = "cmpbrowser.cgi?" . join("&", map $_->[0] . "=" . $_->[1], @hidden);
  my $trackHidden = join("\n",
                         map qq{<input type="hidden" name="$_->[0]" value="$_->[1]">}, @hidden);

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
    print p(a({-href => $URL, -title => "see fitness data", style => "color: black;" },
              i($g->{genus}, $g->{species}), $g->{strain}));

    # svg shows the genes in order
    my $xmin = min(map $_->{begin}, @$genes);
    my $xmax = max(map $_->{end}, @$genes);
    my $xdiff = $xmax - $xmin;
    # min. 1 kb for scale bar
    my $svg_width = 2 * $padding + max($xdiff, 1000) * $kbWidth / 1000.0;

    # SVG has +y axis going down, not up
    my $top = 0;
    my $genetop = $top;
    my $bottom = $trackHeight;
    my $geneymid = ($genetop+$bottom)/2;
    my $right = $padding + $xdiff * $kbWidth/1000.0;
    my $svgHeight = $trackHeight;
    $svgHeight += $barHeight if $bAddScale;

    print qq{<svg width="${svg_width}" height="${svgHeight}">\n};
    print qq{<line x1="$padding" y1="$geneymid" x2="$right" y2="$geneymid" style="stroke:black; stroke-width:1;"/>\n};
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
      my $color = "grey"; #TODO implement colors $genecolor{$track->{orgId}}{$gene->{locusId}} || "white";
      my $poly = qq{<polygon points="$pointstr" style="fill:$color; stroke:black; stroke-width:1;" />};

      my $showId = $gene->{gene} || $gene->{sysName} || $gene->{locusId};
      $showId =~ s/^.*_/_/ if defined $showId;
      my $xlabel = ($xstart+$xstop)/2;
      my $label = qq{<text x="$xlabel" y="$geneymid" fill="black" alignment-baseline="middle" text-anchor="middle">$showId</text>};
      my $URL = encode_entities( "myFitShow.cgi?orgId=$orgId&gene=$gene->{locusId}" );
      print qq{<a xlink:href="$URL">};
      print "<title>" . encode_entities( Utils::gene_link($dbh, $gene, "text") ) . "</title>\n";
      print $poly, $label;
      print qq{</a>\n};
    }
    if ($iTrack == scalar(@tracks) - 1) { # add scale bar for 1 kb at bottom
      my @bary = ($trackHeight + $barHeight * 0.7, $trackHeight + $barHeight * 0.95);
      my $barAt = ($bary[0] + $bary[1])/2;
      my $barright = $padding + $kbWidth;
      print qq{<line x1="$padding" y1="$barAt" x2="$barright" y2="$barAt" style="stroke:black; stroke-width:1"/>\n};
      foreach my $x ($padding, $barright) {
        print qq{<line x1="$x" y1="$bary[0]" x2="$x" y2="$bary[1]" style="stroke:black; stroke-width:1"/>\n};
      }
      my $barcenter = ($padding + $barright)/2;
      my $barAt = $barAt - $barHeight/10;
      print qq{<text x="$barcenter" y="$barAt">1 kb</text>\n};
    }
    print qq{</svg>\n};
  } # end loop over tracks
  #print p(a({-href => $trackURL }, "self"));
  # Form for adding more tracks
  my @orgOptions = ("");
  my %orgLabels = ("" => "Choose an organism:");
  my @orginfo = sort { $a->{genome} cmp $b->{genome} } values(%$orginfo);
  foreach my $hash (@orginfo) {
    my $orgId = $hash->{orgId};
    if ($orgId ne $anchorOrg) {
      push @orgOptions, $orgId;
      $orgLabels{$orgId} = $hash->{genome};
    }
  }
  print start_form(-name => 'input', -method => 'GET', -action => 'cmpbrowser.cgi'),
    $trackHidden,
    p("Add genes by locus tag",
      textfield(-name => 'addRange', -size => 15, -maxlength => 50, -default => '', -override => 1),
      "or find orthologs in",
      popup_menu(-style => 'width: 15em',
                 -name => 'addOrg',
                 -values => \@orgOptions, -labels => \%orgLabels,
                 -default => $orgOptions[0], -override => 1),
      submit(-style => 'text-align: left; float: none; display:inline;',
              -value => "Go", -name => "Go")),
    end_form;
  print q{</div>};
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

