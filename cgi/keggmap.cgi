#!/usr/bin/perl -w

#######################################################
## keggmap.cgi
##
## Copyright (c) 2016 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Show a KEGG map
#
# CGI parameters -- mapId, such as map00010 or 00010
# optional: orgId, to specify which enzymes to mask out ans which to search for members in

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

# dbh, kegg id => map data structure, including
# mapId
# mapdesc
# links -- a list of [coord, url] -- these are the non-enzyme links only (as they do not get shaded)
# ecnums -- a list of hash of ecnum, ecdesc, coord, url, and optionally mask=1 (to grey out)
#	From this routine, mask is never set, and
#	the urls are either to kegg (ww.genome.jp) or they search for the EC #
sub LoadKEGGMap($$);

# map data structure => outputs HTML for the map image and overlays
sub DrawKEGGMap($);

# dbh, orgId, ref. to list of ec numbers => hashref of ecnum to locusId => 1
sub EcToGenes($$$);

my $cgi=CGI->new;
my $mapId = $cgi->param('mapId');
my $orgId = $cgi->param('orgId');
$orgId = "" if !defined $orgId;

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Unknown orgId $orgId" if $orgId && !exists $orginfo->{$orgId};

die "Must specify mapId" unless defined $mapId;
$mapId =~ s/^map//;
die "Invalid mapId" unless $mapId =~ m/^\d+$/;

my $map = &LoadKEGGMap($dbh, $mapId);
my %ecs = map { $_->{ecnum} => 1 } @{ $map->{ecnums} };
my @ecs = keys(%ecs);
if ($orgId) {
    # grey out absent EC#s, and link to actual genes if present
    my $ecGenes = EcToGenes($dbh, $orgId, \@ecs);
    foreach my $row (@{ $map->{ecnums} }) {
        my $ecnum = $row->{ecnum};
        if (exists $ecGenes->{$ecnum}) {
            my @locusSpecs = map "locusId=$_", sort keys %{ $ecGenes->{$ecnum} };
            if (@locusSpecs > 1) {
                $row->{url} = "genesFit.cgi?orgId=$orgId&" . join("&",@locusSpecs);
            } else {
                $row->{url} = "singleFit.cgi?orgId=$orgId&" . $locusSpecs[0];
            }
        } else {
            $row->{mask} = 1;
        }
    }
    # pathway links to remember orgId
    foreach my $row (@{ $map->{links} }) {
        # entry 1 is the url
        $row->[1] .= "&orgId=$orgId" if $row->[1] =~ m/^keggmap.cgi/;
    }
}

my $title = $map->{mapdesc};
my $title2 = $title;
if ($orgId) {
    $title .= " in " . $orginfo->{$orgId}{genome} ;
    $title2 .= " in " . a({href => "org.cgi?orgId=$orgId"}, $orginfo->{$orgId}{genome});
}
my $start = Utils::start_page($title);

print
    header,
    $start, '<div id="ntcontent">',
    h2($title2),
    p( start_form(-name => 'orgselect', -method => 'GET', -action => 'keggmap.cgi'),
       hidden( -name => 'mapId', -value => $mapId, -override => 1),
       "Select organism:",
       Utils::OrgSelector($orgId, $orginfo),
       "<button type='submit'>Go</button>",
       end_form ),
    p(a( {href => "keggmaplist.cgi?orgId=$orgId"}, "Or browse metabolic maps"));

&DrawKEGGMap($map);

print
    p("Enzyme classification numbers are greyed out if $orginfo->{$orgId}{genome} is not predicted to contain this reaction.",
      "These predictions are often incorrect, especially if the EC number is incomplete, as with '2.7.1.-'.",
      "Gene-enzyme associations are based on the last (2011) public release of the",
      a( {href => "http://www.genome.jp/kegg/"}, "Kyoto Encyclopedia of Genes and Genomes"),
      "(KEGG) and also on",
      a( {href => "http://www.jcvi.org/cgi-bin/tigrfams/index.cgi"}, "TIGRFAMS"),
      "and",
      a( {href => "http://theseed.org"}, "SEED") . ".")
    if $orgId ne "";


$dbh->disconnect();
Utils::endHtml($cgi);

# returns a hashref that includes mapId, mapdesc, links, ecnums
sub LoadKEGGMap($$) {
    my ($dbh,$mapid) = @_;
    my ($mapdesc) = $dbh->selectrow_array("SELECT title FROM KEGGMap where mapId = ?",
                                      {}, $mapId);
    die "Unknown map id $mapId" unless defined $mapdesc;
    my $mapobjs = $dbh->selectall_arrayref("SELECT objectId,type,coord,url FROM KEGGConf WHERE mapId = ?",
                                           {}, $mapId);
    my @links = ();
    my @ecnums = ();
    foreach my $row (@$mapobjs) {
        my ($objectId,$type,$coord,$url) = @$row;
        $url = "http://www.genome.jp$url"; # default is, link to KEGG
        if ($type == 1) { # ec number
            my ($ecdesc) = $dbh->selectrow_array("SELECT ecdesc FROM ECInfo WHERE ecnum = ?",
                                             {}, $objectId);
            $url = "myFitShow.cgi?gene=ec:$objectId" unless $objectId =~ m/-/;
            my $ecobj = { ecnum => $objectId, ecdesc => $ecdesc || $objectId,
                          coord => $coord, url => $url };
            push @ecnums, $ecobj;
        } elsif ($type == 0 || $type == 2) { # compounds or maps
            $url = "keggmap.cgi?mapId=$objectId" if $type == 2;
            push @links, [ $coord, $url ];
        } else {
            ; # ignore type 3 or above, i.e. reactions
        }
    }
    return { mapId => $mapId, mapdesc => $mapdesc,
             links => \@links, ecnums => \@ecnums };
}

sub DrawKEGGMap($) {
    my ($map) = @_;
    my $mapId = $map->{mapId};
    print
        qq{<DIV style="position: relative; left:0; top:0;">},
        qq{<IMG src="../kegg/maps/map$mapId.png" usemap="#imagemap" style="position: relative; top: 0; left: 0;">\n};
    print qq{<MAP name="imagemap" id="imagemap" />\n};
    foreach my $row (@{ $map->{links} }) {
        my ($coord, $url) = @$row;
        my ($coordtype, $pos) = split /:/, $coord;
        next unless $coordtype eq "rect" || $coordtype eq "circ";
        print qq{<AREA SHAPE="$coordtype" coords="$pos" href="$url" />\n};
    }
    print qq{</MAP>};

    foreach my $row (@{ $map->{ecnums}}) {
        my $coord = $row->{coord};
        my ($coordtype, $pos) = split /:/, $coord;
        next unless $coordtype eq "rect";
        my ($left,$top,$right,$bottom) = split /,/, $pos;
        die "Invalid coords $coord" unless defined $top;
        my $width = $right-$left+1;
        my $height = $bottom-$top+1;
        my $bg = "rgba(255,255,255,0)";
        if ($row->{mask}) {
            $bg = "rgba(0,0,0,0.3)"; # make it grey
            # $bg = "rgba(0,0,255,0.6)"; # make it blue
        }
        my $alt = $row->{ecdesc};
        my $url = $row->{url};
        print qq{<A style="position:absolute; top:${top}px; left:${left}px; width:${width}px; height:${height}px; background-color: ${bg};" title="$alt" href="$url"></A>\n};
    }
    print "</DIV>\n"; # close the div that contains the image
}

sub EcToGenes($$$) {
    my ($dbh,$orgId,$ecnums) = @_;
    return {} if @$ecnums == 0;
    my %ecGenes = (); # ecnum => locusId => 1
    my @ecspecs = map "'" . $_ . "'", @$ecnums;
    my $ecspec = join(",", map "'" . $_ . "'", @$ecnums);
    # match ECs by TIGRFam, by KEGG ortholog group, and by SEED annotation
    my $hits1 = $dbh->selectall_arrayref(
	qq{ SELECT ec, locusId FROM GeneDomain JOIN Gene USING (orgId,locusId)
            WHERE orgId = ? AND ec IN ( $ecspec ); },
        {}, $orgId);
    my $hits2 = $dbh->selectall_arrayref(
        qq{ SELECT ecnum, locusId FROM KgroupEC JOIN KEGGMember USING (kgroup)
            JOIN BestHitKEGG USING (keggOrg,keggId)
            WHERE orgId = ? AND ecnum IN ( $ecspec ); },
        {}, $orgId);
                                                 
    my $hits3 = $dbh->selectall_arrayref(
        qq{ SELECT num,locusId
		FROM SEEDClass JOIN SEEDAnnotation USING (orgId,locusId)
                WHERE orgId = ? AND num IN ( $ecspec ); },
        {}, $orgId);
    foreach my $hits ($hits1,$hits2,$hits3) {
        foreach my $row (@$hits) {
            my ($ec,$locusId) = @$row;
            $ecGenes{$ec}{$locusId} = 1;
        }
    }
    return \%ecGenes;
}
