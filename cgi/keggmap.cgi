#!/usr/bin/perl -w

#######################################################
## keggmap.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Show a KEGG map
#
# CGI parameters -- mapId, such as map00010 or 00010

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

sub LoadKEGGMap($$);
sub DrawKEGGMap($);

my $cgi=CGI->new;
my $mapId = $cgi->param('mapId');

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);

die "Must specify mapId" unless defined $mapId;
$mapId =~ s/^map//;
die "Invalid mapId" unless $mapId =~ m/^\d+$/;
# The map is a complex data structure including:
# mapId
# links -- a list of [coord, url] -- these are the non-enzyme links only
# ecnums -- a list of hash of ecnum, ecdesc, coord, url, optional color
my $map = &LoadKEGGMap($dbh, $mapId);

my $title = $map->{mapdesc};
my $start = Utils::start_page($title);

print
    header,
    $start, '<div id="ntcontent">',
    h2($title);

&DrawKEGGMap($map);

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
            $ecobj->{color} = "blue" if $objectId eq "2.7.1.41"; # XXX
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
        if (defined $row->{color}) {
            $bg = "rgba(0,0,255,0.6)";
        }
        my $alt = $row->{ecdesc};
        my $url = $row->{url};
        print qq{<A style="position:absolute; top:${top}px; left:${left}px; width:${width}px; height:${height}px; background-color: ${bg};" title="$alt" href="$url"></A>\n};
    }
    print "</DIV>\n"; # close the div that contains the image
}
