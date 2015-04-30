#!/usr/bin/perl -w
#######################################################
## compareExps.cgi -- interactive scatterplot for two experiments
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required parameters: orgId, expName1, expName2
# Optional: tsv -- use tsv=1 to fetch the data instead

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
my $expName1 = $cgi->param('expName1');
my $expName2 = $cgi->param('expName2');
my $tsv = $cgi->param('tsv') ? 1 : 0;
die "Must specify orgId, expName1, expName2" unless defined $orgId && defined $expName1 && defined $expName2;

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
Utils::fail($cgi, "Unknown organism: $orgId") unless exists $orginfo->{$orgId};

my $exp1 = $dbh->selectrow_hashref("SELECT * from Experiment WHERE orgId = ? AND expName = ?",
				   {}, $orgId, $expName1);
die "Unknown experiment: $expName1" unless exists $exp1->{expName};

my $exp2 = $dbh->selectrow_hashref("SELECT * from Experiment WHERE orgId = ? AND expName = ?",
				   {}, $orgId, $expName2);
die "Unknown experiment: $expName2" unless exists $exp2->{expName};

if ($tsv) {
    # fetch the data
    my $genes = $dbh->selectall_arrayref("SELECT * FROM Gene where orgId = ?",
					 { Slice => {} }, $orgId);
    die "No genes" unless @$genes > 0;

    my $fit1 = $dbh->selectall_hashref("SELECT locusId,fit,t from GeneFitness WHERE orgId = ? AND expName = ?",
				       "locusId", {}, $orgId, $expName1);
    Utils::fail("No fitness values for $expName1 in $orgId") unless scalar(keys %$fit1) > 0;
    my $fit2 = $dbh->selectall_hashref("SELECT locusId,fit,t from GeneFitness WHERE orgId = ? AND expName = ?",
				       "locusId", {}, $orgId, $expName2);
    Utils::fail("No fitness values for $expName2 in $orgId") unless scalar(keys %$fit2) > 0;
    print join("\t", qw{locusId sysName gene desc x tx y ty});
    foreach my $gene (@$genes) {
	my $locusId = $gene->{locusId};
	next unless exists $fit1->{$locusId} && exists $fit2->{$locusId};
	print join("\t", $locusId, $gene->{sysName}, $gene->{gene}, $gene->{desc},
		   $fit1->{$locusId}{fit}, $fit1->{$locusId}{t},
		   $fit2->{$locusId}{fit}, $fit2->{$locusId}{t})."\n";
    }
    $dbh->disconnect();
    exit 0;
}
# else

my $style = Utils::get_style();
my $title = "Compare Experiments for $orginfo->{$orgId}{genome}";

print <<END
<!DOCTYPE html>
<head>
<title>$title</title>
<meta name="copyright" content="copyright 2015 UC Berkeley" />

<style>
$style
</style>

<script src="../d3js/d3.min.js"></script>
<body>

<H2>$title</H2>

<H4><A HREF="exp.cgi?orgId=$orgId&expName=$expName1">$expName1</A>: $exp1->{expDescLong}</H3>
<H4><A HREF="exp.cgi?orgId=$orgId&expName=$expName2">$expName2</A>: $exp2->{expDescLong}</H3>

<div id="loading"></div>

<script>
var org = "$orgId";
var xName = "$expName1";
var yName = "$expName2";

var margin = {top: 20, right: 20, bottom: 30, left: 40},
    width = 500 - margin.left - margin.right,
    height = 500 - margin.top - margin.bottom;

var x = d3.scale.linear()
    .range([0, width]);

var y = d3.scale.linear()
    .range([height, 0]);

var color = d3.scale.category10();

var xAxis = d3.svg.axis()
    .scale(x)
    .orient("bottom");

var yAxis = d3.svg.axis()
    .scale(y)
    .orient("left");

var svg = d3.select("body").insert("svg")
    .attr("width",500)
    .attr("height",500)
  .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

d3.select("#loading").html("Loading...");
var tsvUrl = "compareExps.cgi?tsv=1&orgId=" + org + "&expName1=" + xName + "&expName2=" + yName;
d3.tsv(tsvUrl, function(error, data) {
  if (error) {
      d3.select("body").append("div").html("Cannot load data from " + tsvUrl + "<BR>Error: " + error);
      return;
  }
  d3.select("#loading").html("Formatting " + data.length + " genes...");
  data.forEach(function(d) {
    d.x = +d.x;
    d.y = +d.y;
  });

  var extentX = d3.extent(data, function(d) { return d.x; });
  var extentY = d3.extent(data, function(d) { return d.y; });
  var extentXY = d3.extent([ extentX[0], extentX[1], extentY[0], extentY[1] ]);
  x.domain(extentXY).nice();
  y.domain(extentXY).nice();

  svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis)
    .append("text")
      .attr("class", "label")
      .attr("x", width)
      .attr("y", -6)
      .style("text-anchor", "end")
      .text(xName);

  svg.append("g")
      .attr("class", "y axis")
      .call(yAxis)
    .append("text")
      .attr("class", "label")
      .attr("transform", "rotate(-90)")
      .attr("y", 6)
      .attr("dy", ".71em")
      .style("text-anchor", "end")
      .text(yName);

  
  svg.append("line")
       .attr("x1", x(extentXY[0]))
       .attr("x2", x(extentXY[1]))
       .attr("y1", y(extentXY[0]))
       .attr("y2", y(extentXY[1]))
       .style("stroke","darkgrey")
       .style("stroke-width",1);

  svg.selectAll(".dot")
      .data(data)
    .enter().append("circle")
      .attr("class", "dot")
      .attr("r", 2)
      .attr("cx", function(d) { return x(d.x); })
      .attr("cy", function(d) { return y(d.y); })
      .on("click", dotClick);
  // .style("fill", function(d) { return color(d.species); });

  d3.select("#loading").html("");

  /*
  var legend = svg.selectAll(".legend")
      .data(color.domain())
    .enter().append("g")
      .attr("class", "legend")
      .attr("transform", function(d, i) { return "translate(0," + i * 20 + ")"; });
  legend.append("rect").attr("x", width - 18).attr("width", 18).attr("height", 18).style("fill", color);
  legend.append("text").attr("x", width - 24).attr("y", 9).attr("dy", ".35em").style("text-anchor", "end").text(function(d) { return d; });
  */
});

function dotClick(d) {
    d3.select(this).style("fill","green");
    //d3.select("#genesel").append("tr").html("<td>gene</td><td>name</td><td>descasdfsadf</td>");
    columns = [ d.sysName, d.gene, d.desc ];
    var i = d3.select("#genesel").append("tr");
    var showId = d.sysName === "" ? d.locusId : d.sysName;
    var URL = "myFitShow.cgi?orgId=" + org + "&gene=" + d.locusId;
    var beginHref = "<A target='_blank' HREF='" + URL + "'>";
    i.append("td").attr("class","locusId").attr("locusId",d.locusId).html(beginHref + showId + "</A>");
    i.append("td").html(d.gene);
    i.append("td").html(d.desc);
    i.append("td").html("<a href='#' onclick='removeRow(this)'>remove</a>");
}

function removeRow(a) {
    row = a.parentNode.parentNode; // href to td to row
    row.parentNode.removeChild(row);
}

function geneList() {
    var tds = document.getElementsByClassName("locusId");
    if (tds.length > 0) {
	var URL;
	if (tds.length == 1) {
	    URL = "myFitShow.cgi?orgId=" + org + "&gene=" + tds[0].getAttribute("locusId");
	} else {
            var i;
	    URL = "genesFit.cgi?orgId=" + org;
            for (i = 0; i < tds.length; i++ ) {
                URL += "&locusId=" + tds[i].getAttribute("locusId");
            }
	}
        window.open(URL);
   }
}

</script>

<H4>Selected genes:</H4>

<TABLE id="genesel" cellspacing=0 cellpadding=3 >
<tr><th>gene</th><th>name</th><th>description</th><th>&nbsp;</th></tr>
</TABLE>

<P><A href="#" onclick="geneList()">Heatmap for selected genes</A>

</body>
</html>
END
;

