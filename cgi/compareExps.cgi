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
# Key parameeters: orgId, expName1 or query1, expName2 or query2
#	If query1 is set, expName1 is ignored, and similarly for query2
#	Cannot query on both simultaneously, however
# Optional: tsv -- use tsv=1 to fetch the data instead

use strict;
use CGI qw(:standard Vars -nosticky);
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
my $query1 = $cgi->param('query1');
my $query2 = $cgi->param('query2');
my $tsv = $cgi->param('tsv') ? 1 : 0;
die "Must specify orgId" unless defined $orgId && $orgId ne "";
die "Must specify expName1 or query1"
    unless (defined $expName1 || defined $query1)
    && !($expName1 eq "" && $query1 eq "");
die "Must specify expName2 or query2"
    unless (defined $expName2 || defined $query2)
    && !($expName2 eq "" && $query2 eq "");
die "Cannot query both 1 and 2" if defined $query1 && defined $query2 && $query1 ne "" && $query2 ne "";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
Utils::fail($cgi, "Unknown organism: $orgId") unless exists $orginfo->{$orgId};

my @expCand = ();
my $choosing = undef;
if (defined $query1 && $query1 ne "") {
    my @exps = @{ Utils::matching_exps($dbh,$orgId,$query1) };
    Utils::fail($cgi, qq{No experiment matching "$query1"}) if @exps == 0;
    @exps = grep { $_->{expName} ne $expName2 } @exps if @exps > 1;
    if (@exps == 1) {
	$expName1 = $exps[0]{expName};
    } else {
	@expCand = @exps;
	$choosing = 1;
    }
} elsif (defined $query2 && $query2 ne "") {
    my @exps = @{ Utils::matching_exps($dbh,$orgId,$query2) };
    Utils::fail($cgi, qq{No experiment matching "$query2"}) if @exps == 0;
    @exps = grep { $_->{expName} ne $expName1 } @exps if @exps > 1;
    if (@exps == 1) {
	$expName2 = $exps[0]{expName};
    } else {
	@expCand = @exps;
	$choosing = 2;
    }
}

if (scalar(@expCand) > 0) {
    die "Cannot use tsv mode with queries" if $tsv;
    # show table of these experiments
    my $expNameConst = $choosing == 1 ? $expName2 : $expName1;
    my $expConst = $dbh->selectrow_hashref("SELECT * from Experiment WHERE orgId = ? AND expName = ?",
					   {}, $orgId, $expNameConst);
    die "Unknown experiment: $expNameConst" unless exists $expConst->{expName};
    my $notChoosing = $choosing == 1 ? 2 : 1;

    my @trows = ();
    my @headings = qw{&nbsp; name group condition description};
    push @trows, $cgi->Tr({-valign => 'top', -align => 'center'}, $cgi->th(\@headings));
    my $isFirst = 1;
    foreach my $exp (@expCand) {
	my $checked = $isFirst ? "CHECKED" : "";
	push @trows, $cgi->Tr({-valign => 'top', -align => 'left'},
			      $cgi->td([ qq{<input type="radio" name="expName$choosing" value="$exp->{expName}" $checked >},
					 $cgi->a({href => "exp.cgi?orgId=$orgId&$exp->{expName}"}, $exp->{expName}),
					 $exp->{expGroup}, $exp->{condition_1}, $exp->{expDesc} ]));
	$isFirst = 0;
    }
    print
	start_html( -title => "Select experiment to compare to", -style => {-code => $style},
		    -author => 'Morgan Price', -mata => {'copyright'=>'copyright 2015 UC Berkeley'} ),
	h2("Select experiment in $orginfo->{$orgId}{genome}"),
	p("Selected experiment will be compared to "
	  . a( { href => "exp.cgi?orgId=$orgId&expName=$expNameConst" }, $expConst->{expName} )
	  . " : $expConst->{expDescLong}"),
	start_form(-name => 'input', -method => 'GET', -action => 'compareExps.cgi'),
	hidden('orgId', $orgId),
	hidden("expName$notChoosing", $expNameConst),
	table( {cellpadding => 3, cellspacing => 0}, @trows),
	submit('Go'),
	end_form;
    Utils::endHtml($cgi);
}

# else

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

    # Tried adding an index to GeneFitness so we could look up by orgId and expName but it did not really
    # speed things up. Probably the rows are in gene order so the whole part of the table for that orgId has to
    # be scanned anyway. Do one query for both experiments so that it is scanned once not twice.
    my $fit = $dbh->selectall_arrayref("SELECT * FROM GeneFitness WHERE orgId = ? AND expName IN (?,?)",
				       { Slice => {} }, $orgId, $expName1, $expName2);
    my %fit1 = ();
    my %fit2 = ();
    foreach my $row (@$fit) {
	$fit1{$row->{locusId}} = $row if $row->{expName} eq $expName1;
	$fit2{$row->{locusId}} = $row if $row->{expName} eq $expName2;
    }
    Utils::fail($cgi, "No fitness values for $expName1 in $orgId") unless scalar(keys %fit1) > 0;
    Utils::fail($cgi, "No fitness values for $expName2 in $orgId") unless scalar(keys %fit2) > 0;

    print join("\t", qw{locusId sysName gene desc x tx y ty});
    foreach my $gene (@$genes) {
	my $locusId = $gene->{locusId};
	next unless exists $fit1{$locusId} && exists $fit2{$locusId};
	print join("\t", $locusId, $gene->{sysName}, $gene->{gene}, $gene->{desc},
		   $fit1{$locusId}{fit}, $fit1{$locusId}{t},
		   $fit2{$locusId}{fit}, $fit2{$locusId}{t})."\n";
    }
    $dbh->disconnect();
    exit 0;
}
# else

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

<P><A HREF="exp.cgi?orgId=$orgId&expName=$expName2">$expName2</A>: $exp2->{expDescLong}<BR>
vs. <A HREF="exp.cgi?orgId=$orgId&expName=$expName1">$expName1</A>: $exp1->{expDescLong}</H3>

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

d3.select("#loading").html("Fetching data...");
var tsvUrl = "compareExps.cgi?tsv=1&orgId=" + org + "&expName1=" + xName + "&expName2=" + yName;
d3.tsv(tsvUrl, function(error, data) {
  if (error || data.length == 0) {
      d3.select("#loading").html("Cannot load data from " + tsvUrl + "<BR>Error: " + error);
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

<p>
<form method="get" action="compareExps.cgi" enctype="multipart/form-data" name="input">
<input type="hidden" name="orgId" value="$orgId" />
<input type="hidden" name="expName2" value="$expName2" />
Change x axis: <input type="text" name="query1"  size="20" maxlength="100" />
</form>

<form method="get" action="compareExps.cgi" enctype="multipart/form-data" name="input">
<input type="hidden" name="orgId" value="$orgId" />
<input type="hidden" name="expName1" value="$expName1" />
Change y axis: <input type="text" name="query2"  size="20" maxlength="100" />
</form>

<form method="get" action="compareExps.cgi" enctype="multipart/form-data" name="input">
<input type="hidden" name="orgId" value="$orgId" />
<input type="hidden" name="expName1" value="$expName2" />
<input type="hidden" name="expName2" value="$expName1" />
<input type="submit" name="flip" value="Flip axes" />
</form>

</p>

<p>Click on genes to add them to the table:</p>

<TABLE id="genesel" cellspacing=0 cellpadding=3 >
<tr><th>gene</th><th>name</th><th>description</th><th>&nbsp;</th></tr>
</TABLE>

<P><A href="#" onclick="geneList()">View heatmap for selected genes</A>

</body>
</html>
END
;

