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
#	outlier -- list outlying genes (xlow, xhigh, ylow, or yhigh)
#       with minabs -- minimum |abs| on selected axis.

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
my $outlier = $cgi->param('outlier');
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

my $genes; # locusId => genes => attribute, with additional values x, y, tx, ty

if ($tsv || $outlier) {
    # fetch the data
    $genes = $dbh->selectall_hashref("SELECT * FROM Gene where orgId = ?", "locusId", {}, $orgId);
    die "No genes" unless scalar(keys %$genes) > 0;

    # Tried adding an index to GeneFitness so we could look up by orgId and expName but it did not really
    # speed things up. Probably the rows are in gene order so the whole part of the table for that orgId has to
    # be scanned anyway. Do one query for both experiments so that it is scanned once not twice.
    my $fit = $dbh->selectall_arrayref("SELECT * FROM GeneFitness WHERE orgId = ? AND expName IN (?,?)",
				       { Slice => {} }, $orgId, $expName1, $expName2);
    my $found1 = 0;
    my $found2 = 0;
    foreach my $row (@$fit) {
	my $locusId = $row->{locusId};
	die "Unrecognized locus $locusId for org $orgId" unless exists $genes->{$locusId};
	my $gene = $genes->{$locusId};
	if ($row->{expName} eq $expName1) {
	    $gene->{x} = $row->{fit};
	    $gene->{tx} = $row->{t};
	    $found1 = 1;
	}
	if ($row->{expName} eq $expName2) {
	    $gene->{y} = $row->{fit};
	    $gene->{ty} = $row->{t};
	    $found2 = 1;
	}
    }
    Utils::fail($cgi, "No fitness values for $expName1 in $orgId") unless $found1 > 0;
    Utils::fail($cgi, "No fitness values for $expName2 in $orgId") unless $found2 > 0;
}
$dbh->disconnect();

if ($tsv) { # tab delimited values, not a page
    print join("\t", qw{locusId sysName gene desc x tx y ty})."\n";
    while (my ($locusId,$gene) = each %$genes) {
	next unless exists $gene->{x} && exists $gene->{y};
	print join("\t", $locusId, $gene->{sysName}, $gene->{gene}, $gene->{desc},
		   $gene->{x}, $gene->{tx}, $gene->{y}, $gene->{ty})."\n";
    }
    exit 0;
} elsif ($outlier) { # table of outlying genes
    my $minabs = $cgi->param('minabs');
    $minabs = 2 unless defined $minabs && $minabs > 0;

    my $outlierCode = "";
    if ($outlier eq "lowx") {
	$outlierCode = "exp1 &lt; -$minabs and exp2 &gt; -1.0";
    } elsif ($outlier eq "lowy") {
	$outlierCode = "exp1 &gt; -1.0 and exp2 &lt; -$minabs";
    } elsif ($outlier eq "highx") {
	$outlierCode = "exp1 &gt; $minabs and exp2 &lt; 1.0";
    } elsif ($outlier eq "highy") {
	$outlierCode = "exp1 &lt; 1.0 and exp2 &gt; $minabs";
    } else {
	die "Unrecognized code for outlier";
    }
    $outlierCode .= " and &vert;<i>x</i>-<i>y</i>&vert; &gt; 1.0" if $minabs < 2;
    $outlierCode =~ s!exp1!<i>x</i>!g;
    $outlierCode =~ s!exp2!<i>y</i>!g;

    my @genesShow = ();
    while (my ($locusId,$gene) = each %$genes) {
	next unless exists $gene->{x} && exists $gene->{y};
	my $x = $gene->{x};
	my $y = $gene->{y};
	my $diff = abs($x-$y);
	push @genesShow, $gene
	    if $diff > 1
	    && (($outlier eq "lowx" && $x < -$minabs && $y > -1)
		|| ($outlier eq "lowy" && $y < -$minabs && $x > -1)
		|| ($outlier eq "highx" && $x > $minabs && $y < 1)
		|| ($outlier eq "highy" && $y > $minabs && $x < 1));

    }

    if ($outlier eq "lowx") {
	@genesShow = sort { $a->{x} <=> $b->{x} } @genesShow;
    } elsif ($outlier eq "lowy") {
	@genesShow = sort { $a->{y} <=> $b->{y} } @genesShow;
    } elsif ($outlier eq "highx") {
	@genesShow = sort { $b->{x} <=> $a->{x} } @genesShow;
    } elsif ($outlier eq "highy") {
	@genesShow = sort { $b->{y} <=> $a->{y} } @genesShow;
    }

    print
	start_html( -title => "Outlier genes from $orginfo->{$orgId}{genome}", -style => {-code => $style},
		    -author => 'Morgan Price', -mata => {'copyright'=>'copyright 2015 UC Berkeley'} ),
	h2("Outlier genes from $orginfo->{$orgId}{genome}"),
	h3($outlierCode),
	p(qq{<i>x</i> is <A HREF="exp.cgi?orgId=$orgId&expName=$expName1">$expName1</A>: $exp1->{expDescLong} }
	  . "<BR>"
	  . qq{<i>y</i> is <A HREF="exp.cgi?orgId=$orgId&expName=$expName2">$expName2</A>: $exp2->{expDescLong} }),
	p(scalar(@genesShow) . " genes found");
    if (@genesShow > 0) {
	my @trows = ();
	my @headings = qw{gene name description x y};
	push @trows, $cgi->Tr({-align=>'center',-valign=>'top'}, $cgi->th(\@headings));
	foreach my $gene (@genesShow) {
	    push @trows, $cgi->Tr({-align=>'left',-valign=>'top'},
		                  $cgi->td($cgi->a({href => "myFitShow.cgi?orgId=$orgId&gene=$gene->{locusId}",
						    style => "color:rgb(0,0,0)"},
						   $gene->{sysName} || $gene->{locusId})),
				  $cgi->td($gene->{gene}),
				  $cgi->td($gene->{desc}),
				  $cgi->td({ -bgcolor => Utils::fitcolor($gene->{x}) },
					   $cgi->a({title => sprintf("t = %.1f", $gene->{tx})},
						   sprintf("%.1f", $gene->{x}))),
				  $cgi->td({ -bgcolor => Utils::fitcolor($gene->{y}) },
					   $cgi->a({title => sprintf("t = %.1f", $gene->{ty})},
						   sprintf("%.1f", $gene->{y}))) );
	}
	my $limitString = "";
	if (@genesShow > 20) {
	    @genesShow = @genesShow[0..19];
	    $limitString = "top 20";
	}
	my $heatURL = "genesFit.cgi?orgId=$orgId&" . join("&", map { "locusId=" . $_->{locusId} } @genesShow);

	print
	    table({cellpadding=>3, cellspacing=>0}, @trows),
	    p(a({href => $heatURL}, "Heatmap for $limitString genes"));
	    
    }
    print p(a({href => "compareExps.cgi?orgId=$orgId&expName1=$expName1&expName2=$expName2"},"Show scatterplot"));

    Utils::endHtml($cgi);
}
# else interactive scatterplot

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
<body style="padding-left: 1%">

<H2>$title</H2>

<P>
<i>x</i> axis <A HREF="exp.cgi?orgId=$orgId&expName=$expName1">$expName1</A>: $exp1->{expDescLong}</H3>
<BR>
<i>y</i> axis <A HREF="exp.cgi?orgId=$orgId&expName=$expName2">$expName2</A>: $exp2->{expDescLong}

<TABLE width=100% style="border: none;">
<TR>
<TD valign="top" align="left" style="border: none;"><!-- left column -->

<div id="left"><!-- where SVG goes -->
<div id="loading"><!-- where status text goes -->
Please try another browser if this message remains
</div>
</div>
</TD>
<TD valign="top" align="left" style="border: none;"><!-- right column -->
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

<P>Click on genes to add them to the table:

<TABLE id="genesel" cellspacing=0 cellpadding=3 >
<tr><th>gene</th><th>name</th><th>description</th><th>x</th><th>y</th><th>&nbsp;</th></tr>
</TABLE>
</P>

<P>
<form method="get" action="compareExps.cgi" enctype="multipart/form-data" name="input">
<input type="hidden" name="orgId" value="$orgId" />
<input type="hidden" name="expName1" value="$expName1" />
<input type="hidden" name="expName2" value="$expName2" />
Or list outliers with
<select name="outlier">
   <option value="lowx">Low <i>x</i></option>
   <option value="lowy">Low <i>y</i></option>
   <option value="highx">High <i>x</i></option>
   <option value="highy">High <i>y</i></option>
</select>
and |fit| &gt; <select name="minabs" style="width: 60px;">
    <option value="1" selected>1.0 </option>
    <option value="1.5">1.5</option>
    <option value="2">2.0</option>
    <option value="2.5">2.5</option>
    <option value="3">3.0</option>
</select>
<input type="submit" name="submit" value="Go">
</form>

<P>
<A href="#" onclick="geneList()">View heatmap for selected genes</A>
</TD></TR></TABLE>
</P>

<script>
var org = "$orgId";
var xName = "$expName1";
var yName = "$expName2";
var xDesc = "$exp1->{expDesc}";
var yDesc = "$exp2->{expDesc}";

var margin = {top: 20, right: 20, bottom: 50, left: 50},
    width = 500 - margin.left - margin.right,
    height = 500 - margin.top - margin.bottom;

var x = d3.scale.linear()
    .range([0, width]);

var y = d3.scale.linear()
    .range([height, 0]);

//var color = d3.scale.category10();

var xAxis = d3.svg.axis()
    .scale(x)
    .orient("bottom");

var yAxis = d3.svg.axis()
    .scale(y)
    .orient("left");

var iSelected = 0; /* for color coding */
var selectColors = [ 'red', 'green', 'blue', 'magenta', 'brown', 'orange', 'darkturquoise' ];

var svg = d3.select("#left").append("svg")
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
    d.tx = +d.tx;
    d.ty = +d.ty;
  });

  var extentX = d3.extent(data, function(d) { return d.x; });
  var extentY = d3.extent(data, function(d) { return d.y; });
  var extentXY = d3.extent([ extentX[0], extentX[1], extentY[0], extentY[1] ]);
  x.domain(extentXY).nice();
  y.domain(extentXY).nice();

  svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis);

  svg.append("text")
      .attr("class", "label")
      .attr("x", 350)
      .attr("y", 500-25)
      .style("text-anchor", "end")
      .text(xDesc);

  svg.append("g")
      .attr("class", "y axis")
      .call(yAxis);

  svg.append("text")
      .attr("class", "label")
      .attr("transform", "rotate(-90)")
      .attr("x", -80)
      .attr("y", -30)
      .style("text-anchor", "end")
      .text(yDesc);

  
  svg.append("line")
       .attr("x1", x(extentXY[0]))
       .attr("x2", x(extentXY[1]))
       .attr("y1", y(extentXY[0]))
       .attr("y2", y(extentXY[1]))
       .style("stroke","darkgrey")
       .style("stroke-width",1);

  svg.append("line")
       .attr("x1", x(extentXY[0]))
       .attr("x2", x(extentXY[1]))
       .attr("y1", y(0))
       .attr("y2", y(0))
       .style("stroke","darkgrey")
       .style("stroke-width",1);

  svg.append("line")
       .attr("x1", x(0))
       .attr("x2", x(0))
       .attr("y1", y(extentXY[0]))
       .attr("y2", y(extentXY[1]))
       .style("stroke","darkgrey")
       .style("stroke-width",1);

  svg.selectAll(".dot")
      .data(data)
    .enter().append("circle")
      .attr("class", "dot")
      .attr("r", 3)
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
    var col = selectColors[(iSelected++) % selectColors.length];
    d3.select(this).style("fill",col).attr("r",5);
    columns = [ d.sysName, d.gene, d.desc ];
    var tr = d3.select("#genesel").append("tr").attr("valign","top").style("color", col);
    var showId = d.sysName === "" ? d.locusId : d.sysName;
    var URL = "myFitShow.cgi?orgId=" + org + "&gene=" + d.locusId;
    var beginHref = "<A target='_blank' style='color: " + col + "' HREF='" + URL + "'>";
    tr.append("td").attr("class","locusId").attr("locusId",d.locusId).html(beginHref + showId + "</A>");
    tr.append("td").html(d.gene);
    tr.append("td").html(d.desc);
    tr.append("td").html("<A TITLE='t = " + d.tx.toFixed(1) + "'>" + d.x.toFixed(1) + "</A>");
    tr.append("td").html("<A TITLE='t = " + d.ty.toFixed(1) + "'>" + d.y.toFixed(1) + "</A>");
    tr.append("td").html("<button type='button' onclick='removeRow(this)'>remove</button>");
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

<P><A HREF="myFrontPage.cgi">Go back to front page</A>
</body>
</html>
END
;

