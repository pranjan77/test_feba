#!/usr/bin/perl -w
#######################################################
## orthCond.cgi -- compare specific phenotypes for a given Group & Condition_1 across all bugs
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Given a group and condition, identify all specific sick genes,
# group them into ad hoc ortholog groups, and show the groups.
#
# Required CGI parameters: expGroup and condition1 (condition1 can be empty, but requires exact match)
# Optional: help -- 1 if on help/tutorial mode

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;
use URI::Escape;
use Time::HiRes qw(gettimeofday tv_interval);

sub RowForGene($$$$); # gene object, shade or not, which row (0 indexed) => the row
sub summaryRow($$$); # ortholog group, shade, first to include (1 indexed)

my $cgi=CGI->new;
my $style = Utils::get_style();

my $expGroup = $cgi->param('expGroup') || die "no expGroup parameter";
my $condition1 = $cgi->param('condition1');
die "no condition1 parameter" unless defined $condition1;

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my $debug = $cgi->param('debug');
my $start_time = gettimeofday() if $debug;
my $help = $cgi->param('help') || "";

my $speclist = $dbh->selectall_arrayref(qq{SELECT DISTINCT orgId, locusId, expName, expDesc
                                           FROM Experiment JOIN SpecificPhenotype USING (orgId,expName)
					   WHERE Experiment.expGroup = ? AND Experiment.condition_1 = ?},
					   {}, $expGroup, $condition1);
print STDERR sprintf("Specific sicks: %.3f seconds\n", tv_interval([$start_time])) if $debug;

my $title = "Specific Genes for $expGroup Experiments in $condition1 Across Organisms";
$title = "No Specific Genes" if (@$speclist == 0);
my $start = Utils::start_page("$title");

my $js =  q[<script type="text/javascript" src="../images/jquery-1.11.3.min.js"></script>
<script type="text/javascript">$(document).ready(function(){ 
    $('tr.header2').nextUntil('tr.header').hide(); // after each tr.header2, hide rows until the next header
    $('tr.header2').click(function(){
	// click on header2 to hide it and turn on elements until next header
        $(this).toggle();
        $(this).nextUntil('tr.header').css('display', function(i,v){
            return this.style.display === 'table-row' ? 'none' : 'table-row';
        });
    });
    $('tr.header3').click(function(){
        // click on header3 to hide it and turn off elements until next header. show the previous row (header2)
        $(this).toggle();
        $(this).closest('tr').prev().show();
        // $('tr.header2').show();
        $(this).nextUntil('tr.header').css('display', function(i,v){
            return this.style.display === 'table-row' ? 'none' : 'table-row';
        });
    });
    
});
    </script>];

print
    header, $start, $js, '<div id="ntcontent">', #$tabs,
  #   start_html( -title => $title, -style => { -code => $style }, -author => 'Morgan Price',
		# -meta => { 'copyright' => 'copyright 2015 UC Berkeley' }),
    h2($title);
    # div({-style => "float: right; vertical-align: top;"}, a({href => "help.cgi#specific"}, "Help"));

Utils::fail($cgi,"no genes with specific phenotypes were found for this condition") if (@$speclist == 0);

my %genes = (); # orgId => locusId => attributes from a row in the Gene table
# gene object also includes
# 'orth' => orgId => orthId
# 'fit' => expName => hash with fit, t
my %expDesc = (); # orgId => expName => expDesc
foreach my $row (@$speclist) {
	my ($orgId,$locusId,$expName,$expDesc) = @$row;
	die "Invalid $orgId" unless exists $orginfo->{$orgId};
	if (!exists $genes{$orgId}{$locusId}) {
	    my $gene = $dbh->selectrow_hashref(qq{SELECT * from Gene WHERE orgId = ? AND locusId = ?}, {}, $orgId, $locusId);
	    $gene->{orth} = Utils::get_orths($dbh,$orgId,$locusId);
	    $genes{$orgId}{$locusId} = $gene;
	}
	$expDesc{$orgId}{$expName} = $expDesc; # overwriting is OK
}

print STDERR sprintf("Fetch genes and orths: %.3f seconds\n", tv_interval([$start_time])) if $debug;

# Fetch the actual fitness values
while (my ($orgId, $geneHash) = each %genes) {
    my @expNames = keys %{ $expDesc{$orgId} };
    die if @expNames == 0;
    my $expNameSpec = join(",", map { "'" . $_ . "'" } @expNames);
    while (my ($locusId, $gene) = each %$geneHash) {
	$gene->{values} = $dbh->selectrow_arrayref(qq{ SELECT min(fit), max(fit), min(t), max(t) FROM GeneFitness
                                                       WHERE orgId = ? AND locusId = ? AND expName IN ($expNameSpec) },
                                                       {}, $orgId, $locusId);
        die "No fitness data for $orgId $locusId $expGroup $condition1" unless defined $gene->{values}[0];
    }
}

print STDERR sprintf("Fetch fitness: %.3f seconds\n", tv_interval([$start_time])) if $debug;

# Now, sort the genes to show by ortholog groups
# Put first the ones that have the most orthologs in this set.
# Then, put them with whatever is an ortholog of theirs.
#
# The first step is to count the orthologs in the nOrth field of each gene
my @genes = ();
while (my ($orgId, $geneHash) = each %genes) {
    while (my ($locusId, $gene) = each %$geneHash) {
	my $orth = $gene->{orth};
	my $nOrth = 0;
	while (my ($orthTax,$orthObj) = each %$orth) {
            # avoid a subtle bug -- exists $genes{$orthTax}{$orthId} would silently create
            # $genes{$orthTax} and screw up the each %genes loop so that an organism
            # could get processed twice.
	    $nOrth++ if exists $genes{$orthTax} && exists $genes{$orthTax}{$orthObj->{locusId}};
	}
	$gene->{nOrth} = $nOrth;
	push @genes, $gene;
    }
}
@genes = sort { $b->{nOrth} <=> $a->{nOrth} } @genes;

print STDERR sprintf("Fetch orthology: %.3f seconds\n", tv_interval([$start_time])) if $debug;

# And make ortholog groups
my @OGs = ();
my %og = (); # orgId => locusId => which OG if gene is in there
foreach my $gene (@genes) {
    my $orth = $gene->{orth};
    my $iOG = undef;
    # strictly speaking, this should be based on which is most similar, or #hits, in case of conflicts b/w OGs
    # punt on that for now
    while (my ($orthTax,$orthObj) = each %$orth) {
	   $iOG = $og{$orthTax}{$orthObj->{locusId}};
	   last if defined $iOG;
    }
    if (defined $iOG) {
        push @{ $OGs[$iOG] }, $gene;
    } else {
	   $iOG = scalar(@OGs);
	   push @OGs, [ $gene ];
    }
    $og{$gene->{orgId}}{$gene->{locusId}} = $iOG;
}

print STDERR sprintf("Make OGs: %.3f seconds\n", tv_interval([$start_time])) if $debug;

my @headings = ['&nbsp', 'Organism', 'Gene', 'Name', 'Description', 'Fitness (Lower)', 'Fitness (Upper)'];  #Experiment
my @trows = ( $cgi->Tr({ -valign => 'middle', -align => 'center' }, map { th($_) } @headings ) );

my $shade = 0;
my $group = 1;
my $row = 0;
my $singletonHead = 0;

foreach my $og (@OGs) {
    @$og = sort { ($b->{gene} ne "") <=> ($a->{gene} ne "") } @$og; # sort them so the ones with gene names come up first 
    foreach my $gene (@$og) {
        if (scalar(@$og) == 1 and $singletonHead == 0) {
            push @trows, qq[<th colspan="8"><center>Singletons</center></th>];
            $singletonHead = 1;
        }
        if ($row == 3 and scalar(@$og) > 4) {
            # my $color = $shade % 2 ? "#DDDDDD" : "#FFFFFF";
            push @trows, summaryRow($og, $shade, 4);
            push @trows, qq[<tr class="header3"><th colspan="8"><center><span>Collapse -</span></center></th></tr>];
        }
        push @trows, RowForGene($gene, $shade, $singletonHead == 1 ? "" : $group, $row);
        $row++;
    }
    $shade++;
    $group++;
    $row = 0;
}

if ($help == 1) {
        print qq[<BR><BR><div class="helpbox">
        <b><u>About this page:</u></b><BR><ul>
        <li>View genes involved the $expGroup experiments in $condition1 in all organisms. </li>
        <li>To get to this page, search for any experiment and click on the link at the bottom to view specific phenotypes across organisms.</li> 
        <li>Click on the summary rows to expand each section.</li>
        <li>Click on links (including the fitness values) to view more.</li>
        </ul></div>];
    }


print
    p("Genes with",
      a({ -href => "help.cgi#specific" },  "specific phenotypes"),
      "in $expGroup $condition1 are",
      font({ style => "color: rgb(120,120,120); font-weight: bold;"}, "grouped"),
      "by",
      a({ -href => "help.cgi#ortholog"}, "orthology")),
    table({cellpadding => 3, cellspacing => 0}, @trows);

    print "<BR>";

$dbh->disconnect();
Utils::endHtml($cgi);

sub summaryRow($$$) {
    my ($og, $shade, $firstToInclude) = @_; # firstToInclude is 1-based
    my $i = 1;
    my @row;
    my @showOrgs = ();

    foreach my $gene(@$og) {
        my $orgId = $gene->{orgId};
        my $genome = $orginfo->{$orgId}{genome};
        my $genomeShort = $genome;
        $genomeShort =~ s/^(.)\S+/$1./;
        my $locusId = $gene->{locusId};

        push @showOrgs, $cgi->a({href=>"orthFit.cgi?orgId=$orgId&locusId=$locusId"
                                 . "&expGroup=" . uri_escape($expGroup)
                                 . "&condition1=" . uri_escape($condition1),
                                 title=>"$gene->{desc}"},
                                $genomeShort)
            if $i >= $firstToInclude;
        $i++;
    }
    
    my $count = scalar(@$og) - ($firstToInclude-1);
    push @row, $cgi->Tr(
        {-class=>'header2', -valign => 'middle', -align => 'left', -bgcolor => $shade % 2 ? "#DDDDDD" : "#FFFFFF"},
        td(span({class=>"deco"}, a({title=>"Expand"}, '+'))),
        td({colspan=>"6"}, "$count more from " . join(", ", @showOrgs))
    );
    return @row;
}

# $row is 0-indexed; $group is the ortholog group number (1-indexed)
sub RowForGene($$$$) {
    my ($gene,$shade,$group,$row) = @_;
   
    # my $firstForGene = 1;
    my $orgId = $gene->{orgId};
    my $showId = $gene->{sysName} || $gene->{locusId};
    my $genome = $orginfo->{$orgId}{genome};
    my $genomeShort = $genome;
    $genomeShort =~ s/^(.)\S+/$1./;
    my $locusId = $gene->{locusId};

    my ($min,$max,$minT,$maxT) = @{ $gene->{values} };
    my $rowLabel = "";
    my $collapse = "";
    # q[-valign => 'middle', -align => 'left', -bgcolor => $shade % 2 ? "#DDDDDD" : "#FFFFFF"];
    if ($row == 0) {
        $rowLabel = a({style=>"color: #011B47", title=>"Ortholog group $group"},$group);
        $collapse = 'header';
        # q[-class=>'header', -valign => 'middle', -align => 'left', -bgcolor => $shade % 2 ? "#DDDDDD" : "#FFFFFF"];
    }
    my $orthFitURI = "orthFit.cgi?orgId=$orgId&locusId=$locusId"
        . "&expGroup=" . uri_escape($expGroup)
        . "&condition1=" . uri_escape($condition1);
    $orthFitURI .= "&help=1" if $help == 1;
    return $cgi->Tr(
        { -class=> $collapse, -valign => 'middle', -align => 'left', -bgcolor => $shade % 2 ? "#DDDDDD" : "#FFFFFF" },
        td($rowLabel),
        td( a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId},  -title => $genome },  $genomeShort)),
        td( a({ -href => "myFitShow.cgi?orgId=$orgId&gene=$locusId" }, $showId)),
        td( $gene->{gene}),
        td( $gene->{desc}),
        td( { -bgcolor => Utils::fitcolor($min), -style=>'text-align: center;' },
            a( { -title => sprintf("Click to compare (t = %.1f to %.1f)",$minT,$maxT),
                 -style => "color: rgb(0,0,0)",
                 -onMouseOver=>"this.style.color='#CC0024'",
                 -onMouseOut=>"this.style.color='#000000'",
                 -href => $orthFitURI },
               sprintf("%.1f",$min) ) ),
        td( { -bgcolor => Utils::fitcolor($max), -style=>'text-align: center;' },
            a( { -title => sprintf("Click to compare (t = %.1f to %.1f)",$minT,$maxT),
                 -style => "color: rgb(0,0,0)",
                 -onMouseOver=>"this.style.color='#CC0024'",
                 -onMouseOut=>"this.style.color='#000000'",
                 -href => $orthFitURI },
               sprintf("%.1f",$max) ) )
	    );
}
