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

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use List::MoreUtils qw( minmax );

use lib "../lib";
use Utils;
sub RowsForGene($$$$); # gene object, shade or not => list of rows

my $cgi=CGI->new;
my $style = Utils::get_style();

my $expGroup = $cgi->param('expGroup') || die "no expGroup parameter";
my $condition1 = $cgi->param('condition1');
die "no condition1 parameter" unless defined $condition1;

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);

my $speclist = $dbh->selectall_arrayref(qq{SELECT DISTINCT orgId, locusId, expName, expDesc
                                           FROM Experiment JOIN SpecificPhenotype USING (orgId,expName)
					   WHERE Experiment.expGroup = ? AND Experiment.condition_1 = ?},
					   {}, $expGroup, $condition1);

my $title = "Specific Genes for $expGroup Experiments in $condition1 Across Organisms";
$title = "No Specific Genes" if (@$speclist == 0);
my $start = Utils::start_page("$title");
# my $tabs = Utils::tabsExp($dbh,$cgi,$orgId,$expName,$expGroup,$condition1,"specific");
#$(this).find('span').text(function(_, value){return value=='+'?'-':'+'});
#$(this).nextUntil('tr.header').slideToggle(100, function(){});});
#$('.header').click(function(){
        
    #     $(this).nextUntil('tr.header').css('display', function(i,v){
    #         return this.style.display === 'table-row' ? 'none' : 'table-row';
    #     });
    # });

# $('tr').hide().filter(function () {
#         return $(this).find('td[colspan]').length;
#     }).addClass('header').css('display', 'table-row').click(function () {
#         $(this).nextUntil('tr.header').css('display', function (i, v) {
#             return this.style.display === 'table-row' ? 'none' : 'table-row';
#         });
#     });

# $('.header').hide(function(){
#         $(this).find('span').text(function(_, value){return value=='+'?'-':'+'});
#         $(this).nextUntil('tr.header').slideToggle(100);
#     });
#     $('.header').click(function(){
#         $(this).find('span').text(function(_, value){return value=='-'?'+':'-'});
#         $(this).nextUntil('tr.header').slideToggle(100);
#     });

# $(this).find('span').text(function(_, value){return value=='Collapse -'?'Expand +':'Collapse -'});

# tr.header is the first for each gene
# tr.header2 is the expandable row
# tr.header3 is the row that collapses it back down
#
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
        // click on header3 to hide it and turn off elements until next header, but why show *all* header2s ?
        $(this).toggle();
        $('tr.header2').show();
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

# Fetch the actual fitness values
while (my ($orgId, $geneHash) = each %genes) {
    my @expNames = keys %{ $expDesc{$orgId} };
    die if @expNames == 0;
    my $expNameSpec = join(",", map { "'" . $_ . "'" } @expNames);
    while (my ($locusId, $gene) = each %$geneHash) {
	$gene->{fit} = {};
	my $fitrows = $dbh->selectall_arrayref(qq{ SELECT expName,fit,t FROM GeneFitness
                                                   WHERE orgId = ? AND locusId = ? AND expName IN ($expNameSpec) },
					       {}, $orgId, $locusId);
	foreach my $row (@$fitrows) {
	    my ($expName,$fit,$t) = @$row;
	    $gene->{fit}{$expName} = { 'fit' => $fit, 't' => $t };
	}
    }
}

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
	    $nOrth++ if exists $genes{$orthTax}{$orthObj->{locusId}};
	}
	$gene->{nOrth} = $nOrth;
	push @genes, $gene;
    }
}
@genes = sort { $b->{nOrth} <=> $a->{nOrth} } @genes;


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
            push @trows, qq[<tr class="header2"><th colspan="8"><center>Singletons</center></th></tr>];
            $singletonHead = 1;
        }
        if ($row == 3 and scalar(@$og) > 4) {
            # my $color = $shade % 2 ? "#DDDDDD" : "#FFFFFF";
            push @trows, summaryRow($og, $shade);
            push @trows, qq[<tr class="header3"><th colspan="8"><center><span>Collapse -</span></center></th></tr>];
        }
        if ($singletonHead == 1) {
            push @trows, RowsForGene($gene, $shade, "", $row);
        } else {
            push @trows, RowsForGene($gene, $shade, $group, $row);
        }
        $row++;
    }
    $shade++;
    $group++;
    $row = 0;
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

sub summaryRow($$) {
    my ($og, $shade) = @_;
    my $count = scalar(@$og) - 3;
    my $orgs = "$count more from ";
    my @row;
    my $ind = 0;
    my $len = scalar(@$og) - 1;

    foreach my $gene(@$og) {
        my $orgId = $gene->{orgId};
        my $genome = $orginfo->{$orgId}{genome};
        my $genomeShort = $genome;
        $genomeShort =~ s/^(.)\S+/\1./;
        my $locusId = $gene->{locusId};
       
        if ($ind > 2 and $ind < $len) {
            $orgs = $orgs . $cgi->a({href=>"orthFit.cgi?orgId=$orgId&locusId=$locusId&expGroup=$expGroup&condition1=$condition1", title=>"$gene->{desc}"},$genomeShort).  ", ";
        } elsif ($ind == $len) {
             $orgs = $orgs . $cgi->a({href=>"orthFit.cgi?orgId=$orgId&locusId=$locusId&expGroup=$expGroup&condition1=$condition1", title=>"$gene->{desc}"},$genomeShort);
        }
        $ind += 1;
    }
    
    push @row, $cgi->Tr(
        {-class=>'header2', -valign => 'middle', -align => 'left', -bgcolor => $shade % 2 ? "#DDDDDD" : "#FFFFFF"},
        td(span({class=>"deco"}, a({title=>"Expand"}, '+'))),
        td({colspan=>"6"}, $orgs)
    );
    return @row;
}

	
sub RowsForGene($$$$) {
    my ($gene,$shade,$group,$row) = @_;
    my @trows = ();
   
    # my $firstForGene = 1;
    my $orgId = $gene->{orgId};
    my $showId = $gene->{sysName} || $gene->{locusId};
    my $genome = $orginfo->{$orgId}{genome};
    my $genomeShort = $genome;
    $genomeShort =~ s/^(.)\S+/\1./;
    my $locusId = $gene->{locusId};

    my @fitList = ();
    my %fitHash;
    foreach my $expName (sort keys %{ $gene->{fit} }) {
        my $fit = $gene->{fit}{$expName}{fit};
  	    my $t = $gene->{fit}{$expName}{t};
  	    push @fitList, $fit;
          $fitHash{$fit} = $t;
    }
    my ($min, $max) = minmax @fitList;
    my $rowLabel = "";
    my $collapse = "";
    # q[-valign => 'middle', -align => 'left', -bgcolor => $shade % 2 ? "#DDDDDD" : "#FFFFFF"];
    if ($row == 0) {
        $rowLabel = a({style=>"color: #011B47", title=>"Ortholog group $group"},$group);
        $collapse = 'header';
        # q[-class=>'header', -valign => 'middle', -align => 'left', -bgcolor => $shade % 2 ? "#DDDDDD" : "#FFFFFF"];
    }
    push @trows, $cgi->Tr(
        { -class=> $collapse, -valign => 'middle', -align => 'left', -bgcolor => $shade % 2 ? "#DDDDDD" : "#FFFFFF" },
                td($rowLabel),
			      td( 
                  #$firstForGene ? 
                  a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId},  -title => $genome },  $genomeShort)),# : "&nbsp;" ),
			      td( 
                    #$firstForGene ? 
                    a({ -href => "myFitShow.cgi?orgId=$orgId&gene=$locusId" }, $showId)),
				  #: "&nbsp"),
			      td( #$firstForGene ? 
                  $gene->{gene}),#: "&nbsp;" ),
			      td( 
                  #$firstForGene ? 
                  $gene->{desc}), #: "&nbsp;" ),
			      # td( a({ -href => "exp.cgi?orgId=$orgId&expName=$expName" }, $expDesc{$orgId}{$expName}) ),
			      td( { -bgcolor => Utils::fitcolor($min), -style=>'text-align: center;' },
				  a( { -title => sprintf("Click to compare (t = %.1f)",$fitHash{$min}), 
                    -style => "color: rgb(0,0,0)",
                     onMouseOver=>"this.style.color='#CC0024'", onMouseOut=>"this.style.color='#000000'",
				       -href => "orthFit.cgi?orgId=$orgId&locusId=$locusId&expGroup=$expGroup&condition1=$condition1"},
				     sprintf("%.1f",$min) ) ),
                  td( { -bgcolor => Utils::fitcolor($max), -style=>'text-align: center;' },
                  a( { -title => sprintf("Click to compare (t = %.1f)",$fitHash{$max}), 
                    -style => "color: rgb(0,0,0)",
                    onMouseOver=>"this.style.color='#CC0024'", onMouseOut=>"this.style.color='#000000'",
                       -href => "orthFit.cgi?orgId=$orgId&locusId=$locusId&expGroup=$expGroup&condition1=$condition1"},
                     sprintf("%.1f",$max) ) )
	    );
	# $firstForGene = 0;

    return @trows;
}
