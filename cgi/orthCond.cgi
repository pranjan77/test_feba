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

use lib "../lib";
use Utils;
sub RowsForGene($$); # gene object, shade or not => list of rows

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

my $title = "Specific Genes for $expGroup $condition1 Across Organisms";
$title = "No Specific Genes" if (@$speclist == 0);
print
    header,
    start_html( -title => $title, -style => { -code => $style }, -author => 'Morgan Price',
		-meta => { 'copyright' => 'copyright 2015 UC Berkeley' }),
    h2($title),
    div({-style => "float: right; vertical-align: top;"}, a({href => "help.cgi#specific"}, "Help"));

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

my @headings = qw{organism gene name description experiment fitness};
my @trows = ( $cgi->Tr({ -valign => 'top', -align => 'center' }, map { th($_) } @headings ) );

my $shade = 0;

foreach my $og (@OGs) {
    foreach my $gene (@$og) {
	push @trows, RowsForGene($gene, $shade);
    }
    $shade++;
}

print
    p("Genes with",
      a({ -href => "help.cgi#specific" },  "specific phenotypes"),
      "in $expGroup $condition1 are",
      font({ style => "color: rgb(120,120,120); font-weight: bold;"}, "grouped"),
      "by",
      a({ -href => "help.cgi#ortholog"}, "orthology")),
    table({cellpadding => 3, cellspacing => 0}, @trows);

$dbh->disconnect();
Utils::endHtml($cgi);

	
sub RowsForGene($$) {
    my ($gene,$shade) = @_;
    my @trows = ();
   
    my $firstForGene = 1;
    my $orgId = $gene->{orgId};
    my $showId = $gene->{sysName} || $gene->{locusId};
    my $genome = $orginfo->{$orgId}{genome};
    my $genomeShort = $genome;
    $genomeShort =~ s/^(.)\S+/\1./;
    my $locusId = $gene->{locusId};
    foreach my $expName (sort keys %{ $gene->{fit} }) {
	my $fit = $gene->{fit}{$expName}{fit};
	my $t = $gene->{fit}{$expName}{t};
	push @trows, $cgi->Tr({ -valign => 'top', -align => 'left', -bgcolor => $shade % 2 ? "#DDDDDD" : "#FFFFFF" },
			      td( $firstForGene ? a({ -title => $genome }, $genomeShort) : "&nbsp;" ),
			      td( $firstForGene ? a({ -href => "myFitShow.cgi?orgId=$orgId&gene=$locusId" }, $showId)
				  : "&nbsp"),
			      td( $firstForGene ? $gene->{name} : "&nbsp;" ),
			      td( $firstForGene ? $gene->{desc} : "&nbsp;" ),
			      td( a({ -href => "exp.cgi?orgId=$orgId&expName=$expName" }, $expDesc{$orgId}{$expName}) ),
			      td( { -bgcolor => Utils::fitcolor($fit) },
				  a( { -title => sprintf("Click to compare (t = %.1f)",$t), -style => "color: rgb(0,0,0)",
				       -href => "orthFit.cgi?orgId=$orgId&locusId=$locusId&expGroup=$expGroup&condition1=$condition1"},
				     sprintf("%.1f",$fit) ) )
	    );
	$firstForGene = 0;
    }
    return @trows;
}
