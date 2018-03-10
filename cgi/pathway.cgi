#!/usr/bin/perl -w
#######################################################
## pathway.cgi -- view MetaCyc pathway with candidate
##	genes in that organism and fitness data
##
## Copyright (c) 2018 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId -- which organism
# pathwayId -- which MetaCyc pathway
#
# Optional CGI parameters:
# expName -- which experiment(s) to show fitness data from
# addexp -- add experiment(s) which match this query

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;

my $orgId = $cgi->param('orgId') || die "No orgId parameter";
my $pathId = $cgi->param('pathwayId') || die "No pathwayId parameter";
my @expNames = $cgi->param('expName');
my $addexp = $cgi->param('addexp');

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Unknown organism" unless $orgId eq "" || exists $orginfo->{$orgId};

my ($pathwayName) = $dbh->selectrow_array("SELECT pathwayName FROM MetacycPathway WHERE pathwayId = ?",
                                          {}, $pathId);
die "Unknown pathway: $pathId" unless defined $pathwayName;

# Fetch metadata on experiments
my %exps = map { $_->{expName} => $_ } @{ $dbh->selectall_arrayref("SELECT * FROM Experiment WHERE orgId = ?",
                                                                   { Slice => {} }, $orgId) };

if ($addexp) {
  my %expSeen = map { $_ => 1 } @expNames;
  my $addexps = Utils::matching_exps($dbh, $orgId, $addexp);
  foreach my $exp (@$addexps) {
    my $expName = $exp->{expName};
    next if exists $expSeen{$expName};
    $expSeen{$expName} = 1;
    push @expNames, $expName;
  }
}

my @exps = ();
foreach my $expName (@expNames) {
  die "Invalid expName $expName" unless $expName =~ m/^[a-zA-Z][a-zA-Z0-9_-]+$/;
  die "Invalid expName $expName for org $orgId" unless exists $exps{$expName};
  push @exps, $exps{$expName};
}

print
  header,
  Utils::start_page("MetaCyc Pathway: $pathwayName in $orginfo->{$orgId}{genome}"),
  h2("MetaCyc Pathway: ",
     a({ -href => "https://metacyc.org/META/NEW-IMAGE?type=PATHWAY&object=$pathId" }, $pathwayName),
     "in",
     a({ -href => "org.cgi?orgId=$orgId" }, $orginfo->{$orgId}{genome})),
  # the experiment selector
  start_form(-name => 'input', -method => 'GET', -action => 'pathway.cgi'),
  hidden( -name => 'orgId', -value => $orgId, -override => 1),
  hidden( -name => 'pathwayId', -value => $pathId, -override => 1),
  join("\n", map { hidden( -name => 'expName', -value => $_, -override => 1) } @expNames),
  p({-class => "buttons", style=>"align:left; white-space:nowrap; line-height:40px;"}, "Add experiment(s): ",
    textfield(-name => 'addexp', -default => "", -override => 1, -size => 20, -maxLength => 100),
    submit('Add','Add') ),
  end_form,
  br();

my $rxns = $dbh->selectall_arrayref(qq{ SELECT * FROM MetacycPathwayReaction
                                        JOIN MetacycReaction USING (rxnId)
                                        WHERE pathwayId = ? },
                                    { Slice => {} }, $pathId);
die "Empty pathway: $pathId\n" unless @$rxns > 0;
my %rxns = map { $_->{rxnId} => $_ } @$rxns;

my $preds = $dbh->selectall_arrayref("SELECT * FROM MetacycPathwayReactionPredecessor WHERE pathwayId = ?",
                                    { Slice => {} }, $pathId);
# %pred has rxnId => list of predecessor rxnIds, only including those in this pathway
my %pred = map { $_ => [] } keys %rxns;
foreach my $row (@$preds) {
  push @{ $pred{$row->{rxnId}} }, $row->{predecessorId}
    if exists $rxns{ $row->{predecessorId} } && exists $rxns{ $row->{rxnId} };
}

# Order the reactions by precedence
# First, every node with no precedors gets 0
my %score = ();
my $nNoScore = scalar(keys %score);
while (my ($rxnId, $predlist) = each %pred) {
  if (scalar(@$predlist) == 0) {
    $score{$rxnId} = 0;
  }
}

# Then, every successor of a node with a score gets score+1
# until we make no more assignments
for(;;) {
  my $nSet = 0;
  while (my ($succ, $predlist) = each %pred) {
    foreach my $rxnId (@$predlist) {
      if (exists $score{$rxnId} && !exists $score{$succ}) {
        $score{$succ} = $score{$rxnId} + 1;
        $nSet++;
      }
    }
  }
  last if $nSet == 0;
}

my @rxnIds = sort { $score{$a} <=> $score{$b} } keys %rxns;

my $primary = $dbh->selectall_arrayref("SELECT * FROM MetacycPathwayPrimaryCompound WHERE pathwayId = ?",
                                       { Slice => {} }, $pathId);
my %primary = (); # rxnId => compoundId => side => 1
foreach my $row (@$primary) {
  $primary{ $row->{rxnId} }{ $row->{compoundId} }{ $row->{side} } = 1;
}

my @trows = ();
if (@exps > 0) {
  # make a header row
  my @hrow = ( th("Reactions and Genes") );
  foreach my $exp (@exps) {
    push @hrow, th(a( { -href => "exp.cgi?orgId=$orgId&expName=$exp->{expName}",
                       -title => $exp->{expName} },
                     $exp->{expDesc} ));
  }
  push @trows, Tr(@hrow)."\n";
  # make a removal row
  my @rmrow = ( td("") );
  my $base = "pathway.cgi?orgId=$orgId&pathwayId=$pathId";
  foreach my $expId (@expNames) {
    my @expNames2 = grep { $_ ne $expId } @expNames;
    my $URL = join("&", $base, map { "expName=$_" } @expNames2);
    push @rmrow, td(a( { -href => $URL}, "remove"));
  }
  push @trows, Tr(@rmrow)."\n";
}

my $genestyle = "padding-left: 3em;"; # inset gene descriptions

my %primaryShown = (); # compoundId => 1 if shown as primary already
my $ncol = 1 + scalar(@expNames);
foreach my $rxnId (@rxnIds) {
  my $rxn = $rxns{$rxnId};
  my $ecs = $dbh->selectcol_arrayref("SELECT ecnum FROM MetacycReactionEC WHERE rxnId = ?",
                                    {}, $rxnId);
  my $ecdesc = "";
  ($ecdesc) = $dbh->selectrow_array("SELECT ecdesc FROM ECInfo WHERE ecnum = ? LIMIT 1",
                                     {}, $ecs->[0])
    if @$ecs > 0;
  my $rxnName = $rxn->{rxnName};
  $rxnName = $ecdesc || $rxnId if $rxnName eq "" || $rxnName =~ m/^[0-9][.][0-9.,]+$/;
  $rxnName =~ s/[.]$//;
  my $reverse = $rxn->{direction} == -1 ? " (in reverse)" : "";

  my @left = ();
  my @right = ();
  my $cmps = $dbh->selectall_arrayref(qq{ SELECT * FROM MetacycReaction
                                         JOIN MetacycReactionCompound USING (rxnId)
                                         LEFT JOIN MetacycCompound USING (compoundId)
                                         WHERE rxnId = ? },
                                     { Slice => {} }, $rxnId);
  foreach my $i (0..(scalar(@$cmps)-1)) {
    $cmps->[$i]{row} = $i;
  }
  my @sortedcmps = sort { (exists $primary{$rxnId}{$b->{compoundId}}) - (exists $primary{$rxnId}{$a->{compoundId}})
                            || (exists $primaryShown{$b->{compoundId}}) - (exists $primaryShown{$a->{compoundId}})
                              || $a->{row} - $b->{row} } @$cmps;
  foreach my $row (@sortedcmps) {
    my $compoundId = $row->{"compoundId"};
    # Some reactants are actually compound classes -- for those, just use the name.
    my $name = $row->{compoundName} || $compoundId;
    my $color = "DarkGrey";
    if (exists $primary{$rxnId}{$compoundId}) {
      $color = "Brown";
      $primaryShown{$compoundId} = 1;
    }
    $name = span({ style => "color: $color" }, $name);
    my $side = $rxn->{direction} * $row->{"side"};
    my $list = $side == -1 ? \@left : \@right;
    $name = $row->{coefficient} . " " . $name unless $row->{coefficient} eq "" || $row->{coefficient} eq "1";
    $name = $name . Sub("[" . lc($row->{compartment}) . "]") if $row->{compartment} ne "";
    push @$list, $name;
  }
  my $spontaneous = $rxn->{isSpontaneous} ? " (spontaneous)" : "";

  # And the gene(s) for this reaction
  # First, from BestHitMetacyc
  my $bh = $dbh->selectall_arrayref("SELECT * from BestHitMetacyc WHERE rxnId = ? AND orgId = ?",
                                    { Slice => {} }, $rxnId, $orgId);
  my @loci = map { $_->{locusId} } @$bh;

  # And add EC assignments
  my $ecGenes = Utils::EcToGenes($dbh, $orgId, $ecs); # ec => locusId => 1
  my @loci2 = ();
  while (my ($ec, $hash) = each %$ecGenes) {
    foreach my $locusId (keys %$hash) {
      push @loci2, $locusId;
    }
  }
  push @loci, sort @loci2;

  my %lociShown = (); # to ignore duplicates
  my @locirows = ();
  # show list of genes
  foreach my $locusId (@loci) {
    next if exists $lociShown{$locusId};
    $lociShown{$locusId} = 1;
    my $gene = $dbh->selectrow_hashref("SELECT * from Gene WHERE orgId = ? AND locusId = ?",
                                       {}, $orgId, $locusId);
    my @generow = ();
    push @generow, td({ -style => $genestyle }, Utils::gene_link($dbh, $gene, "name", "myFitShow.cgi"));
    my $fitrows = $dbh->selectall_arrayref("SELECT expName,fit,t FROM GeneFitness WHERE orgId=? AND locusId=?",
                                       {}, $orgId, $locusId);
    my %fit = map { $_->[0] => $_ } @$fitrows;
    foreach my $expName (@expNames) {
      my $fitrow = $fit{$expName};
      my ($fit,$t);
      (undef,$fit,$t) = @$fitrow if defined $fitrow;
      push @generow, td({ -bgcolor => Utils::fitcolor($fit) },
                        a({ -href => "strainTable.cgi?orgId=$orgId&locusId=$locusId&expName=$expName",
                            -title => defined $fitrow ? sprintf("t = %.1f", $t) : "No data",
                            -style => "color:rgb(0,0,0)" },
                          defined $fitrow ? sprintf("%.1f", $fit) : "&endash;"));
    }
    push @locirows, Tr(@generow);
  }
  push @trows, Tr(td({ -colspan => $ncol},
                     a({ -href => "https://metacyc.org/META/NEW-IMAGE?type=REACTION&object=$rxnId"}, $rxnName)
                     . $reverse . ":"
                     . br()
                     . join(" + ", @left) . "&rarr;" . join(" + ", @right)
                     . $spontaneous));
  if (@locirows > 0) {
    push @trows, @locirows;
  } else {
    push @trows, Tr(td({ -style => $genestyle, -colspan => $ncol }, "No genes"));
  }
}
print "\n", table( { cellspacing => 0, cellpadding => 3 }, @trows);
$dbh->disconnect();
Utils::endHtml($cgi);

