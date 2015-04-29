#!/usr/bin/perl -w
#######################################################
## genesFit.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Wenjun Shao (wjshao@berkeley.edu) and Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId -- which organism
# genes -- 1 or more locusIds
# Optional CGI parameters:
# addgene -- a gene to add (locusId, name, sysName)
# showAll -- 1 if showing all fitness values instead of just the most extreme ones (not yet implemented)

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use List::Util 'sum';

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();
my $orgId = $cgi->param('orgId') || die "No orgId found";
my @locusIds = $cgi->param('locusId');
die "No locusId found" if @locusIds == 0;

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my $genome = $orginfo->{$orgId}{genome} || die "Invalid orgId $orgId";
my $expinfo = Utils::expinfo($dbh,$orgId);

my $addgene_error = undef;
my $addgene = $cgi->param('addgene');
if (defined $addgene && $addgene ne "") {
    $addgene =~ s/ +$//;
    if ($addgene !~ m/^[A-Za-z90-9_-]*$/) {
	$addgene_error = "Invalid gene to add";
    } else {
	my ($locusId) = $dbh->selectrow_array(qq{ SELECT locusId FROM Gene WHERE orgId = ?
                                                  AND (locusId = ? OR sysName = ? OR gene = ?) LIMIT 1 },
					      {}, $orgId, $addgene, $addgene, $addgene);
	if (defined $locusId) {
	    if (sum(map { $_ eq $locusId } @locusIds) > 0) {
		$addgene_error = qq{Gene "$addgene" (locus $locusId) is already included};
	    } else {
		push @locusIds, $locusId;
	    }
	} else {
	    $addgene_error = qq{Cannot find gene "$addgene"};
	}
    }
}

my @genes = ();
foreach my $locusId (@locusIds) {
    my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId = ? AND locusId = ?",
				       {}, $orgId, $locusId);
    die "No such locus $locusId in org $orgId" if !defined $gene->{locusId};
    # expName => "fit" => fitness value
    $gene->{fit} = $dbh->selectall_hashref(qq{SELECT expName,fit,t FROM GeneFitness
                                               WHERE orgId = ? AND locusId = ?},
					   "expName", {}, $orgId, $locusId);
    foreach my $expName (keys %{ $gene->{fit} }) {
	die "No such experiment: $expName" unless exists $expinfo->{$expName};
    }
    $gene->{nExps} = scalar(keys %{ $gene->{fit} });
    push @genes, $gene;
}

print $cgi->header;
print $cgi->start_html(
    -title =>"Fitness for " . scalar(@genes) . " genes in $genome",
    -style => {-code => $style},
    -author=>'Morgan Price',
    -meta=>{'copyright'=>'copyright 2015 UC Berkeley'},
);

Utils::fail($cgi,"None of the genes has fitness data")
    if (sum(map { $_->{nExps} > 0 ? 1 : 0} @genes) == 0);

my @exps = sort Utils::CompareExperiments (values %$expinfo);

print $cgi->h2("Fitness data for " . scalar(@genes) . " genes in $genome");
print $cgi->h3($addgene_error) if defined $addgene_error;
print $cgi->p("All " . scalar(@exps) . " fitness values, sorted by group and condition");

my @trows = ();
my @headings = qw{group condition};
my @headings2 = ("", "");
foreach my $gene (@genes) {
    push @headings, $cgi->a({ href => "myFitShow.cgi?orgId=$orgId&gene=$gene->{locusId}",
			      title => $gene->{desc} },
			    $gene->{sysName} || $gene->{locusId});
    push @headings2, $gene->{gene};
}
push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'}, $cgi->th(\@headings));
if (sum(map { $_ ne "" } @headings2) > 0) {
    push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'}, $cgi->th(\@headings2));
}
foreach my $exp (@exps) {
    my @values = ();
    my $expName = $exp->{expName};
    push @values, $cgi->td($exp->{expGroup});
    push @values, $cgi->td($cgi->a({ style => "color:rgb(0,0,0)",
				     title => "$expName: $exp->{expDescLong}",
				     href => "exp.cgi?orgId=$orgId&expName=$expName" },
				   $exp->{expDesc}));
    foreach my $gene (@genes) {
	my $showId = $gene->{sysName} || $gene->{locusId};
	my $fit = $gene->{fit}{$expName}{fit};
	my $t = $gene->{fit}{$expName}{t};
	if (defined $fit) {
	    my $fitShow = sprintf("%.1f",$fit);
	    $t = sprintf("%.1f", $t);
	    push @values, $cgi->td({ -bgcolor => Utils::fitcolor($fit) },
				   qq{<div title="$showId: t = $t">$fitShow</div>});
	} else {
	    push @values, $cgi->td({ -bgcolor => Utils::fitcolor($fit) }, "&nbsp;");
	}
    }
    push @trows, $cgi->Tr({-align=>'left',-valign=>'top', -style => "font-size: 70%" }, @values);
}
if (@genes > 0) {
    my @footer = ("","");
    foreach my $gene (@genes) {
	my @others = grep { $_->{locusId} ne $gene->{locusId} } @genes;
	my $url = "genesFit.cgi?orgId=$orgId&" . join("&", map { "locusId=$_->{locusId}" } @others);
	push @footer, $cgi->a( { href => $url, title => $gene->{desc} },
			       "remove<BR>" . ($gene->{sysName} || $gene->{locusId}) );
    }
    push @trows, $cgi->Tr( { -align=>'center', -valign=>'top' }, $cgi->td(\@footer));
}
print $cgi->table( { cellspacing => 0, cellpadding => 3 }, @trows);

print $cgi->start_form(-name => 'input', -method => 'GET', -action => 'genesFit.cgi');
print "<P>Add gene: ";
print $cgi->hidden( 'orgId', $orgId );
foreach my $locusId (@locusIds) { # avoid CGI sticky
    print qq{<input type="hidden" name="locusId" value="$locusId" />\n};
}
print $cgi->textfield( -name => 'addgene', -default => "", -override => 1, -size => 20, -maxLength => 100 );
print $cgi->end_form;

$dbh->disconnect();
Utils::endHtml($cgi);
