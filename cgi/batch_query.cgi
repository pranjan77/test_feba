#!/usr/bin/perl -w
#######################################################
## batch_query.cgi -- show best hits for a gene from a batch comparison
##
## Copyright (c) 2016 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required CGI parameters:
# jobId -- the job identifier
# queryId -- the query protein
# Optional CGI parameter: mode=BLAST

use strict;
use CGI qw(:standard Vars -nosticky);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use HTML::Entities;

use lib "../lib";
use Utils;
use Batch;

my $cgi=CGI->new;

my $jobId = $cgi->param('jobId');
my $queryId = $cgi->param('queryId');
my $mode = $cgi->param('mode');
die "Must specify jobId and queryId" unless $jobId ne "" && $queryId ne "";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my $bdb = Batch::get_batch_dbh($jobId);
my $jobname = Batch::get_job_name($jobId,$bdb);

my $qinfo = $bdb->selectrow_hashref("SELECT * FROM Query WHERE queryId=?", {}, $queryId);
die "Invalid queryId" unless $qinfo->{queryId} eq $queryId;

my $queryURL = undef;
$queryURL = "http://www.microbesonline.org/cgi-bin/fetchLocus.cgi?locus=$1" 
    if $queryId =~ m/^VIMSS(\d+)$/;

print header,
    Utils::start_page("Fitness BLAST - gene $queryId in $jobname"),
    '<div id="ntcontent">',
    h2($mode eq "BLAST" ? "All Hits" : "Best Hits",
       " for ",
       defined $queryURL ? a({-href => $queryURL}, $queryId) : $queryId,
       " from ", a({-href => "batch_overview.cgi?jobId=$jobId"}, $jobname)),
    p(b("Description: ", encode_entities($qinfo->{desc})));

if ($mode eq "BLAST") {
    my $aaseqSpec = "'" . $qinfo->{aaseq} . "'";
    print <<END
<script src="http://fit.genomics.lbl.gov/d3js/d3.min.js"></script>
<script src="http://fit.genomics.lbl.gov/images/fitblast.js"></script>
<div id="fitblast"></div>
<script>
fitblast_load_table("fitblast", "http://fit.genomics.lbl.gov/", $aaseqSpec);
</script>
END
;
    print p("Or see", a({-href=>"batch_query.cgi?jobId=$jobId&queryId=$queryId"}, "best hits"));

} else { # regular mode not FastBLAST mode
    my $besthits = $bdb->selectall_arrayref("SELECT * from BestHit WHERE queryId=? ORDER BY bits DESC",
                                            { Slice => {} }, $queryId);
    if (@$besthits == 0) {
        print p("Sorry, no good hits for this protein");
    } else {
        my @headings = (th("Orth"),th("Identity"),th("Organism"),th("Description"),
                        th({-colspan=>2},
                           a({-title=>"range of values"}, "Fitness")),
                        th("Specific"),
                        th(a({-title=>"Maximum cofitness"}, "Cofit")));
        my @trows = ( $cgi->Tr({ -valign => 'middle', -align => 'center' }, @headings) );
        foreach my $row (@$besthits) {
            my $orgId = $row->{orgId};
            my $locusId = $row->{locusId};
            my $orthSign = "?";
            $orthSign = "N" if ! $row->{isBBH};
            $orthSign = "Y" if $row->{coverage} >= 0.8 && $row->{identity} >= 30 && $row->{isBBH};
            my $orthText = sprintf("Coverage %.0f%% (minimum of both ways), %s",
                                   100 * $row->{coverage},
                                   $row->{bbh} ? "bidirectional best hit" : "not bidirectional");
            my $gene = $dbh->selectrow_hashref(qq"SELECT * from Gene where orgId=? AND locusId=?",
                                               {}, $orgId, $locusId);
            # Try to fail gracefully if $locus is invalid because of database changes
            my $locusShow = $gene->{sysName} || $locusId;
            my $locusDesc = $gene->{desc} || $locusId;
            my $fitLo = td("");
            my $fitHi = td("");
            my $maxCofit = "";
            my $spec_string = "";
            my $fitURL = "#";
            my $cofitURL = "#";
            if (defined $gene->{locusId}) {
                ($maxCofit) = $dbh->selectrow_array(
                    "SELECT cofit FROM Cofit where orgId = ? AND locusId = ? AND rank=1 LIMIT 1",
                    {}, $orgId, $locusId);
                $maxCofit = sprintf("%.2f", $maxCofit) if defined $maxCofit;
                my ($minFit,$maxFit,$minT,$maxT) = $dbh->selectrow_array(
                    qq{ SELECT min(fit), max(fit), min(t), max(t)
                        FROM GeneFitness WHERE orgId = ? AND locusId = ? ; },
                    {}, $orgId, $locusId);
                $fitURL = "singleFit.cgi?orgId=$orgId&locusId=$locusId";
                if (defined $minFit) {
                    $cofitURL = "cofit.cgi?orgId=$orgId&locusId=$locusId";
                    my $fitTitle = sprintf("t = %.1f to %.1f", $minT, $maxT);
                    $fitLo = td({-bgcolor => Utils::fitcolor($minFit), -style=>'text-align: center;' },
                               a( { -title => $fitTitle,
                                    -style => "color: rgb(0,0,0)",
                                    -href =>  $fitURL },
                                  sprintf("%.1f", $minFit)));
                    $fitHi = td({-bgcolor => Utils::fitcolor($maxFit), -style=>'text-align: center;' },
                               a( { -title => $fitTitle,
                                    -style => "color: rgb(0,0,0)",
                                    -href =>  $fitURL },
                                  sprintf("%.1f", $maxFit)));
                    my $spec = $dbh->selectall_arrayref(
                        "SELECT * from SpecOG WHERE orgId=? AND locusId=? ORDER BY minFit",
                        { Slice => {} }, $orgId, $locusId);
                    my @specs = ();
                    foreach my $row (@$spec) {
                        my $col = $row->{minFit} < 0 ? "darkblue" : "darkgoldenrod";
                        my $group = ($row->{expGroup} eq "carbon source" ? " (C)" :
                                     ($row->{expGroup} eq "nitrogen source" ? " (N)" : ""));
                        push @specs, span({-style => "color: $col", 
                                           -title => sprintf("fitness %.1f to %.1f", $row->{minFit}, $row->{maxFit}) },
                                           $row->{condition} . $group);
                    }
                    $spec_string = join(", ", @specs);
                }
            }
            push @trows, $cgi->Tr( td(a({-title => $orthText}, $orthSign)),
                                   td(sprintf("%.1f", $row->{identity})),
                                   td(small($orginfo->{$orgId}{genome})),
                                   td(small(a({-href => $fitURL,
                                               -title => $locusShow},
                                              encode_entities($locusDesc)))),
                                   $fitLo, $fitHi,
                                   td(small($spec_string)),
                                   td(a({-href=>$cofitURL}, $maxCofit)));
        }
        print table({cellpadding => 3, cellspacing => 0}, @trows);
    }
    print p("Or see",
            a({-href=>"batch_query.cgi?jobId=$jobId&queryId=$queryId&mode=BLAST"},
              "all hits for this protein"));
}

Utils::endHtml($cgi);

