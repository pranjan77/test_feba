#!/usr/bin/perl -w
#######################################################
## geneSearch.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo, Wenjun Shao (wjshao@berkeley.edu), and 
## Morgan Price
#######################################################
#
# CGI parameters: none

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);

use lib "../lib";
use Utils;

my $cgi=CGI->new;
# my $style = Utils::get_style();
my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
my @orgOptions = ("");
my @orginfo = sort { $a->{genome} cmp $b->{genome} } values(%$orginfo);
my %orgLabels = ("" => join(" ", "All", scalar(@orginfo), "organisms"));
foreach my $hash (@orginfo) {
    my $orgId = $hash->{orgId};
    push @orgOptions, $orgId;
    $orgLabels{$orgId} = $hash->{genome};
}
my $start = Utils::start_page("Fitness Browser - Sequence Search");

print header, $start,
    div({-id=>"ntcontent"},
        h2("Search by Gene Sequence"),
      Utils::site_intro_text(),
        start_form( -class => 'search', -name    => 'input', -method  => 'POST', -action  => 'mySeqSearch.cgi' ),
        p("Choose query type: ", 
          popup_menu( -name => 'qtype', -values => [ 'protein', 'nucleotide' ], -default => 'protein' )),
        p("Enter sequence in FASTA or Uniprot format: ",
          br(),
          textarea( -name  => 'query', -value => '', -cols  => 70, -rows  => 10 )),
        p("How many hits to show: ",
          popup_menu( -name => 'numHit', -values  => [ 5,10,25,50,100 ], -default => 25 )),
        p({class=>"buttons"}, submit("Search"),
          reset() ),
        end_form,
 
        h6(q(Developed by Victoria Lo, Wenjun Shao, and
         <A HREF="http://morgannprice.org/">Morgan Price</A>.
         Please report any bugs to <A HREF="mailto:funwithwords26@gmail.com">Morgan</A>.)),
    ),
    end_html;

$dbh->disconnect();

# END
