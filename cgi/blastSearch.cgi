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
 #    start_html(
 #      $start,
	# -title =>"Fitness Browser",
	# # -style => {-code => $style},
	# -author=>'wjshaoATberkeley.edu, victorialoATberkeley.edu',
	# -meta=>{'copyright'=>'copyright 2015 UC Berkeley'}),

    
div({-id=>"ntcontent"},
  h2("Search by Gene Sequence"),
  h5(q{Browse <i>unpublished</i> fitness experiments from the
     <A HREF="http://pbd.lbl.gov/scientists/adam-deutschbauer/">Deutschbauer lab</A>,
     the <A HREF="http://genomics.lbl.gov/">Arkin lab</A>,
     and collaborators. Contact <A HREF="mailto:AMDeutschbauer.lbl.gov">Adam Deutschbauer</A> for more information.}),

    # Gene search form

	# start_form( -class => "search", -name    => 'input', -method  => 'GET', -action  => 'myFitShow.cgi' ),
	# # drop down list of species
	# p("Choose organism:",
	#   popup_menu( -name    => 'orgId', -values  => \@orgOptions, -labels  => \%orgLabels, -default => $orgOptions[0])),
	# p("Enter gene name:",
	#   textfield( -name      => 'gene', -size => 20, -maxlength => 100 ),
	#   br(),
	#   small(qq{examples: "Shewana3_0001" or "recA"})),
	# p(submit("Find gene")),
	# end_form,


 #    # Experiment search form
 #    div({-style => "padding-left: 10px; padding-right: 10px; float: right; vertical-align: top; background-color: rgb(240,240,240);"},
	# h3("Find Experiments"),
	# start_form(-name => 'byExp', -method => 'GET', -action => 'exps.cgi'),
	# p("Choose organism:",
	#   popup_menu(-name => 'orgId', -values => \@orgOptions, -labels => \%orgLabels, -default => $orgOptions[0])),
	# p("Enter condition:", textfield(-name => 'query', -size => 20, -maxlength => 100), submit("Go"),
	#   br(),
	#   small(qq{examples: "cisplatin" or "glucose"})),
	# p("Or show",submit("All experiments"),"for one organism"),
	# p("Or see all conditions:", br(),
	#   a({ href => "cond.cgi?expGroup=carbon source" }, "carbon sources," ),
	#   a({ href => "cond.cgi?expGroup=nitrogen source" }, "nitrogen sources," ),
	#   "or",
	#   a({ href => "cond.cgi?expGroup=stress" }, "stresses" )),
 #        end_form),
 #    div({-style => "clear: right"}),


    # h3(qq(Search by Gene Sequence)),
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

 );
    
    end_html;

$dbh->disconnect();

# END
