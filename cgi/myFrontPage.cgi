#!/usr/bin/perl -w
#######################################################
## myFrontPage.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Wenjun Shao (wjshao@berkeley.edu) and Morgan Price
#######################################################
#
# CGI parameters: none

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();

print $cgi->header;
print $cgi->start_html(
    -title =>"Fitness Web Site",
    -style => {-code => $style},
    -author=>'wjshaoATberkeley.edu',
    -meta=>{'copyright'=>'copyright 2015 UC Berkeley'},
#    -BGCOLOR=>'#fffacd'
);

print $cgi->h2("Fitness Web Site");
print $cgi->h5(q(This web site contains <i>unpublished</i> fitness experiments from the Deutschbauer lab, the Arkin lab, and collaborators. Contact <A HREF="mailto:AMDeutschbauer.lbl.gov">Adam Deutschbauer</A> for more information.));

my $dbh = Utils::get_dbh();

# Gene search form
print $cgi->h3("Search by Gene Name");
print $cgi->start_form(
        -name    => 'input',
        -method  => 'GET',
        -action  => 'myFitShow.cgi',
);
print qq{<div style="padding-left: 25px;">};
# drop down list of species
my @orgOptions = ("");
my $orginfo = Utils::orginfo($dbh);
my @orginfo = sort { $a->{genome} cmp $b->{genome} } values(%$orginfo);
my %orgLabels = ("" => join(" ", "All", scalar(@orginfo), "organisms"));
foreach my $hash (@orginfo) {
    my $orgId = $hash->{orgId};
    push @orgOptions, $orgId;
    $orgLabels{$orgId} = $hash->{genome};
}
print "<P>Choose organism: ";
print $cgi->popup_menu(
    -name    => 'orgId',
    -values  => \@orgOptions,
    -labels  => \%orgLabels,
    -default => $orgOptions[0]
);
print q{<P>Enter gene name: };
print $cgi->textfield(
    -name      => 'gene',
    -size      => 20,
    -maxlength => 100,
);
print q{<BR><small>examples: "7022746" or "Shewana3_0001" or "recA"</small>};
print $cgi->p(qq{<INPUT TYPE="submit" VALUE="Find gene">});
print $cgi->end_form;
print qq{</div>};

# Experiment search form
print $cgi->h3("Search by Condition");
print $cgi->start_form(-name => 'input', -method => 'GET', -action => 'exps.cgi');
print qq{<div style="padding-left: 25px;">};
print "<P>Choose organism: ";
print $cgi->popup_menu(-name => 'orgId', -values => \@orgOptions, -labels => \%orgLabels, -default => $orgOptions[0]);
print q{<P>Enter condition or experiment: };
print $cgi->textfield(-name => 'query', -size => 20, -maxlength => 100);
print $cgi->p(qq{<INPUT TYPE="submit" VALUE="Find experiments">});
print $cgi->end_form;
print qq{</div>};

# Gene sequence search form
print $cgi->h3(qq(Search by Gene Sequence));
print qq{<div style="padding-left: 25px;">};
print $cgi->start_form(
        -name    => 'input',
        -method  => 'POST',
        -action  => 'mySeqSearch.cgi',
);
print q(<P>Choose query type: );
my @qtype = ("protein","nucleotide");
print $cgi->popup_menu(
    -name    => 'qtype',
    -values  => \@qtype,
    -default => $qtype[0],
);
print "<P>Enter query sequence:<BR>";
print $cgi->textarea(
    -name  => 'query',
    -value => '',
    -cols  => 70,
    -rows  => 10,
);
print qq(<P>How many hits to show: );
my @num = (5,10,25,50,100);
print $cgi->popup_menu(
    -name    => 'numHit',
    -values  => \@num,
    -default => $num[2],
);
print $cgi->p(qq{<INPUT TYPE="submit" VALUE="Search"> <INPUT TYPE="reset" VALUE="Clear sequence" onClick="input.query.value=''">});
print $cgi->end_form;
print qq{</div>};

print $cgi->h6(q(Developed by Wenjun Shao and Morgan Price. Please report any bugs to <A HREF="mailto:funwithwords26@gmail.com">Morgan</A>.));
print $cgi->end_html;

# END
