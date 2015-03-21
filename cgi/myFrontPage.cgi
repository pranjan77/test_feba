#!/usr/bin/perl -w
#######################################################
## myFrontPage.cgi
##
## Copyright (c) 2015 All Right Reserved by UC Berkeley
##
## Authors:
## Wenjun Shao (wjshao@berkeley.edu) and Morgan Price
#######################################################

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();

print $cgi->header;
print $cgi->start_html(
    -title =>"FitBlast Start",
    -style => {-code => $style},
    -author=>'wjshaoATberkeley.edu',
    -meta=>{'copyright'=>'copyright 2015 UC Berkeley'},
#    -BGCOLOR=>'#fffacd'
);

print $cgi->h2("This is the FitBlast entering page!");

print $cgi->h4("Search by gene name:");
print $cgi->start_form(
        -name    => 'input',
        -method  => 'GET',
        -action  => 'myFitShow.cgi',
);
print $cgi->p("Choose Species Name:");
print $cgi->h6("Note: all species with fitness data are listed here.");

my $dbh = Utils::get_dbh();
my $stmt = qq(SELECT species FROM Organism;);
my $sth = Utils::execute_db($dbh, $stmt);
my @species = ("All");
while(my @row = $sth->fetchrow_array()) {
    my ($species) = @row;
    push @species, $species;
}
$species[0] = "All " . $#species . " genomes";

print $cgi->popup_menu(
    -name    => 'species',
    -values  => \@species,
    -default => $species[0],
);

print $cgi->p(q(Enter Gene Name<font color="red">*</font>:));
print $cgi->h6(q(Example: 7022746 (gene locusId) or Shewana3_0001 (gene sysName) or recA (gene name)));
print $cgi->textfield(
    -name      => 'gene',
    -size      => 20,
    -maxlength => 100,
);
print <<EndOfSubmitOne;
<BR><BR>
<INPUT TYPE="submit" VALUE="Start fitness lookup">
<INPUT TYPE="reset" VALUE="Clear" onClick="input.gene.value=''">
EndOfSubmitOne

print $cgi->end_form;
print "<BR>\n";

print $cgi->h4(qq(Search by gene sequence:));
print $cgi->start_form(
        -name    => 'input',
        -method  => 'POST',
        -action  => 'mySeqSearch.cgi',
);

print $cgi->p(q(Choose Query Type<font color="red">*</font>:));
my @qtype = ("protein","nucleotide");
print $cgi->popup_menu(
    -name    => 'qtype',
    -values  => \@qtype,
    -default => $qtype[0],
);
print $cgi->p(q(Enter Query Sequence<font color="red">*</font>:));
print $cgi->textarea(
    -name  => 'query',
    -value => '',
    -cols  => 70,
    -rows  => 10,
);
print $cgi->p(qq(Enter the number of Hits to show:));
my @num = (5,10,25,50,100);
print $cgi->popup_menu(
    -name    => 'numHit',
    -values  => \@num,
    -default => $num[2],
);

print <<EndOfSubmitTwo;
<BR><BR>
<INPUT TYPE="submit" VALUE="Start sequence search">
<INPUT TYPE="reset" VALUE="Clear sequence" onClick="input.query.value=''">
EndOfSubmitTwo

print $cgi->end_form;

print $cgi->end_html;

# END
