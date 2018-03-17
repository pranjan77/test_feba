#!/usr/bin/perl -w

#######################################################
## seedsubsystem.cgi -- overview of a SEED subsystem
##
## Copyright (c) 2017 University of California
##
## Author: Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId -- which organism
# subsystem -- which subsystem

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;

use lib "../lib";
use Utils;

my $cgi=CGI->new;

my $orgId = $cgi->param('orgId') || die "No orgId parameter";
my $subsystem = $cgi->param('subsystem') || die "No subsystem parameter";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Unknown organism" unless $orgId eq "" || exists $orginfo->{$orgId};

my $nice = $subsystem; $nice =~ s/_/ /g;
print
  header,
  Utils::start_page("SEED Subsystem: $nice in $orginfo->{$orgId}{genome}"),
  h2("SEED Subsystem: $nice in",
     a({ -href => "org.cgi?orgId=$orgId"}, $orginfo->{$orgId}{genome}));

my $roles = $dbh->selectcol_arrayref("SELECT seedrole FROM SEEDRoles WHERE subsystem = ?",
                                     {}, $subsystem);
if (@$roles == 0) {
  print p("No roles -- is this a valid subsystem?");
} else {
  my @th = map th($_), qw{Role Gene Name &nbsp;};
  my @trows = ( \@th );
  foreach my $role (@$roles) {
    my $genes = $dbh->selectall_arrayref(qq{ SELECT * FROM SEEDAnnotationToRoles
                                             JOIN SEEDAnnotation USING (seed_desc)
                                             JOIN Gene USING (orgId, locusId)
                                             WHERE orgId = ? AND seedrole = ? },
                                         { Slice => {} }, $orgId, $role);
    my $showRole = $role;
    # Link to EC number hits from other sources
    $showRole =~ s!EC (\d+[.]\d+[.]\d+[.]\d+)!<A HREF="myFitShow.cgi?orgId=$orgId&gene=ec:$1">EC $1</A>!g;
    if (@$genes == 0) {
      my @trow = map td($_), ( $showRole, "No genes", "", "", "" );
      push @trows, \@trow;
    } else {
      my $first = 1;
      foreach my $gene (@$genes) {
        my @trow = map td($_), ( $first ? $showRole : "&nbsp;",
                                 Utils::gene_link($dbh, $gene, "name", "myFitShow.cgi"),
                                 $gene->{gene} || "",
                                 checkbox('locusId', 1, $gene->{locusId}, '')
                               );
        push @trows, \@trow;
        $first = 0;
      }
    }
  }
  @trows = map Tr( { -valign => "top", -align => "left" }, @$_ ), @trows;
  print start_form(-name => 'input', -method => 'GET', -action => 'genesFit.cgi'),
    hidden('orgId', $orgId),
      table({ cellspacing => 0, style => "text-align: center;", cellpadding => 3}, @trows),
	p(submit(-class=>"heatmap", -name=>"Heatmap of selected genes")),
          end_form;
}

print
  p( start_form(-name => 'orgselect', -method => 'GET', -action => 'seedsubsystem.cgi'),
     hidden( -name => 'subsystem', -value => $subsystem, -override => 1),
     "Or select organism:",
     Utils::OrgSelector($orgId, $orginfo),
     "<button type='submit'>Go</button>",
     end_form ),
  p(a( { -href => "http://pubseed.theseed.org/FIG/seedviewer.cgi?page=Subsystems&subsystem=$subsystem" },
       "Or view this subsystem at the SEED" ),
   "(slow)");

Utils::endHtml($cgi);
exit(0);

