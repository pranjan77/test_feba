#!/usr/bin/perl -w
#######################################################
## myFitShow.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo, Wenjun Shao (wjshao@berkeley.edu) and 
## Morgan Price
#######################################################
#
# Required CGI parameters:
# gene -- a locusId, sysName, or gene name to match on
#	(may show multiple hits)
# Optional CGI parameters:
# orgId -- which organism to search in
# showAll -- 1 if showing all fitness values instead of just the most extreme ones

use strict;

use CGI qw(:standard Vars -nosticky);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use IO::Handle;

use lib "../lib";
use Utils;
sub end;

my $cgi=CGI->new;

my $orgSpec = $cgi->param('orgId') || "";
my $geneSpec = $cgi->param('gene');
my $help = $cgi->param('help') || "";

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
die "Illegal orgId parameter $orgSpec" if $orgSpec ne "" & !exists $orginfo->{$orgSpec};

my $orgTitle = $cgi->param('orgId') ? "in $orginfo->{$orgSpec}{genome}" : "";
my $start = Utils::start_page("Genes matching $geneSpec $orgTitle");

$geneSpec =~ s/[ \t]*$//;
$geneSpec =~ s/^[ \t]*//;
    
# if no gene or locus specified
if (!defined $geneSpec || $geneSpec eq "") {
    print $cgi->header;
    print $start,'<div id="ntcontent">';
    Utils::fail($cgi, "you must enter the gene name or locus tag");
}

# match for exact gene or locus ID
my $hits = Utils::matching_exact($dbh, $orgSpec, $geneSpec);

if (scalar(@$hits) == 1) {
    # just 1 hit, redirect
    my $gene = $hits->[0];
    my $orgId = $gene->{orgId};
    my $locusId = $gene->{locusId};
    my $url = "singleFit.cgi?orgId=$orgId&locusId=$locusId&showAll=0";
    $url = "singleFit.cgi?orgId=$orgId&locusId=$locusId&showAll=0&help=1" if $help == 1;

    print redirect(-url=>"$url");
    exit(0);

}

my $orgDescriptor = exists $orginfo->{$orgSpec}{genome} ? " in " . $orginfo->{$orgSpec}{genome} : "";

#start page
print $cgi->header, $start,'<div id="ntcontent">';

my %used = (); # org -> locusId -> 1 for genes that have been shown already
my $count = 0; 

# make table for exact gene/locus matches
if (scalar(@$hits) > 0) {
    my @trows = ();
    print h3(b("Match by gene/locus for $geneSpec:"));
    push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
                          $cgi->th( [ 'Gene ID','Gene Name','Description','Genome','Fitness' ] ) );

    

    foreach my $gene (@$hits) {
    	if ($count < 100) {
            $count++;
            $used{$gene->{orgId}}{$gene->{locusId}} = 1;
            my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $gene->{orgId}, $gene->{locusId});
			my @trow = map $cgi->td($_), (
                            a( {href => "geneOverview.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}"},
                               $gene->{sysName}||$gene->{locusId}), 
                            $gene->{gene}, 
                            $gene->{desc},
                            $cgi->a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId}}, "$orginfo->{$gene->{orgId}}->{genome}"),
                            a( {href => "myFitShow.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}", title => $fittitle, },
                               $fitstring));
            push @trows, $cgi->Tr(@trow);
        }
    }
    print $cgi->table( { cellspacing=>0, cellpadding=>3 }, @trows);
    print "\n"; # to allow flushing
}

autoflush STDOUT 1; # so preliminary results appear

# Handle queries like ko:K14333 to search for a KEGG ortholog group
if ($geneSpec =~ m/^ko:(K\d+)$/i) {
    my $kgroup = uc($1);
    print p("Searching for members of KEGG ortholog group $kgroup");
    my $hits = Utils::matching_kgroup($dbh, $orgSpec, [ $kgroup ]);
    @$hits = grep { !exists $used{$_->{orgId}}{$_->{locusId}} } @$hits;
    if (@$hits > 0) {
        my @trows = ();
        my ($kgroupDesc) = $dbh->selectrow_array("SELECT desc from KgroupDesc WHERE kgroup = ?",
                                                 {}, $kgroup);
        $kgroupDesc =~ s/ +$//;
        print h3(b("Members of $kgroup ($kgroupDesc)", $orgDescriptor));
        push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
                              $cgi->th( [ 'Gene ID','Gene Name','Description','Genome','Fitness' ] ) );
        foreach my $gene (@$hits) {
            next if $count >= 100;
            $count++;
            my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $gene->{orgId}, $gene->{locusId});
            my @trow = map $cgi->td($_), (
                a( {href => "geneOverview.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}"},
                   $gene->{sysName}||$gene->{locusId} ), 
                $gene->{gene}, 
                $gene->{desc},
                $cgi->a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId}},
                        "$orginfo->{$gene->{orgId}}->{genome}"),
                a( {href => "myFitShow.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}", title => $fittitle, },
                   $fitstring));
            push @trows, $cgi->Tr(@trow);
        }
        print $cgi->table( { cellspacing=>0, cellpadding=>3 }, @trows);
        print "\n";
    }
    &end();
}

# Handle queries like EC:1.4.3.1
if ($geneSpec =~ m/^ec:([0-9A-Za-z.-]+)$/i) {
    my $ecnum = $1;
    my ($ecdesc) = $dbh->selectrow_array("SELECT ecdesc FROM ECInfo WHERE ecnum = ? ;", {}, $ecnum);
    $ecdesc = "unknown" if !defined $ecdesc;
    print p("Searching for Enyzme Commission number $ecnum ($ecdesc) by TIGRFam, by KEGG ortholog group, and then by SEED annotation");
    my $hits1 = Utils::matching_domain_ec($dbh, $orgSpec, $ecnum);
    @$hits1 = grep { !exists $used{ $_->{orgId} }{ $_->{locusId} } } @$hits1;
    if (@$hits1 > 0) {
        print h3(b("Match by TIGRFam's EC number for $ecnum"));
        my @trows = ();
        push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
			      $cgi->th( [ 'Gene ID','Gene Name','Description','Genome', 'Domain ID', 'Fitness' ] ) );
        foreach my $gene (@$hits1) {
            $used{ $gene->{orgId} }{ $gene->{locusId} } = 1;
            next if $count >= 100;
            $count++;
            my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $gene->{orgId}, $gene->{locusId});
            my @trow = map $cgi->td($_), (
                a( {href => "geneOverview.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}"},
                   $gene->{sysName}||$gene->{locusId}), 
                $gene->{gene}, 
                $gene->{desc},
                $cgi->a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId}}, "$orginfo->{$gene->{orgId}}->{genome}"),
                $gene->{domainId},
                # $gene->{domainName}, # always the same as the id for TIGR
                a( {href => "myFitShow.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}", title => $fittitle, },
                   $fitstring));
            push @trows, $cgi->Tr(@trow);
        }
        print $cgi->table( { cellspacing=>0, cellpadding=>3 }, @trows);
        print "\n";
    }
    &end() if $count >= 100;
    my $kgroups = $dbh->selectcol_arrayref("SELECT kgroup FROM KgroupEC WHERE ecnum = ?", {}, $ecnum);
    my $hits2 = Utils::matching_kgroup($dbh, $orgSpec, $kgroups);
    @$hits2 = grep { !exists $used{ $_->{orgId} }{ $_->{locusId} } } @$hits2;
    if (@$hits2 > 0) {
        print h3(b("Match by KEGG's EC number for $ecnum"));
        my @trows = ();
        push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
			  $cgi->th( [ 'Gene ID','Gene Name','KO','KO Description','Genome','Fitness' ] ) );

        my $kgroupSpec = join(",", map {"'".$_."'"} map {$_->{kgroup}} @$hits2);
        my $kgroupDesc = $dbh->selectall_hashref("SELECT * from KgroupDesc WHERE kgroup IN ($kgroupSpec)",
                                                 "kgroup");
        foreach my $gene (@$hits2) {
            $used{ $gene->{orgId} }{ $gene->{locusId} } = 1;
            next if $count >= 100;
            $count++;
            my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $gene->{orgId}, $gene->{locusId});
            my @trow = map $cgi->td($_), (
                a( {href => "geneOverview.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}",
                    title => $gene->{desc}}, $gene->{sysName}||$gene->{locusId}), 
                $gene->{gene},
                a( {href => "http://www.kegg.jp/dbget-bin/www_bget?ko:".$gene->{kgroup} },
                   $gene->{kgroup}),
                $kgroupDesc->{$gene->{kgroup}}{desc},
                $cgi->a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId}}, "$orginfo->{$gene->{orgId}}->{genome}"),
                a( {href => "myFitShow.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}", title => $fittitle, },
                   $fitstring));
            push @trows, $cgi->Tr(@trow);
        }
        print $cgi->table( { cellspacing=>0, cellpadding=>3 }, @trows);
        print "\n";
    }
    &end() if $count >= 100;
    my $seedquery = qq{SELECT orgId,locusId,seed_desc,sysName,gene,desc
                       FROM SEEDClass JOIN SEEDAnnotation USING (orgId,locusId)
                       JOIN Gene USING (orgId,locusId)
                       WHERE num = ? };
    $seedquery .= qq{ AND orgId = "$orgSpec" } if $orgSpec ne "";
    my $hits3 = $dbh->selectall_arrayref($seedquery, { Slice => {} }, $ecnum);
    @$hits3 = grep { !exists $used{ $_->{orgId} }{ $_->{locusId} } } @$hits3;
    if (@$hits3 > 0) {
        print h3(b("Match by SEED's EC number for $ecnum"));
        my @trows = ();
        push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
                              $cgi->th( [ 'Gene ID','Gene Name', 'Seed Annotation', 'Genome','Fitness' ] ) );
        foreach my $gene (@$hits3) {
            $used{ $gene->{orgId} }{ $gene->{locusId} } = 1;
            next if $count >= 100;
            $count++;
            my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $gene->{orgId}, $gene->{locusId});
            my @trow = map td($_), (
                a( {href => "geneOverview.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}",
                    title => $gene->{desc}}, $gene->{sysName} || $gene->{locusId}),
                $gene->{gene},
                $gene->{seed_desc},
                a( {href => "org.cgi?orgId=$gene->{orgId}"}, $orginfo->{$gene->{orgId}}{genome}),
                a( {href => "myFitShow.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}", title => $fittitle},
                   $fitstring));
            push @trows, Tr(@trow);

        }
        print $cgi->table({ cellspacing=>0, cellpadding=>3 }, @trows);
        print "\n";
    }
    &end();
}


# match by description
my $descs;
$descs = Utils::matching_descs($dbh, $orgSpec, $geneSpec) if $count < 100;
@$descs = grep { !exists $used{ $_->{orgId} }{ $_->{locusId} } } @$descs;

if (@$descs >= 1) {
    print h3(b("Match by description for $geneSpec:"));
    my @trows = ();
    push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
                          $cgi->th( [ 'Gene ID','Gene Name','Description','Genome','Fitness' ] ) );
    
    foreach my $gene (@$descs) {
        if ($count < 100) {
            $count++;
            $used{$gene->{orgId}}{$gene->{locusId}} = 1;
            my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $gene->{orgId}, $gene->{locusId});
            my @trow = map $cgi->td($_), (
                a( {href => "geneOverview.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}"},
                   $gene->{sysName}||$gene->{locusId} ), 
                $gene->{gene}, 
                $gene->{desc},
                $cgi->a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId}},
                        "$orginfo->{$gene->{orgId}}->{genome}"),
                a( {href => "myFitShow.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}", title => $fittitle, },
                   $fitstring));
            push @trows, $cgi->Tr(@trow);
        } 
    }
    print $cgi->table( { cellspacing=>0, cellpadding=>3 }, @trows);
    print "\n"; # to allow flush
}

# match by KEGG orthology group's description
my $kegghits;
$kegghits = Utils::matching_kgroup_descs($dbh, $orgSpec, $geneSpec) if $count < 100;
@$kegghits = grep { !exists $used{ $_->{orgId} }{ $_->{locusId} } } @$kegghits;
if (@$kegghits > 0) {
    print h3(b("Match by KEGG ortholog group description for $geneSpec:"));
    my @trows = ();
    push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
			  $cgi->th( [ 'Gene ID','Gene Name','KO','KO Description','Genome','Fitness' ] ) );
    my $kgroupSpec = join(",", map {"'".$_."'"} map {$_->{kgroup}} @$kegghits);
    my $kgroupDesc = $dbh->selectall_hashref("SELECT * from KgroupDesc WHERE kgroup IN ($kgroupSpec)",
                                              "kgroup");
    foreach my $gene (@$kegghits) {
        $used{$gene->{orgId}}{$gene->{locusId}} = 1;
        next if $count >= 100;
        $count++;
        $count++;
        my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $gene->{orgId}, $gene->{locusId});
        my @trow = map $cgi->td($_), (
            a( {href => "geneOverview.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}",
                title => $gene->{desc}}, $gene->{sysName}||$gene->{locusId}), 
            $gene->{gene},
            a( {href => "http://www.kegg.jp/dbget-bin/www_bget?ko:".$gene->{kgroup} },
               $gene->{kgroup}),
            $kgroupDesc->{$gene->{kgroup}}{desc},
            $cgi->a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId}}, "$orginfo->{$gene->{orgId}}->{genome}"),
            a( {href => "myFitShow.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}", title => $fittitle, },
               $fitstring));
        push @trows, $cgi->Tr(@trow);
    }
    print $cgi->table( { cellspacing=>0, cellpadding=>3 }, @trows);
    print "\n"; # to allow flush
}

# match by SEED description
my $seedhits = Utils::matching_seed_descs($dbh, $orgSpec, $geneSpec) if $count < 100;
@$seedhits = grep { !exists $used{$_->{orgId}}{$_->{locusId}} } @$seedhits;
if (@$seedhits > 0) {
    print h3(b("Match by SEED description for $geneSpec:"));
    my @trows = ();
    push @trows, $cgi->Tr({-align =>'CENTER',-valign=>'TOP'},
                          th(['GENE ID','Gene Name','SEED Description','Genome','Fitness']));
    foreach my $gene (@$seedhits) {
        $used{$gene->{orgId}}{$gene->{locusId}} = 1;
        next if $count >= 100;
        $count++;
        my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $gene->{orgId}, $gene->{locusId});
        my @trow = map td($_), (
            a( {href => "geneOverview.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}",
                title => $gene->{desc}}, $gene->{sysName}||$gene->{locusId}),
            $gene->{gene},
            $gene->{seed_desc},
            $cgi->a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId}}, "$orginfo->{$gene->{orgId}}->{genome}"),
            a( {href => "myFitShow.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}", title => $fittitle, },
               $fitstring)
        );
        push @trows, $cgi->Tr(@trow);
    }
    print $cgi->table( { cellspacing=>0, cellpadding=>3 }, @trows);
    print "\n"; # to allow flush
}

# make table for domain matches if total matches < 100 so far, filtering for repeats
my $domains;
$domains = Utils::matching_domains($dbh, $orgSpec, $geneSpec) if $count < 100;
@$domains = grep { !exists $used{ $_->{orgId} }{ $_->{locusId} } } @$domains;
if (@$domains >= 1) {
    print h3(b("Match by domain for $geneSpec:"));
    my @trows = ();
    push @trows, $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
			  $cgi->th( [ 'Gene ID','Gene Name','Description','Genome', 'Domain ID', 'Domain Name', 'Fitness' ] ) );

    foreach my $gene (@$domains) {
        if ($count < 100) {
            next if exists $used{ $gene->{orgId} }{ $gene->{locusId} };
            $count++;
            my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $gene->{orgId}, $gene->{locusId});
            my @trow = map $cgi->td($_), (
                a( {href => "geneOverview.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}"},
                   $gene->{sysName}||$gene->{locusId}), 
                $gene->{gene}, 
                $gene->{desc},
                $cgi->a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId}}, "$orginfo->{$gene->{orgId}}->{genome}"),
                $gene->{domainId},
                $gene->{domainName},
                a( {href => "myFitShow.cgi?orgId=$gene->{orgId}&gene=$gene->{locusId}", title => $fittitle, },
                   $fitstring));
            push @trows, $cgi->Tr(@trow);
            $used{$gene->{locusId}} = 1;
        }
    }
    print $cgi->table( { cellspacing=>0, cellpadding=>3 }, @trows);
}

&end();

sub end {
    if ($count >= 100) {
	print "<BR> Only the first 100 results are shown. <br>";
    } elsif ($count == 0) {
        print $cgi->h3("No gene found for $geneSpec", $orgDescriptor);
    }
    $dbh->disconnect();
    Utils::endHtml($cgi);
    exit(0);
}
