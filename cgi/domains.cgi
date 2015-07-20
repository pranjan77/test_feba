#!/usr/bin/perl -w

#######################################################
## domains.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo and Morgan Price
#######################################################
#
# This page generates the overall results for a single 
# organism, with a link to download the file of all 
# experimental results (generated by createExpData.cgi) 
# in addition to a table with overall data, linking to 
# pages with the individual expGroup results.
#
# Required CGI parameters:
# orgId -- which organism to search for
# locusId -- which locus to search for

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use Bio::SeqIO;

use lib "../lib";
use Utils;

my $cgi=CGI->new;

my $orgId = $cgi->param('orgId') || die "No orgId parameter";
my $locusId = $cgi->param('locusId') || die "No locusId parameter";

# my $style = Utils::get_style();
my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
Utils::fail($cgi, "Unknown organism: $orgId") unless $orgId eq "" || exists $orginfo->{$orgId};


# main table
# gather data and slice it into an array of hashes
my $cond = $dbh->selectall_arrayref(qq{SELECT domainDb, orgId, locusId, domainId, domainName, begin, end, score, evalue, definition FROM GeneDomain WHERE orgId=? AND locusId=? ORDER BY begin;},
    { Slice => {} },
    $orgId, $locusId);
# Utils::fail($cgi, "Unknown locus: $locusId") unless $locusId eq ""; #|| exists $locusId;

# gather number of genes and data
# my $numGenes = $dbh->selectrow_array(qq{SELECT COUNT (DISTINCT locusId) FROM Gene WHERE orgId = ?;}, undef, $orgId);

# my $numData = $dbh->selectrow_array(qq{SELECT COUNT (DISTINCT locusId) FROM GeneFitness WHERE orgId = ?;}, undef, $orgId);

#find length of sequence
my $tmpDir = Utils::tmp_dir();
my $seqFile = "$tmpDir/$orgId+$locusId.fasta";
my $myDB = Utils::blast_db();
my $id = join(":",$orgId,$locusId);
my $fastacmd = '../bin/blast/fastacmd';
system($fastacmd,'-d',$myDB,'-s',$id,'-o',$seqFile)==0 || die "Error running $fastacmd -d $myDB -s $id -o $seqFile -- $!";
my $in = Bio::SeqIO->new(-file => $seqFile,-format => 'fasta');
my $seq = $in->next_seq()->seq;
my $seqLen = length($seq);
unlink($seqFile) || die "Error deleting $seqFile: $!";


# write the title
# my $title = scalar(@$cond) > 0 ? # : "Protein Info";
my $title = "Protein Info for $orginfo->{$orgId}{genome} at Locus $locusId";
my $start = Utils::start_page("$title");
my $tabs = Utils::tabsGene($dbh,$cgi,$orgId,$locusId,0,1,"protein");

my $gene = $dbh->selectall_arrayref(qq(SELECT sysName,locusId,gene,desc FROM Gene where orgid = ? and locusId = ?;),
		{ },
	    $orgId, $locusId);
my $sys;
foreach my $grow(@$gene) {
		my ($sysName, $locusId, $geneName, $desc) = @$grow;
		$sys = $sysName || $locusId;
}

print
    header, $start, $tabs, '<div id="tabcontent">',
    # start_html( -title => $title, -style => { -code => $style }, -author => 'Morgan Price, Victoria Lo',
		# -meta => { 'copyright' => 'copyright 2015 UC Berkeley' }),
	h2("Protein Info for $sys in " . $cgi->a({href => "org.cgi?orgId=$orgId"},$orginfo->{$orgId}{genome})),
	h3("Domains");
    # h3("Gene Domains for ". $cgi->a({href => "org.cgi?orgId=$orgId"}, "$orginfo->{$orgId}{genome}")." at Locus $locusId"); #$title),
    # div({-style => "float: right; vertical-align: top;"}, a({href => "help.cgi#specific"}, "Help"));

#exit if no results
if (@$cond == 0) {
	# Utils::fail($cgi, 
	print "Sorry, no domains for this organism and/or locus."
} else {
	# sysname/locusId (gene): desc => myFitShow.cgi
	foreach my $grow(@$gene) {
		my ($sysName, $locusId, $geneName, $desc) = @$grow;
		my $sys = $sysName || $locusId;
		my $name = $geneName || "";
		my $d = $cgi->a({href=>"myFitShow.cgi?orgId=$orgId&gene=$locusId"}, $desc);
		print $cgi->p($cgi->b("$sys $name: $d"));
	}


	#create table
	my @headings = qw{Family ID Coverage EValue}; # Begin End};
	my @trows = ( Tr({ -valign => 'top', -align => 'center' }, map { th($_) } \@headings) );
	foreach my $row (@$cond) {
		# display result row by row
		my $len = $row->{end}-$row->{begin}; 
		my $begin = $row->{begin};
		my $newBegin = $begin;
		my $newLen = $len;
		my $newSeqLen = $seqLen;
		if ($seqLen > 600) {
			$newBegin = 600*$begin/$seqLen;
			$newLen = 600*$len/$seqLen;
			$newSeqLen = 600;
		}
		if ($row->{domainDb} eq 'PFam') {
	    	push @trows, Tr({ -valign => 'top', -align => 'left' },
	    	td([ $row->{domainName}, #name/description
			 	a({href => "http://pfam.xfam.org/family/$row->{domainId}"},
			 	$row->{domainId}), #ID
			 	# img({src=>"../images/grayHorizLine.png", height=>'7', width='200'})
			 	a({title=>"Amino acids $begin to $row->{end} ($len) of $seqLen"}, div({class=>"line"}, img({src=>"../images/grayHorizLine.png", width=>"$newSeqLen", height=>'7'}), div({class=>"line2", style=>"left:$newBegin".'px'}, img({src=>"../images/darkcyan.png", height=>'7', width=>"$newLen"})))),#$len, # $row->{}, #length diagram: end-begin
			 	# $row->{score}, #score
			 	a({title=>"Score: $row->{score}"},$row->{evalue}), #evalue with hover score
			 	# $row->{begin}, #begin
			 	# $row->{end}, #end
		 	]))
		} elsif ($row->{domainDb} eq 'TIGRFam') {
			push @trows, Tr({ -valign => 'top', -align => 'left' },
			td([ $row->{definition} || $row->{domainName}, #name/description
			 	a({href => "http://www.jcvi.org/cgi-bin/tigrfams/HmmReportPage.cgi?acc=$row->{domainId}"},
			 	$row->{domainId}), #ID
			 	a({title=>"Amino acids $begin to $row->{end} ($len) of $seqLen"}, div({class=>"line"}, img({src=>"../images/grayHorizLine.png", width=>"$newSeqLen", height=>'7'}), div({class=>"line2", style=>"left:$newBegin".'px'}, img({src=>"../images/chocolate.png", height=>'7', width=>"$newLen"})))), # $len, # $row->{}, #length diagram: end-begin
			 	# $row->{score}, #score
			 	a({title=>"Score: $row->{score}"},$row->{evalue}), #evalue with hover score
			 	# $row->{begin}, #begin
			 	# $row->{end}, #end
			 ])),
		}
	}

	print table({cellspacing => 0, cellpadding => 3}, @trows);

}

# %0A encodes "\n" so that it looks like fasta input.
print br(),
    p(a({-href => "http://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?seqinput=>${sys}%0A$seq"},
	(@$cond > 0 ? "Or see" : "See") . " Conserved Domains Database"));

# print sequence
my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId=? AND locusId=?",
				   {}, $orgId, $locusId);
Utils::fail($cgi,"unknown gene") unless defined $gene->{locusId};
Utils::fail($cgi,"sequence information is only available for protein-coding genes.") unless $gene->{type} == 1;

$seq =~ s/(.{60})/$1\n/gs;

# my $showId = $gene->{sysName} || $gene->{locusId};
# my $orginfo = Utils::orginfo($dbh);
# my $title = "Sequence of $showId in $orginfo->{$orgId}{genome}",;
print
  #   start_html( -title => $title,
		# -style => {-code => $style},
		# -author=>'jj326ATberkeley.edu',
		# -meta=>{'copyright'=>'copyright 2015 UC Berkeley'}),
    # h3("Sequence of $showId in " . $cgi->a({href => "org.cgi?orgId=". $orginfo->{$gene->{orgId}}->{orgId}}, "$orginfo->{$gene->{orgId}}->{genome}"),),
    h3("Protein Sequence"),
    qq[<div style="position: relative; margin: 0 auto; width: 59em; font-family: monospace;white-space: pre;"],
    # "<center>",
    pre(">$sys $gene->{desc} ($orginfo->{$orgId}{genome})\n$seq"),
    # "</center>";
    # p(a({href => "myFitShow.cgi?orgId=$orgId&gene=$locusId"}, "Show fitness")),
    # p(a({href => "mySeqSearch.cgi?orgId=$orgId&locusId=$locusId"}, "Check homologs"));
    '</div></div>';

#display number of genes that have data out of total genes
# print $cgi->p("Fitness data for $numData genes of $numGenes genes.");


#display taxonomic information and link
# if ((defined $orginfo->{$orgId}{taxonomyId}) && ($orginfo->{$orgId}{taxonomyId} ne "")) {
	# print $cgi->p($cgi->a({href => "http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=$orginfo->{$orgId}{taxonomyId}"}, "NCBI Taxonomy"));
# }

#file of experimental data - generated by createExpData.cgi
# my $data = `./createExpData.cgi orgId=$orgId`;
# print $cgi->p($cgi->a({href => "download.cgi?orgId=$orgId"}, "Download experimental data"), " - Note: May take a minute or so to load once clicked.");
# print $cgi->p($cgi->a({href => "createExpData.cgi?orgId=$orgId"}, "Download experimental data"), " - Note: May take a few seconds to load once clicked.");


$dbh->disconnect();
Utils::endHtml($cgi);
