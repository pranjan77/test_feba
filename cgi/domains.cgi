#!/usr/bin/perl -w

#######################################################
## domains.cgi -- sequence analysis of a protein
##
## Copyright (c) 2015 University of California
##
## Authors:
## Victoria Lo and Morgan Price
#######################################################
#
# Required CGI parameters:
# orgId -- which organism
# locusId -- which locus

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

my $dbh = Utils::get_dbh();
my $orginfo = Utils::orginfo($dbh);
Utils::fail($cgi, "Unknown organism: $orgId") unless $orgId eq "" || exists $orginfo->{$orgId};
my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE orgId=? AND locusId=?",
				   {}, $orgId, $locusId);

Utils::fail($cgi,"unknown gene") unless defined $gene->{locusId};
Utils::fail($cgi,"sequence information is only available for protein-coding genes.") unless $gene->{type} == 1;

# domains table
# gather data and slice it into an array of hashes
my $cond = $dbh->selectall_arrayref(qq{SELECT domainDb, orgId, locusId, domainId, domainName, begin, end, score, evalue, definition FROM GeneDomain WHERE orgId=? AND locusId=? ORDER BY begin;},
    { Slice => {} },
    $orgId, $locusId);

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
my $title = "Protein Info for $orginfo->{$orgId}{genome} at Locus $locusId";
my $start = Utils::start_page("$title");
my $tabs = Utils::tabsGene($dbh,$cgi,$orgId,$locusId,0,1,"protein");

my $sys = $gene->{sysName} || $gene->{locusId};

print
    header, $start, $tabs, '<div id="tabcontent">',
    h2("Protein Info for $sys in " . $cgi->a({href => "org.cgi?orgId=$orgId"},$orginfo->{$orgId}{genome})),
    p("Description:", $gene->{desc}),
    h3("Domains");

if (@$cond == 0) {
	print "Sorry, no domains for this organism and/or locus."
} else {
	#create domains table
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

# UniProt information, if any
my $bhSprot = $dbh->selectrow_hashref("SELECT * from BestHitSwissProt
                                         JOIN SwissProtDesc USING (sprotAccession,sprotId)
                                         WHERE orgId = ? AND locusId = ?",
                                      {}, $orgId, $locusId);
if (defined $bhSprot->{sprotAccession}) {
    my $acc = $bhSprot->{sprotAccession};
    print
        h3("SwissProt"),
        p(
            sprintf("%.0f%% similar to", $bhSprot->{identity}),
            a({ href => "http://www.uniprot.org/uniprot/$acc",
                title => $bhSprot->{sprotAccession} },
              $bhSprot->{sprotId} ) . ":",
            $bhSprot->{desc},
            $bhSprot->{geneName} ? "($bhSprot->{geneName})" : "",
            br(),
            "from",
            $bhSprot->{organism} );
}
print "\n";

# KEGG information, if any

my $bhKEGG = $dbh->selectrow_hashref("SELECT * from BestHitKEGG where orgId = ? AND locusId = ?",
                                     {}, $orgId, $locusId);
if (defined $bhKEGG->{keggId}) {
    my $keggOrg = $bhKEGG->{keggOrg};
    my $keggId = $bhKEGG->{keggId};
    my $ko = $dbh->selectall_arrayref("SELECT * from KEGGMember JOIN KgroupDesc USING (kgroup)
                                       WHERE keggOrg = ? AND keggId = ?",
                                      { Slice => {} }, $keggOrg, $keggId);
    print h3("KEGG");
    if (@$ko > 0) {
        foreach my $row (@$ko) {
            my $ecs = $dbh->selectcol_arrayref("SELECT ecnum FROM KgroupEC WHERE kgroup = ?",
                                               {}, $row->{kgroup});
            my @ecdesc = ();
            foreach my $ec (@$ecs) {
                if ($ec =~ m/-/) {
                    push @ecdesc, $ec;
                } else {
                    push @ecdesc, a({-href => "http://www.kegg.jp/dbget-bin/www_bget?ec:$ec"}, $ec);
                }
            }
            my $ecdesc = join(" ", @ecdesc);
            $ecdesc = "[EC: $ecdesc]" if $ecdesc ne "";

            print p(a({-href => "http://www.kegg.jp/dbget-bin/www_bget?ko:" . $row->{kgroup}},
                      $row->{kgroup}), " :",
                    $row->{desc} || "(no description)",
                    $ecdesc );
        }
    } else {
        print p("No KEGG ortholog group");
    }
    print p(sprintf("(inferred from %.0f%% similarity to ", $bhKEGG->{identity}),
            a({-href => "http://www.kegg.jp/dbget-bin/www_bget?$keggOrg:$keggId"},
              "$keggOrg:$keggId").")");
}


# print sequence
$seq =~ s/(.{60})/$1\n/gs;

print
    h3("Protein Sequence"),
    qq[<div style="position: relative; margin: 0 auto; width: 59em; font-family: monospace;white-space: pre;"],
    pre(">$sys $gene->{desc} ($orginfo->{$orgId}{genome})\n$seq"),
    '</div></div>';

$dbh->disconnect();
Utils::endHtml($cgi);
