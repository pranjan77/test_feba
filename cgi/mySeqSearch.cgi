#!/usr/bin/perl -w
#######################################################
## mySeqSearch.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Wenjun Shao (wjshao@berkeley.edu) and Morgan Price
#######################################################
#
# Required CGI parameters:
# Either query (the sequence, a.a. by default)
# or orgId and locusId (if coming from a page for that gene)
#
# Optional CGI parameters in query mode:
# qtype -- protein or nucleotide (default is protein)
#
# Optional CGI parameters in either mode:
# numHit -- how many hits to show (default is 20)

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Time::HiRes qw(gettimeofday);
use DBI;
use Bio::SeqIO;

use lib "../lib";
use Utils;
use IO::Handle;

my $cgi=CGI->new;
my $style = Utils::get_style();

# read the input information
my $query = $cgi->param('query') || "";
my $qtype = $cgi->param('qtype') || "protein";
my $numHit = $cgi->param('numHit') || 20;
my $locusSpec = $cgi->param('locusId') || "";
my $orgId = $cgi->param('orgId') || "";
my $dbh = Utils::get_dbh();

# set up $seq and $seqFile
my $procId = $$;
my $timestamp = int (gettimeofday * 1000);
my $filename = $procId . $timestamp;
my $tmpDir = Utils::tmp_dir();
my $seqFile = "$tmpDir/$filename.fasta";
my $blast = '../bin/blast/blastall';
my $myDB = Utils::blast_db();
my $blastOut = "$tmpDir/$filename.blast.out";
my $blastSort = "$tmpDir/$filename.blast.sort";
my $seq;
my $locusShow;

print $cgi->header; # must be printed before using fail()

if ($query =~ m/[A-Za-z]/) {
    # parse and write the input sequence
    $seq = "";
    my @lines = split /[\r\n]+/, $query;
    my $def = "";
    $def = shift @lines if @lines > 0 && $lines[0] =~ m/^>/;
    $def =~ s/^>//;
    $def = "query sequence" if $def eq "";

    foreach (@lines) {
        s/[ \t]//g;
        s/^[0-9]+//; # leading digit/whitespace occurs in UniProt format
        next if $_ eq "//";
        Utils::fail($cgi,"Error: more than one sequence was entered.") if m/^>/;
        Utils::fail($cgi,"Unrecognized characters in $_") unless m/^[a-zA-Z*]*$/;
        s/[*]/X/g;
        $seq .= uc($_);
    }

    my $id = "query";
    open(FAA,">",$seqFile) || die "Cannot write fasta file";
    print FAA Utils::formatFASTA($id,$seq);
    close(FAA) || die "Error writing fasta file";    
} elsif ($locusSpec ne "") {
    Utils::fail($cgi,qq($locusSpec is invalid. Please enter correct gene name!)) unless ($locusSpec =~ m/^[A-Za-z0-9_]*$/);

    # extract sequence for the given gene
    my $gene = $dbh->selectrow_hashref("SELECT * from Gene where orgId=? AND locusId=?", {}, $orgId, $locusSpec);
    Utils::fail($cgi, "no such gene: $locusSpec in $orgId") unless defined $gene->{locusId};
    Utils::fail($cgi, "homology search is only available for protein-coding genes") unless $gene->{type} == 1;
    my $id = join(":",$orgId,$locusSpec);
    $locusShow = $gene->{gene} || $gene->{sysName} || $gene->{locusId};
    my $fastacmd = '../bin/blast/fastacmd';
    system($fastacmd,'-d',$myDB,'-s',$id,'-o',$seqFile)==0 || die "Error running $fastacmd -d $myDB -s $id -o $seqFile -- $!";
    my $in = Bio::SeqIO->new(-file => $seqFile,-format => 'fasta');
    $seq = $in->next_seq()->seq;
} else {
    Utils::fail($cgi, "No sequence or gene specified!");
}

# print start of page
my $orginfo = Utils::orginfo($dbh);
my $title = "";
my $tabs = "";
my $orth;
if (defined $orgId and $locusSpec) {
    $title = "Homologs of $locusShow from $orginfo->{$orgId}{genome}";
    $tabs = Utils::tabsGene($dbh,$cgi,$orgId,$locusSpec,0,1,"homo");
    $orth = $dbh->selectall_hashref("SELECT * FROM Ortholog WHERE orgId1 = ? AND locusId1 = ?",
                             'orgId2',{Slice=>{}}, $orgId, $locusSpec);
} else {
    my $qlen = length($seq);
    my $qchar = $qtype eq "protein" ? "a.a." : "nt.";
    $title = "Blast Results for " . substr($seq,0,20) . "... ($qlen $qchar)";
    $tabs = '<div id = "ntcontent">';
}
my $start = Utils::start_page($title);

print $start, $tabs, $cgi->h2($title), "\n";
autoflush STDOUT 1; # so header shows while blast is being run

# run blast
if ($qtype eq "nucleotide") {
    Utils::fail($cgi,qq($query is invalid. Please enter nucleotide sequence or choose sequence type as protein!)) unless ($seq =~ m/^[ATCGatcg]*$/);
    system($blast,'-p','blastx','-e','1e-2','-d',$myDB,'-i',$seqFile,'-o',$blastOut,'-m','8')==0 || die "Error running blastx: $!";
} elsif ($qtype eq "protein") {
    Utils::fail($cgi,qq($query is invalid. Please enter correct protein sequence!)) unless ($seq =~ m/^[A-Za-z]*$/);
    system($blast,'-p','blastp','-e','1e-2','-d',$myDB,'-i',$seqFile,'-o',$blastOut,'-m','8')==0 || die "Error running blastp: $!";
} else {
    die "Unknown query type $qtype";
}

# parse and report the blast result:
# blast output fields: (1)queryId, (2)subjectId, (3)percIdentity, (4)alnLength, (5)mismatchCount, (6)gapOpenCount, (7)queryStart, (8)queryEnd, (9)subjectStart, (10)subjectEnd, (11)eVal, (12)bitScore
# sort the blast result by bit score, E-value, and percent identity
system('sort','-k1,1','-k12,12gr','-k11,11g','-k3,3gr',$blastOut,'-o',$blastSort)==0 || die "Error running sort: $!";

# output blast result

open(RES,$blastSort) || die "Error reading $blastSort";
my $cnt = 0;
my @hits = ();
while(<RES>) {
    chomp;
    my ($queryId,$subjectId,$percIdentity,$alnLength,$mmCnt,$gapCnt,$queryStart,$queryEnd,$subjectStart,$subjectEnd,$eVal,$bitScore) = split /\t/, $_;
    my ($orgId,$locusId) = split /:/, $subjectId;
    my $cov = sprintf("%.1f", 100*abs($queryEnd - $queryStart + 1)/length($seq));
    $percIdentity = sprintf("%.1f", $percIdentity);

    my ($sys,$geneName,$desc) = $dbh->selectrow_array("SELECT sysName,gene,desc FROM Gene WHERE orgId = ? AND locusId = ?",
                                                      undef, $orgId, $locusId);

    if (!defined $desc) {
	print "Warning! Unknown hit $orgId:$locusId<BR>";
	next;
    }
    

    my ($fitstring, $fittitle) = Utils::gene_fit_string($dbh, $orgId, $locusId);
    my $showId = $sys || $locusId;
    my $seqlen = length($seq);
    my @hit = ($cgi->a({href => "org.cgi?orgId=$orgId"},$orginfo->{$orgId}->{genome}),
               $cgi->a({href => "geneOverview.cgi?orgId=$orgId&gene=$locusId"},$showId),
               $cgi->a({href => "myFitShow.cgi?orgId=$orgId&gene=$locusId"}, $geneName),
               $cgi->a({ -title => Utils::alt_descriptions($dbh,$orgId,$locusId) || "no other information",
                         -href => "domains.cgi?orgId=$orgId&locusId=$locusId"},
                       $desc),
               $cgi->a({href => "myFitShow.cgi?orgId=$orgId&gene=$locusId", title => $fittitle }, $fitstring ),
               $cgi->a({title=>"evalue: $eVal ($bitScore bits)"},$percIdentity),
               $cgi->a({title=>"Query $queryStart..$queryEnd of $seqlen aligns to $showId $subjectStart..$subjectEnd"},$cov));
    
    if (defined $orgId and $locusSpec) {
        # add ortholog indicator
        if (exists $orth->{$orgId} && $orth->{$orgId}{locusId2} eq $locusId) {
            unshift @hit, '<center><a title="ortholog">o</a></center>';
        } else {
            unshift @hit, ' ';
        }
    }
    @hit = map td($_), @hit;
    push @hits, $cgi->Tr({ -align => 'left', -valign => 'top', bgcolor=>'white' }, @hit );
    $cnt++;
}
close(RES) || die "Error reading $blastOut";
$#hits = $numHit-1 if defined $numHit && $numHit > 0 && @hits > $numHit;

if ($cnt > 0) {

    print $cgi->p("Top " . scalar(@hits) . " hits (E < 0.01)");
    my @header = ('Orth?', 'Species','Gene','Name','Description','Fitness', '%Identity', '%Coverage');
    shift @header if $query;
    print $cgi->table(
        { cellspacing=>0, cellpadding=>3 },
        $cgi->Tr({-align=>'left',-valign=>'top'},
		 $cgi->th( \@header )),
            # $cgi->Tr({ -align => 'left', -valign => 'top' }, \@td), 
            @hits
    );

} else {
    print $cgi->p("No hits found!");
}

    print qq[<br><a href="http://www.microbesonline.org/cgi-bin/seqsearch.cgi?qtype=protein&query=$seq">Search for homologs in MicrobesOnline</a><BR><BR>] unless $qtype eq "nucleotide";
$dbh->disconnect();
unlink($seqFile) || die "Error deleting $seqFile: $!";
unlink($blastOut) || die "Error deleting $blastOut: $!";
unlink($blastSort) || die "Error deleting $blastSort: $!";

Utils::endHtml($cgi);

exit 0;

# END

#----------------------------------------

