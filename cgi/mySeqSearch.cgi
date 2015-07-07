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

my $cgi=CGI->new;
my $style = Utils::get_style();

# read the input information

my $query = $cgi->param('query') || "";
my $qtype = $cgi->param('qtype') || "protein";
my $numHit = $cgi->param('numHit') || 20;
my $locusSpec = $cgi->param('locusId') || "";
my $orgId = $cgi->param('orgId') || "";
my $dbh = Utils::get_dbh();

my $start = Utils::start_page("Blast Results");
my $tabs = "";
my $orth;
if (defined $orgId and $locusSpec) {
    $tabs = Utils::tabsGene($dbh,$cgi,$orgId,$locusSpec,0,1,"homo");
    $orth = $dbh->selectall_hashref("SELECT * FROM Ortholog WHERE orgId1 = ? AND locusId1 = ?",
                             'orgId2',{Slice=>{}}, $orgId, $locusSpec);
    my $orginfo = Utils::orginfo($dbh);
    $start = Utils::start_page("Orthologs for locus $locusSpec ($orginfo->{$orgId}{genome})");
} else {
    $tabs = '<div id = "ntcontent">';
}


print $cgi->header, $start, $tabs;
# print $cgi->start_html(
#     -title =>"Blast Results",
#     -style => {-code => $style},
#     -author=>'wjshaoATberkeley.edu',
#     -meta=>{'copyright'=>'copyright 2015 UC Berkeley'},
# #    -BGCOLOR=>'#fffacd'
# );

# check user input

Utils::fail($cgi,qq($locusSpec is invalid. Please enter correct gene name!)) unless ($locusSpec =~ m/^[A-Za-z0-9_]*$/);

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

# blast query sequence from the front page

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
        s/^[0-9]+//;
        Utils::fail($cgi,"Error: more than one sequence was entered.") if m/^>/;
        Utils::fail($cgi,"Unrecognized characters in $_") unless m/^[a-zA-Z*]*$/;
        s/[*]/X/g;
        $seq .= uc($_);
    }

    my $id = "query";
    open(FAA,">",$seqFile) || die "Cannot write fasta file";
    print FAA Utils::formatFASTA($id,$seq);
    close(FAA) || die "Error writing fasta file";    

    # run blast

    if ($qtype eq "nucleotide") {
        Utils::fail($cgi,qq($query is invalid. Please enter nucleotide sequence or choose sequence type as protein!)) unless ($seq =~ m/^[ATCGatcg]*$/);
        system($blast,'-p','blastx','-e','1e-2','-d',$myDB,'-i',$seqFile,'-o',$blastOut,'-m','8')==0 || die "Error running blastx: $!";
    } elsif ($qtype eq "protein") {
        Utils::fail($cgi,qq($query is invalid. Please enter correct protein sequence!)) unless ($seq =~ m/^[A-Za-z]*$/);
        system($blast,'-p','blastp','-e','1e-2','-d',$myDB,'-i',$seqFile,'-o',$blastOut,'-m','8')==0 || die "Error running blastp: $!";
    }

# check homologs

} elsif ($locusSpec ne "") {

    # extract sequence for the given gene

    my $gene = $dbh->selectrow_hashref("SELECT * from Gene where orgId=? AND locusId=?", {}, $orgId, $locusSpec);
    Utils::fail($cgi, "no such gene: $locusSpec in $orgId") unless defined $gene->{locusId};
    Utils::fail($cgi, "homology search is only available for protein-coding genes") unless $gene->{type} == 1;

    my $id = join(":",$orgId,$locusSpec);
    my $fastacmd = '../bin/blast/fastacmd';
    system($fastacmd,'-d',$myDB,'-s',$id,'-o',$seqFile)==0 || die "Error running $fastacmd -d $myDB -s $id -o $seqFile -- $!";

    my $in = Bio::SeqIO->new(-file => $seqFile,-format => 'fasta');
    $seq = $in->next_seq()->seq;

    system($blast,'-p','blastp','-e','1e-2','-d',$myDB,'-i',$seqFile,'-o',$blastOut,'-m','8')==0 || die "Error running blastp: $!";

} else {
    print $cgi->p("No sequence or gene specified!");
}

# parse and report the blast result:
# blast output fields: (1)queryId, (2)subjectId, (3)percIdentity, (4)alnLength, (5)mismatchCount, (6)gapOpenCount, (7)queryStart, (8)queryEnd, (9)subjectStart, (10)subjectEnd, (11)eVal, (12)bitScore
# sort the blast result by bit score, E-value, and percent identity
system('sort','-k1,1','-k12,12gr','-k11,11g','-k3,3gr',$blastOut,'-o',$blastSort)==0 || die "Error running sort: $!";

# output blast result

print $cgi->h2("Blast Result");

my $orginfo = Utils::orginfo($dbh);
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
    my @hit = ($cgi->a({href => "org.cgi?orgId=$orgId"},$orginfo->{$orgId}->{genome}),
        $cgi->a({href => "geneOverview.cgi?orgId=$orgId&gene=$locusId"},$sys || $locusId),
        $cgi->a({href => "myFitShow.cgi?orgId=$orgId&gene=$locusId"}, $geneName),
        $desc,
	    $cgi->a({href => "myFitShow.cgi?orgId=$orgId&gene=$locusId", title => $fittitle }, $fitstring ),
        $cgi->a({title=>"evalue: $eVal ($bitScore bits)"},$percIdentity),
        $cgi->a({title=>"evalue: $eVal ($bitScore bits)"},$cov));
    # @hit = map td($_), @hit;

    my $tr;
    if (defined $orgId and $locusSpec) {
        if (exists $orth->{$orgId}) { #and $orth->{$orgId} == $locusId) {
            unshift @hit, '<a title="ortholog">o</a>';
            @hit = map td($_), @hit;
            $tr = $cgi->Tr({ -align => 'left', -valign => 'top', bgcolor=>'white' }, @hit );
        } else {
            unshift @hit, ' ';
            @hit = map td($_), @hit;
            $tr = $cgi->Tr({ -align => 'left', -valign => 'top' }, @hit );
        }
    } else {
        unshift @hit, ' ';
        @hit = map td($_), @hit;
        $tr = $cgi->Tr({ -align => 'left', -valign => 'top' }, @hit );
    }
        
        push @hits, $tr;
        # die $tr;
    $cnt++;
}
close(RES) || die "Error reading $blastOut";
$#hits = $numHit-1 if defined $numHit && $numHit > 0 && @hits > $numHit;

if ($cnt > 0) {

    print $cgi->p("Top " . scalar(@hits) . " hits (E < 0.01)");

    # my @td = ();
    # foreach my $row (@hits) {
    #     push @td, $cgi->td( $row );
    # }
    print $cgi->table(
        { cellspacing=>0, cellpadding=>3 },
        $cgi->Tr({-align=>'left',-valign=>'top'},
		 $cgi->th( ['', 'Species','Gene ID','Name','Description','Fitness','Identity %','Coverage %'])),#,'E Value','Bit Score' ] ) ),
            # $cgi->Tr({ -align => 'left', -valign => 'top' }, \@td), 
            @hits
    );

} else {
    print $cgi->p("No hit found!");
}

    print qq[<br><a href="http://www.microbesonline.org/cgi-bin/seqsearch.cgi?qtype=protein&query=$seq">Search for homologs in MicrobesOnline</a><BR><BR>];
$dbh->disconnect();
unlink($seqFile) || die "Error deleting $seqFile: $!";
unlink($blastOut) || die "Error deleting $blastOut: $!";
unlink($blastSort) || die "Error deleting $blastSort: $!";

Utils::endHtml($cgi);

exit 0;

# END

#----------------------------------------

