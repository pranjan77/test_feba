#!/usr/bin/perl -w
#######################################################
## mySeqSearch.cgi
##
## Copyright (c) 2015 All Right Reserved by UC Berkeley
##
## Authors:
## Wenjun Shao (wjshao@berkeley.edu) and Morgan Price
#######################################################

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
my $gene = $cgi->param('gene') || "";
my $species = $cgi->param('species') || "";

print $cgi->header;
print $cgi->start_html(
    -title =>"Blast Result",
    -style => {-code => $style},
    -author=>'wjshaoATberkeley.edu',
    -meta=>{'copyright'=>'copyright 2015 UC Berkeley'},
#    -BGCOLOR=>'#fffacd'
);

# check user input

Utils::fail($cgi,qq($gene is invalid. Please enter correct gene name!)) unless ($gene =~ m/^[A-Za-z0-9_]*$/);

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

} elsif ($gene ne "") {

    # extract sequence for the given gene

    my $id = join(":",$species,$gene);
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

# connect to database

my $dbh = Utils::get_dbh();
my ($stmt, $sth, $rv);

# output blast result

print $cgi->h2("Blast Result");

open(RES,$blastSort) || die "Error reading $blastSort";
my $cnt = 0;
my @hits = ();
while(<RES>) {
    last if $cnt >= $numHit;
    chomp;
    my ($queryId,$subjectId,$percIdentity,$alnLength,$mmCnt,$gapCnt,$queryStart,$queryEnd,$subjectStart,$subjectEnd,$eVal,$bitScore) = split /\t/, $_;
    my ($nickname,$locusId) = split /:/, $subjectId;
    my $cov = sprintf("%.1f", 100*abs($queryEnd - $queryStart + 1)/length($seq));
    $percIdentity = sprintf("%.1f", $percIdentity);

    # select gene information from the database
    my $stmt = qq(SELECT genus, species FROM Organism WHERE name = ?;);
    $sth = $dbh->prepare( $stmt );
    $rv = $sth->execute( $nickname ) or die $DBI::errstr;
    print $DBI::errstr if($rv < 0);

    my ($genus,$species) = $sth->fetchrow_array();

    $stmt = qq(SELECT sysName, gene, desc FROM Gene WHERE nickname = ? AND locusId = ?;);
    $sth = $dbh->prepare( $stmt );
    $rv = $sth->execute( $nickname, $locusId ) or die $DBI::errstr;
    print $DBI::errstr if($rv < 0);

    my ($sys,$gene,$desc) = $sth->fetchrow_array();

    #check if any fitness data exists
    $stmt = qq(SELECT expName, fit, t FROM GeneFitness WHERE species = ? AND locusId = ?;);
    $sth = $dbh->prepare( $stmt );
    $rv = $sth->execute( $nickname, $locusId ) or die $DBI::errstr;
    print $DBI::errstr if($rv < 0);

    my $fitness = "no data";
    while(my @row = $sth->fetchrow_array()) {
        my $dest = "myFitShow.cgi?species=$nickname&gene=$locusId";
        $fitness = qq(<a href=$dest>check data</a>);
        last;
    }

    my @hit = ($locusId,$sys,$gene,$desc,$species,$percIdentity,$cov,$eVal,$bitScore,$fitness);
    push @hits, @hit;
    $cnt++;
}
close(RES) || die "Error reading $blastOut";

if ($cnt > 0) {

    print $cgi->p("Top $cnt hits:");
    print $cgi->h5("Note: only significant hits (E-value < 0.01) are considered.");

    my @td = ();
    while ( my @elems = splice @hits, 0, 10 ) {
        push @td, $cgi->td( \@elems );
    }
    print $cgi->table(
        { -border=>1, cellpadding=>3 },
        $cgi->Tr({-align=>'CENTER',-valign=>'TOP'},
            $cgi->th( [ 'geneId','sysName','gene','description','species','identity%','coverage%','eValue','bitScore','fitness' ] ) ),
            $cgi->Tr( \@td )
    );

} else {
    print $cgi->p("No hit found!");
}

$dbh->disconnect();
unlink($seqFile) || die "Error deleting $seqFile: $!";
unlink($blastOut) || die "Error deleting $blastOut: $!";
unlink($blastSort) || die "Error deleting $blastSort: $!";


print $cgi->h4(qq(<a href="myFrontPage.cgi">Go back to front page</a>));

print $cgi->end_html;

exit 0;

# END

#----------------------------------------

