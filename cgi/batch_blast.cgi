#!/usr/bin/perl -w
#######################################################
## batch_blast.cgi -- submit a batch blast job
##
## Copyright (c) 2016 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Optional CGI parameters:
# file -- the fasta file to process

use strict;
use CGI qw(:standard Vars -nosticky);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use IO::Handle; # for autoflush
use URI::Escape;
use Digest::MD5 qw(md5_hex);

use lib "../lib";
use Utils;
use FEBA_Utils;
use Batch;

my $cgi=CGI->new;
my $file = $cgi->param('file');

my $maxBytes = 1e7;
my $maxSeq = 20000;
my $nThreads = 12;

if (!defined $file) {
    print
        header,
        Utils::start_page('Fitness BLAST for Genomes - Submit Sequences'),
        '<div id="ntcontent">',
        h2("Submit Sequences to Fitness BLAST for Genomes"),
        p("Fitness BLAST takes just 30 seconds to run on the typical bacterial genome.",
          "To run it, upload a file of protein sequences in FASTA format.",
          "Please limit your input to at most 20,000 sequences.",
          "Your protein sequence file should be at most 10 MB."),
        start_form( -name => 'input', -method => 'POST', -action => "batch_blast.cgi" ),
        qq{<INPUT type="file" name="file">},
        p(submit( { -style => 'position: absolute; text-align: left; float: none; display:inline;' , -name => "Upload" })),
        p(" &nbsp; "),
        end_form();
    my $exampleId = "27feb3f720f0536e0fbb1408512c8741";
    my $exampleName = "Pseudomonas putida KT2440";
    if (-e "../job_data/$exampleId/job.db") {
        print p(br(), 
                "To see what the output will be like, see the",
                a({-href => "batch_overview.cgi?jobId=$exampleId"}, "results for $exampleName,"),
                "including a list of",
                a({-href => "batch_genes.cgi?jobId=$exampleId&hypo=Hypotheticals&link=Either"},
                  "hypothetical proteins with functional links.")
            );
    }
    print p(small("Fitness BLAST for genomes is powered by",
            a({-href=>"http://www.drive5.com/usearch"}, "usearch.")));
    Utils::endHtml($cgi);
    exit(0);
}
# else process input file
die "Not a file upload as expected" unless ref($file) eq "Fh";

autoflush STDOUT 1; # so preliminary results appear
print
    header,
    Utils::start_page('Fitness BLAST - Processing Your Submission'),
    h2('Fitness BLAST - Processing Your Submission'),
    p("This should take less than 1 minute"),
    p("Checking uploaded file..."), "\n";

my $nSeq = 0;
my $bytes = 0;
my @lines = ();
while(my $line = <$file>) {
    $bytes += length($line);
    $nSeq++ if $line =~ m/^>/;
    if ($nSeq > $maxSeq) {
        print "Sorry, too many sequences. The limit is $maxSeq sequences\n";
        exit(0);
    }
    if ($bytes > $maxBytes) {
        print "Sorry, the input file is too large. The limit is 10 MB\n";
        exit(0);
    }
    push @lines, $line;
}
print p("Read $nSeq sequences, ", sprintf("%.2f", $bytes/1e6), " MB"), "\n";

my $jobId = md5_hex(@lines);
my $dir = "../job_data/$jobId";
if (-d $dir) {
    print p("This file has been processed by Fitness BLAST already!");
    if (! -e "$dir/job.db") {
        print ("But the results are not ready yet. This may be a duplicate submission.",
               "If the link below fails, try it again in 30 seconds or so");
    }
} else {
    my $faa = "$dir/query.faa";
    mkdir("../job_data/$jobId");
    system("chmod","a+w",$dir);
    open(FAA, ">", $faa) || die "Cannot write to $faa";
    foreach my $line (@lines) {
        print FAA $line."\n";
    }
    close(FAA) || die "Error writing to $faa";
    
    my %parse = ReadFastaDesc($faa);
    if (exists $parse{error}) {
        print p("Error parsing the input file: $parse{error}");
        print p("Please check that your file is a valid fasta file.");
        exit(0);
    }
    print p("Parsed input OK"),"\n";
    
    my $desc = $parse{desc};
    foreach my $name (keys %$desc) {
        unless ($name =~ m/^[0-9A-Za-z._-|]+$/) {
            print
                p("Sorry, the identifier " . URI_Escape($name) . " is not valid."),
                p("Identifiers in the fasta file may only contain the characters a-zA-Z0-9._|-");
            exit(0);
        }
    }
    print p("Checked headers"),"\n";
    print p("Running usearch..."), "\n";

    my $usearch = "../bin/usearch";
    die "No executable in $usearch" unless -x $usearch;

    my $udbfile = "../cgi_data/aaseqs.udb";
    die "No such file: $udbfile" unless -e $udbfile;
    my $code = system("nice $usearch -threads 16 -ublast $faa -db $udbfile -maxhits 50 -maxaccepts 50 -blast6out $dir/hits -evalue 0.001 >& $dir/usearch.log");
    if ($code != 0) {
        print p("Sorry, usearch failed. Please check your input file");
        exit(0);
    }
    
    print p("usearch succeeded, now building a database..."),"\n";
    $code = system("nice ../bin/batch_setup_tables.pl -db ../cgi_data/feba.db -dir $dir -aaseqs ../cgi_data/aaseqs >& $dir/setup.log");
    if ($code != 0) {
        print p("Sorry, building the database failed for $jobId.",
                "Please contact the system administrators.",
                "Make sure to include jobId=$jobId in your message");
        exit(0);
    }
    p("Fitness BLAST succeeded.");
}

my $url = url();
$url =~ s!/[a-zA-Z_.]+$!!;
my $rel = "batch_overview.cgi?jobId=$jobId";
$url = $url . "/" . $rel;
print
    p("See results at", br(),
      b(a({-href => $url }, $url))),
    p("Please remember to save the URL. But you can always get back to the results by re-uploading the file.");
Utils::endHtml($cgi);
exit(0);
