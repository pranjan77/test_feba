#!/usr/bin/perl -w
#######################################################
## seqservice.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# Required parameters: seq -- the protein sequence to look for matches to.
# Optional parameters: maxhits -- default is 20. Cannot be raised above 50.
# debug -- write lines starting with # to output to record status
# Must start the ublast service first with bin/start_ublast_service.pl
#
# For an example of a page that uses this service see
# ../images/fitblast_example.html
# ../images/fitblast.js

use strict;
use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use Cwd; # for absolute paths
use Time::HiRes qw(usleep);

use lib "../lib";
use Utils;

my $cgi=CGI->new;
print $cgi->header('text/plain');

my $seq = $cgi->param('seq');
$seq =~ s/[ \r\n]//g;
if (! $seq) {
    print "Error\nNo sequence specified\n";
    exit(0);
}
# Allow * as internal stop codon
unless ($seq =~ m/^[A-Z*]+$/) {
    print "Error\nIllegal sequence\n";
    exit(0);
}

my $maxhits = $cgi->param('maxhits');
if (defined $maxhits && $maxhits > 0) {
    $maxhits = $maxhits + 0; # convert to integer
} else {
    $maxhits = 20;
}

my $debug = $cgi->param('debug');

my $base = "../cgi_data/ublast_service";
die "No such directory: $base" unless -d $base;

my $len = length($seq);
my $dir = "$base/seq$$.$len";
mkdir($dir);
die "Cannot mkdir $dir" unless -d $dir;
system("chmod a+w $dir");

open(FASTA, ">", "$dir/faa") || die "Cannot write to $dir/faa";
print FASTA ">seq\n$seq\n";
close(FASTA) || die "Error writing to $dir/faa";

my $path = getcwd();
# note $dir.q is a file in the main ublast_service directory
open(QFILE, ">", "$dir.q") || die "Cannot write to $dir.q";
print QFILE "$path/$dir/faa $path/$dir\n";
close(QFILE) || die "Error writing to $dir.q";

# poll every 0.05 seconds for up to 15 seconds
for (my $i = 0; $i < 300; $i++) {
    if (-e "$dir/done.q") {
	last;
    }  else {
	usleep(50*1000); # in microseconds
    }
}
if (! -e "$dir/done.q") {
    print "Error\nTimeout\n";
    exit;
}

#else
my $dbh = Utils::get_dbh() || die "Cannot access database";
my @rows = ();
open(OUT, "<", "$dir/ublast.out") || die "No $dir/ublast.out file";
while(<OUT>) {
    chomp;
    my @F = split /\t/, $_;
    push @rows, \@F;
    last if @rows >= $maxhits;
}
close(OUT);

unlink("$dir/ublast.out");
unlink("$dir.x");
unlink("$dir.q");
unlink("$dir/faa");
unlink("$dir/done.q");
rmdir("$dir");

print join("\t", qw{orgId organism locusId sysName name description identity coverage evalue bits minfit maxfit minT maxT maxcofit})."\n";
my $orginfo = Utils::orginfo($dbh);
foreach my $row (@rows) {
    my ($query,$locusspec,$identity,$alnlen,$mm,$gap,$qBeg,$qEnd,$lBeg,$lEnd,$evalue,$bits) = @$row;
    my ($orgId,$locusId) = split /:/, $locusspec;
    die "Invalid locus $locusspec" unless defined $orgId && defined $locusId;
    my ($sysName,$name,$desc) = $dbh->selectrow_array("SELECT sysName, gene, desc FROM Gene WHERE orgId = ? AND locusId = ?",
					  {}, $orgId, $locusId);
    $sysName = "" if !defined $sysName;
    $name = "" if !defined $name;
    my ($minFit,$maxFit,$minT,$maxT) = $dbh->selectrow_array(qq{ SELECT min(fit), max(fit), min(t), max(t)
                                                                 FROM GeneFitness WHERE orgId = ? AND locusId = ? ; },
							     {}, $orgId, $locusId);
    $minFit = "" unless defined $minFit;
    $maxFit = "" unless defined $maxFit;
    $minT = "" unless defined $minT;
    $maxT = "" unless defined $maxT;
    my $maxCofit = "";
    if (defined $minFit) {
	$maxCofit = $dbh->selectrow_array(qq{ SELECT cofit FROM Cofit WHERE orgId = ? AND locusId = ? AND rank = 1 LIMIT 1; },
					  {}, $orgId, $locusId);
    }
    print join("\t",$orgId, $orginfo->{$orgId}{genome},
	       $locusId, $sysName, $name, $desc,
	       $identity,$alnlen/length($seq),$evalue,$bits,
	       $minFit,$maxFit,$minT,$maxT,$maxCofit)."\n";
}

$dbh->disconnect();
exit(0);

