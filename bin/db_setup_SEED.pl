#!/usr/bin/perl -w
# Given a list of organisms, build files for loading the SEEDAnnotation and SEEDClass tables

use strict;
use Getopt::Long;
use FindBin qw($Bin);

{
    my $usage = qq{
Usage: db_setup_SEED.pl [ -gdir g ] [ -out . ] org1 ... orgN

Relies on the gdir/organism/seedanno.tab files, which can be produced from
the amino acid fasta file (aaseq) by using the SEED server, i.e.,
http://servers.nmpdr.org/sapling/server.cgi?code=server_paper_example6.pl

Writes to db.SEEDAnnotation and db.SEEDClass in the out directory
}
    ;
    my $gdir = "g";
    my $outdir = ".";
    GetOptions('gdir=s' => \$gdir,
               'outdir=s' => \$outdir)
        || die $usage;
    my @orgs = @ARGV;
    die "No organisms requested:\n$usage" unless scalar(@orgs) > 0;
    die "No such directory: $gdir" unless -d $gdir;
    die "No such directory: $outdir" unless -d $outdir;
    foreach my $org (@orgs) {
        die "No such directory: $gdir/$org" unless -d "$gdir/$org";
    }

    my $anno_file = "$outdir/db.SEEDAnnotation";
    my $class_file = "$outdir/db.SEEDClass";
    open(ANNO, ">", $anno_file) || die "Cannot write to $anno_file";
    open(CLASS, ">", $class_file) || die "Cannot write to $class_file";
    my %badclass = (); # unparseable strings
    foreach my $org (@orgs) {
        my $in_file = "$gdir/$org/seedanno.tab";
        if (! -e $in_file ) {
            print STDERR "Skipping $org -- no seedanno.tab file\n";
            next;
        }
        open(SEED, "<", $in_file) || die "Cannot read $in_file";
        # avoid duplicate rows: remember locusId => part that is already seen
        my %seen = ();
        while(<SEED>) {
            chomp;
            my ($locusId,$seed_desc) = split /\t/, $_;
            die "Invalid input line\n$_\nfrom $in_file"
                unless defined $locusId && defined $seed_desc;
            next if $seed_desc eq ""; # allow missing annotations
            print ANNO join("\t", $org, $locusId, $seed_desc)."\n";

            # and parse out EC and TC number(s)
            # These will be in the description, of the form (EC 2.4.1.129)
            # or (EC 3.4.-.-) or (TC 3.A.1.4.1)
            my @parts = $seed_desc =~ m/[(][ET]C [0-9A-Za-z.-]+[)]/g;
            foreach my $part (@parts) {
                $part =~ m/^[(](.*)[)]$/ || die;
                $part = $1; # remove parentheses
                my ($type,$num) = split / /, $part;
                if ($type eq "EC") {
                    $badclass{$part} = 1
                        unless $num =~ m/^[0-9][.][0-9-]+[.][0-9-]+[.][0-9a-zA-Z-]+/;
                } elsif ($type eq "TC") {
                    $badclass{$part} = 1
                        unless $num =~ m/[0-9][.][A-Z].[0-9]+[.][0-9]+[.][0-9A-Za-z]+$/;
                }
                print CLASS join("\t", $org, $locusId, ($type eq "EC" ? 1 : 2), $num)."\n"
                    unless exists $badclass{$part} || exists $seen{$locusId}{$part};
                $seen{$locusId}{$part} = 1;
            }
        }
        close(SEED) || die "Error reading $in_file";
    }
    close(ANNO) || die "Error writing $anno_file";
    close(CLASS) || die "Error writing $class_file";
    print STDERR "These unparseable classes were ignored:\n" . join("\n",keys %badclass)."\n"
        if scalar(keys(%badclass)) > 0;
    print STDERR "Wrote $anno_file\n";
    print STDERR "Wrote $class_file\n";
}
