#!/usr/bin/perl -w
# A helper script for db_setup.pl -- populates KEGG related tables from besthit.kegg
use Getopt::Long;
use strict;

my $db;
my $kegghits = "besthit.kegg";
my $dir = ".";
my $usage = <<END
Usage: db_setup.pl -db feba.db [ -kegghits $kegghits ] [ -dir $dir ]

Given the output of KeggBestHit.pl and a sqlite3 database, populates
the KEGG related tables. Will delete any existing data in those
tables, but the tables must already exist.

-dir specifies where temporary files will be stored
END
    ;

die $usage unless (GetOptions('kegghits=s' => \$kegghits,
                              'db=s' => \$db,
                              'dir=s' => \$dir)
                   && @ARGV == 0 && defined $db);
die "No such file: $kegghits" unless -e $kegghits;
die "No such directory: $dir" unless -d $dir;
die "No such file: $db" unless -e $db;

my %keggSeen = (); # kegg org => kegg id => 1 if already output
my %kgroupSeen = (); # kegg group => 1 if already output

open(HITS, "<", $kegghits) || die "Cannot read $kegghits";
open(BH, ">", "$dir/db.BestHitKEGG") || die "Cannot write to $dir/db.BestHitKEGG";
open(MEMBER, ">", "$dir/db.KEGGMember") || die "Cannot write to $dir/db.KEGGMember";
open(DESC, ">", "$dir/db.KgroupDesc") || die "Cannot write to $dir/db.KgroupDesc";
open(EC, ">", "$dir/db.KgroupEC") || die "Cannot write to $dir/db.KgroupEC";

while(<HITS>) {
    chomp;
    my ($locusSpec, $keggSpec, $identity, $kgroups, $kdescs) = split /\t/, $_, -1;
    my ($orgId,$locusId) = split /:/, $locusSpec;
    die unless $orgId ne "" && $locusId ne "";
    my ($keggOrg,$keggId) = split /:/, $keggSpec;
    die unless $keggOrg ne "" && $keggId ne "";
    die "Invalid identity $identity" unless $identity =~ m/^[0-9.]+$/;
    print BH join("\t", $orgId, $locusId, $keggOrg, $keggId, $identity)."\n";
    next if $kgroups eq "" && $kdescs eq "";
    next if exists $keggSeen{$keggOrg}{$keggId};
    $keggSeen{$keggOrg}{$keggId} = 1;
    my @kgroups = split /,/, $kgroups, -1;
    my @kdescs = split /; /, $kdescs, -1;
    @kdescs = ("") if $kdescs eq "";
    die "Non-matching numbers for $locusSpec" unless scalar(@kgroups) == scalar(@kdescs);
    while (@kgroups) {
        my $kgroup = shift @kgroups;
        my $kdesc = shift @kdescs;
        print MEMBER join("\t", $keggOrg, $keggId, $kgroup)."\n";
        next if exists $kgroupSeen{$kgroup};
        $kgroupSeen{$kgroup} = 1;
        my @ec = ();
        if ($kdesc =~ m/\[EC:([0-9. -]+)\]/) {
            @ec = split / /, $1;
            $kdesc =~ s/\[EC:[0-9. -]+\]//;
            foreach my $ec (@ec) {
                die "$locusSpec $ec" unless $ec =~ m/^[0-9.-]+$/;
                print EC join("\t",$kgroup,$ec)."\n";
            }
        }
        print DESC join("\t", $kgroup, $kdesc)."\n";
    }
}

close(HITS) || die "Error reading $kegghits";
close(BH) || die "Error writing to $dir/db.BestHitKEGG";
close(MEMBER) || die "Error writing to $dir/db.KEGGMember";
close(DESC) || die "Error writing to $dir/db.KgroupDesc";
close(EC) || die "Error writing to $dir/db.KgroupEc";

print STDERR "Wrote temporary files in $dir/:\ndb.BestHitKEGG db.KEGGMember db.KgroupDesc db.KgroupEc\n";

open(SQLITE, "|-", "sqlite3", $db) || die "Cannot run sqlite3 on $db";
print SQLITE ".mode tabs\n";

foreach my $table (qw{BestHitKEGG KEGGMember KgroupDesc KgroupEC}) {
    print SQLITE "DELETE from $table;\n";
    print SQLITE ".import $dir/db.$table $table\n";
}
close(SQLITE) || die "Error running sqlite3 on db";

print STDERR "Done, removing temporary files\n";
foreach my $table (qw{BestHitKEGG KEGGMember KgroupDesc KgroupEC}) {
    unlink("db.$table");
}

