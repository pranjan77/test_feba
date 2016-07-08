#!/usr/bin/perl -w
# Convert the metacyc uniprot-seq-ids.seq into fasta format
# Tested with MetaCyc 19.5 -- note this uses only what is in the metacyc specific database
# (not all the organism-specific databases)
use strict;
die "Usage: ParseMetaCycSeq.pl < metacyc_dir/uniprot-seq-ids.seq > output_file\n" unless @ARGV==0;

my $nSeq = 0;
while(my $id = <STDIN>) {
    chomp $id;
    next if $id eq "";
    $id =~ s/^[ (]+//;
    my $seq = undef;
    # in case id and sequence are on the same line
    if ($id =~ m/ /) {
        ($id,$seq) = split / /, $id;
    }
    die "Cannot parse $id" unless $id =~ m/^"([A-Za-z0-9]+)"$/;
    $id = $1;
    $seq = <STDIN> if !defined $seq;
    $seq =~ m/^ *"([A-Z]+)"[)]+$/ || die "Cannot parse sequence line $seq";
    $seq = $1;
    print ">$id\n$seq\n";
    $nSeq++;
}
print STDERR "Parsed $nSeq sequences\n";
