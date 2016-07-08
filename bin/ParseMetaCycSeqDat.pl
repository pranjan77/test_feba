#!/usr/bin/perl -w
# Convert the metacyc uniprot-seq-ids.dat into a simple tab-delimited format of
# rxn_id	ec_number	sequence_id
#
# Tested with MetaCyc 19.5 -- note this uses only what is in the metacyc specific database
use strict;
die "Usage: ParseMetaCycSeqDat.pl < metacyc/uniprot-seq-ids.dat > metacyc_seq_rxn.tab\n"
    unless @ARGV == 0;

my $rxn = undef;
my $ec = undef;
my $nLinks = 0;
while(my $line = <STDIN>) {
    chomp $line;
    next if $line eq "";
    next if $line =~ m/^;/; # comments
    $line =~ s/^ +//;
    my @F = split / /, $line;
    while (@F > 0) {
        my $field = shift @F;
        if ($field =~ m/^[(]/) {
            # start of a rxn
            $rxn = $field;
            $rxn =~ s/^[(]+//;
            # sometimes quoted with | characters
            $rxn = $1 if $rxn =~ m/^[|](.*)[|]$/;
            $rxn =~ s/[()]//g; # not sure why extra parentheses around words appear, but link works if they are removed
            # i.e. compare SULFITE-REDUCTASE-(FERREDOXIN)-RXN and
            # http://metacyc.org/META/NEW-IMAGE?type=REACTION&object=SULFITE-REDUCTASE-FERREDOXIN-RXN
            die "No ec for rxn $rxn" if scalar(@F) == 0;
            $ec = shift @F;
            $ec = $1 if $ec =~ m/^"(.*)"$/;
            $ec = $1 if $ec =~ m/^[|](.*)[|]$/;
            # metacyc does not use - on end -- instead, just is missing the trailing components sometimes
            die "Invalid ec: $ec for $rxn" unless $ec =~ m/^[0-9]+[.][0-9]+[.][0-9]+[.][0-9]+$/
                || $ec =~ m/^[0-9]+[.][0-9]+[.][0-9]+$/
                || $ec =~ m/^[0-9]+[.][0-9]+$/
                || $ec =~ m/^[0-9]+$/;
        } else {
            die "No rxn for $field" unless defined $rxn && defined $ec;
            my $lastid = $field =~ m/[)]+$/ ? 1 : 0;
            $field =~ s/[)]+$//;
            die "Invalid sequence id $field" unless $field =~ m/^"([A-Z][A-Z0-9]+)"$/;
            my $seqid = $1;
            print join("\t", $rxn, $ec, $seqid)."\n";
            $nLinks++;
        }
    }
}
print STDERR "Parsed $nLinks links for ids to reaction and EC number\n";
