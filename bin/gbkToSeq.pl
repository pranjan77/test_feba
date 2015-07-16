#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

{
    die "Usage: gbkToSeq.pl GenbankFile > fastaFile\n"
	unless @ARGV==1;
    my ($file) = @ARGV;

    my $obj = Bio::SeqIO->new(-file => $file, -format => "genbank");
    my $out = Bio::SeqIO->new(-format => "fasta");
    while(my $scaffold = $obj->next_seq) {
	my $sc2 = Bio::Seq->new(-display_id => $scaffold->display_id(), -seq => $scaffold->seq());
	$out->write_seq($sc2); # just the raw identifier (nothing else)
    }
}
