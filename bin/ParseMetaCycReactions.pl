#!/usr/bin/perl -w
# Convert the metacyc reactions.dat file into tab-delimited format
# Tested with MetaCyc 19.5 -- note this uses only what is in the metacyc specific database
# (not all the organism-specific databases)
use strict;

die "Run as a filter: ParseMetaCycReactions.pl < metacyc/reactions.dat > metacyc_reactions.tab\n"
    unless @ARGV == 0;

my $nReaction = 0;
my %attr = ();

my $nReactions = 0;
while(my $line = <STDIN>) {
    chomp $line;
    next if $line =~ m/^#/;
    next if $line =~ m/^\^/;
    if ($line eq "//") {
        my $id = $attr{"UNIQUE-ID"};
        die "No id" unless defined $id && $id ne "";
        my $desc = $attr{"COMMON-NAME"} || $attr{"COMMENT"} || "";
        if ($desc ne "") {
            print join("\t", $id, $desc)."\n";
            $nReactions++;
        }
        %attr = ();
    } elsif ($line =~ m!^/!) {
        next; # skip extensions to comments
    } else {
        $line =~ m/^([a-zA-Z?-]+[0-9]*) - (.*)$/ || die "Cannot parse $line";
        $attr{$1} = $2;
    }
}
    
