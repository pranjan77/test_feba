
package FEBA_Utils;
require Exporter;
use strict;
our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw( ReadTable ReadColumnNames ReadFasta );

# filename and list of required fields => list of hashes, each with field->value
sub ReadTable($*) {
    my ($filename,@required) = @_;
    open(IN, "<", $filename) || die "Cannot read $filename";
    my $headerLine = <IN>;
    $headerLine =~ s/[\r\n]+$//; # for DOS
    my @cols = split /\t/, $headerLine;
    my %cols = map { $cols[$_] => $_ } (0..(scalar(@cols)-1));
    foreach my $field (@required) {
	die "No field $field in $filename" unless exists $cols{$field};
    }
    my @rows = ();
    while(my $line = <IN>) {
	chomp $line;
	my @F = split /\t/, $line, -1;
	die "Wrong number of columns in:\n$line\nin $filename"
	    unless scalar(@F) == scalar(@cols);
	my %row = map { $cols[$_] => $F[$_] } (0..(scalar(@cols)-1));
	push @rows, \%row;
    }
    close(IN) || die "Error reading $filename";
    return @rows;
}

# filename to list of column names
sub ReadColumnNames($) {
    my ($filename) = @_;
    open(IN, "<", $filename) || die "Cannot read $filename";
    my $line = <IN>;
    close(IN) || die "Error reading $filename";

    $line =~ s/[\r\n]+$//; # for DOS
    my @cols = split /\t/, $line;
    return @cols;
}

# returns a reference to a hash of name => sequence
sub ReadFasta ($) {
    my ($filename) = @_;
    open(IN, "<", $filename) || die "Cannot read $filename";
    my %seqs = ();
    my $name = undef;
    while(<IN>) {
	chomp;
	if (m/>(\S+)/) {
	    die "Duplicate entry for $1 in filename" if exists $seqs{$1};
	    $name = $1;
	} else {
	    die "Unexpected line $_" unless defined $name;
	    $seqs{$name} .= $_;
	}
    }
    close(IN) || die "Error reading $filename";
    return(\%seqs);
}

1;
