
package FEBA_Utils;
require Exporter;
use strict;
use File::stat;
our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw( ReadTable ReadColumnNames ReadFasta NewerThan ReadFastaDesc );

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
	$line =~ s/[\r\n]+$//;
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

# Returns a hash containing either "error"
# or hashes of "seq", "desc", and "len"
sub ReadFastaDesc($) {
    my ($file) = @_;
    my %seq = ();
    my %desc = ();
    my $name = undef;
    open(FAA, "<", $file) || return('error' => "Cannot read $file" );
    while(<FAA>) {
        s/[\r\n]+$//;
        if (m/^>(.*)$/) {
            my $header = $1;
            if ($header =~ m/^(\S+)\s+(\S.*)$/) {
                $name = $1;
                $desc{$name} = $2;
            } else {
                return('error' => "bad header for sequence:\n$header\n") unless $header =~ m/^\S+$/;
                $name = $header;
                $desc{$name} = $header;
            }
            return('error' => "Duplicate sequence id:\n$name\n") if exists $seq{$name};
            $seq{$name} = "";
        } else {
            return('error' => "sequence before header:\n$_\n") unless defined $name;
            s/\s//g;
            $seq{$name} .= $_;
        }
    }
    close(FAA) || return('error' => "Error reading $file");
    my %len = ();
    while (my ($name,$seq) = each %seq) {
        $len{$name} = length($seq);
        return('error' => "No sequence for id:\n$name\n") if ($len{$name} == 0);
    }
    return("desc" => \%desc, "seq" => \%seq, "len" => \%len);
}

sub NewerThan($$) {
    my ($file1, $file2) = @_;
    die "Invalid arguments to NewerThan" unless defined $file1 && defined $file2;
    die "No such file: $file2" unless -e $file2;
    return 0 unless -e $file1;
    return stat($file1)->mtime >= stat($file2)->mtime ? 1 : 0;
}

1;
