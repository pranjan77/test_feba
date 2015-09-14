# Compounds.pm -- utilities for reading the media and compounds metadata from (by default)
# feba/metadata/FEBA_COMPOUND_sheet and feba/metadata/media.
#
# It maintains a lot of state in global variables.

package Compounds;
require Exporter;
use strict;
use FEBA_Utils;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(LoadCompounds
             GetCompoundList
             FindCompound GetCompoundCAS GetCompoundMW
             ListValidUnits
             LoadMedia
             GetMedias GetMediaComponents GetUndefMedia GetUnknownComponents WarnReusedComponents);

sub LoadCompounds($); # compound directory, no return value
sub FindCompound($); # compound or synonym => compound or undef
sub GetCompoundCAS($); # compound => CAS or ""
sub GetCompoundMW($); # compound => MW or ""
sub GetCompoundList(); # returns the list of compounds

# From a compound name or synonym to a key to look up.
# Only alphanumeric characters, +, or - are considered -- everything else is removed. Also, case is ignored.
# Ignoring whitespace and case means that many variant names for compounds can be handled without introducing more synonyms
sub SynToKey($);

# Local variables for compound information
my %compounds = (); # compound to list of [ compound name, cas, MW ]; cas or MW are "" (not undef) if unknown
my ($I_COMPOUND,$I_CAS,$I_MW) = 0..2;
my %synonyms = (); # processed synonym (from SynToKey) => compound name


# Local variables for media information
my %media = (); # media => list of [ compound_id, number, units ]
my %unknownComponents = (); # media components with no match in the compounds table
my @undefMedia = (); # media with no definition
my %reuseComponents = (); # component => media => 1 if it is reused

sub LoadCompounds($) {
    my ($compoundsDir) = @_;
    my $compoundsFile = "$compoundsDir/FEBA_COMPOUND_sheet";
    my @compounds = FEBA_Utils::ReadTable($compoundsFile, qw{Synonyms});
    my @headers = FEBA_Utils::ReadColumnNames($compoundsFile);
    die "No rows in $compoundsFile" unless @compounds > 0;
    die "Not enough fields in $compoundsFile" unless @headers >= 5;
    my $keycol = $headers[0];
    my $cascol = $headers[1];
    my $mwcol = $headers[4];

    foreach my $row (@compounds) {
        my $compound = $row->{$keycol};
        die "Duplicate compound id $compound" if exists $compounds{$compound};
        my $mw = $row->{$mwcol};
        $mw = "" if $mw eq "NA";
        die "Invalid weight $mw for compound $compound in $compoundsFile\n"
            unless $mw eq "" || $mw =~ m/^[0-9]+[.]?[0-9]*$/;
        my $cas = $row->{$cascol};
        $cas =~ s/ +$//;
        $cas =~ s/^ +//;
        $cas = "" if $cas eq "NA";
        die "Invalid cas number '$cas' for compound $compound in $compoundsFile\n"
            unless $cas eq "" || $cas =~ m/^\d+[0-9-]*\d$/;
        $compounds{$compound} = [ $compound, $cas, $mw ];

        my @syns = $compound;
        push @syns, split /, /, $row->{Synonyms};
        foreach my $syn (@syns) {
            $syn = SynToKey($syn);
            next if $syn eq "";
            die "Conflicting synonyms for $syn in $compoundsFile\n"
                . " after lower-casing and stripping whitespace and special characters: $compound vs. $synonyms{$syn}"
                if exists $synonyms{$syn} && $synonyms{$syn} ne $compound;
            $synonyms{$syn} = $compound;
        }
    }
}

sub FindCompound($) {
    my ($syn) = @_;
    return $syn if exists $compounds{$syn};
    my $key = SynToKey($syn);
    return $synonyms{$key} if exists $synonyms{$key};
    return undef;
}

sub GetCompoundCAS($) {
    my ($compound) = @_;
    return exists $compounds{$compound} ? $compounds{$compound}[$I_CAS] : "";
}

sub GetCompoundMW($) {
    my ($compound) = @_;
    return exists $compounds{$compound} ? $compounds{$compound}[$I_MW] : "";
}

sub GetCompoundList() {
    return sort keys %compounds;
}

sub SynToKey($) {
    my ($syn) = @_;
    $syn =~ s/[^a-zA-Z0-9+-]//g;
    return(lc($syn));
}

### Media functions

sub LoadMedia($) {
    my ($metadir) = @_;
    my $mediaFile = "$metadir/media";
    open(MEDIA, "<", $mediaFile) || die "Cannot read $mediaFile";
    my ($COMPOUND,$NUMBER,$UNITS) = 0..2;

    my %validUnits = map { $_ => 1 } ListValidUnits();

    my $curMedia = undef;
    while(my $line = <MEDIA>) {
        chomp $line;
        $line =~ s/#.*$//; # strip comments
        $line =~ s/\t+$//; # strip trailing fields that are empty (note this means units *must* be present)
        my @F = split /\t/, $line;
        
        if (scalar(@F) == 0) {
            next;
        } elsif (scalar(@F) == 1) {
            $curMedia = $F[0];
            $curMedia =~ s/ +$//;
            die "Duplicate media entry for $curMedia" if exists $media{$curMedia};
            $media{$curMedia} = [];
        } elsif (scalar(@F) == 3) {
            die "No media id yet at line:\n$line\nin $mediaFile" if !defined $curMedia;
            my ($compoundSyn, $concentration, $units) = @F;

            # check if compound is known, and replace synonyms
            my $compound = FindCompound($compoundSyn);
            if (!defined $compound || $compound eq "") {
                $compound = $compoundSyn;
                $unknownComponents{$compoundSyn} = 1;
            }

            $concentration =~ s/^ +//;
            $concentration =~ s/ +$//;
            $concentration eq "" || $concentration =~ m/^\d+$/
                || $concentration =~ m/^\d+[.]\d*$/
                || $concentration =~ m/^\d+[.]?\d*e[+-]\d+$/
                || die "Invalid concentration $concentration in line\n$line\nfor $curMedia\nin $mediaFile";

            # validate units
            $units =~ s/^ +//;
            $units =~ s/ +$//;
            die $line if $units eq "";
            die "Invalid unit $units for\n$line\nin $curMedia, $mediaFile"
                unless exists $validUnits{$units};

            # check for duplicate entries
            foreach my $row (@{ $media{$curMedia} }) {
                if ($row->[$COMPOUND] eq $compound) {
                    $reuseComponents{$compound}{$curMedia} = 1;
                    last;
                }
            }

            push @{ $media{$curMedia} }, [ $compound, $concentration, $units ];
        } else {
            die "Wrong number of fields in\n$line\n";
        }
    }
    close(MEDIA) || die "Error reading $mediaFile";

    while (my ($media,$components) = each %media) {
        if (scalar(@$components) == 0) {
            push @undefMedia, $media;
        }
    }
}

sub GetMedias() {
    return sort keys %media;
}

# returns undef for unknown media; otherwise, a reference to a list of [compound,concentration,units]
sub GetMediaComponents($) {
    my ($media) = @_;
    return $media{$media};
}

sub GetUnknownComponents() {
    return sort keys %unknownComponents;
}

sub GetUndefMedia() {
    return @undefMedia;
}

sub WarnReusedComponents() {
    foreach my $compound (sort keys %reuseComponents) {
        foreach my $media (keys %{ $reuseComponents{$compound} }) {
            my @components = grep { $_->[0] eq $compound } @{ $media{$media} };
            die "Not dup compound $compound media $media: " . scalar(@components) if scalar(@components) < 2;
            my @show = map { $_->[1] . " " . $_->[2] } @components;
            print "Multiple use of $compound in $media: " . join(", ", @show)."\n";
        }
    }
}

sub ListValidUnits() {
    return qw{g/L mg/L ug/L M mM uM % ml/L};
}

1;
