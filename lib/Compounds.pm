# Compounds.pm -- utilities for reading the media and compounds metadata from (by default)
# feba/metadata/Compounds.tsv, feba/metadata/media, feba/metadata/mixes
#
# It maintains a lot of state in global variables.

package Compounds;
require Exporter;
use strict;
use FEBA_Utils; # for ReadTable(), ReadColumnNames()

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(LoadCompounds
             GetCompoundList
             FindCompound GetCompoundCAS GetCompoundMW
             ListValidUnits
             LoadMedia
             GetMedias GetMediaComponents GetUndefMedia GetUnknownComponents WarnReusedComponents
             SynToKey GetSynonymMap);

sub LoadCompounds($); # metadata directory; no return value
sub FindCompound($); # compound or synonym => compound or undef
sub GetCompoundCAS($); # compound => CAS or ""
sub GetCompoundMW($); # compound => MW or ""
sub GetCompoundList(); # returns the list of compounds

sub LoadMedia($); # metadata directory; no return value
sub GetMedias(); # returns list of media names that have definitions
sub GetUndefMedia(); # returns list of media names that lack definitions
sub GetMediaAttributes($); # media name => hash of attributes such as Description or Minimal
sub GetMediaComponents($); # for each component, a list of [ compound, concentration, units, mix ]
#	where mix indicates which mix it was from (if any)
sub WarnReusedComponents(); # report to STDOUT on components that are in a medium more than once

# From a compound name or synonym to a key to look up.
# Only alphanumeric characters, +, or - are considered -- everything else is removed. Also, case is ignored.
# Ignoring whitespace and case means that many variant names for compounds can be handled without introducing more synonyms
sub SynToKey($);

# Local variables for compound information
my %compounds = (); # compound to list of [ compound name, cas, MW ]; cas or MW are "" (not undef) if unknown
my ($I_COMPOUND,$I_CAS,$I_MW) = 0..2;
my %synonyms = (); # processed synonym (from SynToKey) => compound name
sub GetSynonymMap() { return \%synonyms; }

# Local variables for media information
my %media = (); # media => list of [ compound_id, number, units ]
my %mediaAttr = (); # media => Description or Minimal => value
my %mix = (); # mix => list of [ compound_id, number, units ]
my %mixAttr = (); # mix => Description or X => value
my %unknownComponents = (); # media components with no match in the compounds table
my @undefMedia = (); # media with no definition
my %reuseComponents = (); # component => media => 1 if it is reused

sub LoadCompounds($) {
    my ($compoundsDir) = @_;
    my $compoundsFile = "$compoundsDir/Compounds.tsv";
    my @req = qw{Compound CAS FW Synonyms};
    my @compounds = &ReadTable($compoundsFile, @req);
    my @headers = &ReadColumnNames($compoundsFile);
    die "No rows in $compoundsFile" unless @compounds > 0;

    foreach my $row (@compounds) {
        my $compound = $row->{"Compound"};
        die "Duplicate compound id $compound" if exists $compounds{$compound};
        my $mw = $row->{"FW"};
        $mw = "" if $mw eq "NA";
        die "Invalid weight $mw for compound $compound in $compoundsFile\n"
            unless $mw eq "" || $mw =~ m/^[0-9]+[.]?[0-9]*$/;
        my $cas = $row->{"CAS"};
        $cas =~ s/ +$//;
        $cas =~ s/^ +//;
        $cas = "" if $cas eq "NA";
        die "Invalid cas number '$cas' for compound $compound in $compoundsFile\n"
            unless $cas eq "" || $cas =~ m/^\d+[0-9-]*\d$/;
        $compounds{$compound} = [ $compound, $cas, $mw ];

        my @syns = $compound;
        push @syns, split /, /, $row->{Synonyms};
        foreach my $syn (@syns) {
            my $key = SynToKey($syn);
            next if $key eq "";
            print "Warning: non-unique synonym $syn: $compound or $synonyms{$key}\n"
                if exists $synonyms{$key} && $synonyms{$key} ne $compound;
            $synonyms{$key} = $compound;
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

# Returns two hashes:
# media => list of components
# media => attributes
# Handles both the media file and the mixes file
# Does NOT look up compound synonyms or otherwise check the results
sub LoadMedia($) {
    my ($metadir) = @_;
    my $mediaFile = "$metadir/media";
    my $mixFile = "$metadir/mixes";
    die "No such file: $mediaFile\n" unless -e $mediaFile;
    die "No such file: $mixFile\n" unless -e $mixFile;
    my ($mediaC,$mediaA) = ParseMediaFile($mediaFile);
    my ($mixC,$mixA) = ParseMediaFile($mixFile);
    %media = %$mediaC;
    %mediaAttr = %$mediaA;
    %mix = %$mixC;
    %mixAttr = %$mixA;

    # Validation:
    # Each media must have a Description
    # media should not also be mixes
    while (my ($media, $attr) = each %mediaAttr) {
        die "No Description for media $media" unless exists $attr->{Description} && $attr->{Description} ne "";
        die "Media $media is also a mix" if exists $mix{$media};
    }
    # Each mix must have a Description and a numeric X value
    # mixes should not also be media
    while (my ($mix, $attr) = each %mixAttr) {
        die "No Description for mix $mix" unless exists $attr->{Description} && $attr->{Description} ne "";
        die "Invalid X for mix $mix" unless exists $attr->{X} && $attr->{X} =~ m/^[0-9]+[.]?[0-9]*$/;
        die "Mix $mix is also a media" if exists $media{$mix};
    }
    # Replace compound synonyms with compounds, and record any that are not known or are duplicates
    while (my ($media, $list) = each %media) {
        SetupComponentList($media, $list);
    }
    while (my ($mix, $list) = each %mix) {
        SetupComponentList($mix, $list);
    }
}

sub SetupComponentList($$) {
    my ($media, $list) = @_;
    my $isMedia = exists $media{$media};
    my $isMix = exists $mix{$media};
    die "Unknown $media" unless $isMedia || $isMix;

    my $COMPOUND = 0;
    # transfer synonyms or record that it is unknown
    foreach my $row (@$list) {
        my ($orig,$undef,$units) = @$row;
        if ($isMedia && exists $mix{$orig}) {
            # leave as is, but check that X is specified
            die "Mix must be included with X units" unless $units eq "X";
        } else {
            my $compound = FindCompound($orig);
            if (defined $compound) {
                $row->[$COMPOUND] = $compound;
            } else {
                $unknownComponents{$row->[$COMPOUND]} = 1;
            }
        }
    }
    # record repeat entries
    my %seen = ();
    foreach my $row (@$list) {
        my $compound = $row->[$COMPOUND];
        $reuseComponents{$compound}{$media} = 1 if exists $seen{$compound};
        $seen{$compound} = 1;
    }
}

sub ParseMediaFile($) {
    my ($mediaFile) = @_;
    my %comp = ();
    my %attr = ();
    open(MEDIA, "<", $mediaFile) || die "Cannot read $mediaFile";
    my ($COMPOUND,$NUMBER,$UNITS) = 0..2;

    my %validUnits = map { $_ => 1 } ListValidUnits();
    my %validAttr = map { $_ => 1 } qw{Description Minimal X};

    my $curMedia = undef;
    my $readingCompounds = 0;
    while(my $line = <MEDIA>) {
        $line =~ s/[\r\n]+$//; # handle DOS mode files
        $line =~ s/\t+$//; # strip trailing fields that are empty (note this means units *must* be present)
        my @F = split /\t/, $line;
        
        if (scalar(@F) == 0) {
            $curMedia = undef; # blank lines end media descriptions
        } elsif ($F[0] =~ m/^#/) {
            # skip comment line
            ;
        } elsif (scalar(@F) == 2) {
            my ($attr,$value) = @F;
            if ($attr eq "Media") {
                $curMedia = $F[1];
                $curMedia =~ s/ +$//;
                die "Duplicate media entry for $curMedia" if exists $comp{$curMedia};
                $comp{$curMedia} = [];
                $attr{$curMedia} = {};
                $readingCompounds = 0;
            } elsif (exists $validAttr{$attr}) {
                die "No media id yet at line:\n$line\nin $mediaFile" if !defined $curMedia;
                die "Duplicate attr $attr for media $curMedia" if exists $attr{$curMedia}{$attr};
                $attr{$curMedia}{$attr} = $value;
            } else {
                die "Invalid media attribute $F[0]";
            }
        } elsif (scalar(@F) == 3) {
            die "No media id yet at line:\n$line\nin $mediaFile" if !defined $curMedia;
            if ($F[0] =~ m/^Controlled/ && $F[1] eq "Concentration" && $F[2] eq "Units") {
                $readingCompounds = 1;
            } else {
                die "No compounds header for $curMedia at\n$line\n..." unless $readingCompounds == 1;
                my ($compound, $concentration, $units) = @F;
                $compound =~ s/ +$//; # remove trailing spaces
                $concentration =~ s/^ +//;
                $concentration =~ s/ +$//;
                $concentration eq "" || $concentration =~ m/^\d+$/
                    || $concentration =~ m/^\d+[.]\d*$/
                    || $concentration =~ m/^\d+[.]?\d*[eE][+-]\d+$/
                    || die "Invalid concentration $concentration in line\n$line\nfor $curMedia\nin $mediaFile";
                $units =~ s/^ +//;
                $units =~ s/ +$//;
                die $line if $units eq "";
                die "Invalid unit $units for\n$line\nin $curMedia, $mediaFile"
                    unless exists $validUnits{$units};
                push @{ $comp{$curMedia} }, [ $compound, $concentration, $units ];
            }
        } else {
            die "Wrong number of fields in\n$line\n";
        }
    }
    close(MEDIA) || die "Error reading $mediaFile";
    return (\%comp, \%attr);
}

sub GetMedias() {
    return sort keys %media;
}

# returns undef for unknown media; otherwise, a reference to a list of [compound,concentration,units,mix]
# where mix is empty (not undefined) unless the compound was included indirectly via a mix
sub GetMediaComponents($) {
    my ($media) = @_;
    return undef if !exists $media{$media};
    my $out = [];
    foreach my $row (@{ $media{$media} }) {
        my ($comp,$conc,$units) = @$row;
        if (exists $mix{$comp}) {
            die "Units for mix $comp in media $media are not X" unless $units eq "X";
            die "No X value for mix $comp" unless exists $mixAttr{$comp}{"X"};
            my $rel = $conc / $mixAttr{$comp}{"X"};
            die "Invalid relative X $rel for $comp in $media" unless $rel > 0 && $rel < 1e4;
            foreach my $row2 (@{ $mix{$comp} }) {
                my ($comp2, $conc2, $units2) = @$row2;
                push @$out, [ $comp2, $conc2 * $rel, $units2, $comp ];
            }
        } else {
            push @$out, [ $comp, $conc, $units, "" ];
        }
    }
    return $out;
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
    return qw{g/L mg/L ug/L M mM uM vol% ml/L X};
}

1;
