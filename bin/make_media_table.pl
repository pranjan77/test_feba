#!/usr/bin/perl -w
use strict;
use Getopt::Long;
#use FindBin qw($Bin);
use lib "feba/lib";
use FEBA_Utils;

my @validUnits = qw{g/L mg/L ug/L M mM uM % ml/L};
my $unitString = join(", ", @validUnits);

my $usage = <<END
make_media_table.pl -media media -compounds FEBA_COMPOUND_sheet [ -out . ] [ -db feba.db ]
	compounds is tab-delimited with the fields unique id, CAS, source, catalog #, molecular weight,
	and also contains a Synonyms field.

	media has a single field with the media named followed by lines of the form
	compound name, concentration, concentration units,
	where valid units are $unitString.

        writes db.Compounds and db.MediaComponents to output directory.
        If -db is specified, also loads these into the sqlite3 database.
END
    ;

sub SynToKey($);

{
    my ($mediaFile, $compoundsFile,$db);
    my $outdir = ".";

    die $usage unless GetOptions('media=s' => \$mediaFile,
                                 'compounds=s' => \$compoundsFile,
                                 'out=s' => \$outdir,
                                 'db=s' => \$db)
        && defined $mediaFile && defined $compoundsFile
        && @ARGV == 0;
    die "Not a directory: $outdir" unless -d $outdir;
    die "No such file: $db" if defined $db && ! -e $db;

    my @compounds = FEBA_Utils::ReadTable($compoundsFile, qw{Synonyms});
    my @headers = FEBA_Utils::ReadColumnNames($compoundsFile);
    die "No rows in $compoundsFile" unless @compounds > 0;
    die "Not enough fields in $compoundsFile" unless @headers >= 5;
    my $keycol = $headers[0];
    my $cascol = $headers[1];
    my $mwcol = $headers[4];

    my %synonyms = (); # synonym, in lower case and with various characters removed => controlled vocabulary
    my %compounds = (); # vocabulary => row
    
    my $DbCompoundsFile = "$outdir/db.Compounds";
    open(COMPOUNDS, ">", $DbCompoundsFile) || die "Cannot write to $DbCompoundsFile";
    foreach my $row (@compounds) {
        my $key = $row->{$keycol};
        die "Duplicate compound id $key" if exists $compounds{$key};
        my $mw = $row->{$mwcol};
        $mw = "" if $mw eq "NA";
        die "Invalid weight $mw for compound $key in $compoundsFile\n"
            unless $mw eq "" || $mw =~ m/^[0-9]+[.]?[0-9]*$/;
        my $cas = $row->{$cascol};
        $cas =~ s/ +$//;
        $cas =~ s/^ +//;
        die "Invalid cas number '$cas' for compound $key in $compoundsFile\n"
            unless $cas eq "" || $cas eq "NA" || $cas =~ m/^\d+[0-9-]*\d$/;
        $compounds{$key} = $row;
        print COMPOUNDS join("\t", $key, $mw, $cas)."\n";

        my @syns = $key;
        push @syns, split /, /, $row->{Synonyms};
        foreach my $syn (@syns) {
            $syn = SynToKey($syn);
            next if $syn eq "";
            die "Conflicting synonyms for $syn after lower-casing and stripping whitespace and special characters: $key vs. $synonyms{$syn}"
                if exists $synonyms{$syn} && $synonyms{$syn} ne $key;
            $synonyms{$syn} = $key;
        }
    }
    close(COMPOUNDS) || die "Error writing to $DbCompoundsFile";
    print STDERR "Wrote $DbCompoundsFile\n";

    my %validUnits = map { $_ => 1 } @validUnits;

    open(MEDIA, "<", $mediaFile) || die "Cannot read $mediaFile";
    my %media = (); # media => list of [ compound_id, number, units ]
    my ($COMPOUND,$NUMBER,$UNITS) = 0..2;

    my %noSyn = (); # compounds with no match
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
            my ($compound, $concentration, $units) = @F;

            # check if compound is known, and replace synonyms
            if (!exists $compounds{$compound}) {
                my $syn = SynToKey($compound);
                if (exists $synonyms{$syn}) {
                    $compound = $synonyms{$syn};
                } else {
                    $noSyn{$compound} = 1;
                }
            }

            # validate concentration
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
                    print STDERR "Reuse\t$compound\t$row->[$NUMBER] $row->[$UNITS]\t$concentration $units\tin\t$curMedia\n";
                    last;
                }
            }

            push @{ $media{$curMedia} }, [ $compound, $concentration, $units ];
        } else {
            die "Wrong number of fields in\n$line\n";
        }
    }
    close(MEDIA) || die "Error reading $mediaFile";

    my $nNoSyn = scalar(keys %noSyn);
    if ($nNoSyn > 0) {
        print STDERR "Unrecognized compounds: $nNoSyn\n";
        print STDERR join("\t", sort keys %noSyn)."\n";
    }

    my @undefMedia = ();
    while (my ($media,$components) = each %media) {
        if (scalar(@$components) == 0) {
            push @undefMedia, $media;
        }
    }
    my $nUndefMedia = scalar(@undefMedia);
    if ($nUndefMedia > 0) {
        print STDERR "Media with no definitions: $nUndefMedia\n";
        print STDERR join("\t", sort @undefMedia)."\n";
    }

    # write out all the media definitions
    my $DbComponentsFile = "$outdir/db.MediaComponents";
    open(COMPONENTS, ">", $DbComponentsFile) || die "Cannot write to $DbComponentsFile";
    while (my ($media,$components) = each %media) {
        foreach my $row (@$components) {
            my ($compound,$concentration,$units) = @$row;
            # database schema allows for concentration or units to be missing, but this script does not
            print COMPONENTS join("\t", $media, $compound, $concentration, $units)."\n";
        }
    }
    close(COMPONENTS) || die "Error writing to $DbComponentsFile";
    print STDERR "Wrote $DbComponentsFile\n";

    if (defined $db) {
        print STDERR "Loading tables into $db\n";
        my @commands = (".bail on",
                        ".mode tabs",
                        "DELETE FROM Compounds;",
                        ".import $DbCompoundsFile Compounds",
                        "DELETE FROM MediaComponents;",
                        ".import $DbComponentsFile MediaComponents",
                        ".headers on",
                        "SELECT count(distinct media) nMedia, count(*) nTotalComponents from MediaComponents;",
                        "SELECT count(*) nCompounds from Compounds;"
            );
        open(SQLITE, "|-", "sqlite3", "$db") || die "Cannot run sqlite3 on $db";
        foreach my $command (@commands) {
            print SQLITE "$command\n";
        }
        close(SQLITE) || die "Error running sqlite3";
    }
}

sub SynToKey($) {
    my ($syn) = @_;
    $syn =~ s/[^a-zA-Z0-9+-]//g;
    return(lc($syn));
}
