#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use FEBA_Utils qw{ReadFasta};

my $phobiusDir = "$RealBin/phobius";
my $usage = <<END
Usage: myPhobius.pl -in in.faa -out prefix

Writes to prefix.plp, prefix.R, prefix.tsv, and prefix.svg
The input fasta should contain just one sequence.

Options:
  -dir $phobiusDir -- the phobius directory
   It must contain the decodeanhmm executable,
   phobius.options, and phobius.model
END
;

my ($inFaa, $out);
die $usage unless GetOptions('in=s' => \$inFaa,
                             'out=s' => \$out,
                             'dir=s' => \$phobiusDir)
  && @ARGV == 0
  && defined $inFaa && defined $out;

die "No such directory: $phobiusDir\n" unless -d $phobiusDir;
my $decode = "$RealBin/phobius/decodeanhmm";
die "No such executable: $decode\n" unless -x $decode;
my $optFile = "$RealBin/phobius/phobius.options";
my $modelFile = "$RealBin/phobius/phobius.model";
foreach my $file ($optFile, $modelFile) {
  die "No such file: $file\n" unless -e $file;
}

my $seqs = ReadFasta($inFaa);
die "$inFaa should have exactly one sequence\n" unless keys(%$seqs) == 1;
my ($seq) = values(%$seqs);
my $seqLen = length($seq);
my ($id) = keys(%$seqs);
$id =~ s/ .*//;
$id =~ s/["']//g;

my $tmpFile = "$out.plp.$$.tmp";
system("($decode -f $optFile $modelFile -plp < $inFaa > $tmpFile) >& /dev/null") == 0
  || die "phobius decodeanhmm failed: $!";

# Convert the decodeanhmm output, which will look like
# # BT2157
#         i         o         O         M         n         h         c         C   
# M      0.00002   0.00001   0.00009   0.00000   0.99988   0.00000   0.00000   0.00000
# K      0.00002   0.00001   0.00009   0.00000   0.99988   0.00000   0.00000   0.00000
# K      0.00002   0.00001   0.00009   0.00000   0.99945   0.00043   0.00000   0.00000
# into gnuplot format
# # BT2157
# #pos	aa	i	o	M	S
# 1	M	0.00002	0.0001	0.00000	0.99988
# 2	K	0.00002	0.0001	0.00000	0.99988
# 3	K	0.00002	0.0001	0.00000	0.99988

open(my $fhTmp, "<", $tmpFile) || die "Cannot read $tmpFile";
open(my $fhP, ">", "$out.plp") || die "Cannot write to $out.plp";

my $header1 = <$fhTmp>;
print $fhP $header1;
my $header2 = <$fhTmp>;
die "Invalid line 2 from decodeanhmm: $header2" unless $header2 =~ m/# +i +/;
print $fhP join("\t", qw{pos aa i o M S})."\n";
my $iLine = 0;
while(my $line = <$fhTmp>) {
  $iLine++;
  chomp $line;
  my @F = split / +/, $line;
  die "Invalid data line from decodeanhmm: $line" unless @F == 9;
  my $S = 
  print $fhP join("\t", $iLine, $F[0], $F[1], $F[2]+$F[3], $F[4], $F[5]+$F[6]+$F[7]+$F[8])."\n";
}
close($fhTmp) || die "Error reading $tmpFile";
unlink($tmpFile);
close($fhP) || die "Error writing $out.plp";

open(my $fhR, ">", "$out.R")
  || die "Cannot write to $out.gnuplot";
my $Rcmd = qq{
id = "$id";
out = "$out";
};
$Rcmd .= q{

data = read.delim(paste0(out,".plp"), skip=1);
colors = c("red", "magenta", "blue", "green");
labels = c("transmembrane", "signal peptide", "cytoplasmic", "non-cytoplasmic");
columns = c("M", "S", "i", "o");

svg(paste0(out,".svg"), width=5, height=5);
par(mar=c(5,3,2.5,1), mgp=c(2.5,0.8,0), bty="l", xaxs="i", yaxs="i");
plot(c(1,max(data$pos)), c(0,1), col=NA,
  xlab="",
  ylab="Posterior probability", main=paste("Phobius analysis of\n", id));
for(i in 1:length(labels))
  lines(data$pos, data[[ columns[i] ]], col=colors[i], lwd=2, xpd=T);
legend(0, -4 * strheight("A"), labels, lty=1, col=colors, lwd=2, xpd=T, ncol=2, bty="n");

# arrows for sequential regions
class = apply(data[,columns], 1, which.max);
rle = rle(class);
d = data.frame(constEnd = cumsum(rle$lengths),
               class = rle$values);
d$constBegin = d$constEnd - rle$lengths + 1;
y = -3.5 * strheight("A");
arrows(d$constBegin, y, d$constEnd, y,
  length=0.05, code=3, col=colors[d$class], xpd=T, lwd=2);

write.table(data.frame(begin=d$constBegin, end=d$constEnd, state=labels[d$class]),
            paste0(out,".tsv"),
            sep="\t", row.names=F, quote=F);

invisible( dev.off() );
};

print $fhR $Rcmd;
close($fhR) || die "Error writing to $out.R";

system("Rscript", "$out.R") == 0
  || die "R failed on $out.R -- $!";

