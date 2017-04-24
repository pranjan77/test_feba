#!/usr/bin/perl -w
#######################################################
## help.cgi
##
## Copyright (c) 2015 University of California
##
## Authors:
## Morgan Price
#######################################################
#
# CGI parameters: none

use strict;

use CGI qw(:standard Vars);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);

use lib "../lib";
use Utils;

my $cgi=CGI->new;
my $style = Utils::get_style();
my $dbh = Utils::get_dbh();
my $start = Utils::start_page("Help for the Fitness Browser");

print header,
	$start, 
    h2("Help for the Fitness Browser");

my $helpcontent = <<END

<p>The Fitness Browser was developed by the <A HREF="http://genomics.lbl.gov">Arkin lab</A>.
It displays thousands of genome-wide fitness assays from the
<A HREF="mailto:AMDeutschbauer.lbl.gov">Deutschbauer</A> lab, the Arkin lab,
and collaborators.</p>

<H3><A NAME="technology" style="color: black;">How it works</A></H3>

<P>The fitness data is collected using randomly barcoded transposons (RB-TnSeq). Each fitness experiment is based on a pool of 30,000 to 500,000 mutant strains. Every mutant strain has a transposon inserted at a random location in the genome, and each transposon includes a random barcode that allows us to track the abundance of that strain by using PCR followed by DNA sequencing ("BarSeq"). To link the barcode to the location in the genome, we use a more complicated TnSeq-like protocol.</P>

<P>For each fitness experiment, we compare the abundance of each strain at the end of the experiment to its abundance at the beginning. The beginning sample is also referred to as the "Time0" sample. Typically, we recover the pool of mutants from the freezer in rich media, wash the cells and take Time0 sample(s), and transfer the washed cells into many different tubes or wells. Thus, many different conditions may be compared to the same Time0 sample(s).

<P>For details, see our <A HREF="http://mbio.asm.org/content/6/3/e00306-15.full">methods paper</A> (Wetmore et al, mBio 2015).

<H3><A NAME="fitness" style="color: black;">Gene fitness</A></H3>

<P>Fitness values are log<sub>2</sub> ratios that describe the change
in abundance of mutants in that gene during the experiment.  For most of
the fitness experiments, which are growth experiments, the change reflects
how well the mutants grow.
Fitness = 0 means that mutants in this
gene grew well as other mutants and probably about as well as wild
type strains.  Fitness &lt; 0 means that the gene was important for
fitness and the mutants were less abundant at the end of the
experiment than at the beginning.  For example, fitness = -1 means
that mutants in the gene were half as abundant at the end of the
experiment, compared to the beginning.  Fitness &gt; 0 means that the
gene was detrimental to fitness and that mutants had a growth
advantage.

<P>In general, if -1 &lt; fitness &lt; 1, then the gene has a subtle
phenotype that might be statistically significant (see <A HREF="#t">t
scores</A>) but will probably be difficult to interpret. Fitness &lt;
-2 or fitness &gt; 2 are strong fitness effects.
In the typical experiment, the pool of mutants doubles 4-8 times,
so in principle, a conditionally essential gene should have fitness of
-4 to -8. However, it is not possible to tell the difference between
little or no growth with a pooled assay. (Also, very low fitness
values are more noisy because they are based on a log<sub>2</sub>
ratio with a small numerator &ndash; in the typical experiment, a
fitness value of -1 is reliably different from 0, but -5 is not
reliably different from -4.)

<P>More rigorously, gene fitness is the weighted average of strain fitness,
across strains that has a transposon inserted within that gene. A
strain's fitness is the log<sub>2</sub> ratio of abundance at the end
of the experiment compared to its abundance at the beginning of the
experiment, where we use the number of reads for each strain's barcode
as a proxy for its abundance. The gene fitness is normalized so that
the typical gene has a fitness of zero. For genes on large
chromosomes, the gene fitness values are also normalized for changes
in copy number along the chromosome.

<P>Although most experiments are based on growth, this site also includes assays of motility or survival. For a motility assay, the experimental samples might be the cells that reached the outer ring of an agar plate, or that stayed in the inner ring where the cells were originally placed.
For a survival assay, the cells are stressed or starved for a period of time; then, to distinguish viable cells from dead cells,
all cells are transferred to a rich medium and recovered for a few generations.

<H3><A name="t" style="color: black;"><i>t</i> scores</A></H3>

The t-like test statistic indicates how reliably a gene fitness values is
different from zero. Ideally, they are on the same scale
as <i>z</i> scores or <i>t</i> scores. However, since there are thousands of genes
in each experiment, and there can be hundreds of fitness experiments
for a gene, a higher threshold is needed. We usually ignore any
fitness effects with <i>|t| &lt; 4</i>. In most cases, you can gain
confidence in a fitness effect by comparing the phenotype of a gene in
replicate experiments, or in similar experiments (such as different
concentrations of the same inhibitory compound), or for <A
HREF="#ortholog">orthologous</A> genes.

<H3><A name="cofitness" style="color: black;">Cofitness</A></H3>

<P>Cofitness(gene 1, gene 2) is the linear (Pearson) correlation of their fitness patterns.
Alternatively, if two genes in the same organism have similar fitness patterns, then we say that they are cofit.</P>

<P>If two genes have similar fitness patterns (cofitness &gt; 0.75), and they are among the most cofit genes (rank = 1 or rank = 2), then they are likely to function in the same pathway. For genes with strong fitness patterns, often the most cofit genes are other genes in the same operon, so we look a little farther down the list to find genes that may have related functions.</P>

<P>Conserved cofitness: If two genes have cofitness &gt; 0.6, and their orthologs have cofitness &gt 0.6, then this is stronger evidence of a functional relationship.

<P>If we have relatively little data for an organism, then cofitness results will not be available for any of its genes.

<H3><A name="specific" style="color: black;">Specific phenotypes</A></H3>

We define a gene as having a "specific" phenotype in a condition if the gene has a stronger phenotype in this condition than in most other conditions, and lacks phenotypes in most conditions. More precisely, we require

<UL>
<LI> |fit| > 1
<LI> |t| > 5 
<LI> |fit|<sub>95</sub> < 1, where |fit|<sub>95</sub> is the 95th percentile of |fit| across all experiments for this gene 
<LI> |fit| > |fit<sub>95</sub>| + 0.5

</UL>

If we have relatively little data for an organism, then there may not be any specific phenotypes for any of its experiments. Also, these criteria are stringent and may miss some genes.

<H3><A name="ortholog" style="color: black;">Orthologs</A></H3>

<P>We use "orthologs" to refer to similar proteins in different
organisms that may carry out the same function, without regard to
their evolutionary history. Thus they are putative functional
orthologs, not evolutionary orthologs.  The "orthologs" in this web
site are bidirectional best hits from protein BLAST. We also require
that the BLAST alignment cover 80% of each protein.</P>

<P>Many of these "orthologs" actually have different functions. If either gene has a strong fitness pattern, you may be able to use conserved phenotypes or conserved cofitness to confirm that the genes have conserved functions and are truly functional orthologs.</P>

<H3><A name="protein" style="color: black;">Protein sequence analysis</A></H3>

For each protein, the Fitness Browser includes:

<UL>
<LI><A HREF="http://pfam.xfam.org/">PFam</A> domains, computed with <A HREF="http://hmmer.org/">HMMer3</A>
<LI><A HREF="http://www.jcvi.org/cgi-bin/tigrfams/index.cgi">TIGRFam</A> domains or families, computed with <A HREF="http://hmmer.org/">HMMer3</A>
<LI>The best hit to <A HREF="http://www.genome.jp/kegg/">KEGG</A>, computed with <A HREF="https://github.com/zhaoyanswill/RAPSearch2">RAPSearch2</A> and minimum 80% coverage and 30% identity
<LI>The best hit to <A HREF="http://web.expasy.org/docs/swiss-prot_guideline.html">Swiss-Prot</A> (the curated part of UniProt), computed with <A HREF="https://github.com/zhaoyanswill/RAPSearch2">RAPSearch2</A> and minimum 80% coverage and 30% identity
<LI>The best hit to annotated enzymes in <A HREF="http://metacyc.org/">MetaCyc</A>, computed with <A HREF="https://github.com/zhaoyanswill/RAPSearch2">RAPSearch2</A> and minimum 80% coverage and 30% identity.
<LI>The <A HREF="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3965101/">SEED</A> annotation, computed with the <A HREF="http://www.theseed.org/servers/#mozTocId76305">SEED API</A>
</UL>

<P>Information from <A HREF="http://www.jcvi.org/cgi-bin/tigrfams/index.cgi">TIGRFam</A>, <A HREF="http://www.genome.jp/kegg/">KEGG</A>,
and <A HREF="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3965101/">SEED</A> is used to link proteins to
enzyme commision (EC) numbers and hence to metabolic maps (from the last public release of KEGG).

<P>Fitness Browser includes links to other analysis tools (see the protein page) as well as a homologs page (computed using BLAST).


<H3><A name="links" style="color: black;">Linking to the Fitness Browser</A></H3>

Most of the genomes in this web site were taken
from NCBI (i.e., gene identifiers are locus tags and scaffolds are Genbank accessions)
or from MicrobesOnline (i.e., gene identifiers and scaffold identifiers are numbers).
These identifiers should be stable over
time, so URLs from the web site should continue to work in the long run.
For example, to link to the fitness data for <i>endA</i> from
<i>E. coli</i>, you can use
<pre>
<A HREF="http://fit.genomics.lbl.gov/cgi-bin/singleFit.cgi?orgId=Keio&locusId=17024&showAll=0">http://fit.genomics.lbl.gov/cgi-bin/singleFit.cgi?orgId=Keio&locusId=17024&showAll=0</A>
</pre>

<P>Or, you can use <A HREF="../images/fitblast_example.html">Fitness BLAST</A> to link from any protein sequence to the homologs that have fitness data. You can incorporate this into your web page with just a few lines of <A HREF="../images/fitblast_example.html#code">code</A>.

<P>Or, you can use <A HREF="batch_blast.cgi">Fitness BLAST for genomes</A> to identify orthologs in our data set for an entire genome at once. It takes less than a minute and we plan to store the results indefinitely.</A>

<P>(Both Fitness BLAST and Fitness BLAST for genomes are powered by
<A HREF="http://www.drive5.com/usearch/">usearch</A>,
not BLAST. However, <A HREF="blastSearch.cgi">single sequence search</A> and the homologs page rely on BLAST.)

<H3><A name="code" style="color: black;">About the code</A></H3>

<p>The code for this web site is <A HREF="https://bitbucket.org/berkeleylab/feba/src">freely available at bitbucket.org</A>.
The code was written by <A HREF="http://morgannprice.org/">Morgan Price</A>, Victoria Lo, and Wenjun Shao
in the <A HREF="http://genomics.lbl.gov">Arkin lab</A>.</p>

<H3><A name="sources" style="color: black;">References</A></H3>

<UL>
<LI><A HREF="http://genomics.lbl.gov/supplemental/rbarseq/">Wetmore et al 2015</A> -- carbon source experiments for <i>Escherichia coli</i> BW25113, <i>Shewanella oneidensis</I> MR-1, <I>Shewanella amazonensis</I> SB2B, <i>Phaeobacter inhibens</i> BS107, and Pseudomonas stutzeri</i> RCH2
<LI><A HREF="http://www.pnas.org/content/early/2015/10/26/1519220112.abstract">Rubin et al 2015</A> -- the mutant library for <i>Synechococcus elongatus</i> PCC 7942
</ul>

<p>Most of the data is not published. Contact <A HREF="mailto:AMDeutschbauer.lbl.gov">Adam Deutschbauer</A> for more information about the unpublished data.

<H3><A name="ack" style="color: black;">Funding</A></H3>

<p>This site was developed by <A HREF="http://enigma.lbl.gov">
ENIGMA - Ecosystems and Networks Integrated with Genes
and Molecular Assemblies</A>, a Scientific Focus
Area Program at Lawrence Berkeley National Laboratory, and
supported by the U.S. Department of Energy, Office of Science,
Office of Biological &amp; Environmental Research under contract number
DE-AC02-05CH11231.

END
    ;

print div({-id=>"ntcontent"}, $helpcontent);
Utils::endHtml($cgi);
