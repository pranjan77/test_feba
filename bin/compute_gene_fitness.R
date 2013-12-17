#!/usr/bin/env Rscript

library('optparse')

option_list = list (
                     make_option(c("-p","--pool"), 

source("~/scripts/utils.R");




# Load metadata about experiments and genes
exps = read.delim("g/psRCH2/FEBA_BarSeq.tsv",as.is=T);
expsUsed = exps[grepl("_set[12345678]",exps$SetName) & exps$Description != "",];
expsUsed$name = paste(sub("psRCH2_ML7_","",expsUsed$SetName),expsUsed$Index,sep="");
rules = read.table("~/scripts/desc_short_rules",as.is=T);
expsUsed$short = applyRules(rules, expsUsed$Description);
genes = read.delim("g/psRCH2/genes.tab",as.is=T);

# Load the pool
pool10 = read.delim("g/psRCH2/pool.n10",as.is=T);
pool10g = findWithinGrouped(split(without(pool10[!is.na(pool10$pos),],"scaffold"),pool10$scaffold[!is.na(pool10$pos)]), split(genes, genes$scaffoldId), "pos", "begin", "end");
pool10g$f = (pool10g$pos - pool10g$begin) / (pool10g$end - pool10g$begin);
pool10g2 = pool10g[pool10g$f >= 0.1 & pool10g$f <= 0.9,];

# Load the BarSeq
psRCH2_ML7_set1_JGI = read.delim("g/psRCH2/psRCH2_ML7_set1.poolcount",as.is=T)
psRCH2_ML7_set1_seqUC = read.delim("g/psRCH2/psRCH2_ML7_set1_seqUC.poolcount",as.is=T);
psRCH2_ML7_set2 = read.delim("g/psRCH2/psRCH2_ML7_set2.poolcount",as.is=T)
psRCH2_ML7_set3 = read.delim("g/psRCH2/psRCH2_ML7_set3.poolcount",as.is=T)
psRCH2_ML7_set4 = read.delim("g/psRCH2/psRCH2_ML7_set4.poolcount",as.is=T)
psRCH2_ML7_set5 = read.delim("g/psRCH2/psRCH2_ML7_set5.poolcount",as.is=T)
psRCH2_ML7_set6 = read.delim("g/psRCH2/psRCH2_ML7_set6.poolcount",as.is=T)
psRCH2_ML7_set7 = read.delim("g/psRCH2/psRCH2_ML7_set7.poolcount",as.is=T);
psRCH2_ML7_set8 = read.delim("g/psRCH2/psRCH2_ML7_set8.poolcount",as.is=T);

# Combine the two runs of set1
psRCH2_ML7_set1 = data.frame(psRCH2_ML7_set1_JGI[,(1:5)], apply(psRCH2_ML7_set1_JGI[,-(1:5)],2,na0) + apply(psRCH2_ML7_set1_seqUC[,-(1:5)],2,na0));

# Combine all the lanes
all = data.frame(psRCH2_ML7_set1[,1:5],
	apply(prefixName(psRCH2_ML7_set1[,-(1:5)], "set1"), 2, na0),
	apply(prefixName(psRCH2_ML7_set2[,-(1:5)], "set2"), 2, na0),
	apply(prefixName(psRCH2_ML7_set3[,-(1:5)], "set3"), 2, na0),
	apply(prefixName(psRCH2_ML7_set4[,-(1:5)], "set4"), 2, na0),
	apply(prefixName(psRCH2_ML7_set5[,-(1:5)], "set5"), 2, na0),
	apply(prefixName(psRCH2_ML7_set6[,-(1:5)], "set6"), 2, na0),
	apply(prefixName(psRCH2_ML7_set7[,-(1:5)], "set7"), 2, na0),
	apply(prefixName(psRCH2_ML7_set8[,-(1:5)], "set8"), 2, na0));
# Keep only the columns that are supposed to be there
# (often some few indexes are not used, these were removed from expsUsed by ignoring
# rows with blank Description)
all = cbind(all[,1:5], all[,as.character(expsUsed$name)]);

# gene-relevant counts
all_g2 = merge(pool10g2[,words("barcode rcbarcode strand pos locusId f")], all);

# The actual computation and statistics (still undergoing changes)
# Needs to be run with okLane=TRUE in a few cases where there
# is no matching time0 sample.
fit = FEBA_Fit(expsUsed, all, all_g2, genes);

# FEBA_Save_Tables creates the files in the output directories
# You'd need to change the topdir option
FEBA_Save_Tables(fit, genes, org="psRCH2");


