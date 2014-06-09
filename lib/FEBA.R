# FEBA.R -- analysis scripts for barcode sequencing data
#
# The key routines are:
#
# AvgStrainFitness -- compute fitness values from counts for the post-experiment counts
#     and the Time0 counts
# NormalizeByScaffold -- normalize the fitness values by scaffold and position
# GeneFitness -- combines the above two, and also computes a t-like test statistic ("t").
#	To do this, it also computes fitness values for the 1st and 2nd half of most genes
# FEBA_Fit() -- analyze many fitness experiments with AvgStrainFitness() and NormalizeByScaffold()
#      returns a complex data structure
# FEBA_Save_Tables() -- Save the fitness data structure to tab-delimited files and an R image


# Limitations:
# High memory usage (~10GB to process 200K strains x 500 experiments)
# Genes that wrap around the origin (i.e., begin > end) are ignored (no strain will map within them)


# GeneFitness():
# genes -- must include locusId, scaffoldId, and begin
# strainInfo -- must include locusId and (unless use1 is overridden) f, the fraction of the gene
# 	     that the insertion is at
# countCond and countT0 -- counts for each strain
# strainsUsed & genesUsed -- see AvgStrainFitness()
# genesUsed12 -- ditto, for 1st and 2nd half fitness values
# use1 -- which strains are in 1st half (regardless of whether they are usable or not)
# other arguments are passed on to AvgStrainFitness()
# base_se -- likely amount of error in excess of that given by variation within fitness values
# 	for strains in a gene, due to erorrs in normalization or bias in the estimator
#
# Returns a data frame with a row for each gene in genesUsed. It includes
# locusId,
# fit (unnormalized), fitnorm (normalized),
# fit1 or fit2 for 1st- or 2nd-half (unnormalized, may be NA),
# fitnorm1 or fitnorm2 for normalized versions,
# se (estimated standard error of measurement), and t (the test statistic),
# as well as some other values from AvgStrainFitness(), notably sdNaive,
# which is a different (best-case) estimate of the standard error.
#
GeneFitness = function(genes, strainInfo, countCond, countT0,
	    	    strainsUsed, genesUsed, genesUsed12,
		    use1 = strainInfo$f < 0.5,
		    base_se = 0.1,
		    ...) {
    d = AvgStrainFitness(countCond, countT0, strainInfo$locusId,
      		         strainsUsed=strainsUsed, genesUsed=genesUsed, ...);
    d$fitnorm = NormalizeByScaffold(d$fit, d$locusId, genes);

    d1 = AvgStrainFitness(countCond, countT0, strainInfo$locusId,
      		         strainsUsed=strainsUsed & use1, genesUsed=genesUsed12, ...);
    d2 = AvgStrainFitness(countCond, countT0, strainInfo$locusId,
      		         strainsUsed=strainsUsed & !use1, genesUsed=genesUsed12, ...);
    if(any(as.character(d1$locusId) != as.character(d2$locusId))) stop("Non-matching locusId");

    i = match(d$locusId, d1$locusId);

    d$fit1 = d1$fit[i];
    d$fit2 = d2$fit[i];
    d$fitnorm1 = d$fit1 + (d$fitnorm-d$fit);
    d$fitnorm2 = d$fit2 + (d$fitnorm-d$fit);
    d$tot1 = d1$tot[i];
    d$tot1_0 = d1$tot0[i];
    d$tot2 = d2$tot[i];
    d$tot2_0 = d2$tot0[i];

    # for low n, the estimated variance is driven by the overall variance, which can be estimated
    # from the median difference between 1st and 2nd halves via the assumptions
    # Var(fit) = Var((fit1+fit2)/2) ~= Var(fit1-fit2)/4
    # median abs(normal variable) = qnorm(0.75) * sigma = 0.67 * sigma
    # which leads to Var(fit) = Var(fit1-fit2)/4
    # = sigma12**2/4 = median abs diff**2 / (qnorm(0.75)*2)**2
    # The median difference is used because a few genes may have genuine biological differences
    # between the fitness of the two halves.
    # Furthermore, assuming that genes with more reads are less noisy, this
    # pseudovariance should be rescaled based on sdNaive**2
    #
    pseudovar_std = median(abs(d$fit1-d$fit2),na.rm=T)**2 / (2*qnorm(0.75))**2;
    d$pseudovar = pseudovar_std * (d$sdNaive / median(d$sdNaive[!is.na(d$fit1)]))**2;
    # given the variable weighting in sumsq, it is not intuitive that the degrees of freedom is still n-1
    # however, this is the result given the assumption that the weighting is the inverse of the variance
    est_var = (d$pseudovar + d$sumsq)/d$n;
    d$se = sqrt(est_var);
    d$t = d$fitnorm/sqrt(base_se**2 + pmax(d$sdNaive**2, est_var));
    return(d);
}

# AvgStrainFitness():
# strainCounts -- counts at the end of the experiment condition
# strainT0 -- counts for Time0 for each strain
# strainLocus -- which locus the strain is associated with, or NA
#
# If genesUsed (as a list of locusId) and strainsUsed (as boolean vector) are provided,
# then considers only those strains & genes; minimum requirements.
#
# If provided,
# nACG should include columns nA, nG, and nC with the counts of each nucleotides in the rcbarcodes
#
# debug & even are for testing purposes.
AvgStrainFitness = function(strainCounts, strainT0, strainLocus,
		 minStrainT0 = 4, minGeneT0 = 40,
		 nACG = NULL,
		 genesUsed=NULL, strainsUsed=NULL,
		 even=FALSE,
		 # maxWeight of 100 corresponds to having 100 reads on each side (if perfectly balanced)
		 maxWeight = if(even) 0 else 100,
		 debug=FALSE) {
    if (length(strainCounts) < 1 || length(strainT0) < 1 || length(strainLocus) < 1
        || length(strainCounts) != length(strainT0) || length(strainCounts) != length(strainLocus))
        stop("No or misaligned input data");

    if (is.null(strainsUsed)) strainsUsed = strainT0 >= minStrainT0;
    if (is.null(genesUsed)) {
        geneT0 = aggregate(strainT0[strainsUsed], list(locusId=strainLocus[strainsUsed]), sum);
        genesUsed = geneT0$locusId[geneT0$x >= minGeneT0];
    }
    strainsUsed = strainsUsed & strainLocus %in% genesUsed;
    if (!any(strainsUsed)) stop("No usable strains");

    strainT0 = strainT0[strainsUsed];
    strainCounts = strainCounts[strainsUsed];
    strainLocus = strainLocus[strainsUsed];

    readratio = sum(strainCounts) / sum(strainT0);
    # use sqrt(readratio), or its inverse, instead of 1, so that the expectation
    # is about the same regardless of how well sampled the strain or gene is
    strainFit = mednorm(log2(sqrt(readratio) + strainCounts) - log2(1/sqrt(readratio) + strainT0));
    strainFitAdjust = 0;

    if (!is.null(nACG)) {
        if (nrow(nACG) != length(strainsUsed)) stop("Wrong number of rows in nACG");
    	nACG = nACG[strainsUsed,];
	model = lm(strainFit ~ nACG$nA + nACG$nC + nACG$nG);
	if(debug) print(summary(model));
	strainFitAdjust = mednorm(predict(model));
	if(debug) cat("Median abs adjust: ", median(abs(strainFitAdjust)), "\n");
    }


    # Per-strain "smart" pseudocount to give a less biased per-strain fitness estimate.
    # This is the expected reads ratio, given data for the gene as a whole
    # Arguably, this should be weighted by T0 reads, but right now it isn't.
    # Also, do not do if we have just 1 or 2 strains, as it would just amplify noise

    # note use of as.vector() to remove names -- necessary for speed
    strainLocusF = as.factor(strainLocus);
    nStrains = table(strainLocusF);
    if(!all(names(nStrains)==levels(strainLocusF))) stop("Strain mismatch");
    nStrains = as.vector(nStrains);
    geneFit1 = mednorm(as.vector(tapply(strainFit, strainLocusF, mean)));
    i = as.integer(strainLocusF); # from strain index to gene index
    strainPseudoCount = ifelse(nStrains[i] >= 3, 2**geneFit1[i] * readratio, readratio);

    # And apportion the pseudocount equally (in log space) between condition-count and strain-count
    # to minimize the deviations from pseudocount = 1
    condPseudoCount = sqrt(strainPseudoCount);
    t0PseudoCount = 1/sqrt(strainPseudoCount);
    # (or could do some sort of weighted likelihood-based inference of fitness values, might be better)

    # for each strain: fitness, variance, and weight
    strainFit = log2(condPseudoCount + strainCounts) - log2(t0PseudoCount + strainT0) - strainFitAdjust;
    strainVar = sqrt(1/(1+strainT0) + 1/(1+strainCounts)) / log(2);
    # even weighting is an option for testing purposes
    # use harmonic mean for weighting
    strainWeight = 0.5 + pmin(maxWeight, 2/( 1/(1+strainT0) + 1/(1+strainCounts) ) );

    fitness = lapply(split(1:length(strainT0), list(locusId=strainLocus)),
     	           function(j) {
                       totw = sum(strainWeight[j]);
		       fitRaw = sum(strainWeight[j] * strainFit[j]) / totw;
		       tot = sum(strainCounts[j]);
		       tot0 = sum(strainT0[j]);
		       sd = sqrt(sum(strainWeight[j]**2 * strainVar[j]))/totw;
		       sumsq = sum(strainWeight[j] * (strainFit[j]-fitRaw)**2)/totw;
		       # high-N estimate of the noise in the log2 ratio of fitNaive
		       # But sdNaive is actually pretty accurate for small n -- e.g.
		       # simulations with E=10 on each side gave slightly light tails
		       # (r.m.s.(z) = 0.94).
		       sdNaive = sqrt( 1/(1+tot) + 1/(1+tot0) ) / log(2);
		       n = length(j);
		       nEff = totw/max(strainWeight[j]);
		       c(fitRaw=fitRaw, sd=sd, sumsq=sumsq, sdNaive=sdNaive, n=n, nEff=nEff,
		         tot=tot, tot0=tot0);
		});
    fitness = data.frame(do.call(rbind, fitness));
    fitness$fit = mednorm(fitness$fit);
    fitness$fitNaive = mednorm(log2(1+fitness$tot) - log2(1+fitness$tot0));
    fitness$locusId = row.names(fitness);
    if (is.integer(strainLocus)) fitness$locusId = as.integer(as.character(fitness$locusId));
    return(fitness);
}

# NormalizeByScaffold():
# values -- fitness values (as a vector)
# locusId -- the corresponding locusIds
# genes contains locusId, scaffoldId, and begin
# minForLoess -- how many genes a scaffold needs to contain before we try to use Loess to eliminate bias
#     For those scaffolds, it also estimates the mode and removes that
# minToUse -- if a scaffold has too few genes, cannot correct for possible DNA extraction bias
# 	   so need to remove data for that gene (i.e., returns NA for them).
# returns
# normalized -- data with scaffold and position effects removed
#
NormalizeByScaffold = function(values, locusId, genes, minForLoess=100, minToUse=10, debug=FALSE) {
    i = match(locusId, genes$locusId);
    if(any(is.na(i))) stop("Fitness data for loci not in genes");
    beg = genes$begin[i];

    perScaffoldRows = split(1:length(values), genes$scaffoldId[i]);
    for (scaffoldId in names(perScaffoldRows)) {
        rows = perScaffoldRows[[scaffoldId]];
        if (length(rows) < minToUse) {
            if(debug) cat("Removing ",length(rows)," values for ", scaffoldId, "\n");
	    values[rows] = NA;
        } else {
  	    med = median(values[rows]);
	    if(debug) cat("Subtract median for ", scaffoldId, " ", med, "\n");
	    values[rows] = values[rows] - med;

	    if (length(rows) >= minForLoess) {
                d = lowess(beg[rows], values[rows]);
	        if(debug) cat("Subtract loess for ", scaffoldId, " max effect is ", diff(range(d$y)), "\n");
	        values[rows] = values[rows] - approx(d$x,d$y, xout=beg[rows], rule=2, ties="ordered")$y;
	        d = density(values[rows]);
	        mode = d$x[which.max(d$y)];
	        if (debug) cat("Subtract mode for ", scaffoldId, " which is at ", mode, "\n");
	        values[rows] = values[rows] - mode;
            }
        }
    }
    return(values);
}

# simple log-ratio with pseudocount (of 1) and normalized so each scaffold has a median of 0
# note is *not* normalized except to set the total median to 0
StrainFitness = function(count, countT0, scaffolds) {
    fit = mednorm( log2(1+count) - log2(1+countT0) );
    se = sqrt(1/(1+count) + 1/(1+countT0)) / log(2);
    return(data.frame(fit=fit,se=se));
}

# For each strain, find the closest gene, as a row number -- returns a vector
# If there is no gene on that scaffold, returns NA
StrainClosestGenes = function(strains, genes) {
	genes$index = 1:nrow(genes);
	strainSplit = split(strains, strains$scaffold);
	geneSplit = split(genes, genes$scaffold);
	indexSplit = list();
	for (sc in names(strainSplit)) {
		s = strainSplit[[sc]];
		g = geneSplit[[sc]];
		if (is.null(g)) {
		    indexSplit[[sc]] = rep(NA, nrow(s));
		} else {
		    g$pos = (g$begin + g$end) / 2;
		    g = g[order(g$pos),];
		    # rule 2 means use values from extrema
		    i = round(approx(g$pos, 1:nrow(g), xout = s$pos, rule=2)$y);
		    i = pmax(1, pmin(nrow(g), i));
		    indexSplit[[sc]] = g$index[i];
		}
	}
	unsplit(indexSplit, strains$scaffold);
}

FEBA_Fit = function(expsUsed, all, all_g2, genes,
	   		       genesUsed=NULL, strainsUsed=NULL, genesUsed12=NULL,
	   		       minT0Strain=3, minT0Gene=30,
			       minT0GeneSide=minT0Gene/2,
			       minGenesPerScaffold=10,
			       nACG=NULL, # normalize by ACG content if set
			       pred=CrudeOp(genes),
			       okLane=TRUE, # OK to use t0 from another lane if needed?
	   		       metacol=1:7, # for all_g2
			       metacolAll=1:5,
			       # names of experiments to ignore
			       ignore=NULL,
			       debug=FALSE, computeCofit=TRUE,
			       ...) {

	if(!is.null(ignore)) {
	    cat("Ignoring ",ignore,"\n");
	    expsUsed = expsUsed[!expsUsed$name %in% ignore,];
	    all = all[, !names(all) %in% ignore,];
	    all_g2 = all_g2[, !names(all_g2) %in% ignore,];
	}

	if(!all(expsUsed$name == names(all_g2)[-metacol]))
		stop("names do not match");
	if(!all(expsUsed$name %in% names(all)))
		stop("names missing from all");
	if(is.null(genes$scaffoldId)) stop("No scaffold for genes");
	if(is.null(genes$begin)) stop("No begin for genes");
	if (is.null(genes$GC)) stop("Warning: no GC field in genes");

	expsUsed$name = as.character(expsUsed$name);

	# if Group = Time0, it is a Time0, even if "short" has a different description
	if(!is.null(expsUsed$Group)) expsUsed$short[expsUsed$Group == "Time0"] = "Time0";

    write("Aggregating all_g2",stderr())
    #write("all_g2",stdout())
    #write(paste(head(all_g2)),stdout())
	all_gN = aggregate(all_g2[,-(1:7)], list(locusId=all_g2$locusId), sum);

	expsUsed$t0set = paste(expsUsed$Date_pool_expt_started, expsUsed$Set);
	d = expsUsed[expsUsed$short=="Time0",];
	expsT0 = split(d$name, d$t0set);
	no_t0 = setdiff(expsUsed$t0set, expsUsed$t0set[expsUsed$short=="Time0"]);
	# If there is no matched lane/date t0, use the same date from other lane(s)
	# if okLane is TRUE
	if(length(no_t0) > 0) {
	    if(okLane) {
	        newt0 = unique(expsUsed$Date_pool_expt_started[expsUsed$t0set %in% no_t0]);
	        cat("Using t0 from other lanes instead for ", no_t0,"\n");
	        cat("Experiments affected: ", expsUsed$name[expsUsed$short!="Time0" & expsUsed$t0set %in% no_t0],"\n");
	        expsUsed$t0set[expsUsed$t0set %in% no_t0] = expsUsed$Date_pool_expt_started[expsUsed$t0set %in% no_t0];
		for(n in newt0) {
		    expsT0[[n]] = expsUsed$name[expsUsed$Date_pool_expt_started==n & expsUsed$short=="Time0"];
		}
	    } else {
	        stop("No Time0 for ",no_t0);
	    }
	}

	# note that t0_g2 does not contain any metadata columns, but t0_gN does
	t0_g2 = data.frame(lapply(expsT0, function(x) rowSums(all_g2[,x,drop=F])), check.names=F);
    write("Aggregating t0_g2",stderr())
	t0_gN = aggregate(t0_g2, list(locusId=all_g2$locusId), sum);
	cat("Reads per t0set, in millions:\n");
	print(colSums(t0_g2)/1e6, digits=2);

	if(is.null(strainsUsed)) strainsUsed = rowMeans(t0_g2) >= minT0Strain;
	if(is.null(genesUsed)) {
		t0_gN_used = aggregate(t0_g2[strainsUsed,], list(locusId=all_g2$locusId[strainsUsed]), sum);
		genesUsed = t0_gN_used$locusId[ rowMeans(t0_gN_used[,-1,drop=F]) >= minT0Gene ];
	}
	genesPerScaffold = table(genes$scaffoldId[genes$locusId %in% genesUsed]);
	smallScaffold = names(genesPerScaffold)[genesPerScaffold < minGenesPerScaffold];
	if (length(smallScaffold) > 0) cat("Ignoring genes on small scaffolds ",smallScaffold,"\n");
	genesUsed = genesUsed[!genesUsed %in% genes$locusId[genes$scaffoldId %in% smallScaffold]];

	if(length(strainsUsed) != nrow(all_g2)) stop("Invalid strainsUsed");
	cat("Using ",sum(strainsUsed)," of ",nrow(all_g2)," genic strains\n");
	if(length(genesUsed) < 100 || !all(genesUsed %in% genes$locusId)) stop("Invalid genesUsed");

	cat("Using ",length(genesUsed)," of ",length(unique(all_g2$locusId))," genes with data\n");

	if (is.null(genesUsed12)) {
  	    d1 = aggregate(t0_g2[strainsUsed & all_g2$f < 0.5,],
	     		   list(locusId=all_g2$locusId[strainsUsed & all_g2$f < 0.5]), sum);
	    d2 = aggregate(t0_g2[strainsUsed & all_g2$f >= 0.5,],
	     		   list(locusId=all_g2$locusId[strainsUsed & all_g2$f >= 0.5]), sum);
	    genesUsed12 = intersect(d1$locusId[ MyRowMin(d1[,-1,drop=F]) >= minT0GeneSide],
		      		    d2$locusId[ MyRowMin(d2[,-1,drop=F]) >= minT0GeneSide]);
	}
        cat("For cor12, using ",length(genesUsed12),"genes\n");

	if(!all(expsUsed$t0set %in% names(t0_g2))) stop("Illegal t0set", setdiff(expsUsed$t0set, names(t0_g2)));

	all_fit = list();
	strain_fit = list();
	strain_se = list();
	for(n in names(all_g2)[-metacol]) {
		to_subtract = expsUsed$short[expsUsed$name==n] == "Time0";
		if(debug) cat("GeneFitness() on", n,"t0set",expsUsed$t0set[expsUsed$name==n],
			  			  if(to_subtract) "subtracted" else "", "\n");
                x = all_g2[[n]];
		t0set = as.character(expsUsed$t0set[expsUsed$name==n]);
		t0 = t0_g2[[ t0set ]];
		if(to_subtract) t0 = t0 - x;
		if(any(t0 < 0)) stop("Illegal counts under 0 for ",n);
		if(all(t0 == 0)) {
		    cat("Skipping ",n," which has no control counts\n"); # used to be debug only
		} else {
		    all_fit[[n]] = GeneFitness(genes, all_g2[,words("locusId f")], x, t0,
				    strainsUsed, genesUsed, genesUsed12, nACG=nACG);
                    cntrl = setdiff(expsT0[[ t0set ]], n);
		    if(length(cntrl) < 1) stop("No Time0 experiments for ",n," should not be reachable");
		    if(debug) cat("StrainFitness() on", n, "versus", cntrl,"\n");
		    d = StrainFitness(all[[n]], rowSums(all[,cntrl,drop=F]));
		    strain_fit[[n]] = d$fit;
		    strain_se[[n]] = d$se;
		}
	}
	strain_fit = data.frame(strain_fit);
	strain_se = data.frame(strain_se);

	if(length(all_fit) == 0) stop("All comparisons failed\n");
	if(debug) cat("GeneFitness() succeeded\n");
	fit = list(g = all_fit[[1]]$locusId);
	for(n in setdiff(names(all_fit[[1]]), "locusId"))
	    fit[[n]] = data.frame(lapply(all_fit, function(x) x[[n]]));
	names(fit) = sub("fitnorm","lrn",names(fit));
	names(fit) = sub("fit","lr",names(fit));
	if (debug) cat("Extracted fitness values\n");

	q = expsUsed[expsUsed$name %in% names(all_fit), words("name short t0set")];
	nUse = as.character(q$name);
	if(!all(nUse == names(fit$lrn))) stop("Mismatched names in fit");
	q$nMapped = colSums(all[,nUse,drop=F]);
	q$nPastEnd = colSums(all[all$scaffold=="pastEnd",nUse,drop=F]);
	q$nGenic = colSums(all_g2[,nUse,drop=F]);
	q$nUsed = colSums(fit$tot);
	q$gMed = apply(fit$tot,2,median);
	q$gMedt0 = apply(fit$tot0,2,median);
	q$gMean = apply(fit$tot, 2, mean);

	adj = AdjacentPairs(genes);
	adjDiff = adj[adj$strand1 != adj$strand2,];

	# consistency of 1st and 2nd half; m.a.d. is median absolute difference
	q$cor12 = mapply(function(x,y) cor(x,y,method="s",use="p"), fit$lrn1, fit$lrn2);
	q$mad12 = apply(abs(fit$lrn1-fit$lrn2), 2, median, na.rm=T);

	# consistency of log2 counts for 1st and 2nd half, for sample and for time0
	q$mad12c = apply(abs(log2(1+fit$tot1) - log2(1+fit$tot2)), 2, median, na.rm=T);
	q$mad12c_t0 = apply(abs(log2(1+fit$tot1_0) - log2(1+fit$tot2_0)), 2, median, na.rm=T);

	# correlation of operon or adjacent different-strand pairs
	q$opcor = apply(fit$lrn, 2, function(x) paircor(pred[pred$bOp,], fit$g, x, method="s"));
	q$adjcor = sapply(as.character(q$name), function(x) paircor(adjDiff, fit$g, fit$lrn[[x]], method="s"));
	# GC correlation -- a sign of PCR issues (and often associated with high adjcor)
	# c() to make it be a vector instead of a matrix
	q$gccor = c( cor(fit$lrn, genes$GC[ match(fit$g, genes$locusId) ], use="p") );

	# experiments with very high maximum fitness may not give meaningful results for the typical gene
	q$maxFit = apply(fit$lrn,2,max,na.rm=T);

	fit$q = q;
	fit$genesUsed = genesUsed;
	fit$strainsUsed = strainsUsed;
	fit$genesUsed12 = genesUsed12;
	fit$t0_g2 = t0_g2;
	# these gene totals are based on all strains, not on used strains, and will not match tot or tot0
	fit$gN = all_gN;
	fit$t0_gN = t0_gN;

	# These include all strains, not just those in genes
	i = match(as.character(all$barcode), as.character(all_g2$barcode));
	fit$strains = cbind(all[,metacolAll], cbind(all_g2[,words("locusId f")], used=fit$strainsUsed)[i,]);

	fit$strain_lr = strain_fit;
	fit$strain_se = strain_se;

	# Normalized per-strain values
	strainToGene = StrainClosestGenes(fit$strains, genes[match(fit$g, genes$locusId),]);
	fit$strain_lrn = mapply(function(sfit, gdiff) {
		# Add the relevant gene normalization; or, if NA, normalize the scaffold to a median of 0
		sdiff = gdiff[strainToGene];
		sdiffSc = -ave(sfit, fit$strains$scaffold, FUN=median);
		sdiff = ifelse(is.na(sdiff), sdiffSc, sdiff);
		return(sfit + sdiff);
	}, fit$strain_lr, fit$lrn-fit$lr);
	fit$strain_lrn = data.frame(fit$strain_lrn);

	# Statistics of cofitness on pairs
        status = FEBA_Exp_Status(q, ...);
	u = (status == "OK");
	u[is.na(u)] = FALSE;
	if (computeCofit && sum(u) >= 5) {
		cat("Computing cofitness with ", sum(u), " experiments\n");
		adjDiff$rfit = cor12(adjDiff, fit$g, fit$lrn[,u]);
		pred$rfit = cor12(pred, fit$g, fit$lrn[,u]);
		fit$pairs = list(adjDiff=adjDiff, pred=pred);
		random = data.frame(Gene1 = sample(fit$g, length(fit$g)*2, replace=T),
		       	            Gene2 = sample(fit$g, length(fit$g)*2, replace=T));
		random = random[as.character(random$Gene1) != as.character(random$Gene2),];
		random$rfit = cor12(random, fit$g, fit$lrn[,u]);
		fit$pairs = list(adjDiff=adjDiff, pred=pred, random=random);
	} else {
		cat("Only", sum(u),"experiments passed quality filters!\n");
	}
	return(fit);
}

FEBA_Save_Tables = function(fit, genes, org="?",
		 topdir="public_html/FEBA",
		 dir = paste(topdir,org,sep="/"),
		 writeImage=TRUE,
		 template_file="src/feba/lib/FEBA_template.html",
		 expsU=expsUsed,
		 ... # for FEBA_Quality_Plot
		 ) {
	#if(!file.exists(dir)) stop("No such directory ",dir);
	if(!file.exists(dir)) dir.create(dir);

	for (n in words("q lr lrn lrn1 lrn2 t")) {
	    if (is.null(fit[[n]]) || !is.data.frame(fit[[n]])) {
	        stop("Invalid or missing ",n," entry");
	    }
	}
	if (is.null(fit$genesUsed)) stop("Missing genesUsed");
	if (is.null(fit$g)) stop("Missing g -- versioning issue?");

	if(!all(names(fit$lr) == fit$q$name)) stop("Name mismatch");
	if(!all(names(fit$lrn) == fit$q$name)) stop("Name mismatch");

	nameToPath = function(filename) paste(dir,filename,sep="/");
	wroteName = function(x) cat("Wrote ",nameToPath(x),"\n");

	writeDelim(fit$q, nameToPath("fit_quality.tab"));
	wroteName("fit_quality.tab");

	writeDelim(cbind(genes, used=genes$locusId %in% fit$genesUsed), nameToPath("fit_genes.tab"));
	wroteName("fit_genes.tab");

	d = merge(genes[,c("locusId","sysName","desc")], cbind(locusId=fit$g,fit$lr));
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_logratios_unnormalized.tab"));
	wroteName("fit_logratios_unnormalized.tab");

	d = merge(genes[,c("locusId","sysName","desc")], cbind(locusId=fit$g,fit$lrNaive));
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_logratios_unnormalized_naive.tab"));
	wroteName("fit_logratios_unnormalized_naive.tab");

	d = merge(genes[,c("locusId","sysName","desc")], cbind(locusId=fit$g,fit$lrn));
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_logratios.tab"));
	wroteName("fit_logratios.tab");

        status = FEBA_Exp_Status(fit$q);
	u = (status == "OK");
	if (sum(u) == 0) {
		cat("Warning: FEBA_Exp_Status returned 0 OK entries\n");
	} else {
		d = genes[,c("locusId","sysName","desc")];
		d$comb = paste(d$sysName, d$desc); # for MeV
		d = merge(d, cbind(locusId=fit$g,fit$lrn[,u]));
		names(d)[-(1:4)] = paste(fit$q$name,fit$q$short)[u];
		writeDelim(d, nameToPath("fit_logratios_good.tab"));
		cat("Wrote fitness for ",sum(u), " successful experiments to ", nameToPath("fit_logratios_good.tab"),"\n");
	}

	d = merge(genes[,c("locusId","sysName","desc")], cbind(locusId=fit$g,fit$t));
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_t.tab"));
	wroteName("fit_t.tab");

	d = merge(genes[,c("locusId","sysName","desc")], cbind(locusId=fit$g,fit$se));
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_standard_error_obs.tab"));
	wroteName("fit_standard_error_obs.tab");

	d = merge(genes[,c("locusId","sysName","desc")], cbind(locusId=fit$g,fit$sdNaive));
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_standard_error_naive.tab"));
	wroteName("fit_standard_error_naive.tab");

	writeDelim(cbind(fit$strains,fit$strain_lrn), nameToPath("strain_fit.tab"));
	wroteName("strain_fit.tab");

	FEBA_Quality_Plot(fit$q, nameToPath("fit_quality.pdf"), org, ...);
	wroteName("fit_quality.pdf");

	if(is.null(fit$pairs)) {
		paste("No data for cofitness plot\n");
		unlink(nameToPath("cofitness.pdf"));
	} else {
		FEBA_Cofitness_Plot(fit$pairs, nameToPath("cofitness.pdf"), org);
		wroteName("cofitness.pdf");
	}

	pdf(nameToPath("fit_quality_cor12.pdf"),
		pointsize=10, width=5, height=5,
		title=paste(org,"Fitness Cor12 Plots"));
	for (i in 1:nrow(fit$q)) {
	    n = as.character(fit$q$name[i]);
	    changers = fit$g[abs(fit$t[[n]]) >= 3];
	    plot(fit$lrn1[[n]], fit$lrn2[[n]],
	    		  main=sprintf("Cor12 for %s %s (%.0f %.3f)\n%s",
			  	org, n, fit$q$gMed[i], fit$q$cor12[i], fit$q$short[i]),
	    		  xlab="First Half", ylab="Second Half",
			  col=ifelse(fit$g %in% changers, 2, 1));
	    eqline(); hline(0); vline(0);
	}
	dev.off();
	wroteName("fit_quality_cor12.pdf");

	labelAll = sprintf("%.0f %.2f %s %25.25s",
		      fit$q$gMed, fit$q$cor12, sub("set","",fit$q$name), fit$q$short);
        labelAll = ifelse(fit$q$short=="Time0", paste(labelAll, fit$q$t0set), labelAll);

	use = fit$q$short != "Time0";
	if(sum(use) > 2) {
	    lrClust = hclust(as.dist(1-cor(fit$lrn[,as.character(fit$q$name)[use]])));
	    pdf(nameToPath("fit_cluster_logratios.pdf"),
		pointsize=8, width=0.25*pmax(8,sum(use)), height=8,
		title=paste(org,"Cluster Logratios"));
	    plot(lrClust, labels=labelAll[use], main="");
	    dev.off();
	    wroteName("fit_cluster_logratios.pdf");
	}

	countClust = hclust(as.dist(1-cor(log2(1+fit$gN[fit$gN$locusId %in% fit$genesUsed,-1]))));
	pdf(nameToPath("fit_cluster_logcounts.pdf"),
		pointsize=8, width=0.25*nrow(fit$q), height=8,
		title=paste(org,"Cluster Log Counts"));
	# Some Time0s may be missing from fit$q
	d = match(names(fit$gN)[-1], fit$q$name);
	labelAll2 = ifelse(is.na(d), paste("Time0", sub("set","",names(fit$gN)[-1])), labelAll[d]);
	plot(countClust, labels=labelAll2, main="");
	dev.off();
	wroteName("fit_cluster_logcounts.pdf");

        d = table(genes$scaffoldId[genes$locusId %in% fit$genesUsed]);
	maxSc = names(d)[which.max(d)];
	if (is.null(maxSc)) stop("Invalid scaffoldId?");
	beg = ifelse(fit$g %in% genes$locusId[genes$scaffold==maxSc],
	    genes$begin[match(fit$g, genes$locusId)], NA);

	pdf(nameToPath("fit_chr_bias.pdf"), pointsize=10, width=5, height=5,
	          title=paste(org,"Chromosome Bias"));
	for (i in 1:nrow(fit$q)) {
	    n = as.character(fit$q$name[i]);
	    plot(beg, pmax(-2,pmin(2,fit$lr[[n]])),
	    		  main=sprintf("%s %s (%.0f %.3f)\n%s",
			  	org, n, fit$q$gMed[i], fit$q$cor12[i], fit$q$short[i]),
	    		  xlab="Position on Main Scaffold",
			  ylab="Fitness (Unnormalized)",
			  ylim=c(-2,2), col="darkgrey");
	    o = order(beg);
	    lines(beg[o], (fit$lr[[n]] - fit$lrn[[n]])[o], col="darkgreen", lwd=2);
	    hline(0,lty=1,col=1);
	}
	dev.off();
	wroteName("fit_chr_bias.pdf");

	if (!is.null(expsU)) {
		writeDelim(expsU, nameToPath("expsUsed"));
		wroteName("expsUsed");
	}

	if(writeImage) {
	    img = format(Sys.time(),"fit%Y%b%d.image"); # e.g., fit2013Oct24.image
	    expsUsed = expsU;
	    save(fit, genes, expsUsed, file=nameToPath(img));
	    wroteName(img);
	    unlink(nameToPath("fit.image"));
	    file.symlink(img, nameToPath("fit.image"));
	    cat("Created link for ",nameToPath("fit.image"),"\n");
	}

	if(!is.null(template_file)) {
	    FEBA_Save_HTML(nameToPath("index.html"), template_file, org, nrow(fit$q));
	    wroteName("index.html");
	}
}

FEBA_Quality_Plot = function(q, pdfFile, org,
		min_gMed = 50, max_mad12 = 0.5, min_cor12 = 0.2,
		max_gccor = 0.2, max_adjcor = 0.25) {
        status = FEBA_Exp_Status(q,
				 min_gMed=min_gMed, max_mad12=max_mad12, min_cor12=min_cor12,
				 max_gccor=max_gccor, max_adjcor=max_adjcor);
	print(table(status));
	for(s in c("low_count","high_mad12","low_cor12","high_adj_gc_cor")) {
	    if(sum(status==s) > 0) cat(s, ":", q$name[status==s],"\n");
	}
	qCol = ifelse(q$short=="Time0", "grey",
		ifelse(status=="OK", ifelse(q$maxFit > 5, "blue", "darkgreen"), "red"));
	qLab = sub("set","",q$name);

	if(!is.null(pdfFile)) pdf(pdfFile,
		pointsize=10, width=8, height=8,
		title=paste(org,"Fitness Quality"));
	oldpar = par(mfrow=c(2,2));

	# mad12 vs. gMed
	plotlab(0.1 + q$gMed, pmin(0.75,q$mad12), qLab, col=qCol, ylim=c(0,0.75), log="x",
		            cex=0.8,
			    xlab="Median Reads per Gene",
			    ylab="Median abs. diff.(1st half, 2nd half)",
			    main=paste(org,"mad12"));
	vline(min_gMed); hline(max_mad12);

	# cor12 vs. mad12
	plotlab(pmin(0.75,q$mad12), pmax(0,q$cor12), qLab, col=qCol, xlim=c(0,0.75), ylim=0:1,
			     cex=0.8,
			     xlab="Median abs. diff.(1st half, 2nd half)", ylab="rho(1st half, 2nd half)",
			     main=paste(org,"rho12"));
	vline(max_mad12); hline(min_cor12);

	# opcor vs. adjcor
	plotlab(abs(q$adjcor), pmax(0,q$opcor), qLab, col=qCol, xlim=0:1, ylim=0:1,
			     cex=0.8,
			     xlab="| rho(adj. genes on diff. strands) |",
			     ylab="rho(operon pairs)",
			     main=paste(org,"rho(operons)"));
	eqline();
	vline(max_adjcor);

	# adjcor vs. gccor
	plotlab(abs(q$gccor), abs(q$adjcor), qLab, col=qCol, xlim=0:1, ylim=0:1,
			     cex=0.8,
			     xlab="| cor(gene GC, fitness) |",
			     ylab="rho(adj. genes on diff. strands)",
			     main=paste(org,"GC effects"));
	eqline();
	vline(max_gccor); hline(max_adjcor);

	par(oldpar);
	if(!is.null(pdfFile)) dev.off();
}

FEBA_Cofitness_Plot = function(pairs, pdfFile, org) {
	if(!is.null(pdfFile)) pdf(pdfFile,
		pointsize=10, width=4, height=4,
		title=paste(org,"Cofitness"));
	CompareDensities(list(Operon=withoutNA(pairs$pred$rfit[pairs$pred$bOp]),
			          Adjacent=withoutNA(pairs$adjDiff$rfit),
				  Random=pairs$random$rfit),
			 legendX="topleft",
			 xlim=c(-1,1),
			 xlab="Cofitness", ylab="Density", lwd=c(1,2,2), col=c(3,2,1), lty=c(1,2,4),
			 main=paste("Cofitness in",org));
	if(!is.null(pdfFile)) dev.off();
}

FEBA_Save_HTML = function(outfile, template, org, nexps) {
	lines = readLines(template);
	lines = gsub("ORG", org, lines);
	lines = gsub("DATE", date(), lines);
	lines = gsub("NEXPS", nexps, lines);
	writeLines(lines, outfile);
}

# Utilities

# split a string separated by spaces into a list of word
words = function(s, by=" ", ...) { strsplit(s[1], by, ...)[[1]]; }

# Tab-delimited output
writeDelim = function(table,file,report=FALSE,...) {
	write.table(table,file,sep="\t",quote=FALSE,row.names=FALSE,...);
	if(report) cat("Wrote",nrow(table),"rows","to",file,"\n")
}

# Just like plot.default() but adds on the labels
# By default plotting symbols are off but you can override that
plotlab = function(x,y,labels, cex=1, col=1, pch="", ...) {
	plot.default(x,y,cex=cex,col=col,pch=pch,...);
	text(x,y,labels,cex=cex,col=col);
}

### Graphics utilities
hline <- function(y,col="grey",lty=2,lwd=1) {
	lines(c(-1e20,1e-40,1e20),c(y,y,y),col=col,lty=lty,lwd=lwd);
}

vline <- function(x,col="grey",lty=2,lwd=1) {
	lines(c(x,x,x),c(-1e20,1e-40,1e20),col=col,lty=lty,lwd=lwd);
}

eqline <- function(col="grey",lty=2,lwd=1) {
	x <- 10**(-25:25);
	lines(c(-rev(x),x),c(-rev(x),x),col=col,lty=lty,lwd=lwd);
}

# Crude operon predictions -- pairs of genes that are on the same strand and
# separated by less than the median amount are predicted to be in the same opron
# Input genes is a data frame with locusId, strand, begin, end, with genes in sorted order
# Returns a data frame with Gene1, Gene2, Sep for separation, and bOp (TRUE if predicted operon pair)
CrudeOp = function(genes) {
	d = merge(merge(data.frame(Gene1=genes$locusId[-nrow(genes)],Gene2=genes$locusId[-1]), genes, by.x="Gene1", by.y="locusId"), genes, by.x="Gene2", by.y="locusId",suffixes=1:2);
	d = d[d$strand1==d$strand2,]
	d$Sep = pmin(abs(d$begin1-d$end2),abs(d$end1-d$begin2));
	d$bOp = d$Sep < median(d$Sep);
	return(d);
}

# are values correlated for pairs? Intended for operon pairs but would work with other types too
paircor = function(pairs, genes, values, use="p", method="pearson", names=c("Gene1","Gene2")) {
	d = merge(merge(pairs[,names], data.frame(Gene1=genes, value1=values), by.x=names[1], by.y="Gene1"),
			data.frame(Gene2=genes, value2=values), by.x=names[2], by.y="Gene2");
	return(cor(d$value1, d$value2, use=use, method=method));
}

# median-based normalization
mednorm = function(x) x - median(x);

# replace missing values (NA) with 0s
na0 = function(x) ifelse(is.na(x),0,x);

MyRowMin = function(x) {
	 if (is.vector(x)) return(x);
	 apply(x, 1, min);
}

# For reformatting per-lane pool count tables into the "all" table
prefixName = function(x, prefix) { names(x) = paste(prefix,names(x),sep=""); return(x); }

# For shortening the experiment descriptions
applyRules = function(rules, desc) {
    for (i in 1:nrow(rules)) desc = sub(rules[i,1], rules[i, 2], desc);
    return(desc);
}

# to count number of A, C, and G nucleotides in each barcode
CountACG = function(rcbarcodes) {
	 nA = sapply(rcbarcodes, function(x) sum(strsplit(x,"")[[1]] == "A"));
	 nC = sapply(rcbarcodes, function(x) sum(strsplit(x,"")[[1]] == "C"));
	 nG = sapply(rcbarcodes, function(x) sum(strsplit(x,"")[[1]] == "G"));
	 return(data.frame(nA=nA,nC=nC,nG=nG));
}

# Utilities for mapping items within regions
# e.g. for identifying strains that lie within genes

# by.x is a column name in x and begin and end are column names in y
# only returns the first matching row in y (as sorted by begin) for each row in x
# but if you specify all="y" then it returns the first matching row in x for
# each row in y
# If you specify unique=T then tries to find cases where y is the only match for that row in x
#	but it still returns some duplicates (e.g. it is fooled if row yA is within row yB)
# If minmax=T, allows begin > end (takes pairwise minima and maxima) and writes "left" and "right"
findWithin = function(x, y, by.x, begin, end, all="x", unique=FALSE, minmax=F) {
	if (nrow(x) == 0 || nrow(y) == 0) return(NULL); # one input is empty

	if(is.null(y[[begin]])) stop("no field named",begin," in y argument");
	if (minmax) {
		y$left = pmin(y[[begin]], y[[end]]);
		y$right = pmax(y[[begin]], y[[end]]);
		begin = "left";
		end = "right";
	}
	if (all=="x") {
		y2 = y[order(y[[begin]]),];
		if (nrow(y) == 1) {
			floorI = rep(1,nrow(x));
		} else {
			floorI = floor(approx(y2[[begin]],1:nrow(y2), xout=x[[by.x]], rule=2, ties=min)$y);
		}
		outy = y2[floorI,];
		keepx = x[[by.x]] >= outy[[begin]] & x[[by.x]] <= outy[[end]];
		if (unique) {
			floorI2 = ceiling(approx(y2[[end]],1:nrow(y2), xout=x[[by.x]], rule=2, ties=max)$y);
			keepx = keepx & floorI==floorI2;
			# Is this where it gets fooled?
		}
		return(data.frame(cbind(x[keepx,], outy[keepx,])));
	} else if (all == "y") {
		if(unique) stop("unique not supported with all=y");
		x2 = x[order(x[[by.x]]),];
		floorI = floor(approx(x2[[by.x]],1:nrow(x2), xout=y[[end]], rule=2,
			ties='ordered')$y);
		outx = x2[floorI,];
		keepy = outx[[by.x]] >= y[[begin]] & outx[[by.x]] <= y[[end]];
		return(data.frame(cbind(outx[keepy,], y[keepy,])));
	} else {
		stop("Unknown value of option all in findWithin",all);
	}
}

# given that data is subgrouped by the same markers (such as a scaffold),
# do findWithin on each and merge the results
findWithinGrouped = function(splitx, splity, by.x, begin, end, debug=FALSE, ...) {
	out = NULL;
	for(i in names(splitx)) {
		if (!is.null(splity[[i]])) {
			if(debug) cat("Running group ",i,"\n");
			rows = findWithin(splitx[[i]], splity[[i]], by.x, begin, end, ...);
			if(debug) cat("Ran group ",i,"rows",nrow(rows),"\n");
			if(!is.null(rows) && nrow(rows) > 0) {
				out = if(is.null(out)) rows else rbind(out,rows);
			}
		}
	}
	return(out);
}

without <- function(list,columns=list$without) {
	list2 <- list;
	for (i in columns) { list2[[i]] <- NULL; }
	return(list2);
}

AdjacentPairs = function(genes) {
	adj = data.frame(Gene1 = genes$locusId, Gene2=c(genes$locusId[-1],genes$locusId[1]));
	adj = merge(merge(adj, genes, by.x="Gene1", by.y="locusId"), genes, by.x=c("Gene2","scaffoldId"), by.y=c("locusId","scaffoldId"), suffixes=1:2);
	return(adj);
}

# Given a list of Gene1 Gene2 pairs, and a matrix of data (as genes and data-only matrix),
# compute correlations for each pair or NA
cor12 = function(pairs, genes, data, use="p", method="pearson", names=c("Gene1","Gene2")) {
	i1 = match(pairs[[names[1]]], genes);
	i2 = match(pairs[[names[2]]], genes);
	return(sapply(1:nrow(pairs), function(x) if(is.na(i1[x]) | is.na(i2[x])) NA else
		cor(c(data[i1[x],], recursive=T), c(data[i2[x],], recursive=T), method=method, use=use)));
}

# Returns status of each experiment -- "OK" is a non-Time0 experiment that passes all quality metrics
# Note -- arguably min_cor12 should be based on linear correlation not Spearman.
# 0.1 threshold was chosen based on Marinobacter set5, in which defined media experiments with cor12 = 0.1-0.2
# clearly worked, and Kang Polymyxin B (set1), with cor12 ~= 0.13 and they barely worked.
FEBA_Exp_Status = function(q, min_gMed = 50, max_mad12 = 0.5, min_cor12 = 0.1,
				 max_gccor = 0.2, max_adjcor = 0.25) {
    with(q, ifelse(short=="Time0", "Time0",
	           ifelse(gMed < min_gMed, "low_count",
		   ifelse(mad12 > max_mad12, "high_mad12",
		   ifelse(cor12 < min_cor12, "low_cor12",
		   ifelse(abs(gccor) > max_gccor | abs(adjcor) > max_adjcor, "high_adj_gc_cor", "OK"))))));
}

withoutNA = function(x) x[!is.na(x)];

CompareDensities <- function(list,labels=names(list),xlim=range(unlist(list)),ylim=c(0,3),
		col=1:length(labels), lty=1:length(labels), lwd=rep(1,length(labels)),
		legendX=mean(xlim),legendY=ylim[2],main="",xlab="",ylab="",showCounts=FALSE,
		showLegend=TRUE, bty="o") {

	for (i in 1:length(labels)) {
		x = list[[ names(list)[i] ]];
		d <- density(x,from=xlim[1],to=xlim[2]);
		if(i==1) {
			plot(d$x,d$y,type="l",col=col[i],lty=lty[i],lwd=lwd[i],
				main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim);
		} else {
			lines(d$x,d$y,col=col[i],lty=lty[i],lwd=lwd[i]);
		}
		if (showCounts) labels[i] <- paste(labels[i]," (n=", sum(!is.na(x)), ")", sep="");
	}
	if(showLegend) {
	       if(is.numeric(legendX)) {
	           legend(legendX,legendY,labels,col=col,lty=lty,lwd=lwd, bty=bty);
	       } else { # e.g. "topleft"
	           legend(legendX, labels,col=col,lty=lty,lwd=lwd, bty=bty);
 	       }
	}
}

# Given pairs of adjcent genes, identify the upstream and downstream genes,
# as "up" and "dn"
OpPairUpDn = function(oppairs, genes) {
	oppairs$strand = genes$strand[match(oppairs$Gene1, genes$locusId)];
	if(any(is.na(oppairs$strand))) stop("Unknown genes in Gene1 in OpPairUpDn()");

	# mixing factors is problematic so make sure they are character or integer
	if(is.factor(oppairs$Gene1)) oppairs$Gene1 = as.character(oppairs$Gene1);
	if(is.factor(oppairs$Gene2)) oppairs$Gene2 = as.character(oppairs$Gene2);

	oppairs$up = ifelse(oppairs$strand=="+", oppairs$Gene1, oppairs$Gene2);
	oppairs$dn = ifelse(oppairs$strand=="-", oppairs$Gene1, oppairs$Gene2);
	return(oppairs);
}

# How often upstream-only is sick vs. downstream-only
# oppairs must include "up", "dn" (as in returned values from OpPairUpDn)
OperonPairSick = function(oppairs, g, lrn, t,
			sick = -1, min_diff = 0.75, max_t = -4) {
	oppairs = oppairs[oppairs$up %in% g & oppairs$dn %in% g,];
	if(is.vector(lrn)) {
		i1 = match(oppairs$up, g);
		i2 = match(oppairs$dn, g);
		uponly = sum(lrn[i1] < sick & lrn[i2] > sick & t[i1] < max_t & t[i2] > max_t
			& abs(lrn[i1]-lrn[i2]) > min_diff);
		dnonly = sum(lrn[i2] < sick & lrn[i1] > sick & t[i2] < max_t & t[i1] > max_t
			& abs(lrn[i1]-lrn[i2]) > min_diff);
		both = sum(lrn[i1] < sick & lrn[i2] < sick & t[i1] < sick & t[i2] < sick);
		return(c(uponly=uponly, dnonly=dnonly, both=both));
	}
	#else make a table, where the row names will be the experiment names
	as.data.frame(t(mapply(function(x,y)
	    OperonPairSick(oppairs, g, x, y, sick=sick, min_diff=min_diff, max_t=max_t), lrn, t)));
}
