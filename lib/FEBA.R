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

    d$fit1 = d1$fit[match(d$locusId, d1$locusId)];
    d$fit2 = d2$fit[match(d$locusId, d1$locusId)];
    d$fitnorm1 = d$fit1 + (d$fitnorm-d$fit);
    d$fitnorm2 = d$fit2 + (d$fitnorm-d$fit);

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
		 maxWeight = if(even) 0 else 200,
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
    geneFit1 = aggregate(strainFit, list(locusId=strainLocus), mean);
    geneFit1$x = mednorm(geneFit1$x);
    nStrains = aggregate(strainFit, list(locusId=strainLocus), length);
    i = match(strainLocus, geneFit1$locusId); # from strain index to gene index
    strainPseudoCount = ifelse(nStrains$x[i] >= 3, 2**geneFit1$x[i] * readratio, readratio);

    # And apportion the pseudocount equally (in log space) between condition-count and strain-count
    # to minimize the deviations from pseudocount = 1
    condPseudoCount = sqrt(strainPseudoCount);
    t0PseudoCount = 1/sqrt(strainPseudoCount);
    # (or could do some sort of weighted likelihood-based inference of fitness values, might be better)

    # for each strain: fitness, variance, and weight
    strainFit = log2(condPseudoCount + strainCounts) - log2(t0PseudoCount + strainT0) - strainFitAdjust;
    strainVar = sqrt(1/(1+strainT0) + 1/(1+strainCounts)) / log(2);
    # even weighting is an option for testing purposes
    strainWeight = 0.5 + pmin(maxWeight, strainT0+strainCounts);

    geneFit2 = aggregate(1:length(strainT0), list(locusId=strainLocus),
    	                function(j) {
			    totw = sum(strainWeight[j]);
			    meanFit = sum(strainWeight[j] * strainFit[j]) / totw;
			    c(fit = meanFit,
			    	  sd = sqrt(sum(strainWeight[j]**2 * strainVar[j]))/totw,
				  sumsq = sum(strainWeight[j] * (strainFit[j]-meanFit)**2)/totw,
				  sdNaive = sqrt( 1/(1+sum(strainCounts[j])) + 1/(1+sum(strainT0[j])) ) / log(2),
				  n = length(j),
				  nEff = sum(strainWeight[j])/max(strainWeight[j]),
				  tot = sum(strainCounts[j]),
				  tot0 = sum(strainT0[j]));
                        });
    # a hack to get around the lists within the entries
    geneFit2 = data.frame(locusId = geneFit2$locusId, geneFit2[,2]);

    geneFit2$fitRaw = geneFit2$fit; # without the median normalization
    geneFit2$fit = mednorm(geneFit2$fit);

    geneFit2$fitNaive = mednorm(log2(1+geneFit2$tot) - log2(1+geneFit2$tot0));
    return(geneFit2);
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

FEBA_Fit = function(expsUsed, all, all_g2, genes,
	   		       genesUsed=NULL, strainsUsed=NULL, genesUsed12=NULL,
	   		       minT0Strain=3, minT0Gene=30,
			       minT0GeneSide=minT0Gene/2,
			       minGenesPerScaffold=10,
			       nACG=NULL, # normalize by ACG content if set
			       pred=CrudeOp(genes),
			       okLane=TRUE, # OK to use t0 from another lane if needed?
	   		       metacol=1:7,
			       # names of experiments to ignore
			       ignore=NULL,
			       debug=FALSE) {

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

	expsUsed$name = as.character(expsUsed$name);

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

	t0_g2 = data.frame(lapply(expsT0, function(x) rowSums(all_g2[,x,drop=F])), check.names=F);
    write("Aggregating t0_g2",stderr())
	t0_gN = aggregate(t0_g2, list(locusId=all_g2$locusId), sum);
	cat("Reads per t0set, in millions:\n");
	print(colSums(t0_g2)/1e6, digits=2);

	if(is.null(strainsUsed)) strainsUsed = rowMeans(t0_g2[,-1,drop=F]) >= minT0Strain;
	if(is.null(genesUsed)) {
		d = aggregate(t0_g2[strainsUsed,], list(locusId=all_g2$locusId[strainsUsed]), sum);
		genesUsed = d$locusId[ rowMeans(d[,-1,drop=F]) >= minT0Gene ];
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

	all_fit = lapply(names(all_g2)[-metacol], function(n) {
		to_subtract = expsUsed$short[expsUsed$name==n] == "Time0";
		if(debug) cat("GeneFitness() on", n,"t0set",expsUsed$t0set[expsUsed$name==n],
			  			  if(to_subtract) "subtracted" else "", "\n");
                x = all_g2[[n]];
		t0 = t0_g2[[ expsUsed$t0set[expsUsed$name==n] ]];
		if(to_subtract) t0 = t0 - x;
		if(any(t0 < 0)) stop("Illegal counts under 0 for ",n);
		if(all(t0 == 0)) {
		    if(debug) cat("Skipping ",n," which has no control counts\n");
		    return(NULL);
		}
		GeneFitness(genes, all_g2[,words("locusId f")], x, t0,
				    strainsUsed, genesUsed, genesUsed12, nACG=nACG);
	});
	if(debug) cat("GeneFitness() succeeded\n");
	names(all_fit) = names(all_g2)[-metacol];
	keep = !sapply(all_fit, is.null); # drop comparisons which were "skipped"
	if(!any(keep)) stop("All comparisons failed\n");
	all_fit = all_fit[keep, drop=F];
	fit = list(g = all_fit[[1]]$locusId);
	for(n in names(all_fit[[1]])[-1]) {
	    fit[[n]] = data.frame(lapply(all_fit, function(x) x[[n]]));
	}
	names(fit) = sub("fitnorm","lrn",names(fit));
	names(fit) = sub("fit","lr",names(fit));
	if (debug) cat("Extracted fitness values\n");

	q = expsUsed[keep, words("name short t0set")];
	nUse = as.character(q$name);
	if(!all(nUse == names(fit$lrn))) stop("Mismatched names in fit");
	q$nMapped = colSums(all[,nUse]);
	q$nPastEnd = colSums(all[all$scaffold=="pastEnd",nUse,drop=F]);
	q$nGenic = colSums(all_g2[,nUse,drop=F]);
	q$nUsed = colSums(fit$tot);
	q$gMed = apply(fit$tot,2,median);
	q$gMedt0 = apply(fit$tot0,2,median);
	q$gMean = apply(fit$tot, 2, mean);

	q$cor12 = mapply(function(x,y) cor(x,y,method="s",use="p"), fit$lrn1, fit$lrn2);
	q$opcor = apply(fit$lrn, 2, function(x) paircor(pred[pred$bOp,], fit$g, x, method="s"));
	q$maxFit = apply(fit$lrn,2,max);

	fit$q = q;
	fit$genesUsed = genesUsed;
	fit$strainsUsed = strainsUsed;
	fit$genesUsed12 = genesUsed12;
	fit$t0_g2 = t0_g2;
	# these gene totals are based on all strains, not on used strains, and will not match tot or tot0
	fit$gN = all_gN;
	fit$t0_gN = t0_gN;

	return(fit);
}

FEBA_Save_Tables = function(fit, genes, org="?",
		 topdir="public_html/tmp",
		 dir = paste(topdir,org,sep="/"),
		 writeImage=TRUE) {
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
	wroteName = function(x) cat("Wrote ",dir,"/",x,"\n",sep="");

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

	qCol = ifelse(fit$q$short=="Time0", "grey",
		ifelse(fit$q$gMed < 50, "red",
		ifelse(fit$q$maxFit > 5, "orange", "darkgreen")));
	qLab = sub("set","",fit$q$name);

	pdf(nameToPath("fit_quality.pdf"),
		pointsize=8, width=8, height=4,
		title=paste(org,"Fitness Quality"));
	par(mfrow=c(1,2));
	plotlab(pmax(0,fit$q$opcor), pmax(0,fit$q$cor12), qLab, col=qCol, ylim=0:1, xlim=0:1,
			     cex=0.5,
			     xlab="Operon Correlation", ylab="1st/2nd-half Correlation",
			     main=org);
	hline(0.5); vline(0.4);
	plotlab(0.1+fit$q$gMed, pmax(0,fit$q$cor12), qLab, col=qCol, ylim=0:1, log="x", cex=0.5,
			    xlab="Median Reads per Gene",
			    ylab="1st/2nd-half Correlation");
	hline(0.5); vline(50);
	dev.off();
	wroteName("fit_quality.pdf");

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

	use = fit$q$short != "Time0";
	lrClust = hclust(as.dist(1-cor(fit$lrn[,as.character(fit$q$name)[use]])));
	pdf(nameToPath("fit_cluster_logratios.pdf"),
		pointsize=8, width=0.25*sum(use), height=8,
		title=paste(org,"Cluster Logratios"));
	plot(lrClust, labels=labelAll[use], main="");
	dev.off();
	wroteName("fit_cluster_logratios.pdf");

	countClust = hclust(as.dist(1-cor(log2(1+fit$gN[fit$gN$locusId %in% fit$genesUsed,-1]))));
	pdf(nameToPath("fit_cluster_logcounts.pdf"),
		pointsize=8, width=0.25*nrow(fit$q), height=8,
		title=paste(org,"Cluster Log Counts"));
	# Some Time0s may be missing from fit$q
	d = match(names(fit$gN)[-1], fit$q$name);
	labelAll = sprintf("%.0f %.2f %s %15.15s",
		      fit$q$gMed[d], fit$q$cor12[d], sub("set","",names(fit$gN)[-1]), ifelse(is.na(d),"Time0",fit$q$short[d]));
	plot(countClust, labels=labelAll, main="");
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

	if(writeImage) {
	    img = format(Sys.time(),"fit%Y%b%d.image"); # e.g., fit2013Oct24.image
	    save(fit, genes, file=nameToPath(img));
	    wroteName(img);
	    unlink(nameToPath("fit.image"));
	    file.symlink(img, nameToPath("fit.image"));
	    cat("Created link for ",nameToPath("fit.image"),"\n");
	}
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

