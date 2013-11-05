# FEBA.R -- analysis scripts for barcode sequencing data
#
# The key routines are:
#
# AvgStrainFitness -- compute fitness values from counts for the post-experiment counts
#     and the Time0 counts
# NormalizeByScaffold -- normalize the fitness values by scaffold and position
# FEBA_Fit() -- analyze many fitness experiments with AvgStrainFitness() and NormalizeByScaffold()
#      returns a complex data structure
# FEBA_Save_Tables() -- Save the fitness data structure to tab-delimited files and an R image

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
AvgStrainFitness = function(strainCounts, strainT0, strainLocus,
		 minStrainT0 = 4, minGeneT0 = 40,
		 nACG = NULL,
		 genesUsed=NULL, strainsUsed=NULL, debug=FALSE) {
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

    strainFit = mednorm(log2(1+strainCounts) - log2(1+strainT0));
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
    geneFit1 = aggregate(strainFit, list(locusId=strainLocus), mean);
    geneFit1$x = mednorm(geneFit1$x);
    i = match(strainLocus, geneFit1$locusId); # from strain index to gene index
    readratio = sum(strainCounts) / sum(strainT0);
    strainPseudoCount = 2**geneFit1$x[i] * readratio;

    # And apportion the pseudocount equally (in log space) between condition-count and strain-count
    # to minimize the deviations from pseudocount = 1
    condPseudoCount = sqrt(strainPseudoCount);
    t0PseudoCount = 1/sqrt(strainPseudoCount);
    # (or could do some sort of weighted likelihood-based inference of fitness values, might be better)

    # for each strain: fitness, variance, and weight
    strainFit = log2(condPseudoCount + strainCounts) - log2(t0PseudoCount + strainT0) - strainFitAdjust;
    strainVar = sqrt(1/(1+strainT0) + 1/(1+strainCounts)) / log(2);
    strainWeight = 0.5 + strainT0;

    geneFit2 = aggregate(1:length(strainT0), list(locusId=strainLocus),
    	                function(j) {
			    totw = sum(strainWeight[j]);
			    meanFit = sum(strainWeight[j] * strainFit[j]) / totw;
			    c(fit = meanFit,
			    	  sd = sqrt(sum(strainWeight[j]**2 * strainVar[j]))/totw,
				  sdObs = sqrt(sum(strainWeight[j]**2 * (strainFit[j]-meanFit)**2))/totw,
				  sdNaive = sqrt( 1/(1+sum(strainCounts[j])) + 1/(1+sum(strainT0[j])) ) / log(2),
				  n = length(j),
				  # effective degrees of freedom is n if all weights are even,
				  # otherwise is this much higher
				  nEff = length(j) / sqrt(length(j)) * sqrt(sum(strainWeight**2)) / sum(strainWeight) );
    });
    # a hack to get around the lists within the entries
    geneFit2 = data.frame(locusId = geneFit2$locusId, geneFit2[,2]);

    geneFit2$fitRaw = geneFit2$fit; # without the median normalization
    geneFit2$fit = mednorm(geneFit2$fit);

    geneCounts = aggregate(list(t0=strainT0,counts=strainCounts), list(locusId=strainLocus), sum);
    geneFit2$fitNaive = mednorm(log2(1+geneCounts$counts) - log2(1+geneCounts$t0));
    return(geneFit2);
}

# values -- fitness values (as a vector)
# locusId -- the corresponding locusIds
# genes contains locusId, scaffoldId, and begin
# minForLoess -- how many genes a scaffold needs to contain before we try to use Loess to eliminate bias
#     For those scaffolds, it also estimates the mode and removes that
# minToUse -- if a scaffold has too few genes, cannot correct for possible DNA extraction bias
# 	   so need to remove data for that gene
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

FEBA_Fit = function(expsUsed, all, all_g2, genes, genesUsed=NULL, strainsUsed=NULL,
	   		       minT0Strain=3, minT0Gene=30,
			       minT0GeneSide=minT0Gene/2,
			       minGenesPerScaffold=10,
			       nACG=NULL, # normalize by ACG content if set
			       pred=CrudeOp(genes),
			       okLane=TRUE, # OK to use t0 from another lane if needed?
	   		       metacol=1:7,
			       # names of experiments to ignore
			       ignore=NULL) {

	if(!all(expsUsed$name == names(all_g2)[-metacol]))
		stop("names do not match");
	if(!all(expsUsed$name %in% names(all)))
		stop("names missing from all");
	if(is.null(genes$scaffoldId)) stop("No scaffold for genes");
	if(is.null(genes$begin)) stop("No begin for genes");

	if(!is.null(ignore)) {
	    cat("Ignoring ",ignore,"\n");
	    expsUsed = expsUsed[!expsUsed$name %in% ignore,];
	    all = all[, names(all) != ignore,];
	    all_g2 = all_g2[, names(all) != ignore,];
	}

	expsUsed$name = as.character(expsUsed$name);

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

	t0_g2 = data.frame(lapply(expsT0, function(x) MyRowSums(all_g2[,x])), check.names=F);
	t0_gN = aggregate(t0_g2, list(locusId=all_g2$locusId), sum);
	cat("Reads per t0set, in millions:\n");
	print(colSums(t0_g2)/1e6, digits=2);

	if(is.null(strainsUsed)) strainsUsed = MyRowMeans(t0_g2[,-1]) >= minT0Strain;
	if(is.null(genesUsed)) {
		d = aggregate(t0_g2[strainsUsed,], list(locusId=all_g2$locusId[strainsUsed]), sum);
		genesUsed = d$locusId[ MyRowMeans(d[,-1]) >= minT0Gene ];
	}
	genesPerScaffold = table(genes$scaffoldId[genes$locusId %in% genesUsed]);
	smallScaffold = names(genesPerScaffold)[genesPerScaffold < minGenesPerScaffold];
	if (length(smallScaffold) > 0) cat("Ignoring genes on small scaffolds ",smallScaffold,"\n");
	genesUsed = genesUsed[!genesUsed %in% genes$locusId[genes$scaffoldId %in% smallScaffold]];

	if(length(strainsUsed) != nrow(all_g2)) stop("Invalid strainsUsed");
	cat("Using ",sum(strainsUsed)," of ",nrow(all_g2)," genic strains\n");
	if(length(genesUsed) < 100 || !all(genesUsed %in% genes$locusId)) stop("Invalid genesUsed");

	cat("Using ",length(genesUsed)," of ",length(unique(all_g2$locusId))," genes with data\n");

	d1 = aggregate(t0_g2[strainsUsed & all_g2$f < 0.5,],
	     		list(locusId=all_g2$locusId[strainsUsed & all_g2$f < 0.5]), sum);
	d2 = aggregate(t0_g2[strainsUsed & all_g2$f >= 0.5,],
	     		list(locusId=all_g2$locusId[strainsUsed & all_g2$f >= 0.5]), sum);
	genesUsed12 = intersect(d1$locusId[ MyRowMin(d1[,-1]) >= minT0GeneSide],
		      		d2$locusId[ MyRowMin(d2[,-1]) >= minT0GeneSide]);
        cat("For cor12, using ",length(genesUsed12),"genes\n");

	all_lrs = lapply(names(all_g2)[-metacol], function(n)
		AvgStrainFitness(all_g2[[n]], t0_g2[[ expsUsed$t0set[expsUsed$name==n] ]], all_g2$locusId,
		         genesUsed = genesUsed, strainsUsed = strainsUsed, nACG=nACG));
	names(all_lrs) = names(all_g2)[-metacol];
	all_lr = data.frame(locusId=all_lrs[[1]]$locusId, lapply(all_lrs, function(x) x$fit));
	all_lr_naive = data.frame(locusId=all_lrs[[1]]$locusId, lapply(all_lrs, function(x) x$fitNaive));
	all_sd = data.frame(locusId=all_lrs[[1]]$locusId, lapply(all_lrs, function(x) x$sd));
	all_sdObs = data.frame(locusId=all_lrs[[1]]$locusId, lapply(all_lrs, function(x) x$sdObs));
	all_sdNaive = data.frame(locusId=all_lrs[[1]]$locusId, lapply(all_lrs, function(x) x$sdNaive));
	all_nEff = data.frame(locusId=all_lrs[[1]]$locusId, lapply(all_lrs, function(x) x$nEff));
	nStrain = data.frame(locusId=all_lr$locusId, n=all_lrs[[1]]$n); # should be the same each time

	all_lrn = data.frame(locusId=all_lr$locusId, apply(all_lr[,-1], 2,
		  		NormalizeByScaffold, all_lr$locusId, genes,
				    minToUse=minGenesPerScaffold, debug=T));

	all_lrs1 = lapply(names(all_g2)[-metacol], function(n)
		AvgStrainFitness(all_g2[[n]], t0_g2[[ expsUsed$t0set[expsUsed$name==n] ]], all_g2$locusId,
		         genesUsed = genesUsed12, strainsUsed = strainsUsed & all_g2$f < 0.5));
	names(all_lrs1) = names(all_g2)[-metacol];
	all_lrs2 = lapply(names(all_g2)[-metacol], function(n)
		AvgStrainFitness(all_g2[[n]], t0_g2[[ expsUsed$t0set[expsUsed$name==n] ]], all_g2$locusId,
		         genesUsed = genesUsed12, strainsUsed = strainsUsed & all_g2$f >= 0.5));
	names(all_lrs2) = names(all_g2)[-metacol];

	all_lr1 = data.frame(locusId=all_lrs1[[1]]$locusId, lapply(all_lrs1, function(x) x$fit), check.names=F);
	all_lr2 = data.frame(locusId=all_lrs2[[1]]$locusId, lapply(all_lrs2, function(x) x$fit), check.names=F);
	if(!all(all_lr1$locusId == all_lr2$locusId)) stop("Mismatch of locusId");
	nStrain12 = data.frame(locusId=all_lr1$locusId, nStrain1=all_lrs1[[1]]$n, nStrain2=all_lrs2[[1]]$n);

	q = data.frame(name = expsUsed$name,
			short = expsUsed$short,
			nMapped = colSums(all[,expsUsed$name]));
	q$nPastEnd = colSums(all[all$scaffold=="pastEnd",expsUsed$name]);
	q$nGenic = colSums(all_g2[,-(1:7)]);
	q$nUsed = colSums(all_g2[strainsUsed,-(1:7)]);
	q$gMed = apply(all_gN[,-1],2,median);
	median_t0_gN = if(ncol(t0_gN) > 2) apply(t0_gN[,-1], 2, median) else t0_gN[,2];
	q$gMedt0 = median_t0_gN[ match(expsUsed$t0set, names(t0_gN)[-1]) ];
	q$gMean = apply(all_gN[,-1],2,mean);
	q$cor12 = sapply(as.character(q$name), function(n) cor(rank(all_lr1[[n]]), rank(all_lr2[[n]])));
	q$opcor = apply(all_lrn[,-1], 2, function(x) paircor(pred[pred$bOp,], all_lrn$locusId, x, method="s"));
	q$maxFit = apply(all_lr[,-1],2,max);
	
	return(list(expsUsed = expsUsed,
			gN = all_gN,
			t0_g2 = t0_g2,
			t0_gN = t0_gN,
			lr = all_lr,
			lrNaive = all_lr_naive,
			sd = all_sd,
			sdObs = all_sdObs,
			sdNaive = all_sdNaive,
			nEff = all_nEff,
			nStrain = nStrain,
			lrn = all_lrn,
			lr1 = all_lr1,
			lr2 = all_lr2,
			nStrain12 = nStrain12,
			q = q,
			genesUsed = genesUsed,
			strainsUsed = strainsUsed));
}

FEBA_Save_Tables = function(fit, genes, org="?",
		 topdir="public_html/tmp",
		 dir = paste(topdir,org,sep="/"),
		 writeImage=TRUE) {
	if(!file.exists(dir)) stop("No such directory ",dir);

	for (n in words("q lr lr1 lr2 lrn sd")) {
	    if (is.null(fit[[n]]) || !is.data.frame(fit[[n]])) {
	        stop("Invalid or missing ",n," entry");
	    }
	}
	if (is.null(fit$genesUsed)) stop("Missing genesUsed");

	if(!all(names(fit$lr)[-1] == fit$q$name)) stop("Name mismatch");
	if(!all(names(fit$lrn)[-1] == fit$q$name)) stop("Name mismatch");

	nameToPath = function(filename) paste(dir,filename,sep="/");
	wroteName = function(x) cat("Wrote ",dir,"/",x,"\n",sep="");

	writeDelim(fit$q, nameToPath("fit_quality.tab"));
	wroteName("fit_quality.tab");

	writeDelim(cbind(genes, used=genes$locusId %in% fit$genesUsed), nameToPath("fit_genes.tab"));
	wroteName("fit_genes.tab");

	d = merge(genes[,c("locusId","sysName","desc")], fit$lr);
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_logratios_unnormalized.tab"));
	wroteName("fit_logratios_unnormalized.tab");

	d = merge(genes[,c("locusId","sysName","desc")], fit$lrNaive);
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_logratios_unnormalized.tab"));
	wroteName("fit_logratios_unnormalized_naive.tab");

	d = merge(genes[,c("locusId","sysName","desc")], fit$lrn);
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_logratios.tab"));
	wroteName("fit_logratios_normalized.tab");

	d = merge(genes[,c("locusId","sysName","desc")], fit$sd);
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_standard_error.tab"));
	wroteName("fit_standard_error.tab");

	d = merge(genes[,c("locusId","sysName","desc")], fit$sdObs);
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_standard_error_obs.tab"));
	wroteName("fit_standard_error_obs.tab");

	d = merge(genes[,c("locusId","sysName","desc")], fit$nEff);
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_standard_error_obs.tab"));
	wroteName("fit_nEffective.tab");

	d = merge(genes[,c("locusId","sysName","desc")], fit$sdNaive);
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
	    changers = fit$lrn$locusId[abs(fit$lrn[[n]]/fit$sd[[n]]) >= 3];
	    plot(fit$lr1[[n]], fit$lr2[[n]],
	    		  main=sprintf("Cor12 for %s %s (%.0f %.3f)\n%s",
			  	org, n, fit$q$gMed[i], fit$q$cor12[i], fit$q$short[i]),
	    		  xlab="First Half", ylab="Second Half",
			  col=ifelse(fit$lr1$locusId %in% changers, 2, 1));
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
	labelAll = sprintf("%.0f %.2f %s %15.15s",
		      fit$q$gMed, fit$q$cor12, sub("set","",fit$q$name), fit$q$short);
	plot(countClust, labels=labelAll, main="");
	dev.off();
	wroteName("fit_cluster_logcounts.pdf");

	if(!all(fit$lr$locusId == fit$lrn$locusId)) stop("locusId mismatch");
        d = table(genes$scaffoldId[genes$locusId %in% fit$genesUsed]);
	maxSc = names(d)[which.max(d)];
	if (is.null(maxSc)) stop("Invalid scaffoldId?");
	beg = ifelse(fit$lr$locusId %in% genes$locusId[genes$scaffold==maxSc],
	    genes$begin[match(fit$lr$locusId, genes$locusId)], NA);

	pdf(nameToPath("fit_chr_bias.pdf"), pointsize=10, width=5, height=5,
	          title=paste(org,"Chromosome Bias"));
	for (i in 1:nrow(fit$q)) {
	    n = as.character(fit$q$name[i]);
	    plot(beg, pmax(-2,pmin(2,fit$lr[[n]])),
	    		  main=sprintf("%s %s (%.0f %.3f)\n%s",
			  	org, n, fit$q$gMed[i], fit$q$cor12[i], fit$q$short[i]),
	    		  xlab="Position on Main Scaffold",
			  ylab="Fitness",
			  ylim=c(-2,2), col="darkgrey");
	    o = order(beg);
	    lines(beg[o], (fit$lr[[n]] - fit$lrn[[n]])[o], col="darkgreen", lwd=2);
	    hline(0);
	}
	dev.off();
	wroteName("fit_chr_bias.pdf");

	if(writeImage) {
	    img = format(Sys.time(),"fit%Y%b%d.image"); # e.g., fit2013Oct24.image
	    save(fit, genes, expsUsed, file=nameToPath(img));
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

# A work around for when we select just 1 column and then cannot use rowMeans on it
MyRowMeans = function(x) {
	   if(is.vector(x)) return(x);
	   rowMeans(x);
}
MyRowSums = function(x) {
	   if(is.vector(x)) return(x);
	   rowSums(x);
}
MyRowMin = function(x) {
	 if (is.vector(x)) return(x);
	 apply(x, 1, min);
}

# For reformatting per-lane pool count tables into the "all" table
prefixName = function(x, prefix) { names(x) = paste(prefix,names(x),sep=""); return(x); }

# For shortening the experminet descriptions
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
