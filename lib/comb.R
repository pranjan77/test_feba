# Interactive Utilities for viewing fitness data across organisms:

# get_genes(org, list of gene names)
#         given a list or space-delimited set of names, which may be
#	  locus tags (sysName) or VIMSS ids as well as locusIds, return the 
#	  corresponding locusIds, in the same order (or NA for items that do not match).

# gene_info(org, list of gene names)
#	a short report for each of the specified genes, including whether they have specific phenotypes,
#	whether those are conserved, and their top cofitness hits.
#	The gene names are as for get_genes()

# get_exps(org, pattern)
#	returns a list of sucessful experiment ids ("name", i.e. set1H5) that match the pattern

# exp_info(org, list of experiment names)
#	prints information about those experiments

# show_fit(org, list of gene names)
#	an interactive heatmap of the fitness data for those genes
#	The gene names are as for get_genes()

# compare_genes(org, gene1, gene2)
#	interactive scatterplot of fitness data for the two genes
#	The gene names are as for get_genes()

# volcano_fit(org, gene) -- interactive volcano plot (|t| vs. fitness)

# compare_exps(org, exp1, exp2, minT=4)
#	interactive scatterplot of fitness data for the two experiments
#	each name should be like "setH5", or a list of replicates or near-replicates
#	genes with |t1| < minT are in grey
#

# Programmer Utilities:
# get_genes(), get_exps(), get_fit(), get_t(), all_cofitness()

# Use names(orgs) to list which FEBA nicknames are included.

# Data structures (created by LoadOrgs()):
# orgs -- a separate data structure for each bug
#	genes -- metadata about genes, including VIMSS ids when possible
#	exps -- metadata about experiments (including Time 0s)
#	g -- gene ids in the order that they appear in lrn, t
#	q -- quality metrics for fitness experiments, including success or not (u)
#	lrn -- fitness values, 1 row per g, 1 column per successful experiment
#	t -- similarly for t values (but includes all experiments)
#	specsick -- specific phenotypes [these can be positive, actually, so "sick" is a misnomer]
#	cofit -- most cofit hits for each gene
# specsicks -- the union of the specific phenotypes across all organisms
# css -- specific phenotypes that are conserved for orthologs genes (conserved specific sick)
# ccofit -- cofitness that is conserved for orthologous pairs (conserved cofit)
# 	 only has pairs in the order locus1 < locus2

if(file.access("src/feba/lib/FEBA.R")==0) source("src/feba/lib/FEBA.R");	

# note hard-coding of paths, should be fixed
LoadOrgs = function(orgnames) {
	orgs = list();
	for(n in orgnames) {
	       cat("Loading ",n,"\n",sep="");
		exps = read.delim(paste("public_html/FEBA/",n,"/expsUsed",sep=""),as.is=T);
		genes = read.delim(paste("data/FEBA/g/",n,"/genes.tab",sep=""), as.is=T);
		q = read.delim(paste("public_html/FEBA/",n,"/fit_quality.tab",sep=""), as.is=T);
		if(!all(q$name %in% exps$name)) stop(n, ": Unmatched names in q: ", setdiff(q$name,exps$name));
		lrn = read.delim(paste("public_html/FEBA/",n,"/fit_logratios_good.tab",sep=""),check.names=F,as.is=T);
		names(lrn) = sub(" .*","",names(lrn));
		if(is.null(q$u)) q$u = q$name %in% names(lrn);
		g = lrn$locusId;
		lrn = lrn[,-(1:4)];
		if(!all(names(lrn) %in% q$name)) stop(n, ": Unmatched names in lrn: ", setdiff(names(lrn),q$name));

		tval = read.delim(paste("public_html/FEBA/",n,"/fit_t.tab",sep=""),as.is=T,check.names=F);
		tval = tval[,-(1:3)];
		names(tval) = sub(" .*","",names(tval));
		tval = tval[,q$u];
		if(!all(names(tval)==names(lrn))) stop("Invalid names in fit_t.tab");

		file = paste("public_html/FEBA/",n,"/cofit",sep="");
		if (file.access(file,mode=4) != 0) {
			cat("Warning: cannot read cofitness file for ",n,"\n",sep="");
			cofit = NULL;
		} else {
			cofit = read.delim(file, as.is=T);
		}

		file = paste("public_html/FEBA/",n,"/gffmo",sep="");
		if (file.access(file,mode=4) == 0) {# readable
		    gffmo = read.delim(file,as.is=T);
		    genes = merge(genes, gffmo[,words("locusId VIMSS")], all.x=T);
		}

		specsick = read.delim(paste("public_html/FEBA/",n,"/specsick",sep=""),as.is=T);

		orgs[[n]] = list(genes=genes, exps=exps, g=g, q=q, lrn=lrn, t=tval, specsick=specsick, cofit=cofit);
	}

	d = lapply(orgs, function(x) x$specsick);
	for(i in names(d)) d[[i]]$org = i;
	specsicks <<- do.call(rbind,d);

	css <<- read.delim("data/FEBA/css",as.is=T); # conserved specific sick

	cofitpairs = read.delim("data/FEBA/comb_cofit.pairs",as.is=T);
	# organism name is already present, so have locusIds be simple without the org: prefix
	cofitpairs$locus1 = sub(".*:","",cofitpairs$locus1);
	cofitpairs$locus2 = sub(".*:","",cofitpairs$locus2);
	cofitpairs$orth1 = sub(".*:","",cofitpairs$orth1);
	cofitpairs$orth2 = sub(".*:","",cofitpairs$orth2);
	ccofit <<- cofitpairs;

	para <<- read.delim("data/FEBA/aaseqs.para", as.is=T, header=F, col.names=words("org locusId para bits ratio"));
	para <<- split(para, para$org);
	for(org in names(orgs)) {
	    if (!is.null(para[[org]])) {
	        orgs[[org]]$parahits = para[[org]];
	        orgs[[org]]$para = aggregate(para[[org]][,"ratio",drop=F], para[[org]][,"locusId",drop=F], max);
	    }
	}

	orgs <<- orgs;
}

# turn list (or space-delimited) systematic names or VIMSS ids or locusIds into a list of locusIds
get_genes = function(org, specs) {
	genes = orgs[[org]]$genes;
	if(is.null(genes)) stop("Invalid org: ",org);
	if(is.character(specs) && length(specs)==1 && grepl(" ",specs)) {
		specs = words(specs);
	}
	ids = genes$locusId[match(specs,genes$locusId)]; # fetch locusId instead of reusing values to convert to right kind
	if (!is.null(genes$sysName)) ids = ifelse(is.na(ids), genes$locusId[match(specs, genes$sysName)], ids);
	if (!is.null(genes$name)) ids = ifelse(is.na(ids), genes$locusId[match(specs, genes$name)], ids);
	if (!is.null(genes$VIMSS)) ids = ifelse(is.na(ids), genes$locusId[match(specs, genes$VIMSS)], ids);
	return(ids);
}

get_exps = function(org, pattern, ignore.case=T, perl=T) {
	q = orgs[[org]]$q;
	if(is.null(q)) stop("Invalid organism: ",org);
	return(q$name[q$u & grepl(pattern, q$short, ignore.case=ignore.case, perl=perl)]);
}

gene_info = function(org, loci, n=5) {
	arg = loci;
	loci = get_genes(org, loci);
	loci = loci[!is.na(loci)];
	genes = orgs[[org]]$genes;
	info = genes[match(loci, genes$locusId),];
	info$data = info$locusId %in% orgs[[org]]$g;
	if(is.null(info) || nrow(info) == 0) stop("No such gene in ",org,": ",arg);
	out = info[,words("locusId sysName desc data")];
	if(!is.null(info$VIMSS)) out$VIMSS = info$VIMSS;
	row.names(out) = 1:nrow(out);
	print(out);
	cofit = orgs[[org]]$cofit;
	for(locusId in info$locusId[info$data]) {
		out = specsicks[specsicks$org %in% org & specsicks$locusId %in% locusId,words("name cond lrn t")];
		out$conserved = out$cond %in% css$cond[css$locusId %in% locusId & css$tax %in% org];
		locusShow = locusId;
		if (!is.null(info$VIMSS)) locusShow = paste(locusShow, "VIMSS", info$VIMSS[info$locusId %in% locusId]);
		cat("\nSpecific phenotypes for", locusShow,
			info$sysName[info$locusId %in% locusId],
			info$desc[info$locusId %in% locusId], "\n");
		if(nrow(out) >= 1) print(out,digits=2) else cat("None\n");

		if(!is.null(cofit)) {
		    out = cofit[cofit$locusId %in% locusId & cofit$rank <= n,];
		    out$conserved = out$hitId %in% ccofit$locus2[ccofit$tax == org & ccofit$locus1 %in% locusId] |
					out$hitId %in% ccofit$locus1[ccofit$tax == org & ccofit$locus2 %in% locusId];
		    out$locusId = NULL;
		    out = merge(out, genes, by.x="hitId", by.y="locusId");
		    out = out[order(out$rank),];
		    row.names(out) = 1:nrow(out);
		    out2 = out;
		    out2 = out2[,words("cofit hitId sysName desc conserved")];
		    if(!is.null(out$VIMSS)) out2$VIMSS = out$VIMSS;
		    cat("\nTop cofit hits for", locusShow,
			info$sysName[info$locusId %in% locusId],
			info$desc[info$locusId %in% locusId], "\n");
		    print(out2,digits=2);
		}
	}
}

exp_info = function(org, exps) {
	q = orgs[[org]]$q;
	if(is.null(q)) stop("Invalid org: ",org);
	if(!all(exps %in% q$name)) stop("Invalid exps: ", setdiff(exps,q$name));
	print(q[q$name %in% exps,], digits=2);

}

get_fit = function(org, loci, t=FALSE) {
	loci = get_genes(org, loci);
	g = orgs[[org]]$g;
	q = orgs[[org]]$q;
	indexes = match(loci, g);

	if(t == TRUE) {
		return(t(orgs[[org]]$t[indexes,]));
	} else {
		return(t(orgs[[org]]$lrn[indexes,]));
	}
}

get_t = function(org, loci) get_fit(org, loci, t=TRUE);

identify_exps = function(org, x, y) {
	q = orgs[[org]]$q;
	q = q[q$u,];
	if(is.null(q)) stop("Invalid org: ",org);
	if(length(x) != nrow(q)) stop("Incorrect lengths: x ",length(x), " metadata ",nrow(q));
	if(length(y) != nrow(q)) stop("Incorrect lengths: y ",length(y), " metadata ",nrow(q));
	rows = identify(x, y, q$name);
	q$x = c(x,recursive=T);
	q$y = c(y,recursive=T);
	return(q[rows,words("name short x y")]);
}

compare_genes = function(org, locus1, locus2=NULL, xlab=NULL, ylab=NULL, eq=TRUE, locate=TRUE, ...) {
	if (is.null(locus2)) {
		if(length(locus1) != 2) stop("Invalid length");
		locus2 = locus1[2];
		locus1 = locus1[1];
	}
	if(is.null(xlab)) xlab=locus1;
	if(is.null(ylab)) ylab=locus2;
	locus1 = get_genes(org, locus1);
	locus2 = get_genes(org, locus2);
	if (is.na(locus1)) stop("No data for: ",locus1, " ", xlab);
	if (is.na(locus2)) stop("No data for: ",locus2, " ", ylab);
	mat = get_fit(org, c(locus1,locus2));
	plot(mat[,1], mat[,2], xlab=xlab, ylab=ylab, ...);
	if(eq) eqline();
	if(locate) {
		cat("Click on points, or right click to exit\n");
		identify_exps(org, mat[,1], mat[,2]);
	}
}

compare_exps = function(org, exp1, exp2, xlab=exp1, ylab=exp2, eq=TRUE, locate=TRUE, minT=4, col=NULL, ...) {
	if (is.null(orgs[[org]]$lrn[,exp1])) stop("No experiment ",exp1," in ",org);
	if (is.null(orgs[[org]]$lrn[,exp2])) stop("No experiment ",exp2," in ",org);
	x = rowMeans(orgs[[org]]$lrn[,exp1,drop=F]);
	y = rowMeans(orgs[[org]]$lrn[,exp2,drop=F]);
	tx = rowMeans(orgs[[org]]$t[,exp1,drop=F]);
	if(is.null(col)) col=ifelse(abs(tx) >= minT,"black","darkgrey");
	plot(x, y, xlab=xlab, ylab=ylab, col=col, ...);
	if(eq) eqline();
	if(locate) identify_genes(org, x, y);
}

identify_genes = function(org, x, y, col=2, newcol=2, ...) {
	genes = orgs[[org]]$genes;
	g = orgs[[org]]$g;
	if(length(x) != length(g)) stop("Incorrect lengths");
	if(length(y) != length(g)) stop("Incorrect lengths");
	cat("Click on points, or right click to exit\n");
	n = 1;
	while(TRUE) {
	    i = identify(x, y, as.character(n), n=1, col=col, ...);
	    if(length(i) != 1) break;
	    row = genes[genes$locusId == g[i],];
	    if(nrow(row) != 1) stop("Illegal index ",i);
	    id = as.character(row$locusId);
	    if(!is.null(row$VIMSS)) id = paste(id, row$VIMSS);
	    cat(sprintf("%d: %.3f %.3f %s %s %s\n", n, x[i], y[i], id, row$sysName, row$desc));
	    if(!is.null(newcol)) points(x[i],y[i],col=newcol);
	    n = n+1;
	}
}

volcano_fit = function(org, locus, xlab="Fitness", ylab="|t|", locate=TRUE, ...) {
	locus = get_genes(org,locus);
	if (length(locus) != 1) stop("must input one locus");
	x = get_fit(org, locus);
	if(is.na(x[1])) stop("No data for: ",locus);
	tval = get_t(org, locus);
	plot(x, abs(tval), xlab=xlab, ylab=ylab, ...);
	if(locate) {
		cat("Click on points, or right click to exit\n");
		identify_exps(org, x, abs(tval));
	}
}

eqline <- function(col="grey",lty=2,lwd=1) {
	x <- 10**(-25:25);
	lines(c(-rev(x),x),c(-rev(x),x),col=col,lty=lty,lwd=lwd);
}

andNoNA = function(x) ifelse(is.na(x),FALSE,x);

# around: show the gene neighborhood around the gene of interest
#	Gene order will be flipped if gene is on the - strand
#	Only an option if only one locus specified
# condspec: choose a subset of conditions to show, as a perl-style regular expression against q$short
# scale: by default, colors are for fitness = -2 to +2. To focus on stronger differences use scale=3.
# ylabmax: maximum number of characters for the condition labels
# cex: adjust the font size (e.g. cex=0.5 for tiny text)
#
show_fit = function(org, loci, labels=NULL, locate=TRUE, around=0, condspec=NULL, scale=2, ylabmax=20, cex=1) {
	if (length(loci==1) && grepl(" ",loci)) loci = words(loci);
	if(is.null(labels)) labels = loci;
	loci = get_genes(org, loci);
	if (length(loci)==1 && around > 0) {
	    genes = orgs[[org]]$genes;
	    if(!is.null(genes$scaffold)) genes = genes[genes$scaffold == genes$scaffold[genes$locusId==loci], ];
	    strand = genes$strand[genes$locusId==loci];
	    if(!is.null(genes$begin)) genes = genes[order(genes$begin),];
	    i = match(loci, genes$locusId);
	    i1 = pmax(1, i - around);
	    i2 = pmin(nrow(genes), i + around);
	    indexes = if(strand=="+") i1:i2 else i2:i1;
	    loci = genes$locusId[i1:i2];
	    labels = genes$sysName[i1:i2];
	}
	mat = get_fit(org, loci); # columns as genes, rows as experiments
	tval = get_t(org, loci);

	q = orgs[[org]]$q;
	q = q[q$u,];
	if(!is.null(condspec)) {
	    u = grepl(condspec,q$short,perl=T);
	    if(sum(u) == 0) stop("No conditions matching: ",condspec);
	    q = q[u,];
	    mat = mat[u,];
	    tval = tval[u,];
	}
	if(length(q$short) != nrow(mat)) stop("Wrong number of rows");
	if (!all(q$name == row.names(mat))) stop("Name mismatch");

	o = order(q$short,decreasing=TRUE);
	labRows = q$short[o];
	labRows = sub("[0-9.]+ mM$","",labRows);
	labRows = sub("[0-9.]+ mg/ml$","",labRows);
	labRows = sub("[0-9.]+ vol%$","",labRows);
	# note labRows are shown going up so we want the 1st row to be replaced if a duplicate
	labRows = ifelse(andNoNA(labRows == c(labRows[-1],NA)), ".", labRows);

	# mar is bottom,left,top,right
	ylablen = pmin(ylabmax, max(nchar(labRows)));
	oldpar = par(mgp=c(2,1,0), mar=c(cex*(1+max(nchar(labels))),0.5,0.5,cex*ylablen));
	image(t(mat[o,]), col=myHeatColors(), breaks=breaksUse(scale=scale), useRaster=TRUE,
		xlab="", ylab="", xaxt="n", yaxt="n", bty="n");
	# las = 2 (perpendicular to axis)
	mtext(labRows, las=2, side=4, at=seq(0,1,length.out=length(labRows)), cex=cex); # at right
	mtext(labels, las=2, side=1, at=seq(0,1,length.out=length(labels)), cex=cex); # at bottom

	genes = orgs[[org]]$genes;
	if (locate) {
		cat("Click on cells, or right click to exit\n");
		while(TRUE) {
		    at = locator(1);
		    if(is.null(at)) break;
		    x = at$x[1];
		    y = at$y[1];
		    iGene = pmax(1, pmin(length(labels), 1 + round(x * (length(labels)-1))));
		    g = loci[iGene];
		    geneinfo = genes[genes$locusId == g,];
		    if (nrow(geneinfo) != 1) stop("no metdata for gene: ",g);

		    # iOExp is the sorted order shown, iExp is the order in q
		    iOExp = pmax(1, pmin(length(labRows), 1 + round(y * (length(labRows)-1))));
		    iExp = o[iOExp];
		    locusString = as.character(g);
		    if (!is.null(geneinfo$VIMSS)) locusString = sprintf("%s VIMSS %d",locusString,geneinfo$VIMSS);
		    cat(sprintf("Gene: %s %s %s\nExperiment: %s %s\nFitness: %.3f (t %.3f)\n\n",
				locusString, geneinfo$sysName, geneinfo$desc,
				q$name[iExp], q$short[iExp],
				mat[iExp,iGene], tval[iExp,iGene]));
		}
	}
	par(oldpar);
}

breaksUse = function(scale=2) 
{
    d = seq(-3, 3, 0.25) * scale/3;
    d[1] = -100;
    d[length(d)] = 100;
    return(d);
}

# blue to yellow
myHeatColors = function () {
    c(rgb(0, 0, seq(1, 1/12, -1/12)), rgb(seq(1/12, 1, 1/12), 
        seq(1/12, 1, 1/12), 0));
}

# all_cofitness(org, genes) -- compute cofitness of average profile of genes with all
# 		     other genes that have fitness.
# Returns a table sorted by correlation (descending order).
# 		     
all_cofitness = function(org, genes) {
	vec = unlist(rowMeans(get_fit(org, genes)));
	out = data.frame(locusId = orgs[[org]]$g, r = c(cor(vec, t(orgs[[org]]$lrn)), recursive=T));
	out = merge(orgs[[org]]$genes, out);
	out = out[order(-out$r),];
	return(out);
}
