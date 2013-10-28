autoload("loess","modreg");

### Common operations

countif <- function(x) { return(sum(ifelse(x,1,0))); }

countPerIf <- function(x,condition=NULL) {
	if(!is.null(condition)) { x <- select(x,condition); }
	 return (sum(ifelse(x,1,0))/length(x));
}

without <- function(list,columns=list$without) {
	list2 <- list;
	for (i in columns) { list2[[i]] <- NULL; }
	return(list2);
}

na0 = function(x) ifelse(is.na(x),0,x);

withoutNA <- function(x) { x[!is.na(x)]; }

andNoNA <- function(x) { ifelse(is.na(x),FALSE,x) }
orNA <- function(x) { ifelse(is.na(x),TRUE,x) }

sumNoNA <- function(x) { return(sum(x[!is.na(x)])); }
meanNoNA <- function(x) { return(mean(x[!is.na(x)])); }

rankNA <- function(x) {
	out <- rep(NA,length(x));
	out[!is.na(x)] <- rank(x,na.last=NA); # NA values removed
	return(out);
}

percentileNA <- function(x) {
	rankNA(x)/countif(!is.na(x));
}

findIntervalAvg <- function(x,sorted) {
	index <- findInterval(x,sorted);
	return( ifelse(x < sorted[1] | x > sorted[length(sorted)-1], index,
			rank(sorted)[ifelse(index<1,1,index)] ) );
}

reverse <- function(perm) {
	r <- 1:length(perm);
	r[perm] <- 1:length(perm);
	return(r);
}

# returns a vector of matching indices or -1 in index
findMatch <- function(queries,targets) {
	o <- order(targets); # index -> sorted index
	ro <- reverse(o);
	s <- targets[o]; # sorted values
	ns <- findInterval(queries,s); # query -> sorted index
	n <- ifelse(ns<1,-1,o[ifelse(ns<1,1,ns)]); # query -> target
	return(list(o=o,ro=ro,ns=ns,n=n,index=ifelse(n>=1&queries==targets[ifelse(n<1,1,n)],n,-1)));
}

select <- function(X,Y) { return(split(X,Y)$"TRUE"); }

RenameColumns <- function(data, newnames, keep=FALSE) {
	d2 <- list();
	if (keep) {
		for (n in names(data)) {
			if(!is.null(newnames[[n]])) {
				d2[[ newnames[[n]] ]] = data[[n]];
			} else {
				d2[[ n ]] = data[[n]];
			}
		}
	} else {
		for (n in names(newnames)) {
			if(is.null(data[[n]])) warning("missing name ", n);
			d2[[ newnames[[n]] ]] <- data[[n]];
		}
	}
	return( d2 );
}

skewness <- function(x) {
	x <- x[!is.na(x)];
	dx <- x - mean(x);
	n <- length(x);
	m2 <- mean( dx**2 );
	k2 <- m2*n/(n-1);
	m3 <- mean( dx**3 );
	k3 <- m3 * n**2/((n-1)*(n-2));
	return(k3/k2**1.5);
}

# see Mathematica kurtosis page and k-statistics page
kurtosis <- function(x) {
	x1 <- x[!is.na(x)];
	dx <- x1 - mean(x1);
	n <- length(dx);
	S1 <- sum(dx);
	S2 <- sum(dx^2);
	S3 <- sum(dx^3);
	S4 <- sum(dx^4);
	k2 <- (n*S2-S1^2)/(n*(n-1));
	k4 <- (-6*S1^4 + 12*n*S1^2*S2 - 3*n*(n-1)*S2^2 - 4*n*(n+1)*S1*S3 + n^2*(n+1)*S4)/(n*(n-1)*(n-2)*(n-3));
	return(k4/(k2^2));
}

heavytail <- function(x, q=0.01) {
	x1 <- x[!is.na(x)];
	z <- (x-mean(x))/sd(x);
	z0 <- -qnorm(q);
	countPerIf(abs(z)>z0);
}

pairIndices <- function(x) {
	n <- length(x);
	return(data.frame(i1=1:(n-1),i2=2:n));
}

pairDiff <- function(x) {
	return (x[1:(length(x)-1)] - x[2:length(x)]);
}

combineColumns <- function(x) {
	n <- length(x[,1]);
	ntot <- n*length(names(x));
	out <- data.frame(x=rep(0,ntot),type=factor(rep(names(x)[1],ntot),levels=names(x)));
	for (i in 1:length(names(x))) {
		out$x[(i-1)*n+(1:n)] <- x[,i];
		out$type[(i-1)*n+(1:n)] <- levels(out$type)[i];
	}
	return(out);
}

# average information in a vector of probabilities
avginfo <- function(x) {
	return(mean(-log(ifelse(x<1e-8,1,x)) * x));
}

certainty <- function(x) {
	return(2*mean( abs(x - .5) ));
}

info <- function(x) { y <- x/sum(x); return(sum(-log(ifelse(y<1e-8,1,y))*y)); }

relinfo <- function(x, standard=rep(1,length(x))) {
	x <- x+1;
	x <- x/sum(x);
	standard <- standard+1;
	standard <- standard/sum(standard);
	return(sum(x * log(x/standard)));
}

combineP <- function(p1, p2) {
	top <- p1/(1-p1) * p2/(1-p2);
	pm <- pmax(p1,p2);
	return(ifelse(pm >= .99999, pm, top/(1+top)));
}

logodds <- function(p) { return(log(p/(1-p))); }

logodds2p <- function(x) { return( 1/(1+exp(-x)) ); }

cut4 <- function(x,...) { return(cut(x,breaks=quantile(x),include.lowest=TRUE,...)); }
cutN <- function(x,n=8,...) { return(cut(x,breaks=unique(quantile(x,seq(0,1,1/n))),include.lowest=TRUE,...)); }

ks2 <- function(x) { return(ks.test(x$"TRUE",x$"FALSE")); }
wilcox2 <- function(x) { return(wilcox.test(x$"TRUE",x$"FALSE")); }

ksBy <- function(x, by, restrict=rep(TRUE,length(by))) {
	return(ks.test(x[restrict & by], x[restrict & !by]));
}

mi <- function(x,y,breaks) {
	rangeB <- range(breaks);
	# force to within range
	x1 <- pmin(pmax(x,rangeB[1]),rangeB[2]);
	y1 <- pmin(pmax(y,rangeB[1]),rangeB[2]);
	factorX <- cut(x1,breaks=c(rangeB[1]-1,breaks));
	factorY <- cut(y1,breaks=c(rangeB[1]-1,breaks));
	countX <- table(factorX);
	countY <- table(factorY);
	countXY <- table(countX,countY);
	return( info(as.vector(countX))+info(as.vector(countY))-info(as.vector(countXY)) );
}

# formula should be Y~X1+X2+...+XN
plotN <- function(formula,frame,mf,...) {
	oldpar <- par(no.readonly=TRUE);
	par(mfrow=mf);
	plot.formula(formula,frame,ask=FALSE,...);
	on.exit(par(oldpar));
}

# compute a percent >= vs. %(combined rank) plot, like for the Kolmogorov-Smirnov test
ksRanks <- function(list) {
	join <- c();
	for (name in names(list)) { join <- c(join, list[[name]]); }
	ranksJ <- floor(rank(join)); # need integers for tabulate

	out <- list();
	out$per = (1:length(join))/length(join); # percentile of the union
	index <- 0;
	for (name in names(list)) {
		ranksN <- ranksJ[index + (1:length(list[[name]]))];
		countRanks <- cumsum(tabulate(floor(ranksN),nbins=length(ranksJ)));
		out[[name]] <- countRanks/length(list[[name]]);
		index <- index + length(list[[name]]);
	}
	return(data.frame(out));
}

DratioByCutoff <- function(values, score, default=score>.5,
			minS=10, maxS=min(countif(default),countif(!default)), step=40) {
	HoverL <- countif(default)/countif(!default);
	Dall <- ks.test(select(values,!is.na(values) & default),select(values,!is.na(values) & !default))$statistic;
	sorted <- sort(score);
	Dshort <- vector();
	Dlocal <- vector();
	n <- vector();
	x <- minS;
	while (x + step < maxS && HoverL*(x+step) < maxS) {
		cutoffLow <- sorted[x];
		cutoffHigh <- sorted[length(sorted)-HoverL*x];
		cutoffL2 <- sorted[x+step];
		cutoffH2 <- sorted[length(sorted)-HoverL*(x+step)];
		Dshort[length(Dshort)+1] <- ks.test(select(values,!is.na(values) & score <= cutoffLow),
							select(values,!is.na(values) & score >= cutoffHigh))$statistic;
		Dlocal[length(Dlocal)+1] <- ks.test(select(values,!is.na(values) & score<=cutoffL2 & score>=cutoffLow),
						select(values,!is.na(values) & score<=cutoffHigh & score>=cutoffH2)
						)$statistic;
		n[length(n)+1] <- x;
		x <- x + step;
	}
	Dshort[length(Dshort)+1] <- Dall;
	Dlocal[length(Dlocal)+1] <- Dlocal[length(Dlocal)];
	n[length(n)+1] <- min(countif(default),countif(!default));
	return(data.frame(f=n/max(n),n=n,Dshort=Dshort,Dratio=Dall/Dshort,Dlocal=Dlocal));
}

# mutual information after splitting at the cutoff, versus the pooled set
# not sure this is right
miSplit <- function(x,y,cutoff) {
	xL <- sum(ifelse(x<cutoff,1,0));
	xH <- sum(ifelse(x>=cutoff,1,0));
	yL <- sum(ifelse(y<cutoff,1,0));
	yH <- sum(ifelse(y>=cutoff,1,0));
	return( info(c(xL,xH)+info(c(yL,yH)-info(c(xL+yL,xH+yH)))) );
}

# partial correlation between v1&v2, controlling for v3
# p is two-tailed probability
# using pcorTest(tries=1e4) found correct rates of p<0.05, p<0.01, p<0.001
# Is anova a better test? Not sure
pcor <- function(v1, v2, v3,spearman=FALSE) {
	if (spearman) { v1 <- rank(v1); v2 <- rank(v2); v3 <- rank(v3); }
        c12 <- cor(v1, v2)
        c23 <- cor(v2, v3)
        c13 <- cor(v1, v3)
	r <- (c12-(c13*c23))/(sqrt(1-(c13^2)) * sqrt(1-(c23^2)));
	df <- length(v1)-3;
	t <- r/sqrt((1-r*r)/df);
	return(data.frame(r=r,df=df,t=t,p=2*pt(-abs(t),df=df)));
}

pcorTest <- function(n=100, tries=100, r1=.3, r2=.3, r12=0) {
	out <- NULL;
	for (i in 1:tries) {
		y <- rnorm(n);
		x12 <- rnorm(n);
		x1 <- r1*y + r12*x12 + sqrt(1-r1**2-r12**2)*rnorm(n);
		x2 <- r2*y + r12*x12 + sqrt(1-r2**2-r12**2)*rnorm(n);
		row <- pcor(x1,x2,y,spearman=TRUE);
		if(is.null(out)) { out <- row; } else { out[i,] <- row; }
	}
	return(out);
}

# is x1 more correlated with y than is x2? Test if cor(normalized(x1)-normalized(x2),y) > 0.
# Reported in ppart
# To test if a partial correlation is significant use pcor?
cortest2 <- function(x1, x2, y, spearman=FALSE) {
	if (spearman) { x1 <- rank(x1); x2 <- rank(x2); y <- rank(y); }
	r1 <- cor(x1,y);
	r2 <- cor(x2,y);
	z1 <- .5*(log(1+r1)-log(1-r1))*sqrt(length(x1)-3);
	z2 <- .5*(log(1+r2)-log(1-r2))*sqrt(length(x1)-3);
	z12 <- (z1-z2)/sqrt(2); # variances add - just use sqrt(2) because they are both the same
	df <- length(x1)-2;
	dx1 <- (x1-mean(x1))/sd(x1);
	dx2 <- (x2-mean(x2))/sd(x2);
	dy <- (y-mean(y))/sd(y);
	chisq1 <- sum( (dy-r1*dx1)**2 );
	chisq2 <- sum( (dy-r2*dx2)**2 );
	f <- chisq1/chisq2;
	cortestP <- cor.test(dx1-dx2,dy);
	return(data.frame(r1=r1,r2=r2,r12=cor(rank(x1),rank(x2)),z1=z1,z2=z2,z12=z12,
			p1=2*pnorm(abs(z1),lower=FALSE),
			p2=2*pnorm(abs(z2),lower=FALSE),
			p12=2*pnorm(abs(z12),lower=FALSE),
			df=df,
			f=f,pf=2*pf(ifelse(f>1,f,1/f),df1=df,df2=df,lower=FALSE),
			rpart=cortestP$estimate,ppart=cortestP$p.value));
}

# p12 and ppart both work fine if r12=0
# With high r12, however, p12 is *very* conservative, e.g. try
# ctest2 <- cortest2Test(r1=.175,r2=.175,r12=0.9,n=748,tries=1000)
cortest2Test <- function(n=100, tries=100, r1=.3, r2=.3, r12=0) {
	out <- NULL;
	for (i in 1:tries) {
		y <- rnorm(n);
		x12 <- rnorm(n);
		x1 <- r1*y + r12*x12 + sqrt(1-r1**2-r12**2)*rnorm(n);
		x2 <- r2*y + r12*x12 + sqrt(1-r2**2-r12**2)*rnorm(n);
		row <- cortest2(x1,x2,y,spearman=TRUE);
		if(is.null(out)) { out <- row; } else { out[i,] <- row; }
	}
	return(out);
}

# warning: do NOT use with binary values
spearmans <- function(data,use="pairwise.complete.obs") {
	return(cor(apply(data,2,rank),use=use));
}

corNoNa <- function(x,y, spearman=FALSE) {
	b <- !(is.na(x) | is.na(y));
	if(countif(b) < 2) { return(NA); }
	if (spearman) {
		return(cor( rank( select(x,b) ), rank( select(y,b) ) ));
	}
	#else
	return(cor( select(x,b), select(y,b) ));
}

corNoNaSelect <- function(x,y, cond, spearman=FALSE) {
	b <- cond & !is.na(x) & !is.na(y);
	if (countif(b) < 2) { return(NA); }
	if (spearman) { return( cor(rank(select(x,b)),rank(select(y,b))) ); }
	return (cor(select(x,b),select(y,b)));
}

pairCor <- function(dframe, spearman=FALSE) {
	nVar <- length(names(dframe));
	nVar2 <- nVar * (nVar-1)/2;
	nF <- as.factor(names(dframe));
	out <- data.frame(name1=rep(nF[1],nVar2), name2=rep(nF[1],nVar2), cor=rep(0,nVar2),
				i=rep(0,nVar2), j=rep(0,nVar2));
	square <- squarify(data.frame(nF=nF));
	out$i <- square$i;
	out$j <- square$j;
	out$name1 <- nF[out$i];
	out$name2 <- nF[out$j];
	for (k in 1:length(out$i)) {
		out$cor[k] <- corNoNa(dframe[,out$i[k]],dframe[,out$j[k]],spearman=spearman);
	}
	return(out);
}

squarify <- function(data) {
	nVar <- length(data[,1]);
	nVar2 <- nVar * (nVar-1)/2;
	out <- list();
	out$i = rep(0,nVar2);
	out$j = rep(0,nVar2);
	for (name in names(data)) {
		out[[paste(name,".1",sep="")]] = rep(data[1,name], nVar2);
		out[[paste(name,".2",sep="")]] = rep(data[1,name], nVar2);
	}
	out <- data.frame(out);
	index <- 1;
	for (i in 1:(nVar-1)) {
		out$i[index:(index+nVar-1-i)] = i;
		out$j[index:(index+nVar-1-i)] = (i+1):(i+1+nVar-1-i); # == (i+1):nVar
		index <- index + nVar-i;
	}
	for (name in names(data)) {
		out[,paste(name,".1",sep="")] = data[out$i,name];
		out[,paste(name,".2",sep="")] = data[out$j,name];
	}
	return(out);
}

countSimilarity <- function(square,names1,names2) {
	code <- rep(0, length(square[,1]));
	for (i in 1:length(names1)) {
		code <- ifelse(code == (i-1) & square[,names1[i]]==square[,names2[i]], code+1,code);
	}
	return(code);
}

countSimilarPairs <- function(x,names) {
	n <- length(x[,1]);
	code <- rep(0, n**2);
	for (z in 1:length(names)) {
		name <- names[z];
		for (i in 1:n) {
			index <- (1:n)+(i-1)*n;
			code[index] <- code[index] + ifelse(code[index] == (z-1) & x[i,name] == x[,name], 1,0);
		}
	}
	return(code);
}

# log(likelihood ratio) of normal deviates x and y given correlation r1 over correlation r2
# works on vectors of x,y too
logpRatio <- function(x,y,r1,r2=0) {
	return ( -.5*(x**2-y**2-2*r1*x*y)/(1-r1**2)+.5*(x**2-y**2-2*r2*x*y)/(1-r2**2)
		- .5*log((1-r1**2)/(1-r2**2)) );
}

logpRatioUnnorm <- function(x,y,r1,r2=0) {
	xNa <- select(x,!is.na(x));
	yNa <- select(y,!is.na(y));
	x1 <- (x-mean(xNa))/sd(xNa);
	y1 <- (y-mean(yNa))/sd(yNa);
	return(logpRatio(x1,y1,r1,r2));
}

### Common I/O
writeDelim <- function(table,file,report=FALSE,...) {
	write.table(table,file,sep="\t",quote=FALSE,row.names=FALSE,...);
	if(report) cat("Wrote",nrow(table),"rows","to",file,"\n")
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


defaultBreaks <- function(x,n=NULL) {
	if (is.null(n)) { n <- max(1,round(sqrt(length(x)) )); }
	return( min(x) + (0:n)*(max(x)-min(x))/n );
}

linePairs <- segments;

function(x1,y1,x2,y2,...) {
	x1 <- x1 + 0*(x1+y1+x2+y2);
	y1 <- y1 + 0*(x1+y1+x2+y2);
	x2 <- x2 + 0*(x1+y1+x2+y2);
	y2 <- y2 + 0*(x1+y1+x2+y2);
	i <- 0:(3*length(x1)-1);
	index <- 1+ifelse(i %% 3 == 0, i %/% 3, ifelse(i %% 3 == 1, length(x1)+(i%/%3), NA));
	lines(c(x1,x2)[index],c(y1,y2)[index],...);
}

sortedLines <- function(x, y, plot=FALSE, ...) {
	o <- order(x);
	if (plot) plot(x[o], y[o], type="l", ...) else lines(x[o], y[o], ...);
}

CompareDensitySelects <- function(values, list, # of criteria (1 if included in the list)
	labels=names(list),col=1:length(names(list)),lty=1:length(names(list)),
	xlim=range(values),
	plot=TRUE,
	xlab="Values",ylab="Density",xl=xlim[1],adjust=1,showLegend=TRUE,...) {
	dlist <- list();
	maxd <- 0;
	for (name in names(list)) {
		dlist[[name]] <- density(pmin(xlim[2],pmax(xlim[1],select(values,list[[name]]))),
						from=xlim[1],to=xlim[2],adjust=adjust);
		maxd <- max(maxd,dlist[[name]]$y); # not parallel
	}
	if (plot) {
		plot(1e-6,1e-6,col=0,xlim=xlim,ylim=c(1e-4,maxd),xlab=xlab,ylab=ylab,...);
		for (i in 1:length(names(list))) {
			name <- names(list)[i];
			lines(dlist[[name]]$x,dlist[[name]]$y,col=col[i],lty=lty[i]);
		}
		if(showLegend) legend(xl,maxd,labels,col=col,lty=lty);
	}
	else { return(dlist); }
}

CompareDensityRatio <- function(values, truefalse, plot=TRUE, xlim=range(values),
		xlab="Values", ylab="pTrue", ylim=c(0,1), expand=0,adjust=1,...) {
	d1 <- density(pmin(xlim[2],pmax(xlim[1],select(values,truefalse))),
			from=xlim[1]-expand,to=xlim[2]+expand,adjust=adjust);
	d2 <- density(pmin(xlim[2],pmax(xlim[1],select(values,!truefalse))),
			from=xlim[1]-expand,to=xlim[2]+expand,adjust=adjust);
	if(plot) {
		plot(d1$x,d1$y/(d1$y+d2$y),type="l",xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,...);
		hline(.5);
	} else {
		return(data.frame(x=d1$x,dT=d1$y,dF=d2$y,pTrue=d1$y/(d1$y+d2$y)));
	}
	
}

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

# like CompareHistogramList but shows counts
CompareCountsList <- function(list,breaks,
				xlab="Range",ylab="Count",
				legendX=mean(breaks),legendY=-1,
				labels=names(list),
				log="",
				cumulative = "no", # use "greater" or "less"
				showCounts = FALSE, # add #entries in parens if yes
				xlim=range(breaks),...) {
	maxlen <- 1;
	minB <- min(breaks);
	maxB <- max(breaks);
	for (name in names(list)) {
		x <- list[[name]];
		if (min(x) < minB) { cat("Warning:",name,"truncated on left\n"); }
		if (max(x) > maxB) { cat("Warning:",name,"truncated on right\n"); }
		maxlen <- max(maxlen, length(x));
	}
	if(legendY<0) { legendY <- maxlen; }
	logy <- log=="y" || log=="xy";
	plot(c(minB,maxB),c(ifelse(logy,.5,0),maxlen),pch=0,col=0,
			xlab=xlab,ylab=ylab,xlim=xlim,log=log,...);
	for (i in 1:length(names(list))) {
		name <- names(list)[i];
		x <- pmin(pmax(list[[name]],minB),maxB);
		hx <- hist(x,breaks=breaks,plot=FALSE);
		counts <- hx$counts;
		if (cumulative == "greater") { counts <- sum(counts)-cumsum(counts); }
		if (cumulative == "less") { counts <- cumsum(counts); }
		if(logy) { counts <- ifelse(counts==0,.5,counts); }
		points(hx$mids,counts,pch=i,col=i);
		lines(hx$mids,counts,lty=i,col=i);
		if (showCounts) labels[i] <- paste(labels[i]," (n=", countif(!is.na(x)), ")", sep="");
	}
	legend(legendX,legendY,labels,
		pch=1:length(names(list)),
		col=1:length(names(list)),
		lty=1:length(names(list)));
}

# note -- fails if any of the breaks are empty or overlap
# note -- as currently writen uses 95% confidence intervals not 90% confidence intervals
CutOddsRatios <- function(boolean, by, n=6, breaks=quantile(by,seq(0,1,1/n)), col=1, lty=1, lwd=1, plot=TRUE,
				xlab="", ylab="Odds Ratio",
				angle=60, # use angle=90 for vertical bars on horizontals showing the bin
				error.bars=TRUE, error.color="darkgrey", ...) {
	cutBy <- cut(by,breaks,label=FALSE,include.lowest=TRUE);
	odds <- rep(0,length(breaks)-1);
	boundD <- 0*odds;
	boundU <- 0*boundD;
	for (i in 1:length(odds)) {
		f <- fisher.test(table(boolean,cutBy==i));
		odds[i] <- f$estimate;
		boundD[i] <- f$conf.int[1];
		boundU[i] <- f$conf.int[2];
	}
	if(plot) {
		plot(breaks,c(1,odds),pch=-1,col=0,xlab=xlab,ylab=ylab,...);
		# horizontal bars, showing width of bin
		arrows(breaks[1:(length(breaks)-1)],odds,breaks[2:length(breaks)],odds,
			col=col, lty=lty, lwd=lwd, angle=angle, code=3, length=1/12);
		if (error.bars) {
			# vertical bar showing error
			arrows( (breaks[1:(length(breaks)-1)]+breaks[2:length(breaks)])/2,
					boundD,
					(breaks[1:(length(breaks)-1)]+breaks[2:length(breaks)])/2,
					boundU, angle=90, code=3, length=1/12, col=error.color);
		}
	}
	return(odds);
}

CutPercentages <- function(boolean, by, n=6, breaks=unique(quantile(by,seq(0,1,1/n))),
				col=1, lty=1, lwd=1, plot=TRUE,
				angle=60, # use 90 for vertical bars on ends of ranges of bins
				xlab="", ylab="Proportion", error.bars=TRUE, error.color="darkgrey",
				conf.level=.9, ...) {
	cutBy <- cut(by,breaks,label=FALSE,include.lowest=TRUE);
	percent <- rep(0,length(breaks)-1);
	percentL <- 0*percent;
	percentH <- 0*percent;
	for (i in 1:length(percent)) {
		b <- binom.test(countif(cutBy==i & boolean), countif(cutBy==i));
		percent[i] <- b$estimate;
		percentL[i] <- b$conf.int[1];
		percentH[i] <- b$conf.int[2];
	}
	if (plot) {
		plot(range(breaks),0:1,
			pch=-1,col=0,xlab=xlab,ylab=ylab,...);
		# horizontal bars, showing width of bin
		arrows(breaks[1:(length(breaks)-1)],percent,breaks[2:length(breaks)],percent,
			col=col, lty=lty, lwd=lwd, angle=angle, code=3, length=1/12);
		if (error.bars) {
			# vertical bar showing error
			arrows( (breaks[1:(length(breaks)-1)]+breaks[2:length(breaks)])/2,
					percentL,
					(breaks[1:(length(breaks)-1)]+breaks[2:length(breaks)])/2,
					percentH, angle=90, code=3, length=1/12, col=error.color);
		}
	}
	return(percent);
}

CutAverages <- function(real, by, n=6, breaks=unique(quantile(by,seq(0,1,1/n))),
				col=1, lty=1, lwd=1, plot=TRUE,
				angle=60, # use 90 for vertical bars on ends of ranges of bins
				xlab="", ylab="Mean", error.bars=TRUE, error.color="darkgrey",
				conf.level=.9, ...) {
	cutBy <- cut(by,breaks,label=FALSE,include.lowest=TRUE);
	avg <- rep(0,length(breaks)-1);
	avgL <- 0*avg;
	avgH <- 0*avg;
	for (i in 1:length(avg)) {
		b <- t.test(real[cutBy==i]);
		avg[i] <- b$estimate;
		avgL[i] <- b$conf.int[1];
		avgH[i] <- b$conf.int[2];
	}
	if (plot) {
		plot(range(breaks),range(real),
			pch=-1,col=0,xlab=xlab,ylab=ylab,...);
		# horizontal bars, showing width of bin
		arrows(breaks[1:(length(breaks)-1)],avg,breaks[2:length(breaks)],avg,
			col=col, lty=lty, lwd=lwd, angle=angle, code=3, length=1/12);
		if (error.bars) {
			# vertical bar showing error
			arrows( (breaks[1:(length(breaks)-1)]+breaks[2:length(breaks)])/2,
					avgL,
					(breaks[1:(length(breaks)-1)]+breaks[2:length(breaks)])/2,
					avgH, angle=90, code=3, length=1/12, col=error.color);
		}
	}
	return(avg);
}

CompareCumulativeFractions <- function(list,
				xlab="Range",ylab="Cumulative Fraction",
				labels=names(list),
				log="", minY=ifelse(log=="y" || log=="xy",0.1, 0),
				maxlen=max(sapply(list,length)),
				xlim=range(withoutNA(sapply(list,c))),
				ylim=c(minY/maxlen,1),
				legend=TRUE,
				legendX=if(log=="x" || log=="xy") exp(mean(log(xlim))) else mean(xlim),
				legendY=max(ylim),
				cex=1,
				col=1:length(list),
				lty=1:length(list),
				pch=1:length(list),
				lwd=NULL,
				lesser=TRUE,
				new=TRUE,
				...) {
	if(new) plot(c(1/maxlen,1),c(1/maxlen,1),pch=0,col=0,lwd=1,
			xlab=xlab,ylab=ylab,log=log,xlim=xlim,ylim=ylim,...);
	for (i in 1:length(names(list))) {
		name <- names(list)[i];
		x <- sort(withoutNA(list[[name]]));
		if(length(x)>0) {
			y <- if(lesser) (1:length(x))/length(x) else (length(x):1)/length(x);
			points(x,y,pch=pch[i],col=col[i]);
			lines(x,y,lty=lty[i],col=col[i],lwd=if(is.null(lwd)) 1 else lwd[i]);
		}
	}
	if (legend) legend(legendX,legendY,labels,cex=cex,
					pch=pch, col=col, lty=lty,
					lwd=if(is.null(lwd)) 1 else lwd);
}

CompareHistogramList <- function(list,breaks,
				cumulative="no",
				xlab="Range",
				ylab=ifelse(cumulative=="no","Proportion","Cumulative Proportion"),
				labels=names(list),
				log="", minY=ifelse(log=="y" || log=="xy",0.1, 0),
				maxlen=max(sapply(list,length)),
				ylim=c(minY/maxlen,1),
				xlim=range(breaks),
				legend=TRUE,
				legendX=if(log=="x" || log=="xy") exp(mean(log(xlim))) else mean(xlim),
				returnP=FALSE,plot=TRUE,warn=TRUE,
				legendY=max(ylim),cex=1,
				col=1:length(list),
				lty=1:length(list),
				pch=1:length(list),
				lwd=NULL,
				showCounts = FALSE, # add #entries in parens if yes
				...) {
	minB <- min(breaks);
	maxB <- max(breaks);
	for (name in names(list)) {
		x <- list[[name]];
		x <- x[!is.na(x)];
		if (length(x) > 0 && warn) {
			if (min(x) < minB) { cat("Warning:",name,"truncated on left\n"); }
			if (max(x) > maxB) { cat("Warning:",name,"truncated on right\n"); }
		}
	}
	if (plot) plot(c(1/maxlen,1),c(1/maxlen,1),pch=0,col=0,lwd=1,
			xlab=xlab,ylab=ylab,log=log,xlim=xlim,ylim=ylim,...);
	outlist <- list();
	for (i in 1:length(names(list))) {
		name <- names(list)[i];
		x = list[[name]];
		x <- x[!is.na(x)];
		if (length(x) > 0) {
			x <- pmin(pmax(x,minB),maxB);
			hx <- hist(x,breaks=breaks,plot=FALSE);
			hxAdj <- ifelse(hx$counts==0,minY,hx$counts)/length(x);
			if(sum(hxAdj)>0) {
					if (cumulative == "greater") { hxAdj <- 1 - cumsum(hxAdj)/sum(hxAdj); }
				if (cumulative == "less") { hxAdj <- cumsum(hxAdj)/sum(hxAdj); }
			}
			if (plot) {
				if(!is.null(pch) && !is.na(pch[i]) && pch[i] > 0)
					points(hx$mids,hxAdj,pch=pch[i],col=col[i]);
				lines(hx$mids,hxAdj,lty=lty[i],col=col[i],lwd=if(is.null(lwd)) 1 else lwd[i]);
			}
			if (returnP) {
				if (i==1) outlist[["mids"]] <- hx$mids;
				outlist[[name]] <- hxAdj;
			}
			if (showCounts) labels[i] <- paste(labels[i]," (n=", length(x), ")", sep="");
		}
	}
	if (plot && legend) legend(legendX,legendY,labels,cex=cex,
					pch=pch, col=col, lty=lty,
					lwd=if(is.null(lwd)) 1 else lwd);
	if(returnP) { return(data.frame(outlist)); }
}

CompareTimeSeries <- function(x,list,
				xlab="Time",
				ylab="",
				labels=names(list),
				log="", minY=ifelse(log=="y" || log=="xy",0.1, 0),
				ylim=range(c(list,recursive=TRUE),na.rm = TRUE),
				xlim=range(x),
				legend=TRUE,
				legendX=min(xlim),
				legendY=if(log=="y" || log=="xu") exp(mean(log(ylim))) else mean(ylim),
				cex=1,
				col=1:length(list),
				lty=1:length(list),
				pch=1:length(list),
				lwd=NULL,
				warn=TRUE,
				...) {
	for (name in names(list)) {
		y <- list[[name]];
		y <- y[!is.na(y)];
		if (length(y) > 0 && warn) {
			if (min(y) < ylim[1]) { cat("Warning:",name,"truncated on left\n"); }
			if (max(y) > ylim[2]) { cat("Warning:",name,"truncated on right\n"); }
		}
	}
	plot(xlim,ylim,pch=0,col=0,lwd=1,
		xlab=xlab,ylab=ylab,log=log,xlim=xlim,ylim=ylim,...);
	for (i in 1:length(names(list))) {
		name <- names(list)[i];
		y = list[[name]];
		x2=x[!is.na(y)];
		y2=y[!is.na(y)];
		if(!is.null(pch) && !is.na(pch[i]) && pch[i] > 0)
			points(x2,y2,pch=pch[i],col=col[i]);
		lines(x2,y2,lty=lty[i],col=col[i],lwd=if(is.null(lwd)) 1 else lwd[i]);
	}
	if (legend) legend(legendX,legendY,labels,cex=cex,
					pch=pch, col=col, lty=lty,
					lwd=if(is.null(lwd)) 1 else lwd);
}


CompareHistograms <- function(x,y,breaks,
				distr1="x",distr2="y",
				xlab="Range",ylab="Proportion",
				legendX=mean(breaks),legendY=0.9,...) {
	minB <- min(breaks);
	maxB <- max(breaks);
	if (min(x) < minB) { cat("Warning: first argument truncated on left\n"); }
	if (max(x) > maxB) { cat("Warning: first argument truncated on right\n"); }
	if (min(y) < minB) { cat("Warning: second argument truncated on left\n"); }
	if (max(y) > maxB) { cat("Warning: second argument truncated on right\n"); }
	x <- pmin(pmax(x,minB),maxB);
	y <- pmin(pmax(y,minB),maxB);
	hx <- hist(x,breaks=breaks,plot=FALSE);
	hy <- hist(y,breaks=breaks,plot=FALSE);
	plot(hx$mids,hx$counts/length(x),
		xlim=range(breaks),
		ylim=c(0,max(max(hx$counts)/length(x),max(hy$counts)/length(y))),
		xlab=xlab,ylab=ylab,...);
	lines(hx$mids,hx$counts/length(x));
	points(hy$mids,hy$counts/length(y),pch=2);
	lines(hy$mids,hy$counts/length(y),lty=2);
	legend(legendX,legendY,c(distr1,distr2),pch=c(1,2),lty=c(1,2));
}

CompareLogHistograms <- function(x,y,breaks,
				distr1="x",distr2="y",
				xlab="Range",ylab="Proportion",
				legendX=mean(breaks),legendY=0.9,...) {
	minB <- min(breaks);
	maxB <- max(breaks);
	if (min(x) < minB) { cat("Warning: first argument truncated on left\n"); }
	if (max(x) > maxB) { cat("Warning: first argument truncated on right\n"); }
	if (min(y) < minB) { cat("Warning: second argument truncated on left\n"); }
	if (max(y) > maxB) { cat("Warning: second argument truncated on right\n"); }
	x <- pmin(pmax(x,minB),maxB);
	y <- pmin(pmax(y,minB),maxB);
	hx <- hist(x,breaks=breaks,plot=FALSE);
	hy <- hist(y,breaks=breaks,plot=FALSE);
	plot(hx$mids,ifelse(hx$counts==0,.1,hx$counts)/length(x),
		xlim=range(breaks),
		ylim=c(.1/max(length(x),length(y)),1),
		xlab=xlab,ylab=ylab,log="y",...);
	lines(hx$mids,ifelse(hx$counts==0,.1,hx$counts)/length(x));
	points(hy$mids,ifelse(hy$counts==0,.1,hy$counts)/length(y),pch=2);
	lines(hy$mids,ifelse(hy$counts==0,.1,hy$counts)/length(y),lty=2);
	legend(legendX,legendY,c(distr1,distr2),pch=c(1,2),lty=c(1,2));
}

densityPlot <- function(x, y, nx=NULL, ny=NULL,
			xbreaks=defaultBreaks(x,n=nx), ybreaks=defaultBreaks(y,n=ny),
			nz=30, logZ=FALSE, byBreaks=FALSE,
			xlab="",ylab="") {
	nXY <- xtabs(~x+y, data.frame(x=findInterval(x,xbreaks), y=findInterval(y,ybreaks)));
	rectX <- row(nXY);
	rectY <- col(nXY);
	maxXY <- max(nXY);
	if (logZ) {
		intensity <- log(1+nXY)/log(1+maxXY);
		colors <- grey(1 - intensity);
	} else {
		intensity <- 1 + round(nz*nXY/maxXY);
		colors <- ifelse(nXY==0, rgb(1,1,1), grey(1 - intensity/(nz+1)));
	}
	if (byBreaks) {
		plot(range(row(nXY)),range(col(nXY)),col=0,xlab=xlab,ylab=ylab);
		rect(c(row(nXY)),c(col(nXY)),
			c(row(nXY)+1),c(col(nXY)+1),
		col=c(colors), lty=0 );
	} else {
		plot(range(x),range(y),col=0,xlab=xlab,ylab=ylab);
		rect(xbreaks[c(row(nXY))],ybreaks[c(col(nXY))],
			xbreaks[c(row(nXY)+1)],ybreaks[c(col(nXY)+1)],
		col=c(colors), lty=0 );
	}
	return(nXY);
}

# shows extremes of range only if extrema=TRUE
# never shows individual outliers
# shows 90% confidence interval of medians with a bar instead of with a notch
MediansPlot <- function(groups, range=1.5, extrema=FALSE,
			conf.int=TRUE,
			conf.col="darkgrey", conf.border=conf.col, conf.lty=1, conf.lwd=1,
			conf.width=1/6, box=TRUE, bar=FALSE, ymin=NULL, names=names(groups), ...) {
	n <- length(groups);
	stats <- matrix(0, nrow = 5, ncol = n);
	conf <- matrix(0, nrow = 2, ncol = n);
	count <- rep(0,n);

	for (i in 1:n) {
		statsi <- boxplot.stats(groups[[i]],coef=range);
		if (!extrema) {
			statsi$stats[1] <- statsi$stats[2];
			statsi$stats[5] <- statsi$stats[4];
		}
		stats[,i] <- statsi$stats;
		conf[,i] <- statsi$conf;
		count[i] <- statsi$n;
		if(!is.null(ymin)) conf[1,] <- pmax(ymin,conf[1,]);
		#cat(names(groups)[i],"median",statsi$stats[3],"\n");
	}
	if (bar) {
		barplot2(stats[3,], names=names,
				ci.l=if(conf.int) conf[1,] else NULL,
				ci.u=if(conf.int) conf[2,] else NULL,
				#ci.color=conf.col, ci.lty=conf.border, ci.lwd=conf.lwd,
				plot.ci=conf.int,
				col="darkgrey",
				...);
	} else {
		if (box) {
			bxp(list(stats=stats,n=count,conf=conf,out=c(),group=c(),names=names),notch=FALSE,
				 ...);
		} else {
			plot(1:n, stats[3,], col="white", ...); # invisible
		}
		if(conf.int) {
			rect((1:n)-conf.width,conf[1,],(1:n)+conf.width,conf[2,],
				col=conf.col,border=conf.border,lty=conf.lty,lwd=conf.lwd);
		}
		if (!box) {
			segments((1:n)-conf.width, stats[3,], (1:n)+conf.width, stats[3,]);
		}
	}
}

# barplot with vertical error bars
# groups should be a vector of TRUE and FALSE for each condition (e.g. from split(true/false, factor))
ProportionPlot <- function(groups,
			conf.int=TRUE,conf.level=0.9,
			conf.col="darkgrey", conf.lty=1, conf.lwd=1,
			conf.width=1/6, ylim=0:1, ylab="Proportion",
			plot=TRUE,
			labels=names(groups),
			order=1:length(groups),
			# optionally, put stars between columns that are significantly different
			# according to the fisher exact test
			compare=NULL, # which columns to compare to the next column, in ascending order
			compare.levels=0.05,
			compare.strings=c("*"),
			compare.cex=1.5,
			compare.col="black",
			...) {
	n <- length(groups);
	percent <- rep(0,n);
	percentL <- 0*percent;
	percentH <- 0*percent;
	for (i in 1:n) {
		b <- binom.test(countif(groups[[i]]), length(groups[[i]]), conf.level=conf.level);
		percent[i] <- b$estimate;
		percentL[i] <- b$conf.int[1];
		percentH[i] <- b$conf.int[2];
	}
	names(percent) <- labels;
	if (plot) barplot(percent[order], ylim=ylim, ylab=ylab, ...);

	if(plot && conf.int) {
		arrows( 1.2*(1:n) - .5, percentL[order],
			(1:n)*1.2 - .5, percentH[order], angle=90, code=3, length=conf.width,
				col=conf.col, lty=conf.lty, lwd=conf.lwd);
	}
	if (plot && !is.null(compare)) {
		for (i in compare) {
			p <- fisher.test(matrix(c(countif(groups[[i]]),countif(!groups[[i]]),
							countif(groups[[i+1]]),countif(!groups[[i+1]])),ncol=2))$p.value;
			string <- NULL;
			for (j in 1:length(compare.levels)) {
				if(is.null(string) & p <= compare.levels[j]) string <- compare.strings[j];
			}
			if (!is.null(string)) {
				text(i*1.2+0.1,(countPerIf(groups[[i]])+countPerIf(groups[[i+1]]))/2,
					string, cex=compare.cex, col=compare.col);
			}
		}
	}
	return(data.frame(percent=percent,low=percentL,high=percentH,row.names=names(percent)));
}

### Receiver Operating Curves

# a more accurate variant
makeaoc <- function(s,d=NULL,bigger=TRUE) {
	if (is.null(d)) { d <- s$"FALSE"; s <- s$"TRUE"; }
	if (!bigger) { s <- -s; d <- -d; }
	n <- 0*s;
	for(i in 1:length(s)) { n[i] <- countPerIf(s[i]>d) + .5*countPerIf(s[i]==d); };
	return(mean(n));
}

# default is bigger is better
AOC <- function(valuesT,valuesF=NULL, truth=NULL, bigger=TRUE) {
	if (is.null(valuesF)) {
		valuesF <- select(valuesT,!truth);
		valuesT <- select(valuesT,truth);
	}
	if(!bigger) {
		valuesT <- -valuesT;
		valuesF <- -valuesF;
	}
	cutoff <- sort(unique(valuesT));
	cutoff <- c(min(valuesT,valuesF)-1, cutoff, max(valuesT,valuesF)+1);
	sens <- spec <- vector(length=length(cutoff));
	for (i in 1:length(cutoff)) {
		sens[i] <- countif(valuesT >= cutoff[i]);
		spec[i] <- countif(valuesF < cutoff[i]);
	}
	sens <- sens / length(valuesT);
	spec <- spec / length(valuesF);
	acc <- (sens+spec)/2;
	AOC <- sum( (sens[1:(length(sens)-1)]-sens[2:length(sens)]) *
			(spec[1:(length(spec)-1)]+spec[2:length(spec)])/2 );
	if(!bigger) {
		cutoff <- -cutoff;
	}
	return(list(sens=sens,spec=spec,acc=acc,cutoff=cutoff,AOC=AOC));
}



# assumes higher is better
cutoffs <- function(scoresTest, scoresControl) {
	sortT <- sort(scoresTest);
	sortC <- sort(scoresControl);
	cutoffs <- sortT;
	pRatios <- 0*(1:length(sortT));
	pRatios5 <- 0*(1:length(sortT));

	for(i in 1:length(sortT)) {
		pT <- sum(ifelse(sortT >= cutoffs[i],1,0));
		pC <- sum(ifelse(sortC >= cutoffs[i],1,0));
		pRatios5[i] <- (length(pT)/length(pC)) * (pC+.5)/(pT+.5);
		pRatios[i] <- (length(pT)/length(pC)) * pC/pT;
	}
	return(data.frame(cutoff=cutoffs,pRatio5=pRatios5,pRatio=pRatios));
}

# The really smart thing is probably to do some smoothing and put a cutoff on the
# local ratio of p...
chooseCutoff <- function(scoresTest,scoresControl,p=0.05) {
	cutoffX <- cutoffs(scoresTest,scoresControl);
	indexes <- split(1:length(cutoffX$cutoff),cutoffX$pRatio5 <= p)$"TRUE";
	if (length(indexes) > 0) {
		index <- indexes[1];
	} else {
		index <- length(scoresTest);
	}
	return(data.frame(cutoff=cutoffX$cutoff[index],
		sensitivity=(length(scoresTest)-index+1)/length(scoresTest),
		specificity=sum(ifelse(scoresControl < cutoffX$cutoff[index],1,0))/length(scoresControl)
			));
}

chooseNaiveCutoff <- function(scoresTest,scoresControl,p=0.01) {
	n <- round( length(scoresControl)*(1-p) ); # n = length(scoresControl) OK
	cutoffX <- (sort(scoresControl))[n];
	return(data.frame(cutoff=cutoffX,
			sensitivity=sum(ifelse(scoresTest>cutoffX,1,0))/length(scoresTest),
			specificity=n/length(scoresControl)));
}

### Operon Prediction

# Run operon predictions for a taxon
RunOperons = function(taxon = tail(commandArgs(),1), gnc = "gnc", suffix=".scores") {
	if(is.null(taxon)) {
		stop("No taxon specified in RunOperons");
	}
	file = paste(gnc,taxon,suffix,sep="");
	pairs = read.delim(file);
	if(is.null(pairs) || nrow(pairs)==0) {
		stop("Empty input file: ", file);
	}
	preds = OperonsPredict(pairs);
	cat("Taxon",taxon,"fOperon",preds$fOperon,"Threshold",preds$threshOp,
		"ActualOperons",sum(preds$preds$bOp),sum(!preds$preds$bOp),"\n");
	writeDelim(preds$preds[,c("Gene1","Gene2","bOp","pOp","Sep","MOGScore","GOScore","COGSim","ExprSim")],file=paste(gnc,taxon,".pred",sep=""),report=TRUE);
	#writeDelim(preds$pairs,file=paste(gnc,taxon,".pairs",sep=""),report=TRUE);
}

# The main operon prediction routine takes the output of GeneNeighborScores.pl, as
# a data frame, and returns a data frame for the same-strand pairs called preds
# and for all non-CRISPR pairs called pairs
OperonsPredict <- function(pairs,
				vars=c("MOGScore","GOScore","ExprSim"),
				useCOGSim = TRUE,
				debug=FALSE,
				fOperon = OperonFraction(pairs,debug=debug),
				smooth=TRUE,
				plot=debug) {
	pairs$same = (pairs$Code=="Same");
	# ignore CRISPR-related pairs for now -- will add them back below
	pairs2 = pairs[(!pairs$Type1 %in% c(9,10)) & (!pairs$Type2 %in% c(9,10)),];
	varsFinal = c();
	for (v in vars) {
		if (sum(!is.na(pairs2[[v]])) == 0) {
			if(debug) cat("Skipping",v,"\n");
		} else {
			# Use span=1 so it doesn't blow up on the huge GNScore=0 group
			# Makes sense anyway as we expect roughly linear fits...
			if(debug) cat("BinnedBinaryFit",v,"\n");
			bbf <- BinnedBinaryFit(pairs2[[v]], pairs2$same,smooth=smooth,plot=plot,xlab=v,span=1,debug=debug);
			if(debug) { cat(v,range(bbf$p),"\n"); }
			pairs2[[paste("bbf",v,sep="")]] <- log(bbf$p/(1-bbf$p));
			varsFinal = c(varsFinal,v);
		}
	}
	if (useCOGSim && all(is.na(pairs2$COGSim))) {
		cat("Turning useCOGSim off -- no assignments\n");
		useCOGSim = FALSE;
	}
	if (useCOGSim) {
		if(debug) cat("Computing ClassFit, COGsimClass\n");
		COGsim = unclass(pairs2$COGSim);
		pCOG <- ClassFit(COGsim,pairs2$same);
		pairs2$cfCOG <- log(pCOG$p[COGsim])/(1-pCOG$p[COGsim]);
		if(debug) cat("ClassFit done\n");
	}

	bbfvars <- sapply(varsFinal, function(x) { return(paste("bbf",x,sep="")); } );
	if (useCOGSim) { bbfvars <- c(bbfvars,"cfCOG"); }
	formula = eval(parse(text=paste("same","~",paste(bbfvars,collapse="+"),sep="")));
	gcm <- glm(formula,family=binomial(logit),data=pairs2);
	pairs2$pSame <- fitted(gcm);
	pCoverS <- countif(!pairs2$same)/countif(pairs2$same);

	# Split the same pairs2 into high and low pSame
	# Because pSame is highly skewed, and has a lot of values around a median,
	# use quantiles to get our cutoff for doing bbfDistance
	preds <- pairs2[pairs2$Code=="Same",];
	threshPSame <- quantile(preds$pSame,fOperon,name=FALSE);

	# % of Same and Convergent pairs2 with pSame above the threshold
	# note we use >threshPSame -- this is important if (as is often the case) there
	# is a huge clump of values in the middle with GNScore=0, GNWithin=0 -- this ensures that we
	# count those as non-operon pairs2
	# Note the range of 0.01-0.99 is a hack to avoid getting NaN values below
	fHighS <- pmin(0.99, pmax(0.01, countPerIf(preds$pSame > threshPSame)));
	fHighSC <- pmax(0.01, countPerIf(pairs2$pSame > threshPSame));
	fHighC <- pmax(0.01, countif(pairs2$pSame > threshPSame & !pairs2$same)/countif(!pairs2$same));

	# Now estimate the %operon pairs in the high and low sets, assuming the Convergent distribution
	# matches the operon distribution.
	# H means High, O means operon, g means given:
	# thresholds .95 and .1 are arbitrary - could use Dirichlet priors but with what alpha/beta??
	# also note that we use fHighS not fHighSC because only same pairs are being considered
	# now for same pairs p(!O|H) = p(H|!O)*p(!O)/p(H) = p(H|C)*p(!O)/p(H), so
	fOgH = 1 - (fHighC*(1-fOperon)/fHighS);
	fOgHUse <- pmin(fOgH,.95);
	# similarly, for same pairs p(O) = p(O|H)*p(H) + p(O|L)*p(L), so
	fOgL = (fOperon - fOgH*fHighS)/(1-fHighS);
	fOgLUse <- pmax(fOgL,0.1);

	# for a given separation D, we need to get observed %high pairs fHigh = a+b*p/(c+d*p),
	# where p = P(D|O)/(P(D|O)+P(D|!O)):
	# fHigh = P(H|D) = P(H|D)/(P(H|D)+P(L|D)) = P(D|H)*P(H)/(top + same for L)
	#       = { P(D|O)*P(O|H)*P(H) + P(D|!O)*P(!O|H)*P(H) }/(top + same for L)
	#       = { p*P(O|H)*P(H) + (1-p)*P(!O|H)*P(H) }/(top + same for L)
	#       = { P(!O|H)*P(H) + p*(P(O|H)-P(!O|H))*P(H) }/(top + same for L)
	# and P(H) = fHighS since only same pairs are relevant here...
	a <- (1-fOgHUse)*fHighS;
	b <- (fOgHUse-(1-fOgHUse))*fHighS;
	c <- a + (1-fOgLUse)*(1-fHighS);
	d <- b + (fOgLUse-(1-fOgLUse))*(1-fHighS);
	if (debug) cat("bbfD a b c d fOperon fHighC fHighS fHighSC", a,b,c,d, fOperon,fHighC,fHighS,fHighSC,"\n");

	# Compute pSameD, based on interpolating P(D|O)/(P(D|O)+P(D|!O)) for each group
	# Each group is based on overlapping bins (given by gSize, gSep) determined in rank-space,
	# and the p for that group is determined via Dirichlet priors and a,b,c,d
	# If our samples were large, uncontaminated, and of equal size (fHighS=.5), then p(D) = nH/(nH+nL)

	bbfD <- BinnedBinaryFit(preds$Sep,preds$pSame>threshPSame,a=a,b=b,c=c,d=d,debug=debug,
				plot=plot,xlab="Gene Separation",log="y",main=name,
				smooth=smooth);
	if(debug) cat("Finished BinnedBinaryFit for Sep\n");
	if(countif(is.na(bbfD$p)) > 0 || length(bbfD$p) != length(preds$Sep)) {
		cat("BinnedBinaryFit failed");
		return(NULL);
	}
	preds$pSameD <- bbfD$p;

	# And now make predictions:
	# P(O|D,X) = P(O^D^X)/{P(O^D^X)+P(!O^D^X)} = P(D|O)*P(X|O)*P(O)/(top + same for !O)
	#          = P(D|O)/P(D|!O) * P(X|O)/P(X|!O) * P(O) / (top + P(!O))
	#   where P(D|O)/P(D|!O) = pD/(1-pD) from the above scheme, and we solve
	#   P(X|S) = P(X|O)*P(O) + P(X|!O)*P(!O), plus P(X|!O)=P(X|C), to get
	#   P(X|O)/P(X|!O) = { P(X|S) - P(X|!O)*P(!O) }/{ P(O) * P(X|!O) }
	#                  = { P(X|S)/P(X|C) - P(!O) }/P(O) = { P(C)/P(S) * pSame/(1-pSame) - P(!O) }/P(O)
	# so altogether we have
	# P(O|D,X) = pD/(1-pD) * ( P(C)/P(S) * pSame/(1-pSame) - P(!O) ) / (top + P(!O))

	pLogisticRatio <- pCoverS*preds$pSame/(1-preds$pSame);
	pOpLogistic <- pmax(0.01, (pLogisticRatio - (1-fOperon))); # arbitrary floor to keep it from going -
	pOpLogistic <- pOpLogistic/(pOpLogistic + 1 - fOperon);

	# a variant where we just scale the odds ratio from the regression
	pOpLogisticM <- pLogisticRatio * fOperon/(1-fOperon);
	pOpLogisticM <- pOpLogisticM/(1+pOpLogisticM);

	# a variant where we assume some of the control pairs do look like operon pairs
	# will also reduce the subtraction problem
	# specifically, assume
	# P(X|!O) = (1-f)*P(X|C) + f*P(X|O)
	# P(X|S) = P(O)*P(X|O) + P(!O)*((1-f)*P(X|C) + f*P(X|O))
	# P(X|O) = { P(X|S) - P(!O)*(1-f)*P(X|C) } / { P(O) + f*P(!O) }
	# P(X|C) = { P(X|!O) - f*P(X|O) }/(1-f)
	# P(O|X)/P(!O|X) = P(X|O)/P(X|!O) = P(O)/P(!O) * { P(X|S) - P(!O)*(1-f)*P(X|C) }/ { P(O) + f*P(!O) }
	#					/ { (1-f)*P(X|C) + f*P(X|O) }
	# = P(O)/P(!O) * { P(X|S) - P(!O)*(1-f)*P(X|C) }/ { P(O) + f*P(!O) } / { (1-f)*P(X|C) + f*top }
	#	set indicator = P(X|S)/P(X|C) - P(!O) * (1-f)
	# = P(O)/P(!O) * indicator / ( 1-f + f*indicator )
	#By default, used fControlLikeOp = 0.1
	#indicator <- pmax(.02, pLogisticRatio - (1-fOperon)*(1-fControlLikeOp));
	#pllTop <- indicator / (1-fControlLikeOp + fControlLikeOp*indicator);
	#pOpLogisticLike <- pllTop/(1+pllTop);

	#pOpDistance <- combineP(preds$pSameD, fOperon);
	preds$pOp <- combineP(pOpLogistic,preds$pSameD);

	# bOp gives the "right" #same-operon pairs. Unclear if this would be better
	# than pOp >= .5 as a cutoff or not
	threshOp <- sort(preds$pOp)[round((1-fOperon)*length(preds$pOp))];

	preds$bOp = preds$pOp > threshOp;

	# Add back the CRISPR pairs
	pairs$index = 1:nrow(pairs);
	preds = merge(pairs[pairs$Code=="Same",],preds,all.x=TRUE);
	preds = preds[order(preds$index),];
	pairs$index=NULL;
	preds$bOp = ifelse(is.na(preds$pOp), preds$Type1 %in% c(9,10) & preds$Type2 %in% c(9,10) & preds$Sep < 50, preds$bOp);
	return(list(preds=preds,pairs=pairs,fOperon=fOperon,threshOp=threshOp));
}

# Given f(p) = (a+b*p)/(c+d*p) and a Dirichlet prior on p of p**(alpha-1)*(1-p)**(beta-1),
# and data nT, nF (true/false) based on f(p), reports the
# maximum likelihood estimate of p
pseudocountLL <- function(nT,nF,a,b,c,d,alpha=2,beta=2,plot=FALSE,debug=FALSE) {
	if(debug) cat("pseudocountLL","nT",nT,"nF",nF,"abcd",a,b,c,d,"alphabeta",alpha,beta,"\n");
	f <- (nT+1)/(nF+nT+2);
	naivep <- (c*f-a)/(b-d*f);
	pStart <- if(naivep > .001 && naivep < .999) naivep else (a+b*.5)/(c+d*.5);
	# the function we will minimize = -log likelihood
	func <- function(q) {
		if (q < -1e-8 || q >= 1 || a+b*q < 1e-8 || c+d*q < 1e-8 || (c-a + (d-b)*q) < 1e-8) {
			# nlm is smart enough to backtrack from illegal values
			v  <- 1e10;
			attr(v,"gradient") <- 10;
			attr(v,"hessian") <- 10;
			if(debug) { cat("Illegal",q,"\n"); }
			return(v);
		}
		v <- nT*log(a+b*q) + nF*log((c-a)+(d-b)*q) - (nT+nF)*log(c+d*q) + 
			(alpha-1)*log(q) + (beta-1)*log(1-q);
		gradient <- nT*b/(a+b*q) + nF*(d-b)/(c-a+(d-b)*q) - (nT+nF)*d/(c+d*q) +
				(alpha-1)/q - (beta-1)/(1-q);
		hessian <- -nT*b**2/(a+b*q)**2 - nF*(d-b)**2/(c-a+(d-b)*q)**2 + (nT+nF)*d**2/(c+d*q)**2 -
				(alpha-1)/q**2 - (beta-1)/(1-q)**2;
		# we want to maximize; nlm minimizes
		v <- -v;
		attr(v,"gradient") <- -gradient;
		attr(v,"hessian") <- -hessian;
		return(v);
	}
	opt <- nlm(func,pStart,hessian=TRUE,check.analyticals=FALSE);
	if (plot) {
		cat("Estimate",opt$estimate,"Start",pStart,"MaxLL",-opt$minimum,
			"Gradient",-opt$gradient,"Hessian",-opt$hessian,
			"Iter",opt$iter,"Code",opt$code,"\n");
		oldpar <- par(no.readonly=TRUE);
		on.exit(par(oldpar));
		par(mfrow=c(1,2));
		p=(1:299)/300
		dll <- nT*b/(a+b*p) + nF*(d-b)/(c-a+(d-b)*p) - (nT+nF)*d/(c+d*p) +
			(alpha-1)/p - (beta-1)/(1-p);
		plot(p,dll,ylim=c(-20,20),type="b");
		vline(naivep,lty=1);
		vline(opt$estimate,lty=1,col="red");
		vline(0);
		vline(1);
		hline(0,lty=1);
		plot(p,(a+b*p)/(c+d*p),type="l",ylim=c(0,1),lwd=2);
		vline(0);vline(1);hline(0);hline(1);vline(naivep,lty=1);hline(f,lty=1);eqline();
		vline(opt$estimate,lty=1,col="red");
	} else {
		return(list(p=opt$estimate,code=opt$code,iter=opt$iter,gradient=-opt$gradient));
	}
}


# Estimate fraction of same-strand pairs that are operon pairs
OperonFraction = function(scores, keeptypes = c(1), minPairs = 100, debug=FALSE) {
	scores2 = scores[scores$Type1 %in% keeptypes & scores$Type2 %in% keeptypes,];
	splits = split(cbind(scores2,index=1:nrow(scores2)),scores2$Scaffold);
	scaffolds = names(splits);
	starts = termini = list();
	allForward = allForwardplus = allSumBias = rep(NA,nrow(scores2));
	for (s in names(splits)) {
		scscores = splits[[s]];
		if (nrow(scscores) >= minPairs) {
			ns = 1:nrow(scscores);
			sumBias = cumsum(ifelse(scscores$Strand1=="+",1,-1));
			terminus = (1:nrow(scscores))[sumBias==max(sumBias)][1];
			if (scscores$IsPartial[1]==1) terminus = nrow(scscores);
			rotatedN = ns;
			# only use large scaffolds for strand bias
			rotatedN = if(terminus == nrow(scscores)) ns else c((terminus+1):nrow(scscores),1:terminus);
			rotatedS = scscores$Strand1[rotatedN];
			sumRotBias = cumsum(ifelse(rotatedS=="+",1,-1));
			start = select(rotatedN,sumRotBias==min(sumRotBias))[1];
			forwardplus = if(terminus >= start) (ns > start & ns <= terminus) else
					(ns > start | ns <= terminus);
			forward = (scscores$Strand1==ifelse(forwardplus,"+","-"));

			starts[[s]] = scscores$index[start];
			termini[[s]] = scscores$index[terminus];
			allForward[scscores$index] = forward;
			allForwardplus[scscores$index] = forwardplus;
			allSumBias[scscores$index] = sumBias;
		}
	}
	if(sum(!is.na(allForward)) == 0) return(0.5); # give up on fragmented genome
	scores2$allForward = allForward;
	scores2 = scores2[!is.na(scores2$allForward),];
	fPlus = sum(scores2$allForward)/nrow(scores2);
	fP2M2 = fPlus**2 + (1-fPlus)**2;
	fSame = sum(ifelse(scores2$Code=="Same",1,0))/nrow(scores2);
	fOperonNaive = (1 - fP2M2/fSame)/(1 - fP2M2);
	pSameLeading = countPerIf(scores2$Code=="Same", scores2$allForward);
	pSameLagging = countPerIf(scores2$Code=="Same",!scores2$allForward);
	if(debug) cat("fPlus",fPlus,"pSameLeading",pSameLeading,"pSameLagging",pSameLagging,"\n")
	return( OperonFreqStrandBias(pSameLeading,pSameLagging) );
}

BinnedBinaryFit2 <- function(xT, xF, gSize=NULL, gSep=NULL, xout=c(xT,xF),
				fT=length(xT)/(length(xT)+length(xF)), # to compute abcd
				pActual=.5,
				a=NULL,b=NULL,c=NULL,d=NULL,
				smooth=TRUE,span=NULL,
				...) {
	x <- c(xT,xF);
	y <- ( 1:length(x) <= length(xT) );
	if (!is.null(a)) { fT <- NULL; } # don't override abcd if passed in
	return(BinnedBinaryFit(x,y, gSize=gSize, gSep=gSep, xout=xout,
				fT=fT, pActual=pActual, a=a, b=b, c=c, d=d, smooth=smooth, span=NULL, ...));
}

# x is numeric, y is binary, and a,b,c,d are for pseudocountLL
# returns a groups object with n(Observed F),n(Observed T), and p(Actual|x) for the group,
# with group bounds (min/med/max).
# Also, for each xout, returns  p=p(Actual) and index=which group
#
# More strictly, speaking, expects input x taken by sampling values with y=0 or 1.
# The sampling may be biased, taken according to p(y=T) = fT, as opposed to the actual
# distribution p(y=T) = pActual, and there may be noise in the category assignments,
# leading to the general form P(T observed|x) = a+b*p(T actual|x)/(c + d*p(T actual|x))
#
# Within each bin of x (60-200 overlapping bins, unless you override gSize or gSep),
# it observes nTobserved = n(X ^ y=T) and nFobserved = n(X ^ y=F), and uses pseudocountLL to estimate
# p(Actual|X) from nTobserved and nFobserved
#
# Then it interpolates between bins to get compute p(Actual|xout) for each xout
# Optionally, it uses loess (on rank(x) vs. odds ratios) to smooth this p-value (smooth=T)
#	or a similar loess on the groups (smooth="grouped") (recommended for large data sets)
#
BinnedBinaryFit <- function(x,y,gSize=NULL,gSep=NULL,xout=x,fT=NULL,pActual=.5,a=0,b=1,c=1,d=0,
			alpha=2,beta=2,debug=FALSE,
			plot=FALSE,log="",ylim=NULL,smooth=TRUE,span=NULL,...) {
	xsave <- xout;
	nTrueNA <- countif(is.na(x) & y);
	nFalseNA <- countif(is.na(x) & !y);

	# xsave will handle NAs
	y <- y[!is.na(x)];
	x <- x[!is.na(x)];
	if (!is.null(fT)) {
		# notice this is a no-op if fT==pActual
		ratio <- (pActual/(1-pActual)) * (1-fT)/fT;
		a = 0; b <- 1; c <- ratio; d <- 1-c;
		if(debug) cat("a b c d",a,b,c,d,"\n");
	}
	if (is.null(gSize)) {
		gSize <- 2 * (length(x) %/% 40);
		gSize <- pmin(200,pmax(100,gSize));
		if (length(x) < 100) {
			warning("BinnedBinaryFit on under 100 items, using a single bin");
			nT <- countif(y);
			nF <- if(length(y)>0) countif(!y) else 0;
			p <- pseudocountLL(nT,nF,a,b,c,d,alpha=alpha,beta=beta,debug=debug)$p;
			return(list(xSorted=x,groups=data.frame(min=-Inf,med=0,max=Inf,p=p,nT=nT,nF=nF),
					gSize=length(x),gSep=0,a=a,b=b,c=c,d=d,smooth=smooth,span=span,
					pNA=pseudocountLL(nTrueNA, nFalseNA,a,b,c,d,
							alpha=alpha,beta=beta,debug=debug)$p,
					p=rep(p,length(xsave)),
					index=rep(1,length(xsave))));

		}
		if (length(x) <= 2*gSize) {
			gSize <- length(x) %/% 3;
			warning("BinnedBinaryFit on only ",length(x)," items, reset gSize to ",gSize);
		}
	}
	if (is.null(gSep)) gSep <- gSize %/% 2;

	orderX <- order(x);
	xSorted <- x[orderX];
	nG <- ((length(x)-gSize) %/% gSep);
	gBeg <- 1 + gSep*(0:(nG-1));
	gEnd <- gBeg + gSize-1;
	gEnd[length(gEnd)] <- length(x); # the last group includes up to an extra gSep
	gMin <- xSorted[gBeg];
	gMax <- xSorted[gEnd];
	nT <- nF <- p <- gMed <- p <- 0*(1:length(gBeg)); 
	for (i in 1:length(p)) {
		if(debug) cat(i, gMin[i], gMax[i], gBeg[i], gEnd[i], length(xSorted),"\n");
		# use all values between min,max in case of a large group with the same sep
		gx <- select(x,x >= gMin[i] & x <= gMax[i]);
		gy <- select(y,x >= gMin[i] & x <= gMax[i]);
		gMed[i] <- median(gx);
		nT[i] <- countif(gy);
		nF[i] <- length(gy) - nT[i];
		if(debug) cat("", nT[i], nF[i], gMed[i], "\n");
		pLL <- pseudocountLL(nT[i],nF[i],a,b,c,d,alpha=alpha,beta=beta,debug=debug);
		if (pLL$code > 2) {
			cat("Warning: pseudocount failed","nT",nT[i],"nF",nF[i],
				"a",a,"b",b,"c",c,"d",d,"code",pLL$code,"p",pLL$p,"\n");
		}
		p[i] <- pLL$p;
		if(debug) {
			cat("group",i,"min",gMin[i],"max",gMax[i],
				"n",length(gx),"nT",nT[i],"nF",nF[i],"p",p[i],"\n");
		}
	}

	if (plot) {
		if (is.null(ylim)) { ylim <- (if(log=="y"||log=="xy") c(.001+min(p),0.5+max(c(nT,nF)))
					else c(0,1)); }
		plot(gMed,p,log=log,type="l", ylim=ylim, ...);
		for(i in 1:length(gMin)) {
			lines(c(gMin[i], gMax[i]),c(p[i],p[i]),lty=1,col="grey");
		}
		if (log=="y"||log=="xy") {
			points(gMed,.5+nT,type="l",col=3);
			points(gMed,.5+nF,type="l",col=2);
		}
	}

	bbf <- list(xSorted=xSorted,groups=data.frame(min=gMin,med=gMed,max=gMax,p=p,nT=nT,nF=nF),
			gSize=gSize,gSep=gSep,a=a,b=b,c=c,d=d,smooth=smooth,span=span,
			pNA=pseudocountLL(nTrueNA, nFalseNA,a,b,c,d,alpha=alpha,beta=beta,debug=debug)$p);

	if (length(unique(x)) < 2) smooth = FALSE;

	if (smooth == "groups") {
		bbf$groups$pRough = bbf$groups$p;
		x = rank(bbf$groups$med);
		y = logodds(bbf$groups$pRough);
		odds = loess(y ~ x, bbf$groups)$fitted;
		bbf$groups$p = logodds2p(odds);
	}
	if (smooth == TRUE) {
		out <- BBFPredict(xsave,bbf);
		lo <- logodds(out$p);
		bbf$loess <- loess(lo~rank(xsave),span=if(is.null(span)) 0.5 else span);
	}

	out <- BBFPredict(xsave,bbf);
	bbf$index <- out$index;
	bbf$p <- out$p;
	return(bbf);
}

# interpolate because it should help in fitting separations for small genomes
# (e.g. Buchnera), and makes little difference elsewhere
# rank(median for group i) = 1 + gSep*(i-1) + gSize/2, so
# i = 1 + (rank - 1 - gSize/2)/gSep
BBFPredict <- function(xeval, bbf) {
	xeval2 <- xeval[!is.na(xeval)];
	indexScaled <- pmax(1,pmin(length(bbf$xSorted),findIntervalAvg(xeval2,bbf$xSorted)));
	index <- 1 + ( indexScaled - 1 - bbf$gSize/2 )/bbf$gSep;
	indexMM <- pmin(length(bbf$groups$p),pmax(1,index));
	p <- approx(1:length(bbf$groups$p),bbf$groups$p, xout=indexMM)$y;
	expand <- function(values,default,x) {
		v <- rep(default,length(x));
		v[!is.na(x)] <- values;
		return(v);
	}

	if(!is.null(bbf$loess)) {
		pRough <- p;
		limits <- range(bbf$loess$x);
		indexScaled <- pmin(pmax(indexScaled, limits[1]),limits[2]);
		p <- logodds2p(predict(bbf$loess,indexScaled));
		return(data.frame(index=expand(index,-1,xeval),
					pRough=expand(pRough,bbf$pNA,xeval),
					indexScaled=expand(indexScaled,-1,xeval),
					p=expand(p,bbf$pNA,xeval)));
	} else {
		return(data.frame(index=expand(index,-1,xeval),
					p=expand(p,bbf$pNA,xeval)));
	}
}

CompareBBFGroups <- function(list,xlim=NULL,ylim=NULL,
				labels=names(list),
				col=1:length(list),pch=1:length(list),lty=1:length(list),lwd=rep(1,length(col)),
				ltyH=1:length(list),
				hbar=TRUE,
				xlab="Median of bin", ylab="Log odds",legend=TRUE,
				bins=TRUE,
				...) {
	if (is.null(xlim)) {
		xlim <- c( min(sapply(list, function(x) { min(x$med) ;})),
				max(sapply(list, function(x) { max(x$med) ;})) );
	}
	if (is.null(ylim)) {
		ylim <- c( min(sapply(list, function(x) { min(logodds(x$p)); })),
				max(sapply(list, function(x) { max(logodds(x$p)); })) );
	}
	for (i in 1:length(names(list))) {
		name <- names(list)[i];
		y <- logodds(list[[name]]$p);
		if (i == 1) {
			plot(list[[name]]$med,y,col=if(bins) col[i] else 0,pch=pch[i],lty=lty[i],lwd=lwd[i],
				xlim=xlim,ylim=ylim,type="b",xlab=xlab,ylab=ylab,...);
		} else {
			if(bins) points(list[[name]]$med,y,col=col[i],pch=pch[i],lty=lty[i],lwd=lwd[i],type="b");
		}
		if (hbar) {
			linePairs(list[[name]]$min, y, list[[name]]$max, y, col=col[i],lwd=lwd[i],lty=ltyH[i]);
		}
	}
	if(legend) {
		legend(mean(xlim),max(ylim),yjust=1,labels,pch=pch,col=col,lty=lty,lwd=lwd);
	}
}

# like BinnedBinaryFit but where x is a factor. Returns a data frame of p, nT, nF, levels
ClassFit <- function(x,y,fT=NULL,pActual=.5,a=0,b=1,c=1,d=0,
			alpha=2,beta=2,plot=FALSE,debug=FALSE,
			xlab="",main="",log="") {
	if (length(levels(x)) < 2) stop("ClassFit called with less than two levels, length(x) = ", length(x));
	if (!is.null(fT)) {
		# notice this is a no-op if fT==pActual
		ratio <- (pActual/(1-pActual)) * (1-fT)/fT;
		a = 0; b <- 1; c <- ratio; d <- 1-c;
		if(debug) cat("a b c d",a,b,c,d,"\n");
	}
	nT <- nF <- p <- rep(0, length(levels(x)));
	for (i in 1:length(levels(x))) {
		nT[i] <- sum(unclass(x)==i & y, na.rm=T);
		nF[i] <- sum(unclass(x)==i & !y, na.rm=T);
		if(debug) cat("ClassFit", i, nT[i], nF[i],"\n");
		p[i] <- pseudocountLL(nT[i],nF[i],a,b,c,d,alpha=alpha,beta=beta,debug=debug)$p;
	}
	if (plot) {
		plot(1:length(levels(x)), p,main=main,log=log,xlab=xlab);
	}
	return(data.frame(p=p,nT=nT,nF=nF,levels=levels(x)));
}

BayesianPredict <- function(class,value,vars,data,fT=countPerIf(data[[class]]==value)) {
	bdata <- list();
	bdata[[class]] <- ifelse(data[[class]]==value,1,0);
	formula <- paste(class, sep="");
	for (v in vars) {
		p <- BinnedBinaryFit(data[[v]], data[[class]]==value, fT=fT)$p;
		bdata[[v]] <- log(p/(1-p));
		formula <- paste(formula, ifelse(v==vars[1], "~", "+"),v,sep=" ");
	}
	bdata <- data.frame(bdata);
	formula <- eval(parse(text=formula));
	glm <- glm(formula,family=binomial(logit),data=bdata);
	return(list( glm=glm, data=bdata, cor=cor(fitted(glm),bdata[[class]]) ));
}

# same and non have both Gene1 and Gene2
OperonsKnown <- function(op, same=NULL, non=NULL) {
	knownY <- op$preds$Gene1 %in% same$Gene1 &
			findMatch(op$preds$Gene1,same$Gene1)$index == findMatch(op$preds$Gene2,same$Gene2)$index;
	knownN <- op$preds$Gene1 %in% non$Gene1 &
			findMatch(op$preds$Gene1,non$Gene1)$index == findMatch(op$preds$Gene2,non$Gene2)$index;
	op$preds$known <- ifelse(knownY,1,ifelse(knownN,-1,0));
	return(op);
}

# Takes as input OperonsPredict's output and a list of Gene1's for operon and for non-operon pairs
OperonsAOC <- function(op, same=NULL, non=NULL,
			plot="all", seed=.Random.seed, debug=FALSE, nfold=10) {
	if(is.null(same) || is.null(non)) {
		if (op$preds$Taxon[1] == 83333) {
			same = read.delim("operons/data/ecNewSameOp")$Gene1;
			non = read.delim("operons/data/ecNewAdj")$Gene1;
		} else if (op$preds$Taxon[1] == 1423) {
			same = read.delim("operons/data/bsSameOp")$Gene1;
			non = read.delim("operons/data/bsAdjNotOp")$Gene1;
		}
		else {
			cat("Needs same and non!\n"); return(NULL);
		}
	}
	same <- unique(sort(same));
	non <- unique(sort(non));
	known <- 0*op$preds$pOp;
	known[select(1:length(known),findMatch(op$preds$Gene1,same)$index>0)] <- 1;
	known[select(1:length(known),findMatch(op$preds$Gene1,non)$index>0)] <- -1;

	predsS <- select(op$preds, findMatch(op$preds$Gene1,same)$index > 0); # same operon
	if (length(predsS$Gene1) != length(same)) {
		cat("Warning: not all same-operon ids matched!\n");
	}
	predsD <- select(op$preds, findMatch(op$preds$Gene1,non)$index > 0); # different operon
	if (length(predsD$Gene1) != length(non)) {
		cat("Warning: not all non-operon ids matched!\n");
	}
	runif(1); # force seed creation
	trained <- PredictOperonsTrained(predsS, predsD, fOperon=op$fOperon, seed=seed, debug=debug,
					formula=eval(parse(text=op$formula)),
					vars=op$vars,
					postCAIdxdy=op$postCAIdxdy,
					smooth=op$smooth,
					all=op$preds, nfold=nfold );
	if(debug) { cat("PredictOperonsTrained returned",range(trained$predsS$pTrained),
				range(trained$predsD$pTrained),range(predsS$pOp),range(predsS$GNScore),"\n"); }
	aocTrained <- AOC(trained$predsS$pTrained,trained$predsD$pTrained);
	aocUntrained <- AOC(predsS$pOp,predsD$pOp);
	aocUntrained1 <- AOC(predsS$pOp1,predsD$pOp1);
	if(debug) { cat("3 AOC done"); }

	# Naive is to use GNScore > 0 or smaller sep is better (up to -20)
	naive <- function(preds) {
		return(ifelse(preds$GNScore > 0, -1, ifelse(preds$Sep< -20, 1e4, preds$Sep)));
	}
	aocNaive <- AOC(naive(predsS),naive(predsD),bigger=FALSE);
	aocDistance <- AOC(predsS$pSameD,predsD$pSameD);
	aocPSame <- AOC(predsS$pSame,predsD$pSame);
	if(debug) { cat("6 AOC done\n"); }
	if (plot=="all") {
		old.par <- par(no.readonly=TRUE);
		on.exit(par(old.par));
		par(mfrow=c(2,1));
	}
	if (plot=="all" || plot=="hist") {
		CompareHistogramList(list(Operon=predsS$pOp,Adjacent=predsD$pOp),breaks=c(0:20)/20,
					xlab="p(Operon) for known-operon and adjacent pairs",
					main=op$name);
	}
	if(debug) { cat("2 plots done\n"); }
	if (plot=="all" || plot=="AOC") {
		plot(1-aocUntrained$spec,aocUntrained$sens,type="l",lty=1,col=1,xlab="FP",ylab="TP",
			main=paste(op$name,"--",op$control,"control","--",op$formula,
					if(op$postCAIdxdy) "then CAIdxdy" else "",
					"\nAOC =",format(aocUntrained$AOC,digits=3)));
		if(debug) { cat("count",countPerIf(predsD$bOp),"\n"); }
		points(countPerIf(predsD$bOp), countPerIf(predsS$bOp), col=1, pch=2);
		points(1-countPerIf(predsD$pOp<.9), countPerIf(predsS$pOp>=.9), col=1, pch=3);
		lines(1-aocNaive$spec, aocNaive$sens, lty=4, col=2);
		points(countPerIf(predsD$GNScore>=23 | (predsD$Sep>=-20 & predsD$Sep < 55)),
			countPerIf(predsS$GNScore>=23 | (predsS$Sep>=-20 & predsS$Sep < 55)),
			pch=5, col=3);
		lines(1-aocDistance$spec, aocDistance$sens, lty=6, col=4);
		points(1-countPerIf(predsD$pSameD<.5), countPerIf(predsS$pSameD>=.5), col=4,pch=20);
		lines(1-aocPSame$spec, aocPSame$sens, lty=7, col=5);
		fH <- 1/(1+op$pCoverS);
		points(1-countPerIf(predsD$pSame<fH), countPerIf(predsS$pSame>=fH), col=5,pch=20);
		lines(1-aocTrained$spec, aocTrained$sens, lty=8, col=6, lwd=2);
		points(1-countPerIf(trained$predsD$pTrained<.5), countPerIf(trained$predsS$pTrained>=.5),
			col=6,pch=20);
		vline(0);hline(1);
		legend(1,.8,c("Untrained", "Untrained p>=fOp", "Untrained p>=.9",
				"GNScore > 0 or Sep in [-20,N)",
				"GNScore >= 23 or Sep in [-20,55)",
				"Distance only (after logistic)",
				"pSame only",
				"Trained"),
			pch=c(30,2,3,30,5,30,30,30),lty=c(1,0,0,4,0,6,7,8),col=c(rep(1,3),2:6),
			lwd=c(rep(1,7),2),xjust=1);
		f <- function(x) { return(format(x,digits=3)); }
		cat("Untrained AOC\t", f(aocUntrained$AOC),
					"\tmax acc\t", f(0.5*max(aocUntrained$spec+aocUntrained$sens)),
					"\tdef sens\t", f(countPerIf(predsS$bOp)),
					"\tspec\t",f(countPerIf(!predsD$bOp)), "\n",sep="");
		cat("pSame&Dist AOC\t", f(aocUntrained1$AOC),
					"\tmax acc\t", f(0.5*max(aocUntrained1$spec+aocUntrained1$sens)),
					"\tdef sens\t", f(countPerIf(predsS$bOp1)),
					"\tspec\t",f(countPerIf(!predsD$bOp1)),"\n",sep="");
		cat("Dist-only AOC\t",f(aocDistance$AOC),
					"\tmax acc\t", f(0.5*max(aocDistance$spec+aocDistance$sens)),
					"\tdef sens\t",f(countPerIf(predsS$pSameD>=.5)),
					"\tspec\t",f(countPerIf(predsD$pSameD<.5)),"\n",sep="");
		cat("pSame-only AOC\t",f(aocPSame$AOC),
					"\tmax acc\t", f(0.5*max(aocPSame$spec+aocPSame$sens)),
					"\tdef sens\t",f(countPerIf(predsS$pSame>=fH)),
					"\tspec\t",f(countPerIf(predsD$pSame<fH)),"\n",sep="");
		cat("Trained AOC\t",f(aocTrained$AOC),
					"\tmax acc\t", f(0.5*max(aocTrained$spec+aocTrained$sens)),
					"\tdef sens\t", f(countPerIf(trained$predsS$pTrained>=.5)),
					"\tspec\t", f(countPerIf(trained$predsD$pTrained<.5)),
					"\nGNScore>0 AOC\t",f(aocNaive$AOC), "\n", sep="");
	}
	out <- list(predsS=trained$predsS, predsD=trained$predsD,
		aocUntrained=aocUntrained,aocNaive=aocNaive,aocDistance=aocDistance,aocPSame=aocPSame,
		aocTrained=aocTrained,seed=seed,nfold=trained$nfold,
		gcm=trained$gcm,
		bbfD=trained$bbfD, bbfCAI=trained$bbfCAI,cfCOG=trained$cfCOG,
		defSens=countPerIf(predsS$bOp), defSpec=countPerIf(!predsD$bOp),
		maxAcc=.5*max(aocUntrained$spec+aocUntrained$sens),
		expectedDRatio = (1 - countPerIf(!predsS$bOp) - countPerIf(predsD$bOp)),
		known=known,
		pTrained=trained$pAll, pTrainedD=trained$pAllD, pTrainedLogistic=trained$pAllLogistic,
		pTrainedCAI=trained$pAllCAI,
		trained=data.frame(without(trained,trained$without)));
	for (var in op$vars) {
		bbfvar <- paste("bbf",var,sep="");
		out[[bbfvar]] <- trained[[bbfvar]];
	}
	return(out);
}

# Ten-fold cross-validation...
PredictOperonsTrained <- function(predsS,predsD,fOperon=.5,
          all=NULL, #option to make predictions for all the values
          vars=c("MOGScore","GOScore","ExprSim", "Sep"),
          useCOGsim=TRUE,
          smooth=TRUE,
          nfold=10,seed=.Random.seed,debug=FALSE) {
  # we need to overwrite these...
  predsS$same = NULL;
  predsD$same = NULL;
  save.seed <- seed;
  set.seed(seed);
  set <- order(runif(length(c(predsS[[vars[1]]],predsD[[vars[1]]])))) %% nfold;
  preds <- rbind(cbind(predsS,same=1),cbind(predsD,same=0));
  for (var in vars) {
    preds[[paste("bbfTrained",var,sep="")]] <- 0*preds$same;
  }
  if(useCOGsim) cbind(preds, cfCOGTrained = rep(-1,nrow(preds)));
  if (!is.null(all)) {
    logoddsAll <- logoddsAllD <- logoddsAllLogistic <- 0*all$Gene1;
  }
  out <- NULL;
  bbfvars = sapply(vars,function(v) paste("bbf",v,sep=""));
  if (useCOGsim) bbfvars = c(bbfvars,"cfCOG");
  formula = eval(parse(text=paste("same","~",paste(bbfvars,collapse="+"),sep="")));

  if(debug) { cat("Formula",paste(deparse(formula)),
      "vars",vars,
      "seed start",seed[1],"fOp",fOperon,
      "nS", nrow(predsS), "nD", nrow(predsD), "\n");
      cat("Names",names(predsS),names(predsD),"\n");
  }

  for (i in 0:(nfold-1)) {
    train <- select(preds,set != i);
    test <- select(preds, set == i);
    pCoverS <- countif(train$same==0)/countif(train$same==1);

    # do the bbf thing
    for (var in vars) {
      if(debug) {cat("BBF",var,range(train[[var]]),"\n"); }
      bbf <- BinnedBinaryFit(train[[var]], train$same, fT=1/(1+pCoverS),
            smooth=smooth,span=1,
            plot=debug,xlab=var);
      if(debug) {cat(var,range(bbf$p),"\n");}
      bbfvar <- paste("bbf",var,sep="");
      train[[bbfvar]] <- log(bbf$p/(1-bbf$p));
      testP <- BBFPredict(test[[var]],bbf)$p;
      test[[bbfvar]] <- log(testP/(1-testP));
      out[[bbfvar]] <- bbf;
      if(countif(is.na(test[[bbfvar]])) > 0 | countif(is.na(train[[bbfvar]])) > 0) {
        cat("Bad bbf in PredictOperonsTrained\n");
        if(debug) {
          return(list(bbf=bbf,testVar=test[[var]]));
        } else {
          return(NULL);
        }
      }
      if (!is.null(all)) {
        all[[bbfvar]] <- logodds(BBFPredict(all[[var]],bbf)$p);
      }
    }
    if (useCOGsim) {
      cfCOG <- ClassFit(unclass(train$COGSim),train$same);
      train$cfCOG <- logodds(cfCOG$p[unclass(train$COGSim)]);
      test$cfCOG <- logodds(cfCOG$p[unclass(test$COGSim)]);
    } else {
      cfCOG <- NULL;
    }
    if(debug) cat("nTrain",nrow(preds),"nTrainS",countif(train$same==1),
          "nTrainD",countif(train$same==0),"\n");
    gcm <- glm(formula,family=binomial(logit),data=train);
    if(debug) { s <- summary(gcm); s$family <- NULL; s$call <- paste(deparse(s$call));
        s$terms <- paste(deparse(s$terms)); for(n in names(s)) { cat(n,s[[n]],"\n"); } }
    # now a misnomer...response is so we get p-values not log-odds
    if(debug) { cat("do predict",apply(test,2,length),"\n"); }
    pOp <- predict(gcm,test,type="response");
    if(debug) { cat("did predict",length(pOp),length(test$cfCOG),"\n"); }
    if (length(pOp) != length(test$cfCOG)) {
      if(debug) { cat("Returning bad gcm, test\n"); return(list(gcm=gcm,test=test)); }
      cat("Fatal error in Trained predictions...\n");
      return(NULL);
    }
    if(debug) cat("Names",names(preds),"cfCOGTrained",length(preds$cfCOGTrained),length(test$cfCOG),"\n");
    preds$pOp[set==i] = pOp;
    if (useCOGsim) preds$cfCOGTrained[set==i] = test$cfCOG;
    for (var in vars) {
      preds[[paste("bbfTrained",var,sep="")]][set==i] =
        test[[paste("bbf",var,sep="")]];
    }
    if (!is.null(all)) {
      # note that the bbfvars were modified above to contain values based on new models...
      all$cfCOG <- logodds(cfCOG$p)[unclass(all$COGSim)];
      pAllLogistic <- combineP(pCoverS/(1+pCoverS),predict(gcm,all,type="response"));
      pAllD <- BBFPredict(all$Sep, bbfD)$p;
      logoddsAll <- logoddsAll +
        (logodds(pAllLogistic) + logodds(pAllD))/nfold;
      logoddsAllLogistic <-   logoddsAllLogistic + logodds(pAllLogistic)/nfold;
      logoddsAllD <- logoddsAllD + logodds(pAllD)/nfold;
    }
    if(debug) { cat("Finished round",i,"\n");  }
  }

  out <- c(out, list(predsS=select(preds,preds$same==1),predsD=select(preds,preds$same==0),
        cfCOG=cfCOG, gcm=gcm,
        seed=save.seed, vars=vars ));
  if (!is.null(all)) {
    out$pAll <- exp(logoddsAll)/(1+exp(logoddsAll));
    index <- findMatch(preds$Gene1, all$Gene1)$index;
    out$pAll[index>0] - preds$pTrained[index>0];
  }

  if(debug) { cat("Finished Training\n");  cat("Done\n"); }
  return(c(out, list(without=c("without",names(out)),
      nfold=nfold, formula=paste(deparse(formula)),
      pCoverS=pCoverS
  )));
}

# handles output from OperonsAOC
OperonsCompareAOC <- function(aoc1,aoc2,var="pOp",plot=FALSE) {
	return(CompareAOC(aoc1$predsS[[var]],aoc1$predsD[[var]],aoc2$predsS[[var]],aoc2$predsD[[var]],plot=plot));
}

CompareAOC2 <- function(aoc,var1="pOp",var2="pTrained",plot=FALSE) {
	return(CompareAOC(aoc$predsS[,var1],aoc$predsD[,var1],
			aoc$predsS[,var2],aoc$predsD[,var2],plot=plot));
}

AOCStandardError <- function(s1,d1,nOn=length(s1),nOff=length(d1),area=makeaoc(s1,d1)) {
	Q1 <- area/(2-area);
	Q2 <- 2*area*area/(1+area);
	top <- area*(1-area) + (nOn-1)*(Q1-area*area) + (nOff-1)*(Q2-area*area);
	return( sqrt(top/(nOn*nOff)) );
}

# Uses the method of DeLong et al 1988:
# DeLong ER, DeLong DM, Clarke-Pearson DL.
# Biometrics. 1988 Sep;44(3):837-45.
# Comparing the areas under two or more correlated receiver operating characteristic curves:
# a nonparametric approach.
CompareAOC3 <- function(m1,m2,truth,...) {
	return(CompareAOC(select(m1,truth),select(m1,!truth),select(m2,truth),select(m2,!truth),...));
}

CompareAOC <- function(s1,d1,s2,d2,plot=FALSE, col=1:2, lty=1:2, ...) {
	p1 <- findIntervalAvg(s1, sort(d1))/length(d1); # s1 in d1
	q1 <- 1 - findIntervalAvg(d1, sort(s1))/length(s1); # d1 in s1
	area1 <- mean(p1);

	p2 <- findIntervalAvg(s2, sort(d2))/length(d2); # s2 in d2
	q2 <- 1 - findIntervalAvg(d2, sort(s2))/length(s2); # d2 in s2
	area2 <- mean(p2);

	if(plot) {
		plot((1:length(p1))/length(p1), sort(p1), col=col[1], lty=lty[1], type="l", ... );
		lines((1:length(p2))/length(p2), sort(p2), col=col[2],lty=lty[2], ... );
	}

	m <- length(s1);
	n <- length(d1);

	var1 <- (m*sum(p1**2)-sum(p1)**2)/(m*m*(m-1)) + (n*sum(q1**2)-sum(q1)**2)/(n*n*(n-1)) +
			area1*(1-area1)/(m*n);
	var2 <- (m*sum(p2**2)-sum(p2)**2)/(m*m*(m-1)) + (n*sum(q2**2)-sum(q2)**2)/(n*n*(n-1)) +
			area2*(1-area2)/(m*n);
	cov12 <- sum((p1-mean(p1))*(p2-mean(p2)))/(m*(m-1)) + sum((q1-mean(q1))*(q2-mean(q2)))/(n*(n-1));
	se <- sqrt(var1+var2-2*cov12);
	z <- (area1-area2)/se;
	return(data.frame(area1=area1,area2=area2,
			m=m,n=n,var1=var1,var2=var2,cov12=cov12,se=se,
			varp1=var(p1),varq1=var(q1),
			z=z,p=2*pnorm(abs(z),lower=FALSE)));
}

OperonUarrayCor <- function(op, orfData, aoc=NULL, skip=3, minN=10, plot=TRUE, name="orfId",spearman=FALSE) {
	i1 <- findMatch(op$preds$Gene1,orfData[,name])$index;
	i2 <- findMatch(op$preds$Gene2,orfData[,name])$index;
	preds <- select(cbind(op$preds,i1=i1,i2=i2),i1>0&i2>0);
	highOp <- rank(preds$pOpDistance) > .8*length(preds$pOpDistance) &
			rank(preds$pOpLogistic) > .8*length(preds$pOpLogistic);
	lowOp <- rank(preds$pOpDistance) < .2*length(preds$pOpDistance) &
			rank(preds$pOpLogistic) < .2*length(preds$pOpLogistic);
	cat("kept",length(preds$i1),"high",countif(highOp),"low",countif(lowOp),"\n")
	outNames <- names(orfData)[(skip+1):length(names(orfData))];
	out <- data.frame(name=outNames);
	if (!is.null(aoc)) {
		iS1a <- findMatch(aoc$predsS$Gene1,orfData$orfId)$index;
		iS2a <- findMatch(aoc$predsS$Gene2,orfData$orfId)$index;
		iS1 <- select(iS1a, iS1a>0 & iS2a>0);
		iS2 <- select(iS2a, iS1a>0 & iS2a>0);

		iD1a <- findMatch(aoc$predsD$Gene1,orfData$orfId)$index;
		iD2a <- findMatch(aoc$predsD$Gene2,orfData$orfId)$index;
		iD1 <- select(iD1a, iD1a>0 & iD2a>0);
		iD2 <- select(iD2a, iD1a>0 & iD2a>0);

		cat("knownOp",length(iS1),"adjKnownOp",length(iD1),"\n");
	}
	for (j in (skip+1):length(names(orfData))) {
		x <- orfData[,names(orfData)[j]];
		if (!is.null(aoc)) {
			out$knownOp[j-skip] = corNoNa(x[iS1],x[iS2],spearman=spearman);
			out$adjKnownOp[j-skip] = corNoNa(x[iD1],x[iD2],spearman=spearman);
		}
		x1 <- x[preds$i1];
		x2 <- x[preds$i2];
		if (countif(!is.na(x1) & !is.na(x2) & highOp) >= minN
			&& countif(!is.na(x1) & !is.na(x2) & lowOp) >= minN) {
			out$highOp[j-skip]=corNoNaSelect(x1,x2,highOp,spearman=spearman);
			out$lowOp[j-skip]=corNoNaSelect(x1,x2,lowOp,spearman=spearman);
			out$bOp[j-skip]=corNoNaSelect(x1,x2,preds$bOp,spearman=spearman);
			out$notOp[j-skip]=corNoNaSelect(x1,x2,!preds$bOp,spearman=spearman);
			out$all[j-skip]=corNoNa(x1,x2,spearman=spearman);
		}
	}
	return(out);
}

OperonUarrayAccuracy <- function(op, orfData, dims, plot=TRUE, fixedR=FALSE,name="orfId") {
	acc <- 0*dims;
	i1 <- findMatch(op$preds$Gene1, orfData[,name])$index;
	i2 <- findMatch(op$preds$Gene2, orfData[,name])$index;
	preds <- select(cbind(op$preds,i1=i1,i2=i2),i1>0&i2>0);
	highOp <- rank(preds$pOpDistance) > .8*length(preds$pOpDistance) &
			rank(preds$pOpLogistic) > .8*length(preds$pOpLogistic);
	lowOp <- rank(preds$pOpDistance) < .2*length(preds$pOpDistance) &
			rank(preds$pOpLogistic) < .2*length(preds$pOpLogistic);
	out <- data.frame(name=names(orfData)[dims],dim=dims,corLow=0*dims,corHigh=0*dims,nH=0*dims,nL=0*dims,
				accOp=0*dims,accNotOp=0*dims);
	if(plot) {
		old.par <- par(no.readonly=TRUE);
		on.exit(par(old.par));
		n <- ceiling(sqrt(1*length(dims)));
		par(mfrow=c(n,n));
	}
	for (i in 1:length(dims)) {
		x <- orfData[,dims[i]];

		xH <- select(x[preds$i1],!is.na(x[preds$i1]) & !is.na(x[preds$i2]) & highOp);
		yH <- select(x[preds$i2],!is.na(x[preds$i1]) & !is.na(x[preds$i2]) & highOp);
		xL <- select(x[preds$i1],!is.na(x[preds$i1]) & !is.na(x[preds$i2]) & lowOp);
		yL <- select(x[preds$i2],!is.na(x[preds$i1]) & !is.na(x[preds$i2]) & lowOp);

		if (length(xH) < 25 || length(yH) < 25) {
			cat("Ignoring dim",i," -- length <25\n");
			out$nH[i] <- out$nL[i] <- out$r1[i] <- out$r2[i] <- out$accOp[i] <- out$accNotOp[i] <-
					out$accOp2[i] <- out$accNotOp2[i] <- out$acc3[i] <- out$acc3D[i] <-
					out$accR3[i] <-
					out$acc3Max[i] <- out$ac3D[i] <- out$accR4[i] <- out$acc4[i] <-
					out$p.wilcox[i] <- out$p.wilcoxHL[i] <- -1;
		} else {
			out$corHigh[i] <- cor(xH,yH);
			out$nH[i] <- length(xH);
	
			out$corLow[i] <- cor(xL,yL);
			out$nL[i] <- length(xL);

			xOp <- select(x[preds$i1],!is.na(x[preds$i1]) & !is.na(x[preds$i2]) & preds$bOp);
			yOp <- select(x[preds$i2],!is.na(x[preds$i1]) & !is.na(x[preds$i2]) & preds$bOp);

			xNotOp <- select(x[preds$i1],!is.na(x[preds$i1]) & !is.na(x[preds$i2]) & !preds$bOp);
			yNotOp <- select(x[preds$i2],!is.na(x[preds$i1]) & !is.na(x[preds$i2]) & !preds$bOp);

			if (fixedR) {
				r1 <- .5; r2 <- 0;
			} else {
				r1 <- cor(xH,yH); r2 <- cor(xL,yL);
			}
			out$r1[i] <- r1;
			out$r2[i] <- r2;
			prH <- logpRatioUnnorm(xH,yH,r1,r2);
			prL <- logpRatioUnnorm(xL,yL,r1,r2);
			prOp <- logpRatioUnnorm(xOp,yOp,r1,r2);
			prNotOp <- logpRatioUnnorm(xNotOp,yNotOp,r1,r2);
			h <- CompareHistogramList(list(highOp=prH,lowOp=prL,bOp=prOp,notOp=prNotOp),
							breaks=c(-1,-.5,-.25,-.15,-.05,0,.05,.15,.25,.5,1,1.5,2,2.5),
							ylim=c(0,.3),
							plot=FALSE,returnP=TRUE,main=names(orfData)[dims[i]],
							warn=FALSE);
			model1 <- lm(bOp~lowOp+highOp+0,h,weights=lowOp+highOp+1/(length(prL)+length(prH)));
			model2 <- lm(notOp~lowOp+highOp+0,h,weights=lowOp+highOp+1/(length(prL)+length(prH)));

			out$accOp[i] <- model1$coefficients[2]/sum(model1$coefficients);
			out$accNotOp[i] <- model2$coefficients[1]/sum(model2$coefficients);

			ksR <- ksRanks(list(prH=prH,prL=prL,prOp=prOp,prNotOp=prNotOp));
			ks1 <- lm(prOp~prL+prH+0,ksR);
			ks2 <- lm(prNotOp~prL+prH+0,ksR);
			out$accOp2[i] <- ks1$coefficients[2]/sum(ks1$coefficients);
			out$accNotOp2[i] <- ks2$coefficients[1]/sum(ks2$coefficients);

			ksR <- ksRanks(list(prH=prH,prL=prL,prOp=prOp,prNotOp=prNotOp));

			ks3 <- lm(I(prOp-prNotOp) ~ I(prH-prL)+0, ksR);
			out$acc3[i] <- ks3$coefficients[1];
			out$accR3[i] <- cor(ksR$prOp-ksR$prNotOp, ksR$prH-ksR$prL);
			out$acc3Max[i] <- min(ksR$prOp-ksR$prNotOp)/min(ksR$prH-ksR$prL); # both are negative
			out$acc3D[i] <- ks.test(prOp,prNotOp)$statistic/ks.test(prH,prL)$statistic;
			if (median(prOp) < median(prNotOp) || median(prH) < median(prL) || r1 <= r2) out$acc3D[i] <- 0;

			if (FALSE) {
				plot(ksR$per,ksR$prH,type="l",col=3);
				eqline(col=1);
				lines(ksR$per,ksR$prL,type="l",col=2);
				lines(ksR$per,ksR$prOp,type="l",col=3,lty=2);
				lines(ksR$per,ksR$prNotOp,type="l",col=2,lty=2);
				lines(ksR$per, ksR$prH*out$accOp2[i] + ksR$prL*(1-out$accOp2[i]),col="grey",lty=3);
				lines(ksR$per, ksR$prL*out$accNotOp2[i] + ksR$prH*(1-out$accOp2[i]),col="grey",lty=3);
			}
			if (plot) {
				plot(ksR$per, ksR$prH-ksR$prL, type="l");
				hline(0);
				lines(ksR$per, ksR$prOp-ksR$prNotOp, col=2);
				lines(ksR$per, out$acc3[i]*(ksR$prH-ksR$prL), col="grey", lty=2);
			}
			xBoth <- select(x[preds$i1], !is.na(x[preds$i1]) & !is.na(x[preds$i2]));
			yBoth <- select(x[preds$i2], !is.na(x[preds$i1]) & !is.na(x[preds$i2]));
			pOp <- select(preds$pOp, !is.na(x[preds$i1]) & !is.na(x[preds$i2]));
			prBoth <- logpRatioUnnorm(xBoth,yBoth,r1,r2);
			out$accR4[i] <- cor(rank(pOp),rank(prBoth));
			out$acc4[i] <- lm(I(rank(prBoth)-length(prBoth)/2)~I(rank(pOp)-length(pOp)/2)+0)$coefficients[1];
			if (FALSE) {
				nBreaks <- 20;
				breaks <- (0:nBreaks) * (length(pOp) / nBreaks);
				mids <- (breaks[1:nBreaks]+breaks[2:(nBreaks+1)])/2;
				x <- aggregate(rank(prBoth),list(pOp=cut(rank(pOp),breaks=breaks)),mean);
				plot(mids,x$x,type="l");
				lines(0:length(pOp),length(pOp)/2 + ((0:length(pOp)) - length(pOp)/2)*out$acc4[i],col="grey");
			}
			out$p.wilcox[i] <- wilcox.test(prOp,prNotOp)$p.value;
			out$p.wilcoxHL[i] <- wilcox.test(prH,prL)$p.value;
		}
	}
	return(out);
}

UarrayToSums <- function(gene, orfData, dims, name="orfId", orfs=orfData[,name]) {
	i <- findMatch(gene, orfs)$index;
	out <- rep(0, length(gene));
	out[select(1:length(gene),i>0)] <- apply(orfData[select(i,i>0),dims], 1, sum, na.rm=TRUE);
	return(out);
}

UarrayOrfPairCor <- function(gene1, gene2, orfData, dims, minN=10, orfs=orfData[,"locusId"],debug=FALSE) {
	i1 <- findMatch(gene1, orfs)$index;
	i2 <- findMatch(gene2, orfs)$index;
	out <- vector(length=length(gene1));
	n <- sumD <- sumD1 <- sumDsqr <- diff2 <- diff1 <- 0*(1:length(out));

	for (i in 1:length(out)) {
		if(debug) { cat(i,"\n"); }
		if (i1[i] > 0 && i2[i] > 0) {
			dimsUse <- dims[!is.na(orfData[i1[i],dims]) & !is.na(orfData[i2[i],dims])];
			n[i] <- length(dimsUse);
			if (n[i] >= minN) {
				d1 <- t(orfData[i1[i],dimsUse]);
				d2 <- t(orfData[i2[i],dimsUse]);
				out[i] <- cor(d1,d2);
				sumD[i] <- sum(abs(d1))+sum(abs(d2));
				sumD1[i] <- sum(abs(d1));
				diff1[i] <- sum(abs(d1-d2));
				diff2[i] <- sum((d1-d2)**2);
				sumDsqr[i] <- sum(d1**2)+sum(d2**2);
			}
		} else {
			out[i] <- NA;
		}
	}
	return(data.frame(Gene1=gene1, Gene2=gene2,
				r=out,n=n,sumD=sumD,sumD1=sumD1,sumDsqr=sumDsqr,diff1=diff1,diff2=diff2,
				dist=sqrt(diff2/sumDsqr),dist1=diff1/sumD));
}

UarrayOrfLevelAll <- function(data, treatment, control, dims, minN=10, name="orfId") {
	means <- apply(data[,dims],2,meanNoNA);
	swept <- sweep(data[,dims],2,means);

	out <- data.frame(locusId=data[[name]],
			expd=apply(abs(data[,dims]),1,meanNoNA),
			sd=apply(data[,dims],1,function(x) {
					x <- x[!is.na(x)]; if(length(x)<2) return(NA); return(sd(x));
				}),
			S=log10(apply(control[,dims]+treatment[,dims], 1, meanNoNA)),
			nExp=apply(!is.na(data[,dims]),1,countif)
			);
	return(out[out$nExp>=minN,]);
}

UarrayOrfPairCorAll <- function(preds, data, dims, minN=10, name="orfId",debug=FALSE) {
	means <- apply(data[,dims],2,meanNoNA);
	swept <- sweep(data[,dims],2,means);
	return(UarrayOrfPairCor(preds$Gene1,preds$Gene2,swept,1:length(swept),minN=minN,orfs=data[,name],debug=debug));
}

ReadUarrayOpCors <- function(prefix,op,dir="operons/data",buildOrfPairs=TRUE,minN=10) {
	if (buildOrfPairs) {
		ratios <- read.delim(paste(dir,"/",prefix,"Ratios",sep=""));
		controlDyes <- read.delim(paste(dir,"/",prefix,"ControlDyes",sep=""));
		treatmentDyes <- read.delim(paste(dir,"/",prefix,"TreatmentDyes",sep=""));
		ses <- read.delim(paste(dir,"/",prefix,"SEs",sep=""));

		dims <- select(2:length(ratios),ses$condition != "Genomic");

		S1c <- UarrayToSums(op$gc$Gene1, controlDyes, dims);
		S1t <- UarrayToSums(op$gc$Gene1, treatmentDyes, dims);
		S2c <- UarrayToSums(op$gc$Gene2, controlDyes, dims);
		S2t <- UarrayToSums(op$gc$Gene2, treatmentDyes, dims);
		S4 <- pmax(0,pmin(S1c+S1t,S2c+S2t));
		cors <- cbind(UarrayOrfPairCorAll(op$gc, ratios, dims, minN=minN),
				S4=S4,S1c=S1c,S1t=S1t,S2c=S2c,S2t=S2t);
		cors$orfPairs[select(1:length(cors$n),cors$n < minN)] <- NA;
		writeDelim(cors,file=paste(dir,"/",prefix,"Cors",sep=""));
		return(cors);
	} else {
		cors <- read.delim(paste(dir,"/",prefix,"Cors",sep=""));
		cors$orfPairs[select(1:length(cors$n),cors$n < minN)] <- NA;
		return(cors);
	}
}

SimulateMixtureEstimation <- function(n=c(1000,100,250,250), mean2=.5, sd2=1, from=-2, to=2, tries=100) {
	out <- rep(0,tries);
	for (i in 1:tries) {
		trueOff <- pmax(from,pmin(to, rnorm(n[1]) ));
		trueOn <- pmax(from,pmin(to, rnorm(n[2],mean=mean2,sd=sd2) ));
		test <- pmax(from,pmin(to, c(rnorm(n[3]),rnorm(n[4],mean=mean2,sd=sd2)) ));
		d1 <- density(trueOff,from=from,to=to); d1$y <- d1$y/sum(d1$y);
		d2 <- density(trueOn,from=from,to=to);  d2$y <- d2$y/sum(d2$y);
		d3 <- density(test,from=from,to=to);    d3$y <- d3$y/sum(d3$y);
		out[i] <- lm(I(d3$y-d1$y)~I(d2$y-d1$y)+0)$coefficients[1];
	}
	return(out);
}

AccResampleGenes <- function(op,cors,n=NULL,useDensity=NULL,tries=100) {
	# four classes: 1 for opposing-strand, 2 for pred. not op, 3 for pred op (<.95), 4 for high confidence
	classSame <- cut(op$preds$pOp,c(0,.5,.95,1),labels=FALSE);
	class <- rep(1,nrow(cors));
	class[op$gcSC$same] <- 1+classSame;
	counts <- c(table(class));

	out <- data.frame(fOp=rep(0,tries),accT=rep(0,tries),accF=rep(0,tries));
	for (i in 1:tries) {
		# resample within each class
		index <- rep(0,nrow(cors));
		for (j in 1:4) {
			transfer <- (1:nrow(cors))[class==j];
			index[class==j] <- transfer[ 1+floor(runif(counts[j])*counts[j]) ];
		}			
		F <- fOperonVsUarray(op,cors[index,],n=n,useDensity=useDensity);
		out$fOp[i] <- lm(I(same-controlR)~I(highR-controlR)+0,F)$coefficients[1];
		out$accT[i] <- lm(I(op-controlO)~I(highO-controlO)+0,F)$coefficients[1];
		out$accF[i] <- lm(I(nop-highN)~I(controlN-highN)+0,F)$coefficients[1];
	}
	return(out);
}
		
UarrayOrfPairSubsampleAll <- function(prefix,op,dir="operons/data",minN=10) {
	ratios <- read.delim(paste(dir,"/",prefix,"Ratios",sep=""));
	ses <- read.delim(paste(dir,"/",prefix,"SEs",sep=""));
	out <- NULL;
	for (skip in levels(ses$condition)) {
		if (skip != "Genomic") {
			row <- UarrayOrfPairSubsample(skip, prefix, op, ratios=ratios, ses=ses);
			if (is.null(out)) { out <- row; } else { out <- rbind(out, row); }
		}
	}
	return(out);
}

UarrayOrfPairSubsample <- function(skip,prefix,op,dir="operons/data",minN=10,debug=FALSE,
		ratios=read.delim(paste(dir,"/",prefix,"Ratios",sep="")),
		ses=read.delim(paste(dir,"/",prefix,"SEs",sep=""))) {
	dims <- select(2:length(ratios),ses$condition != "Genomic" & ses$condition != skip);
	cat("Skipping",skip,"remaining arrays:",length(dims),"\n");
	minN <- min(minN,length(dims));
	cors <- UarrayOrfPairCorAll(op$gc, ratios, dims, minN=minN,debug=debug);
	cors$orfPairs[select(1:length(cors$n),cors$n < minN)] <- NA;

	F <- fOperonVsUarray(op,cors);
	fOp <- lm(I(same-controlR)~I(highR-controlR)+0,F)$coefficients[1];
	accT <- lm(I(op-controlO)~I(highO-controlO)+0,F)$coefficients[1];
	accF <- lm(I(nop-highN)~I(controlN-highN)+0,F)$coefficients[1];

	return(data.frame(skip=skip,fOp=fOp,accT=accT,accF=accF));
}

UarrayOrfPairSim <- function(preds, orfData, dim, name="orfId") {
	i1 <- findMatch(preds$Gene1, orfData[,name])$index;
	i2 <- findMatch(preds$Gene2, orfData[,name])$index;
	i12 <- select(cbind(preds,i1=i1,i2=i2), i1>0 & i2>0);
	withData <- select(i12, !is.na(orfData[i12$i1,dim]) & !is.na(orfData[i12$i2,dim]));
	return(cbind(withData,v1=orfData[withData$i1,dim],v2=orfData[withData$i2,dim],
			sim=orfData[withData$i1,dim]*orfData[withData$i2,dim]));
}

UarrayOrfPairTest <- function(preds, orfPairs, aoc=NULL) {
	highOp <- rank(preds$pOpDistance) > .8*length(preds$pOpDistance) &
			rank(preds$pOpLogistic) > .8*length(preds$pOpLogistic);
	lowOp <- rank(preds$pOpDistance) < .2*length(preds$pOpDistance) &
			rank(preds$pOpLogistic) < .2*length(preds$pOpLogistic);

	r9 <- rank(preds$pOp) > .95*length(preds$pOp);
	r1 <- rank(preds$pOp) < .05*length(preds$pOp);

	out <- list();
	out$dHighLow <- ks.test(select(orfPairs, !is.na(orfPairs) & highOp),
				select(orfPairs, !is.na(orfPairs) & lowOp ))$statistic;
	out$d20 <- ks.test(select(orfPairs, !is.na(orfPairs) & r9),
				select(orfPairs, !is.na(orfPairs) & r1 ))$statistic;
	out$dOpNot <- ks.test(select(orfPairs, !is.na(orfPairs) & preds$bOp),
				select(orfPairs, !is.na(orfPairs) & !preds$bOp ))$statistic;
	if (!is.null(aoc)) {
		iS <- findMatch(aoc$predsS$Gene1,preds$Gene1)$index;
		iD <- findMatch(aoc$predsD$Gene1,preds$Gene1)$index;
		out$dExperimental <- ks.test(select(orfPairs[iS], !is.na(orfPairs[iS])),
						select(orfPairs[iD], !is.na(orfPairs[iD])))$statistic;
		out$expectedDRatio <- aoc$expectedDRatio;

	}
	out$accuracyHL <- out$dOpNot/out$dHighLow;
	out$accuracy20 <- out$dOpNot/out$d20;
	out$accuracyExperimental <- out$dOpNot/out$dExperimental;
	return(data.frame(out));
}

# For each single microarray (given by dims), remove pairs where neither gene has changed
# (|Delta| < .3).
SingleUarrayOpTest <- function(preds, orfData, dims, aoc=NULL, plot=TRUE, name="orfId", minD=.3, minN=100) {
	if (plot && length(dims)>1) {
		n <- ceiling(sqrt(length(dims)));
		oldpar <- par(no.readonly=TRUE); on.exit(par(oldpar));
		par(mfrow=c(n,n),mar=c(1,1,1,1)+.15,ps=8);
	}
	known <- if(is.null(aoc)) 0 else aoc$known;

	means <- apply(orfData[,dims],2,function(x) { return(mean(withoutNA(x))); });
	swept <- sweep(orfData[,dims],2,means);

	i1 <- findMatch(preds$Gene1, orfData[,name])$index;
	i2 <- findMatch(preds$Gene2, orfData[,name])$index;
	i12 <- select(cbind(preds,i1=i1,i2=i2,known=known), i1>0 & i2>0);

	out <- NULL;
	for (i in 1:length(dims)) {
		dim <- dims[i];
		withData <- select(i12, !is.na(swept[i12$i1,i]) & !is.na(swept[i12$i2,i])
					& (abs(swept[i12$i1,i]) >= minD | abs(swept[i12$i2,i]) >= minD));
		if (is.null(withData) || nrow(withData) < minN) {
			cat("Skipping dimension",dim,"data",nrow(withData),"\n");
			if (plot && length(dims)>1) plot(0:0,0:0,pch=0);
		} else {
			delta <- cbind(withData,
					v1=swept[withData$i1,i],
					v2=swept[withData$i2,i],
					delta=abs(swept[withData$i1,i]-swept[withData$i2,i]) );
			if (plot) {
				CompareHistogramList(split(delta$delta,delta$bOp),breaks=c(0:10)/5,
							legendX=-100,pch=30,ylim=c(0,.5),lty=1,col=c(2,3));
			}

			ksAll <- ks.test(select(delta$delta,delta$bOp), select(delta$delta,!delta$bOp));

			knownOp <- select(delta$delta,delta$known==1);
			knownNop <- select(delta$delta,delta$known==-1);
			useK <- length(knownOp)>5 & length(knownNop)>5;
			ksD <- NULL;
			if (useK) ksD <- ks.test(select(delta$delta,delta$known==1),
							select(delta$delta,delta$known==-1));
			row <- data.frame(dim=dim,rows=nrow(withData),
					nOp=countif(delta$bOp),
					nNop=countif(!delta$bOp),
					nKnownOp=countif(delta$known==1),
					nKnownNop=countif(delta$known==-1),
					medianOp=median(select(delta$delta,delta$bOp)),
					medianNop=median(select(delta$delta,!delta$bOp)),
					medianKnown=if (useK) median(select(delta$delta,delta$known==1)) else NA,
					medianKnownNop=if (useK) median(select(delta$delta,delta$known==-1)) else NA,
					allD=ksAll$statistic, knownK=if(useK) ksD$statistic else NA,
					allR=cor(rank(delta$delta),rank(delta$pOp)),
					knownR=if(useK) cor(select(delta$delta,delta$known!=0),
								select(delta$known,delta$known!=0)) else NA);
			if(is.null(out)) { out <- row } else { out[nrow(out)+1,] <- row }
		}
	}
	return(out);
}

simulateOperonStrand <- function(fPnop=.7, fOpP=.6, fOpM=.5, n=1e3) {
	opR <- runif(n);
	bOpP <- runif(length(opR)) < fPnop; # used if not an operon
	strand <- op <- opR > 0; # vectors of desired length
	# the model is, we add a new gene, either to an operon or not
	# op[1] is true if 1->2 is an operon
	# same[1] is true if 1 and 2 are on the same strand
	for (i in 2:length(strand)) {
		op[i-1] <- if(strand[i-1]) opR[i] < fOpP else opR[i] < fOpM;
		strand[i] <- if (op[i-1]) strand[i-1] else bOpP[i];
	}
	strandNext <- strand[c(2:length(strand),1)];
	same <- strand == strandNext;
	return(data.frame(same=same,op=op,strand=strand,strandNext=strandNext));
}

summarizeSOPS <- function(sops) {
	return(data.frame(
		FPNO = countPerIf(select(sops$strand,!sops$op)),
		OP = countPerIf(select(sops$op,sops$strand)),
		OM = countPerIf(select(sops$op,!sops$strand)),
		FPlus = countPerIf(sops$strand),
		FPP = countPerIf(select(sops$same,sops$strand)),
		FMM = countPerIf(select(sops$same,!sops$strand)),
		OPP = countPerIf(select(sops$op,sops$same&sops$strand)),
		OMM = countPerIf(select(sops$op,sops$same&!sops$strand)),
		fOperon = OperonFreqStrandBias( countPerIf(select(sops$same,sops$strand)),
						countPerIf(select(sops$same,!sops$strand)) )
		));
}

# given P(+2|+1) and P(-2|-1), returns P(O|+1,+2)== P(O|-1,-2)
# under the assumption that they are equal...
OperonFreqStrandBias <- function(pSamePlus,pSameMinus) {
	a <- pSamePlus/pSameMinus;
	b <- -2*pSamePlus;
	c <- pSamePlus+pSameMinus-1;
	pOgivenMinus <- ifelse(b**2 - 4*a*c < 0, 0, -b - sqrt(b**2 - 4*a*c))/(2*a);
	return( pOgivenMinus/pSameMinus );
}


# Given features v1, v2, and known correlations with the dependent variable
# r1true, r2true, and conditional correlation cond12 (given true or false),
# simulate a data set and run our fancy BBF fit on it.
# Assume half of the data set is T for simplicity...
SimulateTraining <- function(r1true=.3,r2true=.3,cond12=0,n=2000,plot=TRUE,...) {
	truth <- ((1:n)%%2)==0;
	v1r <- r1true*ifelse(truth,1,0)/sd(truth) + rnorm(sd=sqrt(1-r1true**2),n=n);
	v1 <- rank(v1r);
	pv1T <- dnorm(v1r,mean=r1true/sd(truth), sd=sqrt(1-r1true**2));
	pv1F <- dnorm(v1r,mean=0,sd=sqrt(1-r1true**2));
	p1 <- pv1T/(pv1T+pv1F);
	v2r <- r2true*ifelse(truth,1,0)/sd(truth) + cond12*v1/sd(v1) +
			rnorm(sd=sqrt(1-r2true**2-cond12**2),n=n);
	v2 <- rank(v2r);
	pv2T <- dnorm(v2r,mean=r1true/sd(truth), sd=sqrt(1-r2true**2));
	pv2F <- dnorm(v2r,mean=0,sd=sqrt(1-r2true**2));
	p2 <- pv2T/(pv2T+pv2F);

	# put perfect error estimates into bbf
	high <- (v1 > n/2); # split for bbf
	fOgH <- countPerIf(truth,high);
	fOgL <- countPerIf(truth,!high);
	fHigh <- countPerIf(high);

	a <- (1-fOgH)*fHigh;
	b <- (fOgH-(1-fOgH))*fHigh;
	c <- a + (1-fOgL)*(1-fHigh);
	d <- b + (fOgL-(1-fOgL))*(1-fHigh);

	bbf <- BinnedBinaryFit(v2, high, a=a,b=b,c=c,d=d,plot=plot,smooth=FALSE,...);
	p2bbfLogOdds <- logodds(bbf$p);
	loessR <- loess(p2bbfLogOdds~rank(v2));
	p2loess <- logodds2p(predict(loessR));
	lines(v2[order(v2)], p2loess[order(v2)]);
	return(list(bbf=bbf,gc=data.frame(truth=ifelse(truth,1,0),v1=v1,v2=v2,
					v12=rank(r1true*v1+r2true*v2),
					high=ifelse(high,1,0),p1=p1,p2=p2,p12=combineP(p1,p2),
					p2bbf=bbf$p,p12bbf=combineP(p1,bbf$p),
					p2loess=p2loess,p12loess=combineP(p1,p2loess) )));
}

SimulateTrainings <- function(tries=30,...) {
	out <- NULL;
	for(i in 1:tries) {
		d <- SimulateTraining(...);
		row <- cor(d$gc)["truth",c("v1","v2","v12","p1","p2","p12","p2bbf","p12bbf","p2loess","p12loess")];
		row <- as.data.frame(as.list(row));
		if (i==1) { out <- data.frame(row); } else { out <- rbind(out,data.frame(row)); }
	}
	return(out);

}

# adds Log2Sum, Log2Diff, and Log2Norm (of 2/1)
# Use rename to ensure that the input contains Fg1, Bg1, Fg2, Bg2, and either sector or sectorX,sectorY
UarrayNormalize <- function(data, minOverBackground=2.0, rename=NULL, plot=TRUE) {
	if (!is.null(rename)) {
		data <- data.frame(RenameColumns(data, rename));
	}
	if (is.null(data$sector)) { data$sector <- interaction(data$sectorX,data$sectorY); }
	data$Log2Sum <- data$Log2Diff <- data$Log2Norm <- rep(NA,nrow(data));

	keep <- select(cbind(data,row=1:nrow(data)),
			data$Fg1>=minOverBackground*data$Bg1 &
			data$Fg2>=minOverBackground*data$Bg2 &
			!is.na(data$Fg1) & !is.na(data$Bg1) &
			!is.na(data$Fg2) & !is.na(data$Bg2) );
	keep$D1 <- keep$Fg1-keep$Bg1;
	keep$D2 <- keep$Fg2-keep$Bg2;
	keep$Log2Sum <- log(keep$D1)+log(keep$D2);
	keep$Log2Diff <- log(keep$D2)-log(keep$D1);
	if (plot) {
		plot(keep$Log2Sum,keep$Log2Diff,pch=unclass(keep$sector),col=unclass(keep$sector));
	}
	ks <- split(cbind(keep,row=1:nrow(keep)),keep$sector)
	for (s in names(ks)) {
		modelS <- loess(Log2Diff~Log2Sum,ks[[s]]);
		if (plot) lines(sort(ks[[s]]$Log2Sum),predict(modelS)[order(ks[[s]]$Log2Sum)]);
		ks[[s]]$Log2Norm <- ks[[s]]$Log2Diff-predict(modelS)
	}
	for (s in names(ks)) {
		data$Log2Norm[ks[[s]]$row] <- ks[[s]]$Log2Norm;
	}
	data$Log2Sum[keep$row] <- keep$Log2Sum;
	data$Log2Diff[keep$row] <- keep$Log2Diff;
	return(data);
}

# Use sectorDiv=c(14,15) and sectorOff=c(1,1) to subtract x and y by 1 and divide by 14,15
# to compute sectorX,sectorY
UarrayNormalizeDir <- function(dir, rename, pattern=".VIMSS$", suffix=".Norm",
				sectorOff=NULL, sectorDiv=NULL, plot=TRUE,
				minOverBackground=2.0) {
	files <- list.files(dir, pattern);
	if (plot) {
		n <- ceiling(sqrt(length(files))+1e-4);
		old.par <- par(no.readonly=TRUE); on.exit(par(old.par));
		par(mfrow=c(n,n),mar=c(1,1,1,1));
	}
	for (file in files) {
		raw <- read.delim(paste(dir,"/",file,sep=""), comment="");
		if(!is.null(rename)) raw <- data.frame(RenameColumns(raw,rename));
		if (!is.null(sectorDiv)) {
			raw$sectorX <- (raw$x-sectorOff[1]) %/% sectorDiv[1];
			raw$sectorY <- (raw$y-sectorOff[2]) %/% sectorDiv[2];
		} else { # in case no position info available
			if (is.null(raw$sectorX)) raw$sectorX <- 0;
			if (is.null(raw$sectorY)) raw$sectorY <- 0;
		}
		cat("Processing ",dir,"/",file,"\n",sep="");
		norm <- UarrayNormalize(raw, rename=NULL, plot=plot, minOverBackground=minOverBackground);
		norm$Ch1 <- norm$Fg1-norm$Bg1;
		norm$Ch2 <- norm$Fg2-norm$Bg2;
		writeDelim(norm,paste(dir,"/",file,suffix,sep=""));
		cat("Wrote ",dir,"/",file,suffix,"\n",sep="");
	}
}

ShowCorDensities <- function(op,cors,name="orfPairs",aoc=NULL,from=-1,to=1,acc=.9,ylim=c(0,3)) {
	both <- cbind(op$preds,select(cors,op$gcSC$same));
	if(!is.null(aoc)) { both$known <- aoc$known; } else { both$known <- 0; }
	data <- list(predOp=select(both,both$bOp & !is.na(both$orfPairs)),
			predNot=select(both,!both$bOp & !is.na(both$orfPairs)),
			control=select(cbind(op$gcSC,cors),!op$gcSC$same & !is.na(cors$orfPairs)));
	if (!is.null(aoc)) {
		data$knownOp <- select(both,both$known==1 & !is.na(both$orfPairs));
		data$knownNot <- select(both,both$known==-1 & !is.na(both$orfPairs));
	}
	CompareDensities(lapply(data,function(x){x[[name]]}),xlim=c(from,to),ylim=ylim);
	d <- lapply(data,function(x){density(x[[name]],from=from,to=to)});
	if (!is.null(aoc)) {
		acc2 <- 1-acc;
		lines(d$predOp$x,acc*d$knownOp$y+acc2*d$knownNot$y,col="grey",lty=2,lwd=2);
		lines(d$predOp$x,acc2*d$predOp$y+acc*d$control$y,col="grey",lty=2,lwd=2);
	}
	return(data);
}

fOperonVsUarray <- function(op,cors,fOperon=op$fOperon,n=NULL,useDensity=NULL) {
	if(is.null(n)) n <- if(is.null(op$sumD1)) 1 else 4;
	if(is.null(useDensity)) useDensity <- TRUE;

	P2 <- cbind(op$preds,select(cors,op$gcSC$same));
	P2 <- P2[!is.na(P2$orfPairs),];
	control <- select(cors,!op$gcSC$same & !is.na(cors$orfPairs));
	high <- select(P2,P2$pOp>.95);
	breaks <- quantile(select(cors$sumD1,!is.na(cors$orfPairs)),probs=seq(0,1,1/n));

	if (max(breaks)-min(breaks) > 0) {
		cutH <- cut(high$sumD1,breaks,labels=FALSE);
		cutC <- cut(control$sumD1,breaks,labels=FALSE);
		cutP2 <- cut(P2$sumD1,breaks,labels=FALSE);
		countsAll <- table(cut(P2$sumD1,breaks,labels=FALSE));
		countsOp <- table(cut(select(P2$sumD1,P2$bOp),breaks,labels=FALSE));
		countsNop <- table(cut(select(P2$sumD1,!P2$bOp),breaks,labels=FALSE));
	} else {
		cutH <- rep(1,nrow(high));
		cutC <- rep(1,nrow(control));
		cutP2 <- rep(1,nrow(P2));
		countsAll <- nrow(P2);
		countsOp <- countif(P2$bOp);
		countsNop <- countif(!P2$bOp);
	}

	histFunc <- function(x,from=-1,to=1,by=.1) {
		h <- hist(x,breaks=seq(from,to,by),plot=FALSE);
		return(list(x=h$mids,y=h$counts));
	}

	densFunc <- if(useDensity) density else histFunc;
	highModel <- lapply(split(high$orfPairs,cutH),densFunc,from=-1,to=1);
	controlModel <- lapply(split(control$orfPairs,cutC),densFunc,from=-1,to=1);
	sameModel <- densFunc(P2$orfPairs,from=-1,to=1);

	sameModel$y <- sameModel$y/sum(sameModel$y);
	opModel <- densFunc(select(P2$orfPairs,P2$bOp),from=-1,to=1);
	opModel$y <- opModel$y/sum(opModel$y);
	nopModel <- densFunc(select(P2$orfPairs,!P2$bOp),from=-1,to=1);
	nopModel$y <- nopModel$y/sum(nopModel$y);

	x <- highModel[[1]]$x;
	allR <- highR <- controlR <- highO <- controlO <- highN <- controlN <- 0*highModel[[1]]$x;
	for (i in 1:length(highModel)) {
		highModel[[i]]$y <- highModel[[i]]$y/sum(highModel[[i]]$y);
		controlModel[[i]]$y <- controlModel[[i]]$y/sum(controlModel[[i]]$y);

		allR <- allR + (countsAll[i]/sum(countsAll))*
			(fOperon*highModel[[i]]$y + (1-fOperon)*controlModel[[i]]$y);
		highR <- highR + (countsAll[i]/sum(countsAll))*highModel[[i]]$y;
		controlR <- controlR + (countsAll[i]/sum(countsAll))*controlModel[[i]]$y;

		highO <- highO + (countsOp[i]/sum(countsOp))*highModel[[i]]$y;
		controlO <- controlO + (countsOp[i]/sum(countsOp))*controlModel[[i]]$y;

		highN <- highN + (countsNop[i]/sum(countsNop))*highModel[[i]]$y;
		controlN <- controlN + (countsNop[i]/sum(countsNop))*controlModel[[i]]$y;
	}
	data <- list(x=x,mix=allR,same=sameModel$y,highR=highR,controlR=controlR,
			op=opModel$y, nop=nopModel$y,
			highO=highO, controlO=controlO,
			highN=highN, controlN=controlN,
			high1=highModel[[1]]$y,control1=controlModel[[1]]$y);
	if (max(breaks)-min(breaks)>0 && n == 4) {
		data$high2 <- highModel[[2]]$y;
		data$control2 <- controlModel[[2]]$y;
		data$high3 <- highModel[[3]]$y;
		data$control3 <- controlModel[[3]]$y;
    		data$high4 <- highModel[[4]]$y;
		data$control4 <- controlModel[[4]]$y;
	}
	return(data.frame(data));
}

LoessAndCut <- function(x,y,breaks=8, rankY=TRUE, scatter=TRUE,
				scatterCol=if(scatter) "grey" else "white",
				loessCol="red", loessLty=1, loessLwd=1,
				varwidth=FALSE,
				boxwex=10,
				span=.5, xlab="", ylab="", main="", notch=F, ...) {
	x1 <- select(x, !is.na(x) & !is.na(y));
	y1 <- select(y, !is.na(x) & !is.na(y));
	x <- x1; y <- y1;
	if (length(breaks)==1) { breaks <- unique(quantile(x,seq(0,1,1/breaks))); }
	fx <- cut(x,breaks,include.lowest=TRUE);
	model <- if (rankY) loess(rank(y)~rank(x))  else loess(y~rank(x), span=span);
	plot(x,y,col=scatterCol,xlab=xlab,ylab=ylab,main=main,...);
	boxplot(y~fx,add=TRUE,at=sapply(split(x,fx),median),
			width=if(varwidth) sapply(split(x,fx),IQR) else rep(1,nlevels(fx)),
			boxwex=boxwex,
			xaxt="n", yaxt="n", notch=notch);
	sortedLines(x,if(rankY) sort(y)[pmax(1,pmin(length(y),model$fitted))] else model$fitted,
			col=loessCol, lty=loessLty, lwd=loessLwd);
}

# Histogram of best hits by what taxon they are from
# Reads files created by e.g.:
# mysql -h lut -u test -ptest genomics_test -B -e 'select s1.taxonomyId tax1, s1.scaffoldId s1, l1.locusId locusId1, s2.taxonomyId tax2, s2.scaffoldId s2, l2.locusId locusId2, b.score score from BLASTp b, Locus l1, Locus l2, Scaffold s1, Scaffold s2 WHERE l1.scaffoldId=s1.scaffoldId AND l2.scaffoldId=s2.scaffoldId AND s1.isActive=1 AND s2.isActive=1 AND s1.taxonomyId<>s2.taxonomyId AND b.locusId=l1.locusId AND b.subject=l2.locusId AND l1.priority=1 AND l2.priority=1 GROUP BY l1.locusId HAVING score=max(score);' > ~/Jamb/besthits
# mysql -h lut -u test -ptest genomics_test -B -e 'select taxonomyId tax, shortName short, name name from Taxonomy order by taxonomyId' > ~/Jamb/genomeNames
#
BestHitsSummary <- function(list=NULL, besthits = NULL, genomes = NULL, dir="Jamb",minCount=10) {
	if (is.null(besthits)) besthits <- read.delim(paste(dir,"/","besthits",sep=""));
	if (is.null(genomes)) genomes <- read.delim(paste(dir,"/","genomeNames",sep=""));
	if (is.null(list)) list <- genomes$tax;

	bh <- as.data.frame.table(table(besthits$tax1, besthits$tax2));
	names(bh) <- c("tax1","tax2","n");
	bh <- merge(bh,genomes,by.x="tax1",by.y="tax");
	bh <- merge(bh,genomes,by.x="tax2",by.y="tax"); # short.x is source, short.y is target
	bh <- bh[order(bh$place.x,-bh$n),]
	bh <- select(bh,bh$tax1 %in% list & bh$n>=minCount);
	return(data.frame(sourceTax=bh$tax1,sourceName=bh$short.x,Count=bh$n,bestTax=bh$tax2,bestName=bh$short.y));
}

GeneSkewPlot <-  function(taxon=NULL, dir="Jamb", GC=read.delim(paste(dir,"/",taxon,".cont",sep=""))) {
	GC$GeneSkew <- GC$GenesPlus-GC$GenesMinus;
	GC$GCSkewC <- cumsum(GC$GCSkew);
	GC$GeneSkewC <- cumsum(GC$GeneSkew)
	GC$GCProfile <- cumsum(GC$GC-mean(GC$GC));
	plotN(Begin~GCProfile+GeneSkewC+GCSkewC,GC,c(2,2),type="l");
}

CumGCProfile <- function(taxon=207559, dir="Jamb", GC=read.delim(paste(dir,"/",taxon,".cont",sep=""))) {
	GC$GCProfile <- cumsum(GC$GC-mean(GC$GC));
	return(GC);
}

COGFuncPlot <- function(list,
			dir="Jamb",
			l2cft = read.delim(paste(dir,"/","Locus2COGFunTable",sep="")),
			COGFun = read.delim(paste(dir,"/","COGFun",sep="")),
			genomes=read.delim(paste(dir,"/","genomeNames",sep="")),
			pch = 1:length(list), col=1+(1:length(list)),
			lty=1:length(list), lwd=rep(1,length(list)),
			log="", off=ifelse(log=="y",1e-4,0),
			showFun=(log!="y") ) {
	l2cfTable <- matrix(l2cft$Freq, ncol=length(levels(l2cft$Var1)), byrow=TRUE);
	colnames(l2cfTable) <- levels(l2cft$Var1);
	rownames(l2cfTable) <- l2cft$Var2[1 + (0:(nrow(l2cfTable)-1))*ncol(l2cfTable)];
	l2cfTableNorm <- sweep(l2cfTable, 1, apply(l2cfTable,1,sum), "/");
	plot(as.factor(colnames(l2cfTableNorm)[col(l2cfTableNorm)]),off+c(l2cfTableNorm),log=log)
	names <- c();
	for (i in 1:length(list)) {
		genome <- list[i];
		name <- toString( select(genomes$short,genomes$tax==genome) );
		names[length(names)+1] <- name;
		points(1:ncol(l2cfTableNorm), off+l2cfTableNorm[toString(genome),], pch=pch[i], col=col[i]);
		lines(1:ncol(l2cfTableNorm), off+l2cfTableNorm[toString(genome),], lty=lty[i], col=col[i], lwd=lwd[i]);
	}
	legend(ifelse(log=="y",3,24.5), ifelse(log=="y",off,.26), names, col=col, pch=pch, lty=lty, lwd=lwd,
			xjust=ifelse(log=="y",0,1),
			yjust=ifelse(log=="y",0,1));
	if(showFun) {
		old.par <- par(no.readonly=TRUE); on.exit(par(old.par));
		legend(.5,.265,
				levels(COGFun$description)[COGFun$description],
				pch=levels(COGFun$funCode)[COGFun$funCode],
				yjust=1,xjust=0, cex=0.8, bty="n");
	}
}

# returns a list of genome, acc, n, percent, percentile, name
# doesn't report percentiles for count==0
MakeGoPercentiles <- function(dir="Jamb",
			# taxId, acc (GO:NNNNN), n (count), name (string)
			go=read.delim(paste(dir,"/","goCounts",sep="")),
			genomes=read.delim(paste(dir,"/","genomeNames",sep="")))
{
	# get count at top level
	go2 <- merge(go, select(go,go$acc=="GO:0003673"), by.x="taxId", by.y="taxId");
	go2$n <- go2$n.x;
	go2$acc <- go2$acc.x;
	go2$percent <- go2$n.x/go2$n.y;
	go2$percentile <- 0; # fix later
	go2$name <- go2$name.x;
	go2 <- without(go2, c("acc.x","acc.y","n.x","n.y","name.x","name.y"));

	
	for(acc in levels(go2$acc)) {
		percents <- go2$percent[go2$acc==acc];
		if (length(percents) < nrow(genomes)) {
			percents <- c(percents, rep(0, nrow(genomes)-length(percents)));
		}
		go2$percentile[go2$acc==acc] <- findIntervalAvg(go2$percent[go2$acc==acc],sort(percents))/length(percents);
	}
	return(go2);
}

# go2 from MakeGoPercentiles
MakeGoZero <- function(go2, percentile=.1) {
	ngenomes <- length(unique(sort(go2$taxId)));
	rankCutoff <- percentile*ngenomes;
	counts <- table(go2$acc);
	go2High <- go2[counts[unclass(go2$acc)] >= ngenomes*(1-percentile) & counts[unclass(go2$acc)] != ngenomes,];
	go2Split <- split(go2High,go2High$acc);
	countsH <- table(go2High$acc);
	countsH <- countsH[countsH>0];

	indices <- c();
	taxIds <- c();
	names <- c();

	taxa <- unique(sort(go2$taxId));

	for (i in 1:length(countsH)) {
		goAcc <- go2Split[[ names(countsH)[i] ]];
		taxAdd <- setdiff(taxa, goAcc$taxId);
		if (length(taxAdd)>0) {
			taxIds[length(taxIds)+(1:length(taxAdd))] <- taxAdd;
			indices[length(indices)+(1:length(taxAdd))] <- i;
			names[length(names)+(1:length(taxAdd))] <- levels(goAcc$name)[goAcc$name[1]];
			cat("acc",names(countsH)[i],"add",taxAdd,"\n");
		}
	}
	accs <- names(countsH)[indices];
	return(data.frame(taxId=taxIds, n=rep(0,length(taxIds)), acc=accs, percent=rep(0,length(taxIds)),
				percentile=rep(0,length(taxIds)), name=names));
}

MakeGoReport <- function(go2, go2z,
			extraacc=NULL, # list of GOs, as vector of strings NOT as a factor
			taxa=c(211586, 881, 207559,891,28232,243231),
			taxShort=c("Sone", "Dvul","Ddes","Dace","Gmet","Gsul"),
			extrataxa=c(83333, 1423),
			extraShort=c("Ecol","Bsub"),
			percentile=.05,
			go2term=read.delim("~/Jamb/go2term"),
			file="~/Jamb/GoReport.html") {
	# first, the set of accs
	go2Taxa <- go2[go2$taxId %in% taxa,];
	accGo <- levels(go2Taxa$acc)[go2Taxa$acc[go2Taxa$percentile<percentile | go2Taxa$percentile>1-percentile]];
	accGoZ <- levels(go2z$acc)[go2z$acc[go2z$taxId %in% taxa]];
	acc <- union(union(sort(unique(accGo)),sort(unique(accGoZ))), extraacc);
	sumAcc <- sapply(split(go2Taxa$n,go2Taxa$acc),sum);
	accOrdered <- acc[order(-sumAcc[acc])];

	if (!is.null(file)) { sink(file); on.exit(sink(NULL)); }

	alltax <- c(taxa,extrataxa);
	allshort <- c(taxShort,extraShort);
	cat("<HTML><HEAD><TITLE>GO analysis for ",toString(taxShort)," </TITLE></HEAD>\n",sep="");
	cat("<BODY>\n");
	cat("<H2>GO analysis for ",toString(taxShort)," </H2>\n",sep="");
	cat("<P>For each outlier category, and for each genome, shows the count and the percentile of the percent of all\n");
	cat("genes assigned to GO categories that hit this category\n");
	cat("<TABLE cellpadding=2 >\n");
	cat("<TR bgcolor=#A0A0FF>\n");
	for (i in 1:length(alltax)) { cat("  <TD>",allshort[i],"</TD>\n"); }
	cat("  <TD>Category</TD>\n");
	cat("  <TD>Comment</TD>\n");
	cat("</TR>\n");
	for(i in 1:length(accOrdered)) {
		acc <- accOrdered[i];
		term <- go2term$id[go2term$acc==acc];
		go2Acc <- go2[go2$acc==acc,];
		cat("<TR bgcolor=", ifelse(i %% 2  == 0, "#FFFFF", "#F0F0A0"), ">\n", sep="");
		for (i in 1:length(alltax)) {
			count <- go2Acc$n[go2Acc$taxId==alltax[i]];
			percentile <- go2Acc$percentile[go2Acc$taxId==alltax[i]];
			if(length(count)==0) { count <- 0; percentile <- 0; }
			cat("  <TD valign=top align=center><A HREF=\"http://escalante.lbl.gov/cgi-bin/vertiGo.cgi?taxId=",
					alltax[i],"&term=",term,"&getGenes=1&goBrowse=0\">",
				count,"</A><BR><small><em>",format(100*percentile,digits=2),
				"</em></small></TD>\n", sep="");
		}
		cat("  <TD valign=top align=left><A HREF=\"http://escalante.lbl.gov/cgi-bin/vertiGo.cgi?term=",
			term,"&taxId=",paste(alltax,collapse=","),"&goBrowse=0\">",
			levels(go2Acc$name)[go2Acc$name[1]],"</A></TD>\n",sep="");
		cat("  <TD>&nbsp;</TD>"); # comment
		cat("</TR>\n");
	}
	cat("</TR></TABLE></BODY></HTML>\n");
}

# tYes and tAll should be from data.frame(table()), having Var1 and Freq
# allYes and all should be the column that was taken out when doing the table
fisherTests <- function(tYes, tAll, allYes=NULL, all=NULL,
		nAllYes=length(unique(sort(allYes))),
		nAll=length(unique(sort(all))),
		p=0.05, minAll=2, maxAllf=0.2) {
	f <- list();
	f$Var1 <- tYes$Var1;
	f$a <- tYes$Freq;
	f$ac <- tAll$Freq + 0*f$a;
	f$ab <- nAllYes + 0*f$a;
	f$abcd <- nAll + 0*f$a;
	f$b <- f$ab-f$a;
	f$c <- f$ac-f$a;
	f$d <- f$abcd - f$a - f$b - f$c;
	f$p <- 0*f$d;
	f$odds <- 0*f$d;
	f <- data.frame(f);
	f <- f[f$ac >= 2 & f$ac <= maxAllf * f$abcd,];
	for (i in 1:nrow(f)) {
		test <- fisher.test(matrix(c(f$a[i],f$b[i],f$c[i],f$d[i]),nrow=2));
		f$p[i] <- test$p.value;
		f$odds[i] <- test$estimate;
	}
	f <- f[f$p<=p,];
	return(f[order(f$p),]);
}

# given a table that foreach group has a count of items, produce a single table
# The 1st three arguments are the column names in the input table
tableToTable <- function(group,item,count,table) {
	items <- unique(sort(table[[item]]));
	groups <- unique(sort(table[[group]]));
	out <- data.frame(name=items);
	names(out) <- item;
	for (g in groups) {
		colName <- paste("count",g,sep="");
		out[[colName]] <- rep(0,nrow(out));
		merged <- merge(out[,c(item,item)], table[table[[group]]==g,c(item,count)], by=item, all=TRUE, all.y=FALSE);
		out[[colName]] <- ifelse(is.na(merged[[count]]), 0, merged[[count]]);
	}
	return(out);
}

OperonStatus <- function(data, op, name="locusId") {
	if (length(unique(sort(data[[name]]))) != nrow(data)) {
		warning("duplicates not allowed");
		return(NULL);
	}
	reorder <- rank(data[[name]]);
	gc1 <- merge(data.frame(Gene1=data[[name]]), op$gcSC, sort=TRUE);
	gc1 <- gc1[reorder,];
	gc2 <- merge(data.frame(Gene2=data[[name]]), op$gcSC, sort=TRUE);
	gc2 <- gc2[reorder,];
	if (nrow(gc1) != nrow(data)) { warning("invalid locusId"); return(NULL); }
	if (nrow(gc2) != nrow(data)) { warning("invalid locusId"); return(NULL); }
	if (countif(gc1$Gene1 != data[[name]]) > 0) { warning("sort failed"); }
	if (countif(gc2$Gene2 != data[[name]]) > 0) { warning("sort failed"); }
	
	data$sameUp <- ifelse(gc2$Strand2=="+",gc2$same,gc1$same);
	data$sameDown <- ifelse(gc1$Strand1=="+",gc1$same,gc2$same);
	data$singleton <- (!data$sameUp) & (!data$sameDown); # is it in a directon of size 1
	data$SepUp <- ifelse(gc2$Strand2=="+", gc2$Sep, gc1$Sep); # size of upstream region

	data$opSize <- rep(1, nrow(data)); # 1 will mean not in an operon
	lcCont <- data[[name]];
	opCont <- rep(TRUE,nrow(data));
	# first, from Gene1 -> Gene2
	while(sum(ifelse(opCont,1,0)) > 0) {
		preds <- merge(data.frame(Gene1=lcCont,SortBy=data[[name]]), op$preds, sort=TRUE, all.x=TRUE);
		preds <- preds[order(preds$SortBy),];
		preds <- preds[reorder,];
		if (countif(preds$SortBy != data[[name]]) > 0) { warning("sort failed"); }
		lcNext <- preds$Gene2;
		opCont <- opCont & !is.na(lcNext) & preds$bOp;
		data$opSize <- data$opSize + ifelse(opCont,1,0);
		lcCont <- lcNext;
	}
	# then, from Gene2 -> Gene1
	lcCont <- data[[name]];
	opCont <- rep(TRUE,nrow(data));
	while(sum(ifelse(opCont,1,0)) > 0) {
		preds <- merge(data.frame(Gene2=lcCont,SortBy=data[[name]]), op$preds, sort=TRUE, all.x=TRUE);
		preds <- preds[order(preds$SortBy),];
		preds <- preds[reorder,];
		if (countif(preds$SortBy != data[[name]]) > 0) { warning("sort failed"); }
		lcNext <- preds$Gene1;
		opCont <- opCont & !is.na(lcNext) & preds$bOp;
		data$opSize <- data$opSize + ifelse(opCont,1,0);
		lcCont <- lcNext;
	}
	# see if gene is at start of predicted operon (including singletons)
	preds1 <- merge(data.frame(Gene1=data[[name]]), op$preds, all.x=TRUE, sort=TRUE)[reorder,];
	preds2 <- merge(data.frame(Gene2=data[[name]]), op$preds, all.x=TRUE, sort=TRUE)[reorder,];
	data$opFirst <- ifelse(gc1$Strand1=="+", ifelse(gc2$same,!preds2$bOp,TRUE),
						 ifelse(gc1$same,!preds1$bOp,TRUE));
	data$opLast <- ifelse(gc1$Strand1=="+", ifelse(gc1$same,!preds1$bOp,TRUE),
						 ifelse(gc2$same,!preds2$bOp,TRUE));
	return(data);
}

# read operons data including predictions
ReadOp <- function(row,gncDir="operons/out",prefix="gnc",useCAI=FALSE,...) {
	op <- as.list(row);
	op$gcSC <- read.delim(paste(gncDir,"/",prefix,op$genome,".gc",sep=""),comment.char="",quote="");
	op$preds <- read.delim(paste(gncDir,"/",prefix,op$genome,".pred",sep=""),comment.char="",quote="");
	op$bbfD <- list();
	op$bbfD$groups <- read.delim(paste(gncDir,"/",prefix,op$genome,".sep",sep=""),comment.char="",quote="");
	if(useCAI) {
		op$bbfCAI <- list();
		op$bbfCAI$groups <- read.delim(paste(gncDir,"/",prefix,op$genome,".caidxdy",sep=""),
					comment.char="",quote="");
	}
	op$without <- c("without","gcSC","preds","bbfD","bbfCAI","vars");
	op$gcSC$same <- op$gcSC$same != "FALSE";
        op$preds$bOp <- op$preds$bOp != "FALSE";
        op$preds$bOp1 <- op$preds$bOp1 != "FALSE";
	op$preds$forward <- op$preds$forward != "FALSE";
	op$gcSC$forward <- op$gcSC$forward != "FALSE";
	if(!is.null(op$formula)) { op$formula <- levels(op$formula)[op$formula]; }
	op$vars <- c("GNScore","GNMinus","GNWithin","GNAll","HK","MI","CAIdxdyR", "CAIdsqrR");
 	return(op);
}

testGroupVsOperons <- function(experiment, condition, operons, group,
				minAbsScore=0.5, locusCol="VIMSS") {
	op2 <- operons[operons$locusId %in% group,];
	op3 <- merge(op2, operons, by="operon");
	# don't include self, but other group members are OK
	# do not count operon twice if two genes are in the list
	# (but both genes will be counted)
	genes <- unique(op3$locusId.y[op3$locusId.y != op3$locusId.x]);
	test <- merge(data.frame(locusId=genes),
			data.frame(locusId=experiment[[locusCol]],score=experiment[[condition]]));
	test <- test[!is.na(test$score) & abs(test$score) >= minAbsScore,];
	return(list(n=nrow(test), nUp=countif(test$score>0)));
}

testGroupsVsOperons <- function(experiment, condition, operons,
				minAbsScore=0.5, locusCol="VIMSS",
				up=TRUE, nChange=1000,
				changers=experiment[[locusCol]][order(ifelse(up,-1,1)*experiment[[condition]])][1:nChange],
				groupSize=50,
				conf.int=95,
				debug=FALSE) {
	changers <- intersect(changers, experiment[[locusCol]][ experiment[[condition]] * ifelse(up,1,-1) > 0 ]);
	if(debug) cat("condition",condition,"nchangers",length(changers),"\n");
	ngroups <- floor(length(changers)/groupSize);
	out <- data.frame(nFirst=groupSize*(1:ngroups)-groupSize+1,
				nLast=groupSize*(1:ngroups),
				nTest=rep(0,ngroups),
				nConfirm=rep(0,ngroups),
				fConfirm=rep(0,ngroups),
				conf.low=rep(0,ngroups),
				conf.high=rep(0,ngroups),
				p=rep(0,ngroups));

	for (i in 1:ngroups) {
		group <- changers[out$nFirst[i]:out$nLast[i]];
		test <- testGroupVsOperons(experiment, condition, operons, group,
						minAbsScore=minAbsScore, locusCol=locusCol);
		out$nTest[i] <- test$n;
		out$nConfirm[i] <- if (up) test$nUp else test$n - test$nUp;
		bTest <- binom.test(out$nConfirm[i], out$nTest[i]);
		out$fConfirm[i] <- bTest$estimate;
		out$conf.low[i] <- bTest$conf.int[1];
		out$conf.high[i] <- bTest$conf.int[2];
		out$p[i] <- bTest$p.value;
	}
	return(out);
}

plotChangersVsOperons <- function(experiment, # all the data
				 conditions, # vector of conditions (column names in experiment)
				operons, # data frame of operon and locusId
				showByScore=FALSE,
				# for other genes in operon, not for changers themselves
				minAbsScore=rep(0.5,length(conditions)),
				locusCol="VIMSS", # in experiment
				nChange=1000,
				groupSize=50,
				conf.int=95,
				xlab=if(showByScore) "Z-Score" else "# Changers",
				ylab="% True Changers",
				main="Operon-Based Estimates of Local Accuracy",
				showConf=FALSE,
				randomControl=FALSE,
				labels=conditions,
				show="both",
				down=(show=="both" || show=="down"),
				up=(show=="both" || show=="up"),
				nplots=length(conditions) * if(down && up && !showByScore) 2 else 1,
				pch=1:nplots, lty=1:nplots, col=1:nplots,
				ylim=c(0,1), 
				xlim=NULL,
				legend = TRUE,
				mfrow=NULL, # set to e.g. c(3,4) to get 1 panel per condition
				mainSep="\n", # separator b/w main and condition if mfrow
				debug=FALSE
		) {
	if (!is.null(mfrow)) {
		oldpar <- par(no.readonly=TRUE);
		par(mfrow=mfrow);
		on.exit(par(oldpar));
	}
	if (is.null(xlim)) {
		if (showByScore) {
			values <- sapply(conditions, function(x) {
				s <- sort(withoutNA(experiment[,x]));
				return(c(s[groupSize],s[length(s)-groupSize])) });
			xlim <- c(-3,3);#c(-max(range(values)), max(range(values)));
		} else {
			xlim <- (0:1)*nChange;
		}
	}

	for (i in 1:length(conditions)) {
		if (i==1 || !is.null(mfrow)) {
			plot( xlim, ylim, col="white", xlab=xlab, ylab=ylab,
				main = if(is.null(mfrow)) main else paste(main,mainSep,labels[i],sep=""),
				ylim=ylim);
			legends <- c();
		}

		if (showByScore) legends <- c(legends, labels[i]);
		if (up) {
			if (!showByScore) legends <- c(legends, paste(labels[i], "(up)"));
			nplot <- length(legends);
			if(debug) cat("Running up for",conditions[i],"\n");
			testUp <- testGroupsVsOperons(experiment, conditions[i], operons,
							minAbsScore=minAbsScore[i], locusCol=locusCol,
							up=TRUE,
							nChange=nChange, groupSize=groupSize,
							conf.int=conf.int,debug=debug);
			if(debug) cat("Ran\n");
			x <- testUp$nLast;
			if (showByScore) {
				sorted <- -sort(withoutNA(-experiment[[conditions[i]]]));
				x <- sorted[testUp$nLast];
			}
			points(x, 2*testUp$fConfirm-1, pch=pch[nplot], col=col[nplot]);
			lines(x, 2*testUp$fConfirm-1, lty=lty[nplot], col=col[nplot]);
			if (showConf) arrows( x, 2*testUp$conf.low-1,
						x, 2*testUp$conf.high-1,
						angle=90, code=3, length=1/12, col="darkgrey");
		}
		if (down) {
			if (!showByScore) legends <- c(legends, paste(labels[i], "(down)"));
			nplot <- length(legends);
			testDown <- testGroupsVsOperons(experiment, conditions[i], operons,
							minAbsScore=minAbsScore[i], locusCol=locusCol,
							up=FALSE,
							nChange=nChange, groupSize=groupSize,
							conf.int=conf.int,debug=debug);
			x <- testDown$nLast - if (showConf) i*2 else 0;
			if (showByScore) {
				sorted <- sort(withoutNA(experiment[[conditions[i]]]));
				x <- sorted[testDown$nLast];
			}
			points(x, 2*testDown$fConfirm-1, pch=pch[nplot], col=col[nplot]);
			lines(x, 2*testDown$fConfirm-1, lty=lty[nplot], col=col[nplot]);
			if (showConf) arrows( x, 2*testDown$conf.low-1,
						x, 2*testDown$conf.high-1,
						angle=90, code=3, length=1/12, col="darkgrey");
		}
		if (i == length(conditions) || !is.null(mfrow)) {
			if(legend) legend(min(xlim), min(ylim), legends,
					col=col[1:length(legends)],
					lty=lty[1:length(legends)],
					pch=pch[1:length(legends)], yjust=0);
			if(randomControl) {
				random <- sample(experiment[[locusCol]], 8*groupSize);
				testRandom <- testGroupVsOperons(experiment, conditions[1], operons, random,
								minAbsScore=minAbsScore[i], locusCol=locusCol);
				hline(abs(2*testRandom$nUp/testRandom$n-1));
			}
		}						
	}
}

# returns a data frame of nPairs, nGenes, statistic, opAcc (the cumulative accuracy estimate), and other stats
# nPairs can be > nGenes because of genes in large operons;
# nPairs can be < nGenes b/c of genes not in operons
# both Gene1 (the gene in the list) and Gene2 (the gene whose data is being tested -- this may
# not be the same order as in the input) are also included.
#
# input (locusId, statistic) need not be sorted
changersVsOperonPairs <- function(locusId, # vector
				statistic, # must be signed (for up vs. down), e.g., zScore
				# must include Gene1, Gene2, and pOp; including only pOp>.5 is recommended
				pairs,
				flipPairs=TRUE, # need to add reverse of pairs (Gene2,Gene1, same pOp)?
				window=100, # but first set of window is cumulative average
				useP=FALSE, # use P(Operon)?
				up=TRUE,
				both=FALSE,
				minAbs=0) { # only for the 2nd gene in the pair, but this can be biasing...
	if (both) {
		out1 <- changersVsOperonPairs(locusId,statistic,pairs,flipPairs=flipPairs,useP=useP,up=TRUE,minAbs=minAbs);
		out2 <- changersVsOperonPairs(locusId,statistic,pairs,flipPairs=flipPairs,useP=useP,up=FALSE,minAbs=minAbs);
		return(rbind(out1,out2));
	}
	if(flipPairs) {
		pairs2 <- rbind(pairs[,c("Gene1","Gene2","pOp")],
				data.frame(Gene1=pairs$Gene2,Gene2=pairs$Gene1,pOp=pairs$pOp));
	} else {
		pairs2 <- pairs[,c("Gene1","Gene2","pOp")];
	}
	# get statistic for both pairs, remove NA and small (unconfident) values, and sort by stat1
	pairs2 <- merge(pairs2,data.frame(Gene1=locusId,stat1=statistic),by="Gene1");
	pairs2 <- merge(pairs2,data.frame(Gene2=locusId,stat2=statistic),by="Gene2");
	pairs2 <- pairs2[!is.na(pairs2$stat1) & !is.na(pairs2$stat2),];
	pairs2 <- pairs2[abs(pairs2$stat2) >= minAbs,];
	sign <- if(up) 1 else -1;
	pairs2 <- pairs2[order(-sign*pairs2$stat1,pairs2$Gene1),];
	n <- countif(sign*pairs2$stat1>0);
	pairs2 <- pairs2[1:n,];

	out <- data.frame(nPairs=1:n);

	# to compute nGenes, need to sort all of them, whether or not in pairs
	genes <- merge(data.frame(Gene1=locusId,stat1=statistic),pairs2[,c("Gene1","Gene2")],
			by="Gene1", all.x=TRUE,all.y=FALSE);
	genes <- genes[!is.na(genes$stat1),];
	genes <- genes[order(-sign*genes$stat1,genes$Gene1),];
	genes <- genes[sign*genes$stat1>0,];
	# counting different values
	nGenes <- c(1,1+cumsum(ifelse(genes$Gene1[1:(nrow(genes)-1)] != genes$Gene1[2:nrow(genes)],1,0)));
	out$nGenes <- nGenes[!is.na(genes$Gene2)]; # only keep rows that were in pairs2

	out$statistic <- pairs2$stat1;

	out$bAgree <- ifelse(sign*pairs2$stat2>0,1,0);
	out$pOp <- pairs2$pOp;
	out$raw <- ifelse(out$bAgree,1,-1)/(if(useP) out$pOp else rep(1,n));
	# smooth raw over window size to give opAcc
	if(n > window) {
		out$opAcc <- cumsum(out$raw);
		out$opAcc[(window+1):n] <- out$opAcc[(window+1):n] - out$opAcc[1:(n-window)];
		out$opAcc <- out$opAcc/pmin(1:n,window);
	} else {
		out$opAcc <- cumsum(out$raw)/(1:n);
	}
	out$cumAcc <- cumsum(out$raw)/(1:n);
	out$Gene1 <- pairs2$Gene1;
	out$Gene2 <- pairs2$Gene2;
	return(out);
}

# returns s0, sigma, effect, and d as described by SAM, where
# d =(data1-data2)/(sigma + s0)
# and s0 is chosen to minimize the effect of sigma on d
# data1 and data2 must have matching rows (e.g. sets of estimated expression levels)
# if data2 is NULL, does a 1-sample computation
#
# Does NOT do normalization or imputation of missing values
flattenSigma <- function(data1, data2=NULL, nPercentileS0Try=100, nQuantileS0Test=20) {
	if (is.null(data2)) {
		effect <- apply(data1,1,mean,na.rm=TRUE);
		sigma <- apply(data1,1,sd,na.rm=TRUE);
	} else {
		effect <- apply(data1,1,mean,na.rm=TRUE)-apply(data2,1,mean,na.rm=TRUE);
		sigma <- sapply(1:nrow(data1), function(i) {
				x1 <- withoutNA(data1[i,]);
				x2 <- withoutNA(data2[i,]);
				# variance of x1-x2 with assumption of equal variance
				sqrt( ( (1/length(x1) + 1/length(x2)) / (length(x1)+length(x2)-2) ) *
					(sum((x1-mean(x1))**2) + sum((x2-mean(x2))**2) ));
			});
	}
	use <- is.finite(sigma) & is.finite(effect);
	sets <- cutN(sigma[use],n = nQuantileS0Test, label=FALSE);
	s0cand <- quantile(sigma[use],(0:nPercentileS0Try)/nPercentileS0Try);
	s0score <- sapply(s0cand, function(s0) {
		d <- (effect/(sigma+s0))[use];
		groups <- split(d,sets);
		mads <- sapply(groups,mad); # median absolute deviation, scaled to give sd if normal
		sd(mads)/mean(mads);
	});
	s0 <- s0cand[s0score==min(s0score)][1];
	return(data.frame(effect=effect,sigma=sigma,s0=rep(s0,length(effect)),d=effect/(sigma+s0)));
}

genes12 <- function(x) { unique(c(t(x[,c("Gene1","Gene2")]))); }
in12 = function(x, list) x$Gene1 %in% list | x$Gene2 %in% list;
both12 = function(x, list) x$Gene1 %in% list & x$Gene2 %in% list;
sub12 = function(x, list) x[x$Gene1 %in% list | x$Gene2 %in% list,];

genes12Groups <- function(x) { rbind(data.frame(Group=1:nrow(x),locusId=x$Gene1),
					data.frame(Group=1:nrow(x),locusId=x$Gene2)); }
# Mantains counts of genes on both sides. Ignores ordering
genes12Controls <- function(x, n=1e4, excludeSet=NULL, excludeSelf=TRUE, excludeInput=TRUE) {
	excludeSet <- unique(excludeSet[,c("Gene1","Gene2")]);
	genes <- c(x$Gene1,x$Gene2);
	out <- data.frame(Gene1=sample(genes,n,replace=TRUE),Gene2=sample(genes,n,replace=TRUE));
	# exclude self pairs
	if(excludeSelf) out <- out[out$Gene1 != out$Gene2,];
	# excludeSet
	if(!is.null(excludeSet)) {
		out <- merge(out,cbind(excludeSet[,c("Gene1","Gene2")],exclude=rep(1,nrow(excludeSet))),
			all.x=TRUE,all.y=FALSE);
		out <- out[is.na(out$exclude), c("Gene1","Gene2")];
		out <- merge(out,data.frame(Gene1=excludeSet$Gene2,
						Gene2=excludeSet$Gene1,
						exclude=rep(1,nrow(excludeSet))),
			all.x=TRUE,all.y=FALSE);
		out <- out[is.na(out$exclude), c("Gene1","Gene2")];
	}
	# excludeInput
	if(excludeInput) {
		out <- merge(out,cbind(x[,c("Gene1","Gene2")],exclude=rep(1,nrow(x))),
			all.x=TRUE,all.y=FALSE);
		out <- out[is.na(out$exclude), c("Gene1","Gene2")];
		out <- merge(out,data.frame(Gene1=x$Gene2,Gene2=x$Gene1,exclude=rep(1,nrow(x))),
			all.x=TRUE,all.y=FALSE);
		out <- out[is.na(out$exclude), c("Gene1","Gene2")];
	}
	return(out);
}

SimulateUarrayVsOperons <- function(locusIds, opPreds, # with Gene1, Gene2, and pOp, or NULL
				# standard model taken from Shewanella salt data
				fit=NULL,
				a = if(is.null(fit)) 0.3 else fit$a,
				bc = if(is.null(fit)) 0.1 else fit$bc,
				c = if(is.null(fit)) 0.28 else fit$c,
				v = if(is.null(fit)) 3.25 else fit$v,
				b = 1/(1/bc-1/c),
				n = 10, # number of replicates (must be constant)
				control=FALSE, # true if have separate control chips
				n2 = n, # number of replicates for control data (if not on same chip)
				pNA = 0.05, # proportion of missing measurements

				# lots of alternatives, defaulting to standard model
				# sd(mu) = (changer? s : s2) + sqrt(sFactor)*sigma
				pChange = 1.0, # by default all genes change
				s = 0.0, # standard deviation of true mu
				s2 = 0.0, # standard deviation of true mu for non-changers
				sFactor = 1/b, # variance ratio for changers
				sFactorNC = sFactor, # variance ratio for non-changers

				# s.d. of replicate measurements = sigmaFactor/sqrt(rgamma(shape,scale))
				sigmaFactor = 1,
				shape=(v+1)/2, scale=2/a,
				sigma = sigmaFactor/sqrt(rgamma(length(locusIds),shape=shape,scale=scale)),
				sdBias = 0, # systematic error (sigma is between replicates)
				scaledBias = 1/c, # systematic error relative to sigma^2 (like sFactor)
				matchSigma=TRUE, # force sigma/theta to match across operon

				sweep=FALSE, # sweep means to be zero

				debug=FALSE,
				seed=NULL,
				analyze=TRUE, # run AnalyzeUarrayVsOperons
				estimate=TRUE, # if false, analyze with true parameters
				... # for AnalyzeUarrayVsOperons
				) {
	if (is.null(seed)) seed <- floor(runif(1)*1e9);
	set.seed(seed);
	changer <- runif(length(locusIds)) < pChange;
	mu <- rnorm(length(locusIds)) * (ifelse(changer, s, s2) + sigma*sqrt(ifelse(changer,sFactor,sFactorNC)));

	# make mu and perhaps sigma agree for all operon pairs
	if (!is.null(opPreds)) {
		oppairs <- merge(opPreds[,c("Gene1","Gene2","pOp")],
					data.frame(Gene1=locusIds,index1=1:length(locusIds)),by="Gene1");
		oppairs <- merge(oppairs,data.frame(Gene2=locusIds,index2=1:length(locusIds)),by="Gene2");
		oppairs$sOp <- runif(nrow(oppairs)) < oppairs$pOp;
		oppairs <- oppairs[order(oppairs$index1,oppairs$index2),];
		oppairsT <- oppairs[oppairs$sOp,];

		mu[oppairsT$index1] <- mu[oppairsT$index2];
		if (matchSigma) sigma[oppairsT$index1] <- sigma[oppairsT$index2];
		while (countif(badpair <- (mu[oppairsT$index1] != mu[oppairsT$index2])) > 0) {
			if(debug) cat("Syncing operon pairs, mismatch left", countif(badpair),"\n");
			mu[oppairsT$index2] <- mu[oppairsT$index1];
			if(matchSigma) sigma[oppairsT$index2] <- sigma[oppairsT$index1];
		}
	}

	bias <- rnorm(length(locusIds)) * (sdBias + sigma*sqrt(scaledBias));

	ss <- function(x) sum( withoutNA(x - mean(x,na.rm=TRUE))^2 );
	if (control) {
		
		data1 <- t(sapply(1:length(locusIds),
			function(x) ifelse(runif(n)<pNA,NA,rnorm(n,mean=mu[x]+bias[x],sd=sigma[x]))))
		data2 <- t(sapply(1:length(locusIds),
			function(x) ifelse(runif(n2)<pNA,NA,rnorm(n2,mean=0,sd=sigma[x]))));
		muObs <- apply(data1,1,mean,na.rm=TRUE) - apply(data2,1,mean,na.rm=TRUE);
		# estimated standard deviation of mu: Var(SumX/N - SumY/N) = (Var(X)+Var(Y))/N)
		sumsq <- apply(data1,1,ss) + apply(data2,1,ss);
	} else {
		data1 <- t(sapply(1:length(locusIds),
			function(x) ifelse(runif(n)<pNA,NA,rnorm(n,mean=mu[x]+bias[x],sd=sigma[x]))));
		data2 <- NULL;
		muObs <- apply(data1,1,mean,na.rm=TRUE);
		sumsq <- apply(data1,1,ss);
	}
	if(sweep) muObs <- muObs - mean(muObs);
	if (!is.null(opPreds)) oppairs$agree <- (muObs[oppairs$index1]>0) == (muObs[oppairs$index2]>0);
	n1 <- apply(data1, 1, function(x) countif(!is.na(x)))
	n2 <- if(is.null(data2)) NULL else apply(data2, 1, function(x) countif(!is.na(x)));
	sim <- list(mu=mu,sigma=sigma,bias=bias,muObs=muObs,sumsq=sumsq,
			n1=n1,n2=n2,
			changer=changer,
			oppairs=if(is.null(opPreds)) NULL else oppairs,data1=data1,data2=data2,
			control=control,
			locusId=locusIds,
			# summary statistics
			medianSS=median(sumsq,na.rm=TRUE),
			sdMuObs=sd(muObs,na.rm=TRUE),
			seed=seed,
			trueFit=list(a=a,bc=bc,b=b,c=c,v=v),
			pOpAgree=if(is.null(opPreds)) NULL else countPerIf(withoutNA(oppairs$agree[oppairs$pOp>.5])));
	if(!estimate) sim$fit <- sim$trueFit;
	if(analyze) return(AnalyzeUarrayVsOperons(sim,debug=debug,...));
	return(sim);
}

# sim must contain locusId, muObs, sumsq, n1, n2, and oppairs
# with index1, index2, and either pOp (if useP) or sOp (if !useP)
# If not provided, sim$fit is estimated with FitUarrayData
# Then out is computed (from TestUarrayVsModelOpWise)
# Then changers is computed (from changersVsOperonPairs) and also aggregated stats changersCMeans, changersPAgrees(Adj)
AnalyzeUarrayVsOperons <- function(sim, useP=TRUE,useBias=TRUE,bcMethod="ll",debug=FALSE,nConf=8) {
	if(is.null(sim$fit)) {
		sim$fit <- FitUarrayData(sim$oppairs,sim$muObs,sim$sumsq,sim$n1,sim$n2,
					useP=useP,useBias=useBias,bcMethod=bcMethod,debug=debug);
	}
	sim$out <- TestUarrayVsModelOpWise(sim,sim$fit,useP=useP);
	sim$changers <- changersVsOperonPairs(sim$locusId,logodds(sim$out$p),sim$oppairs,both=TRUE);
	sim$changers$confidence <- logodds2p(abs(sim$changers$statistic));
	sim$changers$cfactor <- cutN(sim$changers$confidence,n=nConf);
	sim$changersCMeans <- sapply(split(sim$changers$confidence,sim$changers$cfactor),mean);
	sim$changersPAgrees <- sapply(split(sim$changers$bAgree,sim$changers$cfactor),mean);
	sim$changersPAgreesConf <- t(sapply(split(sim$changers$bAgree,sim$changers$cfactor),
					function(x) t.test(x)$conf.int));
	# note this uses +-1 not 1,0
	sim$changersPAgreesAdj <- sapply(split(sim$changers$raw/sim$changers$pOp,sim$changers$cfactor),mean);
	sim$changersPAgreesAdjConf <- t(sapply(split(sim$changers$raw/sim$changers$pOp,sim$changers$cfactor),
					function(x) t.test(x)$conf.int));
	return(sim);
	
}

# take real data and produce a list for input to AnalyzeUarrayVsOperons
# data and (if genomic control) data2 should be data frames or matrices with 1 row per gene, 1 column per replicate,
# and NA for missing values
#
# For genomic control, muObs will be mean(data1)-mean(data2)
#
# If sweep is TRUE then each column is swept (using genes present in all columns)
# If loess is non-null then a loess normalization is performed on the mean vs. A (which needs
#	to be provided for direct comparisons but not for genomic control)
#		In this case,
# If sector is non-null then a sector normalization is performed (after the loess), using
#	sector as a factor over rows
# If spots is non-null then the data is over spots not over locusId (so that locusId is non-unique)
SummarizeUarrayVsOperons <- function(opPred, locusId, data1, data2=NULL, sweep=FALSE, loess=FALSE,
					A=NULL, sector=NULL, spots=FALSE, name=NULL) {
	locusInput <- locusId;
	n1 <- apply(data1,1,function(x) countif(is.finite(x)));
	n2 <- if(!is.null(data2)) apply(data2,1,function(x) countif(is.finite(x))) else NULL;

	# sweep before averaging or normalization or normalization
	if(sweep) {
		good <- n1==ncol(data1);
		if(!is.null(data2)) good <- good & n2==ncol(data2);
		for (i in 1:ncol(data1)) { data1[,i] <- data1[,i] - mean(data1[good,i]); }
		if(!is.null(data2)) for (i in 1:ncol(data2)) { data2[,i] <- data2[,i] - mean(data2[good,i]); }
	}

	muObs <- apply(data1,1,mean,na.rm=TRUE);
	if(!is.null(data2)) {
		mu1 <- muObs;
		mu2 <- apply(data2,1,mean,na.rm=TRUE);
		muObs <- mu1 - mu2;
	}

	if(loess) {
		use <- !is.na(muObs);
		if(is.null(data2)) use <- use & !is.null(A);
		if(is.null(data2)) {
			l <- loess(muObs[use]~A[use])$fitted;
		} else {
			l <- loess(muObs[use]~I(mu1[use]+mu2[use]))$fitted;
		}
		muObs[use] <- muObs[use] - l;
	}
	if(!is.null(sector)) {
		sectorMedians <- sapply(split(muObs,sector),median,na.rm=TRUE);
		muObs <- muObs - sectorMedians[unclass(sector)];
	}

	ss <- function(x) sum( withoutNA(x - mean(x,na.rm=TRUE))^2 );

	if(!spots) {
		sumsq <- apply(data1,1,ss);
		if (!is.null(data2)) sumsq <- sumsq + apply(data2,1,ss);
	} else {
		# aggregate everything over spots
		n1 <- aggregate(n1,list(locusId),sum)$x;
		muObs <- aggregate(muObs,list(locusId),mean,na.rm=TRUE)$x;
		sumsq <- aggregate(c(data1,recursive=TRUE),list(rep(locusId,ncol(data1))),ss)$x;
		if(!is.null(n2)) {
			n2 <- aggregate(n2,list(locusId),sum)$x;
			sumsq <- sumsq + aggregate(c(data2,recursive=TRUE),list(rep(locusId,ncol(data2))),ss)$x;
			mu1 <- aggregate(mu1,list(locusId),mean,na.rm=TRUE)$x;
			mu2 <- aggregate(mu2,list(locusId),mean,na.rm=TRUE)$x;
		}
		locusId <- unique(sort(locusId));
	}

	oppairs <- merge(opPred[opPred$pOp>.5,],data.frame(index1=1:length(locusId),Gene1=locusId));
	oppairs <- merge(oppairs,data.frame(index2=1:length(locusId),Gene2=locusId));
	data <- list(muObs=muObs,sumsq=sumsq,n1=n1,n2=n2,oppairs=oppairs,locusId=locusId,
			locusInput=locusInput,data1=data1,data2=data2,
			sweep=sweep,loess=loess,sector=sector,spots=spots,A=A);
	if(!is.null(n2)) { data$mu1 <- mu1; data$mu2 <- mu2; }
	if(!is.null(name)) { data$name <- name; }
	return(data);
}

# returns p(true mu > 0)
TestUarrayVsModel <- function(m, sumsq, fit=list(a=0.3, b=0.1,c=Inf,bc=0.1,v=3.25), n1=10, n2=NULL) {
	ntot <- if(is.null(n2)) n1 else n1 + n2 - 1;
	nn <- if(is.null(n2)) n1 else n1 * n2/(n1 + n2);
	nprime <- 1/(1/nn + 1/fit$c);
	V <- (fit$a + m^2 * nprime * fit$b/(fit$b + nprime) + sumsq)/((fit$b + nprime)*(fit$v + ntot + 1));
	mu0 <- m * nprime/(nprime + fit$b);
	p <- pt(-mu0/sqrt(V), df=fit$v + ntot + 1, lower=FALSE);
	return(ifelse(ntot<1 | andNoNA(nn==0),0.5,p));
}

# Tests data for a SINGLE group of genes believed to have matching expression patterns
# All arguments except for fit should be vectors
TestUarrayVsModelN <- function(m, sumsq, n1, n2=NULL,fit=list(a=0.3, b=0.1555555,c=.28,bc=0.1,v=3.25)) {
	ntot <- if(is.null(n2)) n1 else n1 + n2 - 1;
	nn <- if(is.null(n2)) n1 else 1/(1/n1 + 1/n2); # = 0 if either n1 or n2 is 0 -> no blow ups

	good <- (ntot >= 0) & is.finite(m);
	if(!any(good)) return(0.5); # no data

	m <- m[good]; sumsq <- sumsq[good]; ntot <- ntot[good]; nn <- nn[good];

	nprime <- 1/(1/nn + 1/fit$c); # n1c, n2c, etc. in text
	V <- (fit$a + sum(sumsq) + sum(nprime*m^2) - (sum(nprime*m))^2/(fit$b+sum(nprime))) /
		((fit$v + sum(ntot) + 1) * (fit$b + sum(nprime)));
	mu0 <-sum(m * nprime)/(sum(nprime) + fit$b);
	p <- pt(-mu0/sqrt(V), df=fit$v + sum(ntot) + 1, lower=FALSE);
	return(p);
}

TestUarrayVsModelPairs <- function(pair,mu,sumsq,n1,n2=NULL,useP=TRUE,useBias=TRUE,
				fit=FitUarrayData(pair,mu,sumsq,n1,n2=n2,useP=useP,useBias=useBias)) {
	sapply(1:nrow(pair),function(i) { j <- c(pair$index1[i],pair$index2[i]);
					TestUarrayVsModelN(mu[j],sumsq[j],n1[j],
								if(is.null(n2)) NULL else n2[j], fit=fit); });
}

# requires sim$oppair, sim$muObs, sim$sumsq, sim$n1, and sim$n2 but not other components
# oppair needs index1, index2, and either pOp or sOp
TestUarrayVsModelOpWise <- function(sim,fit,useP=TRUE) {
	p <- TestUarrayVsModel(sim$muObs,sim$sumsq,sim$n1,sim$n2,fit=fit);
	pairs <- sim$oppair[if(useP) sim$oppair$pOp>.5 else sim$oppair$sOp,];
	# compute posterior probability for each operon
	if(useP) {
		if (is.null(sim$n2)) {
			ntot <- sim$n1;
			nn <- sim$n1;
		} else { # genomic control
			ntot <- sim$n1 + sim$n2 - 1;
			nn <- sim$n1*sim$n2/(sim$n1+sim$n2);
		}
		pairs$pOpPost <- logodds2p(logodds(pairs$pOp)
					+ log(relLik(sim$muObs[pairs$index1],
						sim$sumsq[pairs$index1],
						sim$muObs[pairs$index2],
						sim$sumsq[pairs$index2],
						a=sim$fit$a,
						bc=sim$fit$bc,
						c=sim$fit$c,
						n1tot=ntot[pairs$index1],
						nn1=nn[pairs$index1],
						n2tot=ntot[pairs$index2],
						nn2=nn[pairs$index2])));
		pairs$pOpPost[is.na(pairs$pOpPost)] <- pairs$pOp[is.na(pairs$pOpPost)];
	} else {
		pairs$pOpPost <- rep(1,nrow(pairs));
	}
	pairs$row <- 1:nrow(pairs);
	ppair <- TestUarrayVsModelPairs(pairs,sim$muObs,sim$sumsq,sim$n1,sim$n2,fit=fit);

	# 1.Gene2=2.Gene1=center
	# index1.x,index2.x,index2.y is the sequence of indices
	triples <- merge(pairs,pairs,by.x="Gene2",by.y="Gene1");
	if(nrow(triples) > 0) {
		triples$row <- 1:nrow(triples);
		ptriple <- apply(triples[,c("index1.x","index2.x","index2.y")], 1,
			function(j) TestUarrayVsModelN(sim$muObs[j],sim$sumsq[j],
							sim$n1[j],if(is.null(sim$n2)) NULL else sim$n2[j], fit=fit) );
	} else {
		ptriple <- 0.5;
	}
	pTot <- p;
	inPair1 <- merge(data.frame(index1=1:length(sim$muObs)),pairs,all.x=TRUE,all.y=FALSE);
	inPair2 <- merge(data.frame(index2=1:length(sim$muObs)),pairs,all.x=TRUE,all.y=FALSE);
	if (nrow(triples)>0) {
		inTriple <- merge(data.frame(index2.x=1:length(sim$muObs)),triples,all.x=TRUE,all.y=FALSE);
	} else {
		inTriple <- data.frame(index2.x=1:length(sim$muObs),row=rep(NA,length(sim$muObs)));
	}

	# for genes not in pairs, have no change
	# for genes in a pair but not in a triple, use P(Op)*ppair + (1-P(Op))*p
	# for genes in a triple, with P(OpL) and P(OpR), use all three:
	# pTot = P(OpL)*P(OpR)*ptriple + (1-P(OpL))*P(OpR)*ppairR + P(OpL)*(1-P(OpR))*ppairL
	#	+ (1-P(OpL))*(1-P(OpR))*p
	# (If useP is false than all P(Op)=1)

	out <- data.frame(locusId=sim$locusId,
				index=1:length(sim$muObs),
				p=p,
				pOp1=pairs$pOp[inPair1$row],
				pOpPost1=pairs$pOpPost[inPair1$row],
				ppair1=ppair[inPair1$row],
				pOp2=pairs$pOp[inPair2$row],
				pOpPost2=pairs$pOpPost[inPair2$row],
				ppair2=ppair[inPair2$row]);
	out$ptriple <- ptriple[inTriple$row];
	out$pOp1[is.na(out$pOp1)] <- 0;
	out$pOpPost1[is.na(out$pOpPost1)] <- 0;
	out$ppair1[is.na(out$ppair1)] <- 0.5;
	out$pOp2[is.na(out$pOp2)] <- 0;
	out$pOpPost2[is.na(out$pOpPost2)] <- 0;
	out$ppair2[is.na(out$ppair2)] <- 0.5;
	out$ptriple[is.na(out$ptriple)] <- 0.5;

	# associated pOp will be 0 if value is not relevant
	out$pTot <- (out$pOpPost1*out$pOpPost2*out$ptriple
			+ out$pOpPost1*(1-out$pOpPost2)*out$ppair1
			+ (1-out$pOpPost1)*out$pOpPost2*out$ppair2
			+ (1-out$pOpPost1)*(1-out$pOpPost2)*out$p);
	return(out);
}

# m1 is observed mean, ss1 is sum of squared variation from mean, n1 is #data points,
# and similarly for gene 2
relLik <- function(m1,ss1,m2,ss2,a=0.3,bc=.1,c=0.5,v=3.25,
			n1tot=10,n2tot=10, nn1=n1tot, nn2=n2tot) {
	if(c-bc < 1e-6) return( 0 );
	b <- 1/(1/bc-1/c);
	n1c <- 1/(1/nn1+1/c);
	n2c <- 1/(1/nn2+1/c);

	return( (a/2)^((v+1)/2-(v+1))
		* sqrt(b/(b+n1c+n2c)) * 1/sqrt((1+nn1/c)*(1+nn2/c)) * sqrt((bc+nn1)*(bc+nn2))/bc
		* gamma((v+n1tot+n2tot+1)/2)*gamma((v+1)/2)/(gamma((v+n1tot+1)/2)*gamma((v+n2tot+1)/2))
		* ((a + ss1 + m1^2*nn1*bc/(bc+nn1))/2)^((v+n1tot+1)/2) 
		* ((a + ss2 + m2^2*nn2*bc/(bc+nn2))/2)^((v+n2tot+1)/2) 
		/ ((a + ss1 + ss2 + n1c*m1^2 + n2c*m2^2 - (m1*n1c + m2*n2c)^2/(b+n1c+n2c))/2)^((v+n1tot+n2tot+1)/2) );
}

# derivative of the logarithm of relLik() with respect to c
relLikLogDeriv <- function(m1,ss1,m2,ss2,a=0.3,bc=.1,c=0.5,v=3.25,
			n1tot=10,n2tot=10, nn1=n1tot, nn2=n2tot) {
	if(c-bc < 1e-6) return( NA );

	b <- 1/(1/bc-1/c);
	n1c <- 1/(1/nn1+1/c);
	n2c <- 1/(1/nn2+1/c);

	# derivatives
	dbdc <- -(c/bc - 1)^-2;
	dn1cdc <- (1 + c/nn1)^-2;
	dn2cdc <- (1 + c/nn2)^-2;

	dcombdc <- (dn1cdc*m1^2 + dn2cdc*m2^2
			- 2*(m1*n1c + m2*n2c)*(m1*dn1cdc + m2*dn2cdc)/(b+n1c+n2c)
			- -(m1*n1c + m2*n2c)^2/(b+n1c+n2c)^2 * (dbdc+dn1cdc+dn2cdc) )/2;

	return(0.5*dbdc/b - 0.5*(dbdc+dn1cdc+dn2cdc)/(b+n1c+n2c)
		- 0.5*(-nn1/c^2)/(1+nn1/c) - 0.5*(-nn2/c^2)/(1+nn2/c)
		- (v+n1tot+n2tot+1)/2 * dcombdc
			/ ((a + ss1 + ss2 + n1c*m1^2 + n2c*m2^2 - (m1*n1c+m2*n2c)^2/(b+n1c+n2c))/2) );
}

# pair must include pOp, index1, and index2 (or binary sOp instead of pOp if useP is FALSE)
# mu and sumsq are observations
#
# Methods for estimating bc include "ll" (maximum log likelihood, the default), "moments" (from
# mean of mu^2), and mad (from mad(mu^2) -- doesn't work)
# "ll" is slightly more accurate in simulations and somewhat more robust against
# heavy tails in the data than "moments"
FitUarrayData <- function(pair,mu,sumsq, n1, n2=NULL, useP=TRUE, useBias=TRUE, bcMethod="ll", debug=FALSE) {
	if (is.null(n2)) {
		ntot <- n1;
		nn <- n1;
	} else { # genomic control
		ntot <- n1 + n2 - 1; # powers of sqrt(theta) in likelihood function for each gene
		nn <- n1*n2/(n1+n2); # multiplier of (muObs - mu)^2 in likelihood function for each gene
	}
	fit <- fitFDist(sumsq/(ntot-1),ntot-1); # sumsq depends on ntot, not on nn
	a <- fit$scale*fit$df2;
	v <- fit$df2-1;
	good <- andNoNA(is.finite(sumsq/(ntot-1)) & nn >= 1);
	meanM2 <- if(bcMethod=="mad") mad(mu[good])^2 else mean(mu[good]^2);
	bc <- 1/( (v-1)/a * meanM2 - mean(1/nn[good]));
	if (bcMethod=="ll") {
		if(bc < 1e-6) bc <- 0.01;
		f <- function(x) { if(x < 1e-6) return(1e20);
					v <- OpUarrayLLNoBias(mu,sumsq,n1,n2,a,x,v);
					v <- -v;
					attr(v,"gradient") <- -attr(v,"gradient");
					return(v); };
		nlmOut <- nlm(f, bc, check.analyticals=debug);
		if(debug) cat("nlm bc","min",nlmOut$minimum,"est(bc)",nlmOut$estimate,"grad",nlmOut$gradient,
					"code",nlmOut$code,"iter",nlmOut$iterations,"\n");
		bc <- nlmOut$estimate;
	}
	fit <- list(a=a,bc=bc,b=bc,c=Inf,v=v);
	if(!is.null(pair) && bc > 1e-5 && useBias) {
		p <- if(useP) pair[pair$pOp>.5,] else pair[pair$sOp,];
		if(!useP) p$pOp <- rep(1,nrow(p));
		p <- p[good[p$index1] & good[p$index2],];
		if(debug) cat("Using",nrow(p),"good operon pairs to estimate bias\n")
		# - sign to maximize, 4*bc a legal and plausible starting point
		# use 1/c as the parameter because nlm behaves better
		#	(numerical derivatives issue?)
		#
		# For gradient computation: f is the log likelihood if an operon
		# d(log(1-p+p*exp(f(c))))/dc = p*exp(f(c))*df/dc/(1-p+p*exp(f(c)))
		# dg/d(1/c) = dg/dc / (d(1/c)/dc) = -c^2 * df/dc = -(1/c)^-2 * df/dc
		relLikCInv <- function(x) {
					if(length(x) != 1) error("More than 1 value");
					liks <- relLik(mu[p$index1],sumsq[p$index1],
						mu[p$index2],sumsq[p$index2],
						a=a,bc=bc,c=1/x,v=v,
						n1tot=ntot[p$index1],
						n2tot=ntot[p$index2],
						nn1=nn[p$index1],
						nn2=nn[p$index2]);
					dllsdc <- relLikLogDeriv(mu[p$index1],sumsq[p$index1],
							mu[p$index2],sumsq[p$index2],
							a=a,bc=bc,c=1/x,v=v,
							n1tot=ntot[p$index1],
							n2tot=ntot[p$index2],
							nn1=nn[p$index1],
							nn2=nn[p$index2]);
					ratios <- 1-p$pOp + p$pOp*liks;
					out <- -sum(log(ratios));
					attr(out,"gradient") <- x^-2*sum(p$pOp*liks*dllsdc/ratios);
					return(out);
		}
		# c=4*bc (3:1 signal:bias) is a reasonable starting point
		nlmOut <- nlm(relLikCInv, 1/(4*bc), check.analyticals=debug);
		if(debug) cat("nlm c","min",nlmOut$minimum,"est(1/c)",nlmOut$estimate,"grad",nlmOut$gradient,
					"code",nlmOut$code,"iter",nlmOut$iterations,"\n");
		fit$c <- 1/nlmOut$estimate;
		fit$b <- 1/(1/fit$bc-1/fit$c);
		# method of moments
		fit$c2 <- mean(p$pOp)/(mean((mu[p$index1]-mu[p$index2])^2)*(v-1)/(2*a)
					- mean(1/nn[c(p$index1,p$index2)]) - mean(1-p$pOp)/bc);
		# for max. likelihood ratio test (df=1)
		fit$chisqC <- 2*((-nlmOut$minimum) - (-relLikCInv(0))); # c=Inf or 1/c=0
		fit$chisqCP <- 2*pchisq(fit$chisqC,df=1,lower=FALSE);
	}
	return(fit);
}

### Taken from ebayes.R in LIMMA, by Gordon Smyth
### Reference: G.K. Smyth (2004) "Linear Models and Empirical Bayes Methods for Assessing
### 	Differential Expression in Microarray Experiments." Statistical Applications in
###	Genetics and Molecular Biology 3(1) Article 3.
###
fitFDist <- function(x,df1) {
#	Moment estimation of the parameters of a scaled F-distribution
#	The first degrees of freedom is given
#	Gordon Smyth
#	8 Sept 2002.  Last revised 12 Apr 2003.

#	Remove missing or infinite values and zero degrees of freedom
	o <- is.finite(x) & is.finite(df1) & (x > 0) & (df1 > 0)
	if(any(!o)) {
		x <- x[o]
		df1 <- df1[o]
	}
	n <- length(x)
	if(n==0) return(list(scale=NA,df2=NA))
	
#	Better to work on with log(F)
	z <- log(x)
	e <- z-digamma(df1/2)+log(df1/2)
	emean <- mean(e)
	evar <- mean(n/(n-1)*(e-emean)^2-trigamma(df1/2))
	if(evar > 0) {
		df2 <- 2*trigammaInverse(evar)
		s20 <- exp(emean+digamma(df2/2)-log(df2/2))
	} else {
		df2 <- Inf
		s20 <- exp(emean)
	}
	list(scale=s20,df2=df2)
}

trigammaInverse <- function(x) {
#	Solve trigamma(y) = x for y
#	Gordon Smyth
#	8 Sept 2002.  Last revised 12 March 2004.

#	Non-numeric or zero length input
	if(!is.numeric(x)) stop("Non-numeric argument to mathematical function")
	if(length(x)==0) return(numeric(0))

#	Treat out-of-range values as special cases
	omit <- is.na(x)
	if(any(omit)) {
		y <- x
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x < 0)
	if(any(omit)) {
		y <- x
		y[omit] <- NaN
		warning("NaNs produced")
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x > 1e7)
	if(any(omit)) {
		y <- x
		y[omit] <- 1/sqrt(x[omit])
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x < 1e-6)
	if(any(omit)) {
		y <- x
		y[omit] <- 1/x[omit]
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}

#	Newton's method
#	1/trigamma(y) is convex, nearly linear and strictly > y-0.5,
#	so iteration to solve 1/x = 1/trigamma is monotonically convergent
	Rversion <- as.numeric(version$major)+as.numeric(version$minor)/10
	if(Rversion > 1.81)
		tetraGamma <- function(x) psigamma(x,deriv=2)
	else
		tetraGamma <- tetragamma
	y <- 0.5+1/x
	iter <- 0
	repeat {
		iter <- iter+1
		tri <- trigamma(y)
		dif <- tri*(1-tri/x)/tetraGamma(y)
		y <- y+dif
		if(max(-dif/y) < 1e-8) break
		if(iter > 50) {
			warning("Iteration limit exceeded")
			break
		}
	}
	y
}

TestPValues <- function(p, truth,rank=FALSE,plot=TRUE,worst=FALSE,minE=20) {
	o <- order(p);
	p <- p[o];
	truth <- truth[o];
	nTrue <- cumsum(truth);
	nExpected <- cumsum(p);
	Variance <- cumsum(p*(1-p));
	if(plot && !worst) {
		plot(if(rank) rank(p) else nExpected,
			(nTrue-nExpected)/sqrt(Variance),
			type="l",
			ylab="Z score",
			xlab=if(rank) "Rank(p)" else "# Expected",
			log=if(rank) "" else "x"); hline(0);
	} else {
		z <- (nTrue-nExpected)/sqrt(Variance);
		if(worst) return(max(abs(z[nExpected>minE])));
		return(z[nExpected>minE]);
	}
}

IntegrateThMu <- function(fThMu) {
	integrate(function(x) {sapply(x,function(xx) IntegrateMu(xx,fThMu))}, 0, Inf);
}

IntegrateMu <- function(y,fThMu) { integrate(function(z) fThMu(y,z),-Inf,Inf)$value }

## obsolete
tcomb <- function(x,x1=1,x2=1,V1=1,V2=1,d1=2,d2=2) {
	return( dt((x-x1)/sqrt(V1),df=d1) * dt((x-x2)/sqrt(V2),df=d2) );
}

## obsolete
# tests if true mean > 0 given two independent sources of evidence with different
# estimates, variances, and dfs...
CombineTTest <- function(x1,x2,V1,V2,d1,d2) {
	sapply(1:length(x1), function(i) CombineTTest1(x1[i],x2[i],V1[i],V2[i],d1[i],d2[i]));
}

## obsolete
CombineTTest1 <- function(x1,x2,V1,V2,d1,d2) {
	return( integrate(tcomb,0,Inf,x1=x1,x2=x2,V1=V1,V2=V2,d1=d1,d2=d2)$value
		/ integrate(tcomb,-Inf,Inf,x1=x1,x2=x2,V1=V1,V2=V2,d1=d1,d2=d2)$value );
}

# returns a string of 0 and 1, starting with the 0th power
IntToBits <- function(x) {
	if(x==0) { return(0); }
	n <- 1 + floor(log(x)/log(2));
	bits <- sapply(0:n,function(i) floor(1e-9 + x/2^i)%%2);
	return(bits);
}

BitsToInt <- function(bits) {
	on <- (0:(length(bits)-1))[bits==1];
	return(sum(2^on));
}

# argument is a list of analyzes sims from AnalyzeUarrayVsOperons
CompareOpAgreements <- function(list, adjust=TRUE,
				col=1:length(list), pch=1:length(list), lty=1:length(list), lwd=rep(1,length(list)),
				labels=names(list),
				error.bars=c(1) # which elements to show error bars for
			)
{
	plot(c(0.5,1),c(0.5,if(adjust) 1.2 else 1.0),col="white",xlab="Confidence within Group", ylab=if(adjust) "Adjusted Agreement with Operon Predictions" else "Agreement with Operon Predictions");
	for (i in 1:length(names(list))) {
		s <- list[[names(list)[i]]];
		y <- if(adjust) s$changersPAgreesAdj else s$changersPAgrees;
		points(s$changersCMeans, y, col=col[i], pch=pch[i]);
		lines(s$changersCMeans, y, col=col[i], lty=lty[i], lwd=lwd[i]);
		if (i %in% error.bars) {
			conf.int <- if(adjust) s$changersPAgreesAdjConf else s$changersPAgreesConf;
			arrows( s$changersCMeans, conf.int[,1], s$changersCMeans, conf.int[,2],
					angle=90, code=3, length=1/12, col="darkgrey" );
		}
	}
	eqline();
	legend(.5,1,labels,col=col,pch=pch,lty=lty,lwd=lwd);
}

AverageOpAgreements <- function(sims) {
	changersCMeans <- sapply(sims,function(x) x$changersCMeans); # column for each sim
	changersPAgrees <- sapply(sims,function(x) x$changersPAgrees);
	changersPAgreesAdj <- sapply(sims,function(x) x$changersPAgreesAdj);
	return(list(changersCMeans=apply(changersCMeans,1,mean),
			changersPAgrees=apply(changersPAgrees,1,mean),
			changersPAgreesAdj=apply(changersPAgreesAdj,1,mean)));
}

OpUarrayLLNoBias <- function(m,sumsq,n1,n2,a,b,v) {
	ntot <- if(is.null(n2)) n1 else n1 + n2 - 1;
	nn <- if(is.null(n2)) n1 else 1/(1/n1 + 1/n2);
	good <- (ntot >= 0) & is.finite(m);
	if(!any(good)) return(0); # no data
	m <- m[good]; sumsq <- sumsq[good]; ntot <- ntot[good]; nn <- nn[good];

	ll <- sum( log(a)*(v+1)/2 + log(b)/2 - lgamma((v+1)/2) - log(2)*(v+1)/2
			- (ntot+1)/2 * log(2*pi)
			+ (log(2*pi) - log(b+nn))/2 + lgamma((v+ntot+1)/2)
			- ((v+ntot+1)/2) * log((a + sumsq + m^2*nn*b/(b+nn))/2) );
	attr(ll,"gradient") <- 0.5*sum(1/b - 1/(b+nn) - (v+ntot+1)*m^2*nn^2/((a+sumsq)*(b+nn)^2 + m^2*nn*b*(b+nn)));
	return(ll);

}

OpUarrayLLNoBiasDerivative <- function(m,sumsq,n1,n2,a,b,v) {
	ntot <- if(is.null(n2)) n1 else n1 + n2 - 1;
	nn <- if(is.null(n2)) n1 else 1/(1/n1 + 1/n2);
	good <- (ntot >= 0) & is.finite(m);
	if(!any(good)) return(0); # no data
	m <- m[good]; sumsq <- sumsq[good]; ntot <- ntot[good]; nn <- nn[good];

}

# Given a normalized data set with raw data1, data2, get back the raw
# data aggregated by spots. (Just returns data1 & data2 if !data$spots)
OpUarrayCombineSpots <- function(data,renormalize=FALSE) {
	if (!data$spots) {
		data1 <- data$data1;
		data2 <- data$data2;
	} else {
		nspots <- max(table(data$locusInput));
		n1 <- ncol(data$data1)*nspots;
		n2 <- if(is.null(data$data2)) 0 else ncol(data$data2)*nspots;
		data1 <- t(sapply(split(c(data$data1,recursive=TRUE),rep(data$locusInput,ncol(data$data1))),
					function(x) { length(x) <- n1; x}));
		if (n2 > 0) {
			data2 <- t(sapply(split(c(data$data2,recursive=TRUE),rep(data$locusInput,ncol(data$data2))),
						function(x) { length(x) <- n2; x}));
		}
		dataIndex <- order(data$locusId);
		data1 <- data1[dataIndex,];
		if(n2 > 0) data2 <- data2[dataIndex,];
	}
	if (renormalize) {
		muExpected <- apply(data1,1,mean,na.rm=TRUE);
		if (!is.null(data$data2)) muExpected <- muExpected - apply(data2,1,mean,na.rm=TRUE);
		data1 <- data1 + data$muObs-muExpected;
	}
	return(list(data1=data1,data2=data2));
}

# Computes non-parametric p-values for each gene:
#
# To do this, we perform all balanced permutations or sign-shifts of the data
# (maintaining the original data) and count how often the single-gene p-values for each
# are more extreme than in the given data object.
#
OpUarrayPermuteTest <- function(data,debug=FALSE,permutations=NULL,useBias=FALSE,useP=TRUE,returnP=FALSE) {
	nspots <- 1;
	if(data$spots) { nspots <- max(table(data$locusInput)); }
	n1 <- ncol(data$data1)*nspots;
	n2 <- if(is.null(data$data2)) 0 else ncol(data$data2)*nspots;
	# flipping all==flipping none
	# for +/- flips, flipping top is equivalent to flipping all others
	nPerm <- if(n2>0) 2^(n1+n2)-1 else 2^(n1-1);
	nPermUse <- rep(0,length(data$locusId));
	permutedN <- rep(0,length(data$locusId));

	fit <- data$fit;
	if (!useBias) { fit$c <- Inf; fit$b <- fit$bc; }

	# make a simplified data set -- will use the normalization offset implied by data$muObs
	r <- OpUarrayCombineSpots(data,renormalize=TRUE);
	data1 <- r$data1;
	if (n2 > 0) data2 <- r$data2;

	conf <- NULL; # will be set during iteration 0
	for (perm in 0:(nPerm-1)) {
		bits <- IntToBits(perm);
		flip <- rep(FALSE,n1+n2);
		if(length(bits) > length(flip)) length(bits) <- length(flip);
		if(perm > 0) flip[1:length(bits)][bits==1] <- TRUE;
		use <- if(is.null(permutations) || perm==0) TRUE else perm %in% permutations;
		if (n2 > 0) {
			flip1 <- flip[1:n1];
			flip2 <- flip[(n1+1):(n1+n2)];
			use <- use && countif(flip1) == countif(flip2);
		}
		if(use) {
			if(debug) cat("Using permutation",perm,"bits",bits,"flip",flip,"\n");
			if (n2 > 0) {
				data1p <- cbind(data1[,!flip1],data2[,flip2]);
				data2p <- cbind(data1[,flip1],data2[,!flip2]);
			} else {
				data1p <- data1;
				for (j in 1:n1) if(flip[j]) data1p[,j] <- -data1p[,j];
				data2p <- NULL;
			}
			ss <- function(x) sum( withoutNA(x - mean(x,na.rm=TRUE))^2 );

			muPerm <- apply(data1p,1,mean,na.rm=TRUE);
			sumsqPerm <- apply(data1p,1,ss);
			if (n2 > 0) {
				muPerm <- muPerm - apply(data2p,1,mean,na.rm=TRUE);
				sumsqPerm <- sumsqPerm + apply(data2p,1,ss);
			}
			n1Perm <- apply(data1p,1,function(x) countif(!is.na(x)));
			n2Perm <- if (n2 > 0) apply(data2p,1,function(x) countif(!is.na(x))) else NULL;
			permP <- TestUarrayVsModel(muPerm, sumsqPerm,fit=fit,n1=n1Perm,n2=n2Perm);
			if(returnP) return(permP); # for debugging
			confPerm <- 2*abs(permP - 0.5);
			if (perm==0) conf <- confPerm; # compare to non-permuted confidence values
			usableData <- n1Perm >= 1;
			if (n2 > 0) usableData <- usableData & n2Perm >= 1;
			nPermUse <- nPermUse + ifelse(usableData,1,0);
			permutedN <- permutedN + ifelse(confPerm >= conf & usableData, 1, 0);
			if(debug) cat("Confidence >= unpermuted for", countif(confPerm >= conf & usableData),
					"of", countif(nPermUse),"loci\n");
		}
	}
	if(debug) cat("ran",max(nPermUse),"permutations\n");
	return(permutedN/nPermUse);
}

SimulateUarrayVsOperonsMultiple <- function(locusIds,opPreds,tries=100,...) {
	sims <- list();
	mus <- NULL;
	fits <- NULL;
	for(i in 1:tries) {
		s <- SimulateUarrayVsOperons(locusIds,opPreds,...);
		sims[[i]] <- s;
		pTrueFit <- TestUarrayVsModelOpWise(s,s$trueFit);

		if(is.null(s$n2)) s$n2 <- 0*s$n1;
		mu <- data.frame(mu=s$mu,sigma=s$sigma,muObs=s$muObs,sumsq=s$sumsq,
					n1=s$n1,n2=s$n2,locusId=s$locusId);
		mus <- if(is.null(mus)) mu else rbind(mus,mu);
		s$fit$pOpAgree <- s$pOpAgree;

		# T is short for true, I is short for ideal
		model1T <- summary(glm(s$mu>0 ~ logodds(s$out$p)+0,family=binomial(logit)));
		s$fit$slope1T <- model1T$coefficients[,1];
		s$fit$slope1Terr <- model1T$coefficients[,2];

		inpair <- s$out$pOp1 > 0 | s$out$pOp2 > 0;

		modelOpT <- summary(glm(s$mu[inpair]>0 ~ logodds(s$out$pTot[inpair])+0,family=binomial(logit)));
		s$fit$slopeOpT <- modelOpT$coefficients[,1];
		s$fit$slopeOpTerr <- modelOpT$coefficients[,2];

		model1I <- summary(lm(logodds(s$out$p)~logodds(pTrueFit$p)+0));
		s$fit$slope1I <- model1I$coefficients[,1];
		s$fit$slope1Ierr <- model1I$coefficients[,2];

		modelOpI <- summary(lm(logodds(s$out$pTot[inpair])~logodds(pTrueFit$pTot[inpair])+0));
		s$fit$slopeOpI <- modelOpI$coefficients[,1];
		s$fit$slopeOpIerr <- modelOpI$coefficients[,2];

		fits <- if(is.null(fits)) as.data.frame(s$fit) else rbind(fits,as.data.frame(s$fit));
	}
	return(list(mus=mus,fits=fits,sims=sims));
}

# see Storey 2001
# use a simple lambda cutoff instead of the cubic spline
PositiveFDR <- function(p, lambda=0.5, debug=FALSE) {
	index <- order(p);
	o <- reverse(index);
	p0 <- min(countPerIf(p>lambda)/(1-lambda), 1);
	if(debug) cat("p0",p0,"\n");
	q <- p0 * p * length(p)/ o;
	for (i in (length(q)-1):1) q[index[i]] <- min(q[index[i]],q[index[i+1]]);
	return(pmin(q,1));
}

# by Marc Schwartz
# see
# http://www.r-project.org/nocvs/mail/r-devel/2002/1374.html

barplot2 <- function(height, ...) UseMethod("barplot2")

barplot2.default <-
function(height, width = 1, space = NULL, names.arg = NULL,
       legend.text = NULL, beside = FALSE, horiz = FALSE,
       density = NULL, angle = 45,
       col = heat.colors(NR), prcol = NULL, border = par("fg"),
       main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
       xlim = NULL, ylim = NULL, xpd = TRUE, log = "",
       axes = TRUE, axisnames = TRUE,
       cex.axis = par("cex.axis"), cex.names = par("cex.axis"),
       inside = TRUE, plot = TRUE,
       plot.ci = FALSE, ci.l = NULL, ci.u = NULL,
       ci.color = "black", ci.lty = "solid", ci.lwd = 1,
       plot.grid = FALSE, grid.inc = 5,
       grid.lty = "dotted", grid.lwd = 1, grid.col = "black",
       ...)
{
    if (!missing(inside)) .NotYetUsed("inside", error = FALSE)
    if (!missing(border)) .NotYetUsed("border", error = FALSE)

    if (missing(space))
      space <- if (is.matrix(height) && beside) c(0, 1) else 0.2

    space <- space * mean(width)

    if (plot && axisnames && missing(names.arg))
      names.arg <- if(is.matrix(height)) colnames(height) else names(height)

    if (is.vector(height))
    {
      height <- cbind(height)
      beside <- TRUE
    }
    else if (is.array(height) && (length(dim(height)) == 1))
    {
      height <- rbind(height)
      beside <- TRUE
    }
    else if (!is.matrix(height))
      stop("`height' must be a vector or a matrix")

    if(is.logical(legend.text))
    {
      if(legend.text && is.matrix(height))
        legend.text <- rownames(height)
      else
        legend.text <- NULL
    }

    # Check for log scales
    logx <- FALSE
    logy <- FALSE

    if (log != "")
    {
      if (any(grep("x", log)))
        logx <- TRUE

      if (any(grep("y", log)))
        logy <- TRUE
    }

    # Cannot "hatch" with rect() when log scales used
    if ((logx || logy) && !is.null(density))
      stop("Cannot use shading lines in bars when log scale is used")

    NR <- nrow(height)
    NC <- ncol(height)

    if (beside)
    {
      if (length(space) == 2)
        space <- rep(c(space[2], rep(space[1], NR - 1)), NC)

      width <- rep(width, length = NR * NC)
    }
    else
      width <- rep(width, length = NC)

    delta <- width / 2
    w.r <- cumsum(space + width)
    w.m <- w.r - delta
    w.l <- w.m - delta

    #if graphic will be stacked bars, do not plot ci
    if (!beside && (NR > 1) && plot.ci)
      plot.ci = FALSE

    # error check ci arguments
    if (plot && plot.ci)
    {
      if ((missing(ci.l)) || (missing(ci.u)))
        stop("confidence interval values are missing")

      if (is.vector(ci.l))
        ci.l <- cbind(ci.l)
      else if (is.array(ci.l) && (length(dim(ci.l)) == 1))
        ci.l <- rbind(ci.l)
      else if (!is.matrix(ci.l))
        stop("`ci.l' must be a vector or a matrix")

      if (is.vector(ci.u))
        ci.u <- cbind(ci.u)
      else if (is.array(ci.u) && (length(dim(ci.u)) == 1))
        ci.u <- rbind(ci.u)
      else if (!is.matrix(ci.u))
        stop("`ci.u' must be a vector or a matrix")

      if (dim(height) != dim(ci.u))
        stop("'height' and 'ci.u' must have the same dimensions.")

      if (dim(height) != dim(ci.l))
        stop("'height' and 'ci.l' must have the same dimensions.")
    }

    if (beside)
      w.m <- matrix(w.m, nc = NC)

    # check height/ci.l if using log scale to prevent log(<=0) error
    # adjust appropriate ranges and bar base values
    if ((logx && horiz) || (logy && !horiz))
    {
      if (min(height) <= 0)
        stop("log scale error: at least one 'height' value <= 0")

      if (plot.ci && (min(ci.l) <= 0))
        stop("log scale error: at least one lower c.i. value <= 0")

      if (logx && !is.null(xlim) && (xlim[1] <= 0))
        stop("log scale error: 'xlim[1]' <= 0")

      if (logy && !is.null(ylim) && (ylim[1] <= 0))
        stop("'log scale error: 'ylim[1]' <= 0")

      # arbitrary adjustment to display some of bar for min(height) or min(ci.l)
      if (plot.ci)
        rectbase <- rangeadj <- (0.9 * min(c(height, ci.l)))
      else
        rectbase <- rangeadj <- (0.9 * min(height))

      # if axis limit is set to < above, adjust bar base value
      # to draw a full bar
      if (logy && !is.null(ylim))
        rectbase <- ylim[1]
      else if (logx && !is.null(xlim))
        rectbase <- xlim[1]

      # if stacked bar, set up base/cumsum levels, adjusting for log scale
      if (!beside)
        height <- rbind(rectbase, apply(height, 2, cumsum))

      # if plot.ci, be sure that appropriate axis limits are set to include range(ci)
      lim <-
        if (plot.ci)
          c(height, ci.l, ci.u)
        else
          height
    }
    else
    {
      # Use original bar base value
      rectbase <- 0

      # if stacked bar, set up base/cumsum levels
      if (!beside)
        height <- rbind(rectbase, apply(height, 2, cumsum))

      # if plot.ci, be sure that appropriate axis limits are set to include range(ci)
      lim <-
        if (plot.ci)
          c(height, ci.l, ci.u)
        else
          height

      # use original range adjustment factor
      rangeadj <- (-0.01 * lim)
    }

    # define xlim and ylim, adjusting for log-scale if needed
    if (horiz)
    {
      if (missing(xlim)) xlim <- range(rangeadj, lim, na.rm=TRUE)
      if (missing(ylim)) ylim <- c(min(w.l), max(w.r))
    }
    else
    {
      if (missing(xlim)) xlim <- c(min(w.l), max(w.r))
      if (missing(ylim)) ylim <- range(rangeadj, lim, na.rm=TRUE)
    }

    if(plot) ##-------- Plotting :
    {
      opar <-
        if (horiz)
          par(xaxs = "i", xpd = xpd)
        else
          par(yaxs = "i", xpd = xpd)

      on.exit(par(opar))

      plot.new()
      plot.window(xlim, ylim, log = log, ...)

      # if prcol specified, set plot region color
      if (!missing(prcol))
      {
        usr <- par("usr")

        # adjust par("usr") values if log scale(s) used
        if (logx)
        {
          usr[1] <- 10 ^ usr[1]
          usr[2] <- 10 ^ usr[2]
        }

        if (logy)
        {
          usr[3] <- 10 ^ usr[3]
          usr[4] <- 10 ^ usr[4]
        }

        rect(usr[1], usr[3], usr[2], usr[4], col = prcol)
      }

      # if plot.grid, draw major y-axis lines if vertical or x axis if horizontal
      # Due to issues with xaxp and yaxp when using log scale, use pretty() to set
      # grid and axis increments for for both linear and log scales when plotting grid lines
      if (plot.grid)
      {
        par(xpd = FALSE)

        if (horiz)
        {
          grid = pretty(xlim, n = grid.inc)
          abline(v = grid, lty = grid.lty, lwd = grid.lwd, col = grid.col)
        }
        else
        {
          grid = pretty(ylim, n = grid.inc)
          abline(h = grid, lty = grid.lty, lwd = grid.lwd, col = grid.col)
        }

         par(xpd = xpd)
      }

      # Beware : angle and density are passed using R scoping rules
      xyrect <- function(x1,y1, x2,y2, horizontal = TRUE, ...)
      {
        if(horizontal)
          rect(x1,y1, x2,y2, angle = angle, density = density, ...)
        else
          rect(y1,x1, y2,x2, angle = angle, density = density, ...)
      }

      if (beside)
        xyrect(rectbase, w.l, c(height), w.r, horizontal=horiz, col = col)
      else
      {
        for (i in 1:NC)
          xyrect(height[1:NR, i], w.l[i], height[-1, i], w.r[i], horizontal=horiz, col = col)
      }

      if (plot.ci)
      {
        # CI plot width = barwidth / 2
        ci.width = width / 4

        if (horiz)
        {
          segments(ci.l, w.m, ci.u, w.m, col = ci.color, lty = ci.lty, lwd = ci.lwd)
          segments(ci.l, w.m - ci.width, ci.l, w.m + ci.width, col = ci.color, lty = ci.lty, lwd = ci.lwd)
          segments(ci.u, w.m - ci.width, ci.u, w.m + ci.width, col = ci.color, lty = ci.lty, lwd = ci.lwd)
        }
        else
        {
          segments(w.m, ci.l, w.m, ci.u, col = ci.color, lty = ci.lty, lwd = ci.lwd)
          segments(w.m - ci.width, ci.l, w.m + ci.width, ci.l, col = ci.color, lty = ci.lty, lwd = ci.lwd)
          segments(w.m - ci.width, ci.u, w.m + ci.width, ci.u, col = ci.color, lty = ci.lty, lwd = ci.lwd)
        }
      }

      if (axisnames && !is.null(names.arg)) # specified or from {col}names
      {
        at.l <-
        if (length(names.arg) != length(w.m))
        {
          if (length(names.arg) == NC) # i.e. beside (!)
            colMeans(w.m)
          else
            stop("incorrect number of names")
        }
        else w.m

        axis(if(horiz) 2 else 1, at = at.l, labels = names.arg, lty = 0, cex.axis = cex.names, ...)
      }

      if(!is.null(legend.text))
      {
        legend.col <- rep(col, length = length(legend.text))

        if((horiz & beside) || (!horiz & !beside))
        {
          legend.text <- rev(legend.text)
          legend.col <- rev(legend.col)
          density <- rev(density)
          angle <- rev(angle)
        }

            xy <- par("usr")
            legend(xy[2] - xinch(0.1), xy[4] - yinch(0.1),
             legend = legend.text, angle = angle, density = density,
             fill = legend.col, xjust = 1, yjust = 1)
      }

      title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...)

      # if axis is to be plotted, adjust for grid "at" values
      if (axes)
      {
        if(plot.grid)
          axis(if(horiz) 1 else 2, at = grid, cex.axis = cex.axis, ...)
        else
          axis(if(horiz) 1 else 2, cex.axis = cex.axis, ...)
      }

      invisible(w.m)

    }

    else w.m
}


# pairs has columns Gene1 and Gene2. Returns a dataframe of pairs with added columns func & size,
# where func is the smallest function (and size is the #genes with that func) shared by the two,
# or size=0 if none is shared (with func=NA), or both=NA if one or other gene is uncharacterized
# Note return may be reordered
PairFuncSimilarity <- function(pairs, genes, functions) {
	fCount <- as.data.frame(table(functions));
	names(fCount) <- c("func","size");
	assign <- merge(data.frame(locusId=genes,func=functions),fCount);
	expanded <- merge(merge(pairs[,c("Gene1","Gene2")],assign,by.x="Gene1",by.y="locusId"),
			assign,
			by.x=c("Gene2","func","size"),
			by.y=c("locusId","func","size"));
	expanded <- expanded[order(expanded$size),];
	# retain only first func shared by each gene
	row1 <- aggregate(1:nrow(expanded),expanded[,c("Gene1","Gene2")],min); # Gene1,Gene2,x
	expanded <- merge(cbind(expanded,x=1:nrow(expanded)),row1);
	out <- merge(pairs,expanded[,c("Gene1","Gene2","func","size")],all.x=TRUE,all.y=FALSE);
	# and handle unannotated genes: first set unshared pairs to 0
	out$size[is.na(out$size)] <- 0;
	out$size[out$size==0 & !(out$Gene1 %in% genes & out$Gene2 %in% genes)] <- NA;
	return(out);
}

FisherTables <- function(table1,table2,ncol=2) {
	return(fisher.test(matrix(c(table1,table2),ncol=ncol)));
}

TestLevels <- function(f,by,level=levels(f),min=0.05) {
	out <- NULL;
	for (x in level) {
		if(countif(f==x)>=2) {
			test <- fisher.test(f==x,by);
			row <- data.frame(level=factor(x,levels=levels(f)),odds=test$estimate,p=test$p.value);
			if (row$p <= min) {
				if(is.null(out)) { out <- row; } else { out[nrow(out)+1,] <- row; }
			}
		}
	}
	row.names(out) <- 1:nrow(out);
	out <- out[order(out$p),];
}

LogLevels <- function(data,minN=4) { # first column of data must be locusIds; control should be first
	n <- apply(!is.na(data[,2:length(data)]),1,countif);
	nmax <- length(data)-1;
	swept <- sweep(data[,2:length(data)],2,apply(data[n==nmax,2:length(data)],2,mean));
	control <- apply(swept[,1:(nmax/2)],1,mean,na.rm=TRUE);
	treat <- apply(swept[,(nmax/2 + 1):nmax],1,mean,na.rm=TRUE);
	out <- data.frame(locusId=data[,1],llControl=control,llTreatment=treat);
	out <- out[n>=minN,];
	out <- out[order(out$locusId),];
	return(out);
}

# compare changers plots. Assumes outP is already set
OpWiseRepVsNoRep <- function(data,nConf=8) {
	changersP <- changersVsOperonPairs(data$locusId,logodds(data$outP$p0),data$oppairs,both=TRUE);
	confidence <- logodds2p(abs(changersP$statistic));
	cfactor <- cutN(confidence,n=nConf);
	CMeans <- sapply(split(confidence,cfactor),mean);
	PAgrees <- sapply(split(changersP$bAgree,cfactor),mean);
	PAgreesConf <- t(sapply(split(changersP$bAgree,cfactor), function(x) t.test(x)$conf.int));
	PAgreesAdj <- sapply(split(changersP$raw/changersP$pOp,cfactor),mean);
	PAgreesAdjConf <- t(sapply(split(changersP$raw/changersP$pOp,cfactor), function(x) t.test(x)$conf.int));
	OperonAgreementVsConfidencePlot(data);
	points(CMeans,PAgreesAdj,col=2);
	lines(CMeans,PAgreesAdj,col=2);
}

ShowPubmedId <- function(id,browser="mozilla") {
	system(paste(browser, " http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve\\&DB=pubmed\\&dopt=abstract\\&list_uids=",as.character(id),sep=""));
}

# only if all on one line
ReadPhylip = function(file) {
	d = read.table(file,header=F,skip=1);
	namesd = d[,1];
	d=d[,2:ncol(d)];
	names(d)=namesd;
	row.names(d) = namesd;
	return(d);
}

# given a list of values, returns the Chao2 estimator
Chao2 = function(x, counts = as.data.frame(table(table(x))), warn=TRUE) {
	if(sum(counts$Var1==1) != 1 || sum(counts$Var1==2) != 1) {
		if(warn) cat("Warning: Not enough values with distribution of 1 or 2 in Chao2\n");
		return(NULL);
	}
	f1 = counts$Freq[counts$Var1==1];
	f2 = counts$Freq[counts$Var1==2];
	tot = sum(counts$Freq);
	return( tot + f1**2/(2*f2) );
}

# Given a sorted list of values, and a vector of lengths to
# use to subsample the data, reports how the Chao2 estimator varies with the
# length. Returns a data frame
Chao2Curve = function(x, lengths=unique(floor(length(x)*seq(0.01,1.0,0.01)))) {
	out = NULL;
	for(len in lengths) {
		estimate = Chao2(x[1:len],warn=FALSE);
		if(!is.null(estimate)) {
			row = data.frame(len=len,nuniq=length(unique(x[1:len])),Chao2=estimate);
			out = if(is.null(out)) as.data.frame(row) else rbind(out,row);
		}
	}
	return(out);
}

# Not sure if this is still used
COGSimToCode <- function(x) { ifelse(x=="","U",ifelse(x=="-","N","Y")); }

# returns single character codes if identical, empty strings if missing or unknown categories,
# or "-" if different
# Assumes factors as input
COGSimilarity <- function(x1,x2) {
	out <- rep("",length(x1));
	x1 <- as.character(x1);
	x2 <- as.character(x2);
	x1[is.na(x1)]<-"";
	x1[x1 %in% c("R","S","-")]<-"";
	x2[is.na(x2)]<-"";
	x2[x2 %in% c("R","S","-")]<-"";
	v1 <- strsplit(x1,"");
	v2 <- strsplit(x2,"");

	for (i in 1:length(x1)) { out[i] <- paste(intersect(v1[[i]],v2[[i]]),collapse=""); }
	out[out=="" & x1!="" & x2!=""] <- "-";
	return(out);
}

# for each element of x, if it is near an element in y, report the nearest one
# returns a merged data set or, if merge=F, a list of pairs of indices
findNearby = function(x, y, by.x, by.y=by.x, within=1) {
	if(is.null(y[[by.y]])) stop("no field named",by.y," in y argument");
	y2 = y[order(y[[by.y]]),];

	if (nrow(y) == 1) {
		closestI = 1;
	} else {
		closestI = round(approx(y2[[by.y]],1:nrow(y2), xout=x[[by.x]], rule=2,
			ties='ordered')$y);
	}
	closestV = y2[closestI,by.y];
	keep = andNoNA(abs(closestV-x[[by.x]]) <= within);
	# the call to data.frame forces the names to be unique
	return(data.frame(cbind(x[keep,], y2[closestI[keep],])));
}

# given data that is subgrouped by the same markers, do findNearby on each and merge the results
findNearbyGrouped = function(splitx, splity, by.x, by.y=by.x, within=1) {
	out = NULL;
	for(i in names(splitx)) {
		if (!is.null(splity[[i]])) {
			rows = findNearby(splitx[[i]], splity[[i]], by.x, by.y, within);
			if(nrow(rows) > 0) {
				out = if(is.null(out)) rows else rbind(out,rows);
			}
		}
	}
	return(out);
}

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

# match rows in x to the rows in y such that y[[by.y]] >= x[[by.x]], maximizing x[[by.x]] for each y
# note requires >1 column in each
# note x may match to more than one column in y
# If sign is -1, then calls for less, not more
findAfter = function(x, y, by.x, by.y, sign=1) {
	x = x[order(x[[by.x]]*sign),];
	if (nrow(x) == 1) {
		# Cannot use approx
		iy = rep(1, nrow(y));
	} else {
		iy = floor(approx(x[[by.x]]*sign, 1:nrow(x), xout=y[[by.y]]*sign, rule=2, ties='ordered')$y);
	}
	use = y[[by.y]]*sign >= x[[by.x]][iy] * sign;
	return(data.frame(x[iy,], y, check.names=T)[use,]);
}

findAfterSS = function(x, y, by.x, by.y) {
	xs = split(x, x[,c("scaffoldId","strand")]);
	ys = split(without(y, c("scaffoldId","strand")), y[,c("scaffoldId","strand")]);
	out = NULL;
	for (i in names(xs)) {
		if (!is.null(ys[[i]])) {
			sign = ifelse(xs[[i]]$strand[1] == "+", 1, -1);
			rows = findAfter(xs[[i]], ys[[i]], by.x, by.y, sign=sign);
			if (nrow(rows) > 0) out = if(is.null(out)) rows else rbind(out,rows);
		}
	}
	return(out);
}

# given that data is subgrouped by the same markers, do findWithin on each and merge the results
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

# like a filter() -- a fast estimate to local correlation of data to rep(0, sz), rep(1, sz)
# not sure why not exactly the same (numerical problem?), but very close in practice
fastlocalR = function(data, sz=50, circular=T) {
	f = c(rep(0,sz),rep(1,sz));
	f = (f-mean(f))/sd(f);
	xy = filter(data,f,circular=circular);
	meanx = filter(data,rep(1,2*sz)/(2*sz),circular=circular);
	meanxx = filter(data**2,rep(1,2*sz)/(2*sz),circular=circular)
	return( -xy/(sqrt(meanxx - meanx**2)*2*sz) );
}

# local peak-finding; use -argument to get minima
peaks = function(series,span=3)
{
	if (span %% 2 != 1) stop("span must be odd in peaks() -- ", span);
	z <- embed(series, span);
	s <- span%/%2;
	v<- max.col(z) == 1 + s;
	return(c(rep(FALSE,s),v,rep(FALSE,s)));
}

# given boolean input, aggregate regions with T or F; returns dataframe with indices ifirst, ilast
trueRegions = function(b) {
	len = length(b);
	ifirst = (2:len)[b[2:len] & !b[1:(len-1)]];
	if (b[1]) { ifirst = c(1, ifirst); }
	ilast = (1:(len-1))[b[1:(len-1)] & !b[2:len]];
	if (b[len]) { ilast = c(ilast, len); }
	if (length(ifirst) < 1 || length(ilast) < 1) return(NULL);
	return(data.frame(ifirst=ifirst[1:length(ilast)], ilast=ilast));
}

# given a frame with begin, end, normll, for a single strand and scaffold,
# smooth by the indicated number of adjacent probes and then
# identify regions that are above the threshold
onRegions = function(d, name, thresh=median(d[[name]]), smooth=40) {
	d = d[!is.na(d[[name]]),];
	d = d[order(d$end),];
	d$smooth = filter(d[[name]], rep(1,smooth)/smooth, circular=T);
	d$mid = (d$begin+d$end)/2;
	r = trueRegions(d$smooth > thresh);
	# and compute the midpoint and the average for each region
	r$begin = d$mid[r$ifirst];
	r$end = d$mid[r$ilast];
	r$avg = sapply(1:nrow(r), function(x) mean(d[r$ifirst[x]:r$ilast[x],name]));
	return(r[,c("begin","end","avg")]);
}

# run onRegions separately for each strand and scaffold and combine the results
onRegionsSS = function(d, name="normll", ...) {
	s = split(d, d[,c("strand","scaffoldId")]);
	out = NULL;
	for (i in names(s)) {
		rows = onRegions(s[[i]], name, ...);
		rows = cbind(scaffoldId = s[[i]]$scaffoldId[1], strand = s[[i]]$strand[1], rows);
		out = if(is.null(out)) rows else rbind(out, rows);
	}
	return(out);
}

# span2 is to avoid reporting multiple local maxima
peaksByStrandScaffold = function(data, name="normll", thresh=0.6, sz=50, span2=21, circular=T) {
	if(is.null(data[[name]])) fail("bad name");
	data$index = 1:nrow(data);
	data$up = FALSE;
	data$down = FALSE;
	data$r = NA;
	d = split(data,data[,c("strand","scaffoldId")]);
	if(is.null(data)) fail("bad data");
	
	for(i in names(d)) {
		v = d[[i]];
		v = v[order(v$end),];
		v = v[!is.na(v[[name]]),];
		r = fastlocalR(v[[name]],circular=circular,sz=sz);
		cat(nrow(v), length(r),"\n");

		data$r[v$index] = r;
		sign = if(v$strand[1] == "+") 1 else -1;
		data$up[v$index] = peaks(sign*r, span=span2) & sign*r >= thresh;
		data$down[v$index] = peaks(-sign*r, span=span2) & -sign*r >= thresh;
	}
	return(data[data$up | data$down,]);
}

# Both data and normalizeby should already be log-transformed
# Both should be vectors
# After subtracting out loess effects it then forces the median to be 0
# Reports NA for values with normalizeby whose percentile is below the argument
NormalizeBy = function(data, normalizeby, scaffold=NULL, percentile=0.02) {
	if (length(data) != length(normalizeby)) stop("Data lengths do not match in NormalizeBy");
	limit = quantile(normalizeby, percentile);
	l = lowess(normalizeby[normalizeby >= limit], data[normalizeby >= limit]);
	out = data - approx(l$x, l$y, xout=normalizeby, rule=2, ties="ordered")$y;
	out[normalizeby < limit] = NA;
	if (is.null(scaffold)) {
		out = out - median(out, na.rm=T);
	} else {
		# set median separately for each scaffold
		if (length(data) != length(scaffold)) stop("Scaffold length incorrect in NormalizeBy");
		out = unsplit(lapply(split(out, scaffold), function(x) { x - median(x, na.rm=T) }), scaffold);
	}
	return(out);
}

# data has strand, scaffold, mid, name
# ranges has strand, scaffold, begin, end, and perhaps other fields
RangeAverages = function(data, range, name) {
	splitd = split(data,data[,c("strand","scaffoldId")]);
	range$density = NA;
	range$avg = NA;
	splitr = split(range,range[,c("strand","scaffoldId")]);
	for (n in names(splitd)) {
		d = splitd[[n]];
		d = d[order(d$mid),];
		r = splitr[[n]];
		if (!is.null(r)) {
			beginI = floor(0.5 + approx(d$mid, 1:nrow(d), xout = r$begin, rule=2)$y);
			endI = floor(0.5 + approx(d$mid, 1:nrow(d), xout = r$end, rule=2)$y);
			splitr[[n]]$density = (endI-beginI+1)/(r$end-r$begin+1);
			splitr[[n]]$avg = sapply(1:nrow(r), function(x) mean(d[beginI[x]:endI[x], name], na.rm=TRUE));
		}
	}
	return(unsplit(splitr, range[,c("strand","scaffoldId")]));
}

# returns norm and quantile for the same rows as x & by
# by default, drops the lowest fraction
# use keep to keep only a subset of rows before normalizing
QuantileNormalize = function(x, by, drop=0.01, keep=TRUE) {
	d = data.frame(index = 1:length(x), x=x, by=by, quantile=NA);
	d = d[keep,];
	d$index2 = 1:nrow(d);

	d$quantile = approx(x=sort(d$x), y=(1:nrow(d))/nrow(d), xout=d$x, rule=2)$y;
	d$norm = approx(x=(1:length(by))/length(by), y=sort(by), xout=d$quantile, rule=2)$y;
	d$quantile[d$quantile < drop] = NA;
	d$norm[is.na(d$quantile)] = NA;

	out = data.frame(index = 1:length(x), norm=NA, quantile=NA);
	out$norm[d$index] = d$norm;
	out$quantile[d$index] = d$quantile;
	return(out[,c("norm","quantile")]);
}

# averages, or alternate function, on the columns of x that have matching names in names
ByRowByName = function(x, names, func, ...) {
	if (ncol(x) != length(names)) stop("Columns and names do not match in length")
	out = list();
	u = unique(names);
	for (n in u) {
		out[[n]] = apply(cbind(x[,names == n]), 1, func, ...);
	}
	out = data.frame(out);
	names(out) = u;
	return(out);
}

# range of correlations for columns of x that have matching names in names
CorByName = function(x, names, use="pairwise", skip=T) {
	if (ncol(x) != length(names)) stop("Columns and names do not match in length")
	d = split(1:ncol(x), names);
	out = lapply(d,
		function(col) { if(length(col) == 1) return(c(n=1, mincor=NA, medcor=NA, maxcor=NA));
				cors = cor(x[,col], use=use);
				diag(cors)=NA;
				r = range(cors,na.rm=T);
				return(c(n=length(col), mincor=r[1], medcor=median(r,na.rm=T), maxcor=r[2])); });
	out = data.frame(name=names(d), t(as.data.frame(out)));
	row.names(out) = 1:nrow(out);
	if (skip) return(out[out$n>1,]);
	return(out);
}

nunique = function(x) { if (is.data.frame(x)) nrow(unique(x)) else length(unique(x)); }

# Just like plot.default() but adds on the labels
# By default plotting symbols are off but you can override that
plotlab = function(x,y,labels, cex=1, col=1, pch="", ...) {
	plot.default(x,y,cex=cex,col=col,pch=pch,...);
	text(x,y,labels,cex=cex,col=col);
}

# get exact matches to all of query froma  file
read_blat_matches = function(file) {
	d = read.delim(file, skip=5, col.names=c("match","mm","repmatch","Ns","nQgap","Qgapsz","nTgap","Tgapsz","strand","queryId","queryLen","pepbeg","pepend","locusId","locusLen","locusBeg","locusEnd","nBlock","bSizes","qStarts","tStarts"),stringsAsFactors=F);
	d = d[d$mm==0 & d$match==d$queryLen & d$nQgap==0 & d$nTgap==0,
		c("queryId","queryLen","locusId","locusLen","locusBeg","locusEnd")];
	d$locusBeg = d$locusBeg + 1; # off by one oddity of BLAT
	return(d);
}

# requires a "sorted" data frame with scaffoldId, strand, end, norm (or other data name),
# and index that matches the values in indexes
# sorted must be sorted by scaffold, strand, and end
SideAverages = function(sorted, indexes, off1, off2, minpoints=10, nt=250, name="norm") {
	sorted$index2 = 1:nrow(sorted);
	sorted2 = sorted[sorted$index %in% indexes,];
	return( sapply(indexes, function(index) {
		row = sorted2[sorted2$index == index,];
		if (is.null(row$scaffoldId)) stop("Invalid index ",index);
		i = (row$index2+off1):(row$index2+off2);
		i = i[i >= 1 & i <= nrow(sorted)];
		i = i[sorted$scaffoldId[i] == row$scaffoldId
			& sorted$strand[i] == row$strand
			& abs(sorted$end[i] - row$end) <= nt
			& !is.na(sorted[i,name])];
		if (length(i) < minpoints) return(NA);
		return(mean(sorted[i,name]));
	}) );
}

# Given two data frames x and y that contain objects associated with scaffolds and strands,
# use findWithin() to find pairs of objects such that
# x$by.x >= y$begin & x$by.x <= y$end
# Like merge(), it returns a data frame that contains both sets of columns
findWithinByLoc = function(x, y, by.x, begin, end, scaffoldCol="scaffoldId", ...) {
	return( findWithinGrouped(split(x, x[,c(scaffoldCol,"strand")]),
			split(without(y,c(scaffoldCol,"strand")), y[,c(scaffoldCol,"strand")]),
			by.x=by.x, begin=begin, end=end, ...));
}

findNearbyByLoc = function(x, y, ...) {
	return( findNearbyGrouped(split(x, x[,c("scaffoldId","strand")]),
			split(without(y,c("scaffoldId","strand")), y[,c("scaffoldId","strand")]), ...));
}

t2 = function(x) t.test(x$"TRUE",x$"FALSE",paired=F);

# Assumes that where and value line up and are sorted on where
# at is a list of locations (approximately matching values in where) to compute a slope near
# offsets are index offsets, e.g. 1:20
slopeAround = function(where, value, at, offsets, minpoints=length(offsets)/2) {
	if (length(at) == 0) return(NULL);
	out = data.frame(index = findInterval(at, where), r=NA, slope=NA, pslope=NA, avg=NA, wmin=NA, wmax=NA);
	for (i in 1:nrow(out)) {
		off = offsets + out$index[i];
		off = off[off >= 1 & off <= length(value)];
		off = off[!is.na(value[off])];
		if (out$index[i] > 0 && length(off) >= minpoints) {
			x = where[off];
			y = value[off];
			out$avg[i] = mean(y);
			test = cor.test(x,y);
			out$r[i] = test$estimate;
			out$pslope[i] = test$p.value;
			out$slope[i] = lm(y~x)$coefficients[2];
			out$wmin[i] = min(x);
			out$wmax[i] = max(x);
		}
	}
	return(out);
}

# data is not assumed to be sorted, with wherecol and valuecol being the names of the
# relevant columns and scaffoldId/strand being presnt in data and also in loc
# atcol is the column name in loc to try to match up with wherecol
slopeAroundSS = function(data, loc, offsets, wherecol, valuecol, atcol, minpoints=length(offsets)/2) {
	d = split(data.frame(where=data[[wherecol]], value=data[[valuecol]]), data[,c("scaffoldId","strand")])
	a = split(loc, loc[,c("scaffoldId","strand")]);
	out = unsplit(lapply(names(a),  function(n) {
		d2 = d[[n]];
		d2 = d2[!is.na(d2$value),];
		d2 = d2[order(d2$where),];
		sign = ifelse(a[[n]]$strand[1] == "+", 1, -1);
		cbind(a[[n]], slopeAround(d2$where, d2$value, a[[n]][[atcol]], sign*offsets, minpoints=minpoints));
	}), loc[,c("scaffoldId", "strand")]);
	out$slope = ifelse(out$strand=="+",1,-1)*out$slope;
	out$r = ifelse(out$strand=="+",1,-1)*out$r;
	return(out);
}

# Find rows of data where the selected column is above or below the others by at least thresh
# Returns a boolean array
specificTo = function(data, col, thresh=0.5) {
	if (is.character(col)) col = match(col, names(data));
	if (any(!is.integer(col)) || any(col < 1) || any(col > ncol(data))) stop("Invalid column");
	mins = apply(data[,-col], 1, min, na.rm=T);
	maxs = apply(data[,-col], 1, max, na.rm=T);
	if (length(col) == 1) values  = data[,col] else values = apply(data[,col], 1, mean );
	return(andNoNA(values < mins-thresh | values > maxs+thresh));
}

# input is usually a matrix with matching row and column names such as produced by cor(x),
# but any numeric matrix should work (regardless of matching of names)
# requires ggplot2 (and R-2.11)
# ... is for opts, e.g. use title to set the title
CorPlot = function(inp, names=c("Var1","Var2","r"),
	breaks=seq(-1,1,0.1), limits=range(breaks), midpoint=mean(limits),
	low="red", mid="black", high="green",
	...) {
	m = melt(inp);
	names(m) = names;
	m[,1] = as.character(m[,1]);
	m[,2] = as.character(m[,2]);
	(ggplot(m, aes_string(x=names[1],y=names[2])) +
		geom_tile(aes_string(fill=names[3])) +
		scale_fill_gradient2(names[3], limits=limits, breaks=breaks, low=low, mid=mid, high=high, midpoint=midpoint) +
		opts(axis.text.x=theme_text(angle=90,hjust=1), ...));
}

RedToBlackToGreenColors = function() {
	c(rgb(seq(1,1/12,-1/12), 0, 0),  rgb(0, seq(0/12,1,1/12), 0));
}

# variant of CorPlot() based on heatmap
CorPlot2 = function(inp,
		col=RedToBlackToGreenColors(),
		breaks=seq(-1,1,2/(length(col))),
		Rowv=NA, Colv=NA,
		 ...) {
	heatmap(inp, Rowv=Rowv, Colv=Colv, scale="none", col=col, breaks=breaks, na.rm=F, ...);
}

# Given a list of Gene1 Gene2 pairs, and a matrix of data (as genes and data-only matrix),
# compute correlations for each pair or NA
cor12 = function(pairs, genes, data, use="p", method="pearson", names=c("Gene1","Gene2")) {
	i1 = match(pairs[[names[1]]], genes);
	i2 = match(pairs[[names[2]]], genes);
	return(sapply(1:nrow(pairs), function(x) if(is.na(i1[x]) | is.na(i2[x])) NA else
		cor(c(data[i1[x],], recursive=T), c(data[i2[x],], recursive=T), method=method, use=use)));
}

# Given a list of Gene1 Gene2 pairs, and a matrix of data (as genes and data-only matrix),
# compute a function for each pair or NA
func12 = function(pairs, genes, data, func="cor", names=c("Gene1","Gene2"), ...) {
	i1 = match(pairs[[names[1]]], genes);
	i2 = match(pairs[[names[2]]], genes);
	return(sapply(1:nrow(pairs), function(x) if(is.na(i1[x]) | is.na(i2[x])) NA else
		func(c(data[i1[x],], recursive=T), c(data[i2[x],], recursive=T), ...)));
}

# are values correlated for pairs? Intended for operon pairs but would work with other types too
paircor = function(pairs, genes, values, use="p", method="pearson", names=c("Gene1","Gene2")) {
	d = merge(merge(pairs[,names], data.frame(Gene1=genes, value1=values), by.x=names[1], by.y="Gene1"),
			data.frame(Gene2=genes, value2=values), by.x=names[2], by.y="Gene2");
	return(cor(d$value1, d$value2, use=use, method=method));
}

# are values correlated for pairs? Intended for operon pairs but would work with other types too
pairfunc = function(pairs, genes, values, func=cor, names=c("Gene1","Gene2")) {
	d = merge(merge(pairs[,names], data.frame(Gene1=genes, value1=values), by.x=names[1], by.y="Gene1"),
			data.frame(Gene2=genes, value2=values), by.x=names[2], by.y="Gene2");
	return(func(d$value1, d$value2));
}

# strains should include strain and pool.
straincor = function(strains, values, use="p", method="pearson") {
	d = data.frame(row=1:nrow(strains), strain=strains$strain, pool=strains$pool);
	d = merge(d[d$pool=="up",], d[d$pool=="dn",], by="strain", suffixes=c(".up",".dn"));
	return(cor(values[d$row.up], values[d$row.dn], use=use, method=method));
}

# assumes that the data table includes scaffoldId, strand, and mid, is sorted by those
# assumes that the ranges table includes scaffoldId, strand, begin, end
fetchValuesForRanges = function(data, ranges,
		FUN = mean,
		name="norm", begin="begin", end="end",
		debug = FALSE,
		 ...) {
	if (is.null(data[[name]])) stop("No data column named",norm);
	if (is.null(data$mid)) stop("No data column named mid");
	r = ranges;
	r$i = 1:nrow(r);
	r$b = pmin(r[[begin]],r[[end]]);
	r$e = pmax(r[[begin]],r[[end]]);
	rSplit = split(r, r[,c("scaffoldId","strand")], drop=T);
	indexes = split(1:nrow(data), data[,c("scaffoldId","strand")]);

	ranges$i1 = NA;
	ranges$i2 = NA;
	ranges$value = NA;

	iSplit = list();
	for (n in names(rSplit)) {
		r2 = rSplit[[n]];
		dataindex = indexes[[n]];
		if(is.null(dataindex)) stop("No data for",n);
		bi = findInterval(r2$b, data$mid[dataindex], all.inside=T);
		ei = findInterval(r2$e, data$mid[dataindex], all.inside=T);
		# note drop above to ensure this does not get invoked
		for (row in 1:nrow(r2)) {
			i = r2$i[row];
			ranges$i1[i] = dataindex[bi[row]];
			ranges$i2[i] = dataindex[ei[row]];
			if(debug) cat("Name",n,"Range",r2$b[row],r2$e[row],"becomes",ranges$i1[i],ranges$i2[i],"name",name,"\n");
			ranges$value[i] = FUN(data[dataindex[bi[row]:ei[row]],name], ...);
		}
	}
	return(ranges);
}

# returns a matrix with columns as -span to span and rows as
# the centers of the adjacent items
embedCircular = function(series, span=5) {
	n = length(series);
	if (span < 0) stop("embedCircular requires nonnegative span");
	if (span >= n-1) stop("span too high in embedCircular");
	d = sapply(-span:span, function(off) series[if(off>0) c((off+1):n, 1:off) else if (off<0) c((n+off+1):n,1:(n+off)) else 1:n]);
	return(d);
}

# Make an area under the receiver operating characteristic curve (AOC) plot
# Fast because it uses interpolation, so it is NOT exactly accurate in the case of ties,
#	but the error is small if you have >200 points
# If you ask for the values, it returns the sensitivities [fraction of trues above threshold]
#	and the false positive rate [fraction of falses above threshold]
AOCPlot = function(trues, falses=NULL, n=NULL, value=F, new=T,
		xlab="False Positives", ylab="True Positives", xlim=0:1, ylim=0:1, ...) {
	if (is.null(falses)) { falses = trues$"FALSE"; trues = trues$"TRUE"; }
	trues = withoutNA(trues);
	falses = withoutNA(falses);
	if (!(length(trues) > 1) || !(length(falses) > 1)) stop("too few values in AOCPlot");
	if(is.null(n)) { n = floor(1.5*(length(trues) + length(falses))); if(n>2000) n=2000; }
	frac = seq(0,1,length.out=n);
	at = quantile(c(trues,falses), frac);
	# sensitivity is fraction of true values above those thresholds, or 
	# 1 minus the quantile of true values, which is quickly estimated with approx()
	# Use ties=max in case of a large clump of values at (say) the minimum threshold
	sens = 1 - approx(quantile(trues,frac), frac, xout=at,rule=2,ties=max)$y;
	# And similarly for the false positive rate, or fraction of false values above thresholds
	fps = 1 - approx(quantile(falses,frac), frac, xout=at,rule=2,ties=max)$y;
	if(new) plot(fps,sens,type="l",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...) else lines(fps,sens,...);
	if(value) return(data.frame(thresh=at,sens=sens,fps=fps));
}

# Compare mulitple AOC curves in one plot, with a data frame of multiple predictors (data) and
# a binary indicator of the true prediction (by).
# If showLine, adds an x=y line for non-informative predictions.
AOCPlots = function(data, by, vars=names(data), labels=vars,
			col=1:length(vars), lty=1:length(vars), lwd=rep(1,length(vars)),
			xlim=0:1, ylim=0:1, showLegend=TRUE,
			xlab="False Positives", ylab="True Positives",
			showLine=TRUE,
			legendX=mean(xlim), legendY=ylim[1]+diff(ylim)/10, xjust=0, yjust=0, ...) {
	for(i in 1:length(vars)) {
		AOCPlot(split(data[[i]], by), col=col[i], lty=lty[i], lwd=lwd[i],
				xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, new=(i==1), ...);
	}
	if(showLine) eqline();
	if(showLegend)	legend(legendX, legendY, labels,col=col,lty=lty,lwd=lwd,pch="", xjust=xjust, yjust=yjust);
}

# Uncentered correlation, i.e. mean(xy)/sqrt(mean(x**2) * mean(y**2)) is not centered like the
# Pearson correlation is
uncenterCor = function(x,y) {
	if (length(x) != length(y) || length(x) < 3 || length(y) < 3) stop("Wrong lengths in uncenterCor:",length(x)," ",length(y));
	v1 = sum(x**2)/(length(x));
	v2 = sum(y**2)/(length(x));
	return( mean(x*y) / sqrt(v1*v2) );
}

# split a string separated by spaces into a list of word
words = function(s, by=" ", ...) { strsplit(s[1], by, ...)[[1]]; }

# returns TRUE or FALSE indicating whether the value is repeated
# input should be a list of values
is.unique = function(x) {
	counts = as.data.frame.table(table(x), dnn="x");
	return(x %in% counts$x[counts$Freq==1]);
}

# Test a subset of genes for enrichment or depletion of various types of proteins
# Handles various types of functional classifications at once
# Returns a data frame with one row per functional classification
# geneset should be a list of locusIds
# protinfo should include locusId and a value for each column in by
#	usually do not include non-protein-coding genes
# COG is handled specially -- COGfunctions is used to map these values to function codes
#	COGfunctions must include COG and funCode
# The values are in the output are:
#	class -- classification type (TIGR, COG, or conserv)
#	enrich -- TRUE or FALSE (FALSE means depleted)
#	propin -- proportion of genes in the set that are in this category
#	propout -- proportion of genes not in the set that are in this category
#	n -- number of genes in the category, regardless of being in the set or not
#	fset -- proportion of genes in the category that are in the set
#	fset.low & fset.high -- confidence interval for fset (95% confidence by default)
#	p -- P-value from Fisher exact test. By default, these are adjusted by fdr and
#		rows with adjusted p-values above 0.05 are ignored
#	cat -- string value of the category
#
testGenes = function(geneset, protinfo, COGfunctions = cogfuns, by=intersect(names(protinfo),c("TIGR","COG","conserv")),
			p=0.05, adjust=TRUE, adjust.method="fdr", debug=FALSE,
			conf.level=0.95, filter=TRUE) {
	if (length(by)==0) stop("No columns to test on in testGenes() -- check protinfo");
	if (any(!by %in% names(protinfo))) stop("Invalid columns to test on in testGenes");
	geneset = intersect(geneset,protinfo$locusId);
	if (length(geneset)==0) return(NULL);

	out = NULL;
	nMax = 0;
	for (class in by) {
		if (class=="COG") {
			data = merge(protinfo[,c("locusId","COG")], COGfunctions, by="COG");
			names(data)[names(data)=="funCode"] = "cat";
		} else {
			data = data.frame(locusId=protinfo$locusId, cat = protinfo[[class]]);
		}
		data = data[!is.na(data$cat),];
		data$inset = data$locusId %in% geneset;
		if(debug) cat("By",class,"Rows",nrow(data),"geneset rep",sum(data$inset),"\n");
		if (any(data$inset)) {
			tables = as.data.frame(table(data$inset,data$cat)); names(tables) = c("inset","cat","n");
			tables = merge(tables[tables$inset=="TRUE",], tables[tables$inset=="FALSE",], suffixes=c("In","Out"), by="cat");
			nInTot = length(unique(data$locusId[data$inset]));
			nTot = length(unique(data$locusId));
			nOutTot = nTot - nInTot;
			tables$propin = tables$nIn/nInTot;
			tables$propout = tables$nOut/nOutTot;

			nBoth = tables$nIn + tables$nOut;
			tables$n = nBoth;
			tables$fset = tables$nIn/nBoth;
			tables$fset.low = sapply(1:nrow(tables), function(x)
				if(nBoth[x]==0) NA else binom.test(tables$nIn[x], tables$nIn[x]+tables$nOut[x],
								conf.level=conf.level)$conf[1]);
			tables$fset.high = sapply(1:nrow(tables), function(x)
				if(nBoth[x]==0) NA else binom.test(tables$nIn[x], tables$nIn[x]+tables$nOut[x],
								conf.level=conf.level)$conf[2]);

			tables$p = apply(tables[,c("nIn","nOut")], 1, function(d) {
				fisher.test(matrix(c(d[1], nInTot-d[1], d[2], nOutTot-d[2]),ncol=2))$p.value; } );

			if(debug) cat("p values",tables$p,"\n");
			tables$p = p.adjust(tables$p, method=adjust.method);
			if(filter) tables = tables[tables$p < p,];
			if(nrow(tables) > 0) {
				if(debug) cat("By",class,"Sig",nrow(tables),"\n");
				tables$enrich = tables$propin > tables$propout;
				tables$class = as.character(class);
				row.names(tables) = nMax+(1:nrow(tables));
				if(debug) cat("By",class,"rownames",paste(row.names(tables),collapse=" "),"\n");
				nMax = nMax + nrow(tables);
				tables$cat = as.character(tables$cat);
				d = tables[,c("class","enrich","propin","propout",
						"n", "fset","fset.low","fset.high","p","cat")];
				out = if(is.null(out)) d else rbind(out,d);
			}
		}
	}
	if(is.null(out)) { cat("No significant categories\n"); return(NULL); }
	out = out[order(out$enrich,out$class,out$p),];
	return(out);
}

# Plot the proportion of genes in each category that are in the set, with confidence intervals (95% by default)
# The input should be the result of testGenes, usually with filter=FALSE and COG/conserv only
# and perhaps with smaller categories removed
plotTestGenes = function(tests,
		ylab="Gene Category", xlab="Proportion", xlim=c(0,max(tests$fset.high)),
		col=1,
		pch=ifelse(tests$p<0.05, 8, 1),
		arrowLength=0.03,
		...) {
	plot(tests$fset, 1:nrow(tests), xlab=xlab, ylab=ylab, xlim=xlim, yaxt="n", col=col, pch=pch, ...);
	arrows(tests$fset.low, 1:nrow(tests), tests$fset.high, 1:nrow(tests),
			code=3, angle=90, length=arrowLength, col=col);
	# side=2 means at right; las=2 means perpendicular (i.e., labels are horizontal)
	if(length(col)==1) {
		axis(2, at=1:nrow(tests), labels=tests$cat, las=2, col.axis=col);
	} else {
		d = split(1:nrow(tests), col);
		for (i in 1:length(d)) {
			axis(2, at=d[[i]], labels=tests$cat[ d[[i]] ], las=2, col.axis=col[ d[[i]][1] ]);
		}
	}
}

testCorFisher = function(n, r, rho=0, conf.level=0.95) {
	# Fisher transform and divide by sigma = 1/sqrt(n-3)
	z = 0.5 * log((1+r)/(1-r));
	mu = 0.5 * log((1+rho)/(1-rho));
	sigma = 1/sqrt(n-3);
	zStandard = (z-mu)/sigma;
	p = 2.0 * pnorm(-abs(zStandard));
	q = (1 - conf.level)/2;
	confZ = z + sigma * qnorm(c(q,1-q)); # offset by z to get confidence interval around the observation
	# transform back to correlation space
	confR = (exp(2*confZ)-1)/(exp(2*confZ)+1);
	return(list(p=p,conf.int=confR));
}


# at each position p, the correlation between x and y over the range of positions p+offset1 to p+offset2
#	(or p-offset1 to p-offset2 if on the minus strand),
# or NA if not enough data points are available.
# On either strand, a rise will give a positive value and a drop will give a negative one
LocalCorrelationAt = function(at, x, y, strand, offset1, offset2, minN=10) {
	LocalFuncAt(at, x, y, strand, offset1, offset2, function(x,y,ssign) ssign*cor(x,y), minN=minN);
}

# at each position p, returns the average value of y for those with x in p+offset1 to p+offset2
#	(or p-offset1 to p-offset2 if on the minus strand)
# or NA if not enough data points are available
LocalLevelAt = function(at, x, y, strand, offset1, offset2, minN=5) {
	LocalFuncAt(at, x, y, strand, offset1, offset2, function(x,y,ssign) mean(y), minN=minN);
}

# at each position p, compute the function 
# func should take in arguments subsetX, subsetY, and ssign, where ssign is 1 for the + strand and -1 for the minus strand
# Note the rounding of the indices implies that the last point chosen could be 1 off in the series
LocalFuncAt = function(at, x, y, strand, offset1, offset2, func, minN=2) {
	ssign = if (strand=="+") 1 else -1;

	iPre = round(approx(x, 1:length(x), xout=at + ssign*offset1, rule=2, ties='ordered')$y);
	iPost = round(approx(x, 1:length(x), xout=at + ssign*offset2, rule=2, ties='ordered')$y);

	len = length(x);
	iPre = pmin(len, pmax(1, iPre));
	iPost = pmin(len, pmax(1, iPost));

	sapply(1:length(at), function(i) {
		i1 = iPre[i];
		i2 = iPost[i];
		# if have 3 points, then i1-i2 = +-2
		if (abs(i1-i2) < minN-1) return(NA);
		return(do.call(func, list(x[i1:i2],y[i1:i2],ssign)))
	});
}

CountACG = function(rcbarcodes) {
	 nA = sapply(rcbarcodes, function(x) sum(strsplit(x,"")[[1]] == "A"));
	 nC = sapply(rcbarcodes, function(x) sum(strsplit(x,"")[[1]] == "C"));
	 nG = sapply(rcbarcodes, function(x) sum(strsplit(x,"")[[1]] == "G"));
	 return(data.frame(nA=nA,nC=nC,nG=nG));
}

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
    if (!any(strainsUsed)) error("No usable strains");

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

    geneFit1 = aggregate(strainFit, list(locusId=strainLocus), mean);
    geneFit1$x = mednorm(geneFit1$x);

    i = match(strainLocus, geneFit1$locusId); # from strain index to gene index
    strainFitPseudo = mednorm(log2(2**geneFit1$x[i] + strainCounts) - log2(1 + strainT0));
    strainFitPseudo = strainFitPseudo - strainFitAdjust; # correct for nt content

    readratio = sum(strainCounts) / sum(strainT0);
    geneFit2 = aggregate(1:length(strainT0), list(locusId=strainLocus),
    	                function(j) {
		 	    # for each strain: fitness, variance, and weight
			    fit = log2(2**strainFitPseudo[j] * readratio + strainCounts[j]) - log2(1 + strainT0[j]);
			    var = sqrt(1/(1+strainT0[j]) + 1/(1+strainCounts[j]))/log(2);
			    weight = 0.5 + strainT0[j];
			    meanFit = sum(weight*fit)/sum(weight);
			    c(fit = meanFit,
			    	  sd = sqrt(sum(weight**2 * var))/sum(weight),
				  sdObs = sqrt(sum(weight**2 * (fit-meanFit)**2))/sum(weight),
				  n = length(j),
				  # effective degrees of freedom is n if all weights are even,
				  # otherwise is this much higher
				  wbias = sqrt(length(j)) * sqrt(sum(weight**2)) / sum(weight) );
    });
    # a hack to get around the lists within the entries
    geneFit2 = data.frame(locusId = geneFit2$locusId, geneFit2[,2]);
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
    if(any(is.na(i))) error("Fitness data for loci not in genes");
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

prefixName = function(x, prefix) { names(x) = paste(prefix,names(x),sep=""); return(x); }

applyRules = function(rules, desc) {
    for (i in 1:nrow(rules)) desc = sub(rules[i,1], rules[i, 2], desc);
    return(desc);
}

mednorm = function(x) x - median(x);

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

FEBA_Fit = function(expsUsed, all, all_g2, genes, genesUsed=NULL, strainsUsed=NULL,
	   		       minT0Strain=3, minT0Gene=30,
			       minT0GeneSide=minT0Gene/2,
			       minGenesPerScaffold=10,
			       nACG=NULL,
			       pred=CrudeOp(genes),
			       okLane=FALSE, # OK to use t0 from another lane if needed?
	   		       metacol=1:7) {
	if(!all(expsUsed$name == names(all_g2)[-metacol]))
		stop("names do not match");
	expsUsed$name = as.character(expsUsed$name);
	if(is.null(genes$scaffoldId)) stop("No scaffold for genes");
	if(is.null(genes$begin)) stop("No begin for genes");

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
	q$cor12 = sapply(as.character(q$name), function(n) cor(all_lr1[[n]], all_lr2[[n]]));
	q$opcor = apply(all_lrn[,-1], 2, function(x) paircor(pred[pred$bOp,], all_lrn$locusId, x));
	q$maxFit = apply(all_lr[,-1],2,max);
	
	return(list(expsUsed = expsUsed,
			gN = all_gN,
			t0_g2 = t0_g2,
			t0_gN = t0_gN,
			lr = all_lr,
			sd = all_sd,
			nStrain = nStrain,
			lrn = all_lrn,
			lr1 = all_lr1,
			lr2 = all_lr2,
			nStrain12 = nStrain12,
			q = q,
			genesUsed = genesUsed,
			strainsUsed = strainsUsed));
}

CrudeOp = function(genes) {
	d = merge(merge(data.frame(Gene1=genes$locusId[-nrow(genes)],Gene2=genes$locusId[-1]), genes, by.x="Gene1", by.y="locusId"), genes, by.x="Gene2", by.y="locusId",suffixes=1:2);
	d = d[d$strand1==d$strand2,]
	d$Sep = pmin(abs(d$begin1-d$end2),abs(d$end1-d$begin2));
	d$bOp = d$Sep < median(d$Sep);
	return(d);
}

FEBA_Save_Tables = function(fit, genes, org="?", topdir="public_html/tmp",
		dir = paste(topdir,org,sep="/")) {
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

	d = merge(genes[,c("locusId","sysName","desc")], fit$lrn);
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_logratios.tab"));
	wroteName("fit_logratios_normalized.tab");

	d = merge(genes[,c("locusId","sysName","desc")], fit$sd);
	names(d)[-(1:3)] = paste(fit$q$name,fit$q$short);
	writeDelim(d, nameToPath("fit_standard_error.tab"));
	wroteName("fit_standard_error.tab");

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
		      sampleq$gMed, sampleq$cor12, sub("set","",sampleq$name), sampleq$short);

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
		      sampleq$gMed, sampleq$cor12, sub("set","",sampleq$name), sampleq$short);
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
}
