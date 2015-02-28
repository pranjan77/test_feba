#!/usr/bin/Rscript

usage = paste("",
              "Usage: StrainData.R netcdf4File locusId > output.table",
	      "    Given a netcdf4File that contains per-strain fitness values",
	      "    makes a tab-delimited table of the strains within that gene",
	      "    and their fitness values. Includes strains near the edges",
	      "    of the gene or with low coverage that were not used to",
	      "    estimate the gene's fitness (used=FALSE).",
	      "", sep="\n");

library("ncdf4");

StrainData = function(args = commandArgs(trailingOnly=TRUE)) {
	source(file.path(GetPath(), "NetCDFInterface_Functions.R"));

	if (length(args) != 2) stop(usage);
	filename = args[1];
	locusId = args[2];

	nc_file = nc_open(filename);
	strainCols = dimensionsForVariable(nc_file,"fit/strains")[[3]]$vals;
	strains = ncvar_get(nc_file,"fit/strains");
	locusIdCol = which(strainCols == "locusId");
	# values are sometimes padded with spaces
	strainLocus = sub("^ +", "", strains[,locusIdCol], perl=T);
	strainLocus = sub(" +$", "", strainLocus, perl=T);
	lrnCols = dimensionsForVariable(nc_file, "fit/strain_lrn")[[2]]$vals;
	outCols = c(strainCols, lrnCols);
	iRows = which(strainLocus == locusId);
	if (length(iRows) == 0) {
		# empty table
		write(outCols, stdout(), sep="\t", ncol=length(outCols));
	} else {
		values = sapply(iRows, function(i) ncvar_get(nc_file,"fit/strain_lrn",start=c(i,1),count=c(1,-1)));
		values = t(values) # now, 1 column per fitness experiment, one row per strain
		d = cbind(strains[iRows,], values);
		colnames(d) = outCols;
		for (i in 1:length(strainCols)) {
		    d[,i] = sub("^ +", "", d[,i], perl=T);
		    d[,i] = sub(" +$", "", d[,i], perl=T);
		}
		write.table(d, stdout(), quote=F, row.names=F, sep="\t");
	}
}

GetPath = function() {
    argv = commandArgs(trailingOnly = FALSE)
    dirname(substring(argv[grep("--file=", argv)], 8));
}

# Actually do the work
if(!interactive()) {
	StrainData();
	quit();
}
