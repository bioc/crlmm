setAs("SnpSuperSet", "CNSet",
      function(from, to){
	      CA <- CB <- matrix(NA, nrow(from), ncol(from))
	      dimnames(CA) <- dimnames(CB) <- list(featureNames(from), sampleNames(from))		  
	      new("CNSet",
		  call=calls(from),
		  callProbability=assayData(from)[["callProbability"]],  ##confs(from) returns 1-exp(-x/1000)
		  alleleA=A(from),
		  alleleB=B(from),
		  CA=CA,
		  CB=CB,
		  phenoData=phenoData(from),
		  experimentData=experimentData(from),
		  annotation=annotation(from),
		  protocolData=protocolData(from),
		  featureData=featureData(from))
      })

setMethod("computeCopynumber", "CNSet", function(object, cnOptions){
	## to do the bias adjustment, initial estimates of the parameters are needed
	##  The initial estimates are gotten by running computeCopynumber with cnOptions[["bias.adj"]]=FALSE
	bias.adj <- cnOptions[["bias.adj"]]
	if(bias.adj & all(is.na(CA(object)))){
		cnOptions[["bias.adj"]] <- FALSE
	}
	object <- computeCopynumber.CNSet(object, cnOptions)				
	if(bias.adj & !cnOptions[["bias.adj"]]){
		## Do a second iteration with bias adjustment
		cnOptions[["bias.adj"]] <- TRUE
		object <- computeCopynumber.CNSet(object, cnOptions)
	}
	object
})

setMethod("computeCopynumber", "character", function(object, cnOptions){
	crlmmFile <- object
	isCNSet <- length(grep("cnSet", crlmmFile[1])) > 0
	for(i in seq(along=crlmmFile)){
		cat("Processing ", crlmmFile[i], "...\n")
		load(crlmmFile[i])
		if(isCNSet){
			object <- get("cnSet")
		} else {
			object <- get("callSetPlus")
		}
		CHR <- unique(chromosome(object))
		##if(length(CHR) > 1) stop("More than one chromosome in the object. This method requires one chromosome at a time.")		
		if(all(CHR==24)){
			message("skipping chromosome 24")
			next()
		}
		cat("----------------------------------------------------------------------------\n")
		cat("-        Estimating copy number for chromosome", CHR, "\n")
		cat("----------------------------------------------------------------------------\n")
		cnSet <- computeCopynumber(object, cnOptions)
		save(cnSet, file=file.path(dirname(crlmmFile), paste("cnSet_", CHR, ".rda", sep="")))
		if(!isCNSet) if(cnOptions[["unlink"]]) unlink(crlmmFile[i])
		rm(object, cnSet); gc();
	}	
})





##setMethod("computeHmm", "SnpSuperSet", function(object, hmmOptions){
##	cnSet <- computeCopynumber(object, hmmOptions)
##	computeHmm(cnSet, hmmOptions)
##})

## Genotype everything to get callSetPlus objects
## Go from callSets to Segments sets, writing only the segment set to file
## Safe, but very inefficient. Writes the quantile normalized data to file several times...
##setMethod("computeHmm", "character", function(object, hmmOptions){
##	outdir <- cnOptions[["outdir"]]
##	hmmOptions <- hmmOptions[["hmmOpts"]]
##	filenames <- object
##	for(i in seq(along=filenames)){
##		chrom <- gsub(".rda", "", strsplit(filenames[i], "_")[[1]][[2]])
##		if(hmmOptions[["verbose"]])
##			message("Fitting HMM to chromosome ", chrom)
##		if(file.exists(filenames[i])){
##			message("Loading ", filenames[i])
##			load(filenames[i])
##			cnSet <- get("cnSet")
##		} else {
##			stop("File ", filenames[i], " does not exist.")
##		}
##		hmmOptions$emission <- computeEmission(filenames[i], hmmOptions)
##		cnSet <- computeHmm(cnSet, hmmOptions)
##		##MIN.MARKERS <- hmmOptions[["MIN.MARKERS"]]
##		##segmentSet <- segments[segments$nprobes >= MIN.MARKERS, ]
##		message("Saving ", file.path(outdir, paste("cnSet_", chrom, ".rda", sep="")))
##		save(cnSet,
##		     file=file.path(outdir, paste("cnSet_", chrom, ".rda", sep="")))
##		unlink(file.path(outdir, paste("cnSet_", chrom, ".rda", sep="")))
##	}
##	fns <- list.files(outdir, pattern="cnSet", full.names=TRUE)
##	return(fns)	
##})

setMethod("copyNumber", "CNSet", function(object){
	I <- isSnp(object)
	CA <- CA(object)
	CB <- CB(object)
	CN <- CA + CB
	##For nonpolymorphic probes, CA is the total copy number
	CN[!I, ] <- CA(object)[!I, ]
	CN
})


setMethod("ellipse", "CNSet", function(x, copynumber, batch, ...){
	ellipse.CNSet(x, copynumber, batch, ...)
})

##setMethod("ellipse", "CNSet", function(x, copynumber, ...){
ellipse.CNSet <- function(x, copynumber, batch, ...){
	if(nrow(x) > 1) stop("only 1 snp at a time")
	##batch <- unique(x$batch)
	if(missing(batch)){
		stop("must specify batch")
	}
	if(length(batch) > 1) stop("batch variable not unique")
	nuA <- getParam(x, "nuA", batch)
	nuB <- getParam(x, "nuB", batch)
	phiA <- getParam(x, "phiA", batch)
	phiB <- getParam(x, "phiB", batch)
	tau2A <- getParam(x, "tau2A", batch)
	tau2B <- getParam(x, "tau2B", batch)
	sig2A <- getParam(x, "sig2A", batch)
	sig2B <- getParam(x, "sig2B", batch)
	corrA.BB <- getParam(x, "corrA.BB", batch)
	corrB.AA <- getParam(x, "corrB.AA", batch)
	corr <- getParam(x, "corr", batch)
	for(CN in copynumber){
		for(CA in 0:CN){
			CB <- CN-CA
			A.scale <- sqrt(tau2A*(CA==0) + sig2A*(CA > 0))
			B.scale <- sqrt(tau2B*(CB==0) + sig2B*(CB > 0))
			scale <- c(A.scale, B.scale)
			if(CA == 0 & CB > 0) rho <- corrA.BB
			if(CA > 0 & CB == 0) rho <- corrB.AA
			if(CA > 0 & CB > 0) rho <- corr
			if(CA == 0 & CB == 0) rho <- 0
			lines(ellipse(x=rho, centre=c(log2(nuA+CA*phiA),
					     log2(nuB+CB*phiB)),
				      scale=scale), ...)
		}
	}
}






setMethod("show", "CNSet", function(object){
	callNextMethod(object)
##	cat("emissionPr\n")
##	cat("   array:", nrow(object), "features,", ncol(object), "samples,", dim(emissionPr(object))[3], "states\n")
##	cat("   hidden states:\n")
##	cat("      ", dimnames(emissionPr(object))[[3]], "\n")
##	cat("   Missing values:", sum(is.na(emissionPr(object))), "\n")
##	if(!all(is.na(emissionPr(object)))){
##		cat("   minimum value:", min(emissionPr(object), na.rm=TRUE), "\n")
##	} else  cat("   minimum value: NA (all missing)\n")
##	cat("rangedData:  ")
##	cat("    ", show(rangedData(object)), "\n")
})
##setMethod("rangedData", "CNSet", function(object) segmentData(object))
##setReplaceMethod("rangedData", c("CNSet", "RangedData"), function(object, value){
##	segmentData(object) <- value
##})
##
##setMethod("segmentData", "CNSet", function(object) object@segmentData)
##setReplaceMethod("segmentData", c("CNSet", "RangedData"), function(object, value){
##	object@segmentData <- value
##	object
##})
##
##setMethod("emissionPr", "CNSet", function(object) object@emissionPr)
##setReplaceMethod("emissionPr", c("CNSet", "array"), function(object, value){
##	object@emissionPr <- value
##	object
##})
##setMethod("start", "CNSet", function(x, ...) start(segmentData(x), ...))
##setMethod("end", "CNSet", function(x, ...) end(segmentData(x), ...))
##setMethod("width", "CNSet", function(x) width(segmentData(x)))
