setMethod("initialize", "CNSet",
          function(.Object,
		   assayData,
                   phenoData,
		   featureData,
		   annotation,
		   experimentData,
		   protocolData,
                   calls=new("matrix"),
                   callsConfidence=new("matrix"),
                   senseThetaA=new("matrix"),
                   senseThetaB=new("matrix"),
		   CA=new("matrix"),
		   CB=new("matrix"),
		   segmentData=new("RangedData"),
		   emissionPr=new("array"), ... ){
		  if(missing(assayData)){
			  assayData <- assayDataNew("lockedEnvironment",
						    calls=calls,
						    callsConfidence=callsConfidence,
						    senseThetaA=senseThetaA,
						    senseThetaB=senseThetaB,
						    CA=CA,
						    CB=CB)
		  } 
		  assayData(.Object) <- assayData
		  if (missing(phenoData)){
			  phenoData(.Object) <- annotatedDataFrameFrom(calls, byrow=FALSE)
		  } else phenoData(.Object) <- phenoData
		  if (missing(featureData)){
			  featureData(.Object) <- annotatedDataFrameFrom(calls, byrow=TRUE)
		  } else featureData(.Object) <- featureData
		  if(!missing(annotation)) annotation(.Object) <- annotation
		  if(!missing(experimentData)) experimentData(.Object) <- experimentData
		  if(!missing(protocolData)) protocolData(.Object) <- protocolData
		  if(!missing(emissionPr)) .Object@emissionPr <- emissionPr
		  segmentData(.Object) <- segmentData
		  .Object	    
          })

setAs("SnpCallSetPlus", "CNSet",
      function(from, to){
	      CA <- CB <- matrix(NA, nrow(from), ncol(from))
	      dimnames(CA) <- dimnames(CB) <- list(featureNames(from), sampleNames(from))		  
	      new("CNSet",
		  calls=calls(from),
		  callsConfidence=callsConfidence(from),
		  senseThetaA=A(from),
		  senseThetaB=B(from),
		  CA=CA,
		  CB=CB,
		  phenoData=phenoData(from),
		  experimentData=experimentData(from),
		  annotation=annotation(from),
		  protocolData=protocolData(from),
		  featureData=featureData(from))
      })

setValidity("CNSet", function(object) {
	assayDataValidMembers(assayData(object), c("CA", "CB", "call", "callProbability", "senseThetaA", "senseThetaB"))
})



setMethod("computeCopynumber", "CNSet", function(object, cnOptions){
	computeCopynumber.CNSet(object, cnOptions)
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

setMethod("pr", signature(object="CNSet",
			  name="character",
			  batch="ANY",
			  value="numeric"), 
	  function(object, name, batch, value){
		  label <- paste(name, batch, sep="_")
		  colindex <- match(label, fvarLabels(object))
		  if(length(colindex) == 1){
			  fData(object)[, colindex] <- value
		  } 
		  if(is.na(colindex)){
			  stop(paste(label, " not found in object"))
		  }
		  if(length(colindex) > 1){
			  stop(paste(label, " not unique"))
		  }
		  object
	  })

setMethod("computeEmission", "CNSet", function(object, hmmOptions){
	computeEmission.CNSet(object, hmmOptions)
})

setMethod("computeEmission", "character", function(object, hmmOptions){
	filename <- object
	chrom <- gsub(".rda", "", strsplit(filename, "_")[[1]][[2]])
	if(hmmOptions[["verbose"]])
		message("Compute emission probabilities for chromosome ", chrom)
	if(file.exists(filename)){
		load(filename)
		cnSet <- get("cnSet")
	} else {
		stop("File ", filename, " does not exist.")
	}
	emission <- computeEmission(cnSet, hmmOptions)
	message("Saving ", file.path(dirname(filename), paste("emission_", chrom, ".rda", sep="")))
	if(hmmOptions[["save.it"]]){
		save(emission,
		     file=file.path(dirname(filename), paste("emission_", chrom, ".rda", sep="")))
	}
	return(emission)
})

setMethod("computeHmm", "CNSet", function(object, hmmOptions){
	computeHmm.CNSet(object, hmmOptions)
})

##setMethod("computeHmm", "SnpCallSetPlus", function(object, hmmOptions){
##	cnSet <- computeCopynumber(object, hmmOptions)
##	computeHmm(cnSet, hmmOptions)
##})

## Genotype everything to get callSetPlus objects
## Go from callSets to Segments sets, writing only the segment set to file
## Safe, but very inefficient. Writes the quantile normalized data to file several times...
setMethod("computeHmm", "character", function(object, hmmOptions){
	outdir <- cnOptions[["outdir"]]
	hmmOptions <- hmmOptions[["hmmOpts"]]
	filenames <- object
	for(i in seq(along=filenames)){
		chrom <- gsub(".rda", "", strsplit(filenames[i], "_")[[1]][[2]])
		if(hmmOptions[["verbose"]])
			message("Fitting HMM to chromosome ", chrom)
		if(file.exists(filenames[i])){
			message("Loading ", filenames[i])
			load(filenames[i])
			cnSet <- get("cnSet")
		} else {
			stop("File ", filenames[i], " does not exist.")
		}
		hmmOptions$emission <- computeEmission(filenames[i], hmmOptions)
		cnSet <- computeHmm(cnSet, hmmOptions)
		##MIN.MARKERS <- hmmOptions[["MIN.MARKERS"]]
		##segmentSet <- segments[segments$nprobes >= MIN.MARKERS, ]
		message("Saving ", file.path(outdir, paste("cnSet_", chrom, ".rda", sep="")))
		save(cnSet,
		     file=file.path(outdir, paste("cnSet_", chrom, ".rda", sep="")))
		unlink(file.path(outdir, paste("cnSet_", chrom, ".rda", sep="")))
	}
	fns <- list.files(outdir, pattern="cnSet", full.names=TRUE)
	return(fns)	
})



setValidity("CNSet", function(object) {
	msg <- validMsg(NULL, assayDataValidMembers(assayData(object), c("CA", "CB")))
	if (is.null(msg)) TRUE else msg
})
##may want to allow thresholding here (... arg)
setMethod("CA", "CNSet", function(object) assayData(object)[["CA"]]/100)
setMethod("CB", "CNSet", function(object) assayData(object)[["CB"]]/100)

setReplaceMethod("CA", signature(object="CNSet", value="matrix"),
		 function(object, value){
			 assayDataElementReplace(object, "CA", value)		
		 })

setReplaceMethod("CB", signature(object="CNSet", value="matrix"),
		 function(object, value){
			 assayDataElementReplace(object, "CB", value)
		 })

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
##	nuA <- as.numeric(fData(x)[, match(paste("nuA", batch, sep="_"), fvarLabels(x))])
##	nuB <- as.numeric(fData(x)[, match(paste("nuB", batch, sep="_"), fvarLabels(x))])	
##	phiA <- as.numeric(fData(x)[, match(paste("phiA", batch, sep="_"), fvarLabels(x))])
##	phiB <- as.numeric(fData(x)[, match(paste("phiB", batch, sep="_"), fvarLabels(x))])	
##	tau2A <- as.numeric(fData(x)[, match(paste("tau2A", batch, sep="_"), fvarLabels(x))])
##	tau2B <- as.numeric(fData(x)[, match(paste("tau2B", batch, sep="_"), fvarLabels(x))])
##	sig2A <- as.numeric(fData(x)[, match(paste("sig2A", batch, sep="_"), fvarLabels(x))])
##	sig2B <- as.numeric(fData(x)[, match(paste("sig2B", batch, sep="_"), fvarLabels(x))])
##	corrA.BB <- as.numeric(fData(x)[, match(paste("corrA.BB", batch, sep="_"), fvarLabels(x))])
##	corrB.AA <- as.numeric(fData(x)[, match(paste("corrB.AA", batch, sep="_"), fvarLabels(x))])
##	corr <- as.numeric(fData(x)[, match(paste("corr", batch, sep="_"), fvarLabels(x))])
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

setMethod("segmentData", "CNSet", function(object) object@segmentData)
setReplaceMethod("segmentData", c("CNSet", "RangedData"), function(object, value){
	object@segmentData <- value
	object
})

setMethod("emissionPr", "CNSet", function(object) object@emissionPr)
setReplaceMethod("emissionPr", c("CNSet", "array"), function(object, value){
	object@emissionPr <- value
	object
})

setMethod("show", "CNSet", function(object){
	callNextMethod(object)
	cat("emissionPr\n")
	cat("   array:", nrow(object), "features,", ncol(object), "samples,", dim(emissionPr(object))[3], "states\n")
	cat("   hidden states:\n")
	cat("      ", dimnames(emissionPr(object))[[3]], "\n")
	cat("   Missing values:", sum(is.na(emissionPr(object))), "\n")
	if(!all(is.na(emissionPr(object)))){
		cat("   minimum value:", min(emissionPr(object), na.rm=TRUE), "\n")
	} else  cat("   minimum value: NA (all missing)\n")
	cat("segmentData:  ")
	cat("    ", show(segmentData(object)), "\n")
##	cat("   ", nrow(segmentData(object)), "segments\n")
##	cat("    column names:", paste(colnames(segmentData(object)), collapse=", "), "\n")
#	cat("    mean # markers per segment:", mean(segmentData(object)$nprobes))
})







