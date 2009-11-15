setMethod("initialize", "CrlmmSet",
          function(.Object,
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
		   CB=new("matrix"), ... ){
		  ad <- assayDataNew("lockedEnvironment",
				     calls=calls,
				     callsConfidence=callsConfidence,
				     senseThetaA=senseThetaA,
				     senseThetaB=senseThetaB,
				     CA=CA,
				     CB=CB)
		  assayData(.Object) <- ad
		  if (missing(phenoData)){
			  phenoData(.Object) <- annotatedDataFrameFrom(calls, byrow=FALSE)
		  } else phenoData(.Object) <- phenoData
		  if (missing(featureData)){
			  featureData(.Object) <- annotatedDataFrameFrom(calls, byrow=TRUE)
		  } else featureData(.Object) <- featureData
		  if(!missing(annotation)) annotation(.Object) <- annotation
		  if(!missing(experimentData)) experimentData(.Object) <- experimentData
		  if(!missing(protocolData)) protocolData(.Object) <- protocolData
		  .Object	    
          })

setAs("SnpCallSetPlus", "CrlmmSet",
      function(from, to){
	      CA <- CB <- matrix(NA, nrow(from), ncol(from))
	      dimnames(CA) <- dimnames(CB) <- list(featureNames(from), sampleNames(from))		  
	      new("CrlmmSet",
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

setValidity("CrlmmSet", function(object) {
	assayDataValidMembers(assayData(object), c("CA", "CB", "calls", "callsConfidence", "senseThetaA", "senseThetaB"))
})



setMethod("computeCopynumber", "CrlmmSet", function(object, cnOptions){
	computeCopynumber.CrlmmSet(object, cnOptions)
})

setMethod("computeCopynumber", "character", function(object, cnOptions){
	crlmmFile <- object
	isCrlmmSet <- length(grep("crlmmSet", crlmmFile[1])) > 0
	for(i in seq(along=crlmmFile)){
		cat("Processing ", crlmmFile[i], "...\n")
		load(crlmmFile[i])
		if(isCrlmmSet){
			object <- get("crlmmSet")
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
		crlmmSet <- computeCopynumber(object, cnOptions)
		save(crlmmSet, file=file.path(dirname(crlmmFile), paste("crlmmSet_", CHR, ".rda", sep="")))
		if(!isCrlmmSet) if(cnOptions[["unlink"]]) unlink(crlmmFile[i])
		rm(object, crlmmSet); gc();
	}	
})

setMethod("pr", signature(object="CrlmmSet",
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

setMethod("computeEmission", "CrlmmSet", function(object, hmmOptions){
	computeEmission.CrlmmSet(object, hmmOptions)
})

setMethod("computeEmission", "character", function(object, hmmOptions){
	filename <- object
	chrom <- gsub(".rda", "", strsplit(filename, "_")[[1]][[2]])
	if(hmmOptions[["verbose"]])
		message("Compute emission probabilities for chromosome ", chrom)
	if(file.exists(filename)){
		load(filename)
		crlmmSet <- get("crlmmSet")
	} else {
		stop("File ", filename, " does not exist.")
	}
	emission <- computeEmission(crlmmSet, hmmOptions)
	message("Saving ", file.path(dirname(filename), paste("emission_", chrom, ".rda", sep="")))
	if(hmmOptions[["save.it"]]){
		save(emission,
		     file=file.path(dirname(filename), paste("emission_", chrom, ".rda", sep="")))
	}
	return(emission)
})

setMethod("computeHmm", "CrlmmSet", function(object, hmmOptions){
	computeHmm.CrlmmSet(object, hmmOptions)
})

##setMethod("computeHmm", "SnpCallSetPlus", function(object, hmmOptions){
##	crlmmSet <- computeCopynumber(object, hmmOptions)
##	computeHmm(crlmmSet, hmmOptions)
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
			crlmmSet <- get("crlmmSet")
		} else {
			stop("File ", filenames[i], " does not exist.")
		}
		hmmOptions$emission <- computeEmission(filenames[i], hmmOptions)
		segmentSet <- computeHmm(crlmmSet, hmmOptions)
		##MIN.MARKERS <- hmmOptions[["MIN.MARKERS"]]
		##segmentSet <- segments[segments$nprobes >= MIN.MARKERS, ]
		message("Saving ", file.path(outdir, paste("segmentSet_", chrom, ".rda", sep="")))
		save(segmentSet,
		     file=file.path(outdir, paste("segmentSet_", chrom, ".rda", sep="")))
		unlink(file.path(outdir, paste("crlmmSet_", chrom, ".rda", sep="")))
	}
	fns <- list.files(outdir, pattern="segmentSet", full.names=TRUE)
	return(fns)	
})

##initializeCrlmmSet <- function(annotation, sns, backingfile="./"){
##	require(bigmemory)
##	message("Initializing CrlmmSet file--be patient...")
##	path <- system.file("extdata", package=paste(annotation, "Crlmm", sep=""))
##	load(file.path(path, "snpProbes.rda"))
##	load(file.path(path, "cnProbes.rda"))
##	nr <- nrow(snpProbes)+nrow(cnProbes)
##	pos <- c(snpProbes[, "position"], cnProbes[, "position"])
##	chr <- c(snpProbes[, "chrom"], cnProbes[, "chrom"])
##	index <- order(chr, pos)
##	fns <- c(rownames(snpProbes), rownames(cnProbes))[index]
####	fns <- fns[1:20]
##	nr <- length(fns)
##	nc <- length(sns)
####	if(!file.exists(file.path(backingfile, "A.bin"))){
##		AA <- filebacked.big.matrix(nr, nc, 
##					    type="integer",
##					    init=0,
##					    backingpath=backingfile,
##					    backingfile="A.bin",
##					    descriptorfile="A.desc",
##					    dimnames=list(fns, sns))
####	} else {
####		AA <- attach.big.matrix("A.desc", backingpath=backingfile)
####	}
####	if(!file.exists(file.path(backingfile, "B.bin"))){
##		BB <- filebacked.big.matrix(nr, nc, 
##					    type="integer",
##					    init=0,
##					    backingpath=backingfile,
##					    backingfile="B.bin",
##					    descriptorfile="B.desc",
##					    dimnames=list(fns, sns))			   
####	}  else {
####		BB <- attach.big.matrix("B.desc", backingpath=backingfile)
####	}
####	if(!file.exists(file.path(backingfile, "GT.bin"))){
##		gt <- filebacked.big.matrix(nr, nc, 
##					    type="integer",
##					    init=0,
##					    backingpath=backingfile,
##					    backingfile="GT.bin",
##					    descriptorfile="GT.desc",
##					    dimnames=list(fns, sns))
##		confs <- filebacked.big.matrix(nr, nc, 
##					       type="integer",
##					       init=0,
##					       backingpath=backingfile,
##					       backingfile="confs.bin",
##					       descriptorfile="confs.desc",
##					       dimnames=list(fns, sns))	
####	} else {
####		gt <- attach.big.matrix("GT.desc", backingpath=backingfile)
####	}
####	if(!file.exists(file.path(backingfile, "CA.bin"))){
##		ca <- filebacked.big.matrix(nr, nc, 
##					    type="integer",
##					    init=0,
##					    backingpath=backingfile,
##					    backingfile="CA.bin",
##					    descriptorfile="CA.desc",
##					    dimnames=list(fns, sns))
####	}  else {
####		ca <- attach.big.matrix("CA.desc", backingpath=backingfile)
####	}
####	if(!file.exists(file.path(backingfile, "CB.bin"))){
##		cb <- filebacked.big.matrix(nr, nc, 
##					    type="integer",
##					    init=0,
##					    backingpath=backingfile,
##					    backingfile="CB.bin",
##					    descriptorfile="CB.desc",
##					    dimnames=list(fns, sns))
####	} else {
####		cb <- attach.big.matrix("CB.desc", backingpath=backingfile)
####	}
##	return(new("CrlmmSet", A=AA, B=BB, CA=ca, CB=cb, GT=gt, confs=confs))
##}


##setMethod("A", signature(object="CrlmmSet"),
##          function(object) assayDataElement(object,"A"))
##setMethod("B", signature(object="CrlmmSet"),
##          function(object) assayDataElement(object,"A"))
##setMethod("CA", signature(object="CrlmmSet"),
##          function(object) assayDataElement(object,"A"))
##setMethod("CB", signature(object="CrlmmSet"),
##          function(object) assayDataElement(object,"A"))
##setMethod("GT", signature(object="CrlmmSet"),
##          function(object) assayDataElement(object,"GT"))

##setReplaceMethod("CA", signature(object="CrlmmSetList", value="big.matrix"),
##		 function(object, value){
##			 assayDataElementReplace(object, "A", value)
##		 })
##setReplaceMethod("CB", signature(object="CrlmmSetList", value="big.matrix"),
##		 function(object, value){
##			 CB(object[[3]]) <- value
##			 object
##			 })
##
##setReplaceMethod("A", signature(object="CrlmmSet", value="big.matrix"),
##		 function(object, value){
##			 assayDataElementReplace(object, "A", value)
##		 })

##setReplaceMethod("GT", signature(object="CrlmmSet", value="big.matrix"),
##		 function(object, value){
##			 assayDataElementReplace(object, "GT", value)
##		 })
##setReplaceMethod("confs", signature(object="CrlmmSet", value="big.matrix"),
##		 function(object, value){
##			 assayDataElementReplace(object, "callProbability", value)
##		 })
##setReplaceMethod("confs", signature(object="CrlmmSet", value="matrix"),
##		 function(object, value){
##			 assayDataElementReplace(object, "callProbability", value)
##		 })
##setReplaceMethod("GT", signature(object="CrlmmSet", value="matrix"),
##		 function(object, value){
##			 assayDataElementReplace(object, "GT", value)
##		 })
##setReplaceMethod("B", signature(object="CrlmmSetList", value="big.matrix"),
##		 function(object, value){
##			 B(object[[1]]) <- value
##			 object
##		 })
