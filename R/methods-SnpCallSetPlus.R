##How to make the initialization platform-specific?
setMethod("initialize", "SnpCallSetPlus",
          function(.Object,
                   phenoData, featureData,
                   calls=new("matrix"),
                   callsConfidence=new("matrix"),
                   senseThetaA=new("matrix"),
                   senseThetaB=new("matrix"),
		   annotation,
		   experimentData,
		   protocolData, ... ){
		  ad <- assayDataNew("lockedEnvironment",
				     calls=calls,
				     callsConfidence=callsConfidence,
				     senseThetaA=senseThetaA,
				     senseThetaB=senseThetaB)
		  assayData(.Object) <- ad
		  if (missing(phenoData)){
			  phenoData(.Object) <- annotatedDataFrameFrom(calls, byrow=FALSE)
		  } else{
			  phenoData(.Object) <- phenoData
		  }
		  if (missing(featureData)){
			  featureData(.Object) <- annotatedDataFrameFrom(calls, byrow=TRUE)
		  } else{
			  featureData(.Object) <- featureData
		  }
		  if(!missing(annotation)) annotation(.Object) <- annotation
		  if(!missing(experimentData)) experimentData(.Object) <- experimentData
		  if(!missing(protocolData)) protocolData(.Object) <- protocolData
		  .Object
          })
getParam.SnpCallSetPlus <- function(object, name, batch){
		  label <- paste(name, batch, sep="_")
		  colindex <- grep(label, fvarLabels(object))
		  if(length(colindex) == 1){
			  param <- fData(object)[, colindex]
		  }
		  if(length(colindex) < 1){
			  param <- NULL
		  }
		  if(is.na(colindex)){
			  stop(paste(label, " not found in object"))
		  }
		  if(length(colindex) > 1){
			  stop(paste(label, " not unique"))
		  }
		  return(param)
	  }

setMethod("getParam", signature(object="SnpCallSetPlus",
				name="character",
				batch="ANY"),
	  function(object, name, batch){
		  if(length(batch) > 1){
			  warning("batch argument to getParam should have length 1.  Only using the first")
			  batch <- batch[1]
		  }
		  getParam.SnpCallSetPlus(object, name, batch)
})

setMethod("splitByChromosome", "SnpCallSetPlus", function(object, cnOptions){
	tmpdir <- cnOptions[["tmpdir"]]
	outdir <- cnOptions[["outdir"]]	
	save.it <- cnOptions[["save.it"]]
	path <- system.file("extdata", package=paste(annotation(object), "Crlmm", sep=""))	
	load(file.path(path, "snpProbes.rda"))
	snpProbes <- get("snpProbes")
	load(file.path(path, "cnProbes.rda"))
	cnProbes <- get("cnProbes")	
	k <- grep("chr", colnames(snpProbes))
	if(length(k) < 1) stop("chr or chromosome not in colnames(snpProbes)")
	for(CHR in 1:24){
		cat("Chromosome ", CHR, "\n")
		snps <- rownames(snpProbes)[snpProbes[, k] == CHR]
		cnps <- rownames(cnProbes)[cnProbes[, k] == CHR]
		index <- c(match(snps, featureNames(object)),
			   match(cnps, featureNames(object)))
		index <- index[!is.na(index)]
		callSetPlus <- object[index, ]
		if(CHR != 24){
			crlmmSet <- computeCopynumber(callSetPlus, cnOptions)
			
		} else{
			message("Copy number estimates not available for chromosome Y.  Saving only the 'callSetPlus' object for this chromosome")
			save(callSetPlus, file=file.path(outdir, paste("callSetPlus_", CHR, ".rda", sep="")))
		}
		if(cnOptions[["hiddenMarkovModel"]] & CHR != 24){
			segmentSet <- computeHmm(crlmmSet, cnOptions)
			save(segmentSet, file=file.path(outdir, paste("segmentSet_", CHR, ".rda", sep="")))
			saved.objects <- list.files(outdir, pattern="segmentSet", full.names=TRUE)
		} else{ ## save crlmmSet to outdir
			save(crlmmSet, file=file.path(outdir, paste("crlmmSet_", CHR, ".rda", sep="")))
			saved.objects <- list.files(outdir, pattern="crlmmSet", full.names=TRUE)			
		}		
	}
	saved.objects
})
setMethod("addFeatureAnnotation", "SnpCallSetPlus", function(object, ...){
	addFeatureAnnotation.SnpCallSetPlus(object, ...)
})
setMethod("computeCopynumber", "SnpCallSetPlus",
	  function(object, cnOptions){
		  computeCopynumber.SnpCallSetPlus(object, cnOptions)
	  })
setMethod("chromosome", "SnpCallSetPlus", function(object) fData(object)$chromosome)
setMethod("position", "SnpCallSetPlus", function(object) fData(object)$position)
gtConfidence <- function(object) 1-exp(-callsConfidence(object)/1000)
	
