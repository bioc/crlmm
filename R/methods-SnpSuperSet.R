##How to make the initialization platform-specific?

setMethod("initialize", "SnpSuperSet",
          function(.Object,
		   assayData,
                   call=new("matrix"),
                   callProbability=new("matrix"),
                   alleleA=new("matrix"),
                   alleleB=new("matrix"),
		   featureData,
		   annotation,
		   ...){
		  if(!missing(assayData)){
			  .Object <- callNextMethod(.Object, assayData=assayData,...)
		  } else{
			  ad <- assayDataNew("lockedEnvironment",
					     call=call,
					     callProbability=callProbability,
					     alleleA=alleleA,
					     alleleB=alleleB)
			  .Object <- callNextMethod(.Object,
						    assayData=ad, ...)
		  }		  
		  if(missing(annotation)){
			  stop("must specify annotation")
		  } else{
			  stopifnot(isValidCdfName(annotation))
			  .Object@annotation <- annotation
		  }		  
		  if (missing(featureData)){
			  featureData(.Object) <- annotatedDataFrameFrom(call, byrow=TRUE)
		  } else{
			  featureData(.Object) <- featureData
		  }
		  ## Do after annotation has been assigned
		  if(!(all(c("chromosome", "position", "isSnp")  %in% colnames(.Object@featureData)))){
			  ##update the featureData
			  .Object@featureData <- addFeatureAnnotation.SnpSuperSet(.Object)
		  }
		  .Object
          })

setMethod("addFeatureAnnotation", "SnpSuperSet", function(object, ...){
	addFeatureAnnotation.SnpSuperSet(object, ...)
})

getParam.SnpSuperSet <- function(object, name, batch){
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



setMethod("splitByChromosome", "SnpSuperSet", function(object, cnOptions){
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
			cnSet <- computeCopynumber(callSetPlus, cnOptions)
			
		} else{
			message("Copy number estimates not available for chromosome Y.  Saving only the 'callSetPlus' object for this chromosome")
			save(callSetPlus, file=file.path(outdir, paste("callSetPlus_", CHR, ".rda", sep="")))
		}
		if(cnOptions[["hiddenMarkovModel"]] & CHR != 24){
			cnSet <- computeHmm(cnSet, cnOptions)
		}
		save(cnSet, file=file.path(outdir, paste("cnSet_", CHR, ".rda", sep="")))
		saved.objects <- list.files(outdir, pattern="cnSet", full.names=TRUE)
##		} else{ ## save crlmmSet to outdir
##			save(cnSet, file=file.path(outdir, paste("cnSet_", CHR, ".rda", sep="")))
##			saved.objects <- list.files(outdir, pattern="cnSet", full.names=TRUE)			
##		}		
	}
	saved.objects
})

setMethod("computeCopynumber", "SnpSuperSet",
	  function(object, cnOptions){
		  computeCopynumber.SnpSuperSet(object, cnOptions)
	  })

##gtConfidence <- function(object) 1-exp(-confs(object)/1000)
	
