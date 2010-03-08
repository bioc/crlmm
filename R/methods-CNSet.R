setMethod("initialize", "CNSet", function(.Object, options=new("CrlmmOptions"), ...){
	.Object <- callNextMethod(.Object, ...)
	.Object@options <- options
	storageMode(.Object) <- "environment"
	return(.Object)
})

setMethod("show", "CNSet", function(object){
	callNextMethod(object)
	cat("lM: ", length(lM(object)), " elements \n")
	print(names(lM(object)))
})

setMethod("[", "CNSet", function(x, i, j, ..., drop=FALSE){
	x <- callNextMethod(x, i, j, ..., drop=drop)
	if(!missing(i)){
		if(class(lM(x)) == "ffdf"){
			lM(x) <- lapply(physical(lM(x)), function(x, i){open(x); x[i, ]}, i=i)
		} else {
			lM(x) <- lapply(lM(x), function(x, i) x[i, , drop=FALSE], i=i)
		}
	}
	x
})
.getCA <- function(object){
	value <- assayDataElement(object, "CA")
	if(class(value) == "ff_matrix"){
		return(value)
	} else {
		return(value/100)
	}
}
.getCB <- function(object){
	value <- assayDataElement(object, "CB")
	if(class(value) == "ff_matrix"){
		return(value)
	} else {
		return(value/100)
	}
}
.replaceCA <- function(object, value){
	value <- matrix(as.integer(value*100), nrow(value), ncol(value), dimnames=dimnames(value))
	assayDataElementReplace(object, "CA", value)
}
.replaceCB <- function(object, value){
	value <- matrix(as.integer(value*100), nrow(value), ncol(value), dimnames=dimnames(value))			 
	assayDataElementReplace(object, "CB", value)
}
##setMethod("CA", "AffymetrixCNSet", .getCA)
##setMethod("CA", "IlluminaCNSet", .getCA)
##setMethod("CB", "AffymetrixCNSet", .getCB)
##setMethod("CB", "IlluminaCNSet", .getCB)
setMethod("CA", "CNSet", .getCA)
setMethod("CB", "CNSet", .getCB)
##setReplaceMethod("CA", signature(object="AffymetrixCNSet", value="matrix"),
##		 function(object, value) .replaceCA(object, value))
##setReplaceMethod("CA", signature(object="IlluminaCNSet", value="matrix"),
##		 function(object, value) .replaceCA(object,value))
setReplaceMethod("CA", signature(object="CNSet", value="matrix"),
		 function(object, value) .replaceCB(object,value))
setReplaceMethod("CB", signature(object="CNSet", value="matrix"),
		 function(object, value) .replaceCB(object,value))

##setAs("SnpSuperSet", "CNSet",
##      function(from, to){
##	      CA <- CB <- matrix(NA, nrow(from), ncol(from))
##	      dimnames(CA) <- dimnames(CB) <- list(featureNames(from), sampleNames(from))		  
##	      new("CNSet",
##		  call=calls(from),
##		  callProbability=assayData(from)[["callProbability"]],  ##confs(from) returns 1-exp(-x/1000)
##		  alleleA=A(from),
##		  alleleB=B(from),
##		  CA=CA,
##		  CB=CB,
##		  phenoData=phenoData(from),
##		  experimentData=experimentData(from),
##		  annotation=annotation(from),
##		  protocolData=protocolData(from),
##		  featureData=featureData(from))
##      })



##setMethod("computeCopynumber", "character", function(object, cnOptions){
##	crlmmFile <- object
##	isCNSet <- length(grep("cnSet", crlmmFile[1])) > 0
##	for(i in seq(along=crlmmFile)){
##		cat("Processing ", crlmmFile[i], "...\n")
##		load(crlmmFile[i])
##		if(isCNSet){
##			object <- get("cnSet")
##		} else {
##			object <- get("callSetPlus")
##		}
##		CHR <- unique(chromosome(object))
##		##if(length(CHR) > 1) stop("More than one chromosome in the object. This method requires one chromosome at a time.")		
##		if(all(CHR==24)){
##			message("skipping chromosome 24")
##			next()
##		}
##		cat("----------------------------------------------------------------------------\n")
##		cat("-        Estimating copy number for chromosome", CHR, "\n")
##		cat("----------------------------------------------------------------------------\n")
##		cnSet <- computeCopynumber(object, cnOptions)
##		save(cnSet, file=file.path(dirname(crlmmFile), paste("cnSet_", CHR, ".rda", sep="")))
##		if(!isCNSet) if(cnOptions[["unlink"]]) unlink(crlmmFile[i])
##		rm(object, cnSet); gc();
##	}	
##})
setMethod("copyNumber", "CNSet", function(object) .copyNumber(object))
setMethod("ellipse", "CNSet", function(x, copynumber, batch, ...){
	ellipse.CNSet(x, copynumber, batch, ...)
})
##setMethod("[", "AffymetrixCNSet", function(x, i, j, ..., drop=FALSE){
##	x <- callNextMethod(x, i, j, ..., drop=drop)
##	if(!missing(i)){
##		if(class(linearModelParam(x)) == "ffdf"){
##			linearModelParam(x) <- lapply(physical(linearModelParam(x)), function(x, i){open(x); x[i, ]}, i=i)
##		} else {
##			linearModelParam(x) <- lapply(linearModelParam(x), function(x, i) x[i, , drop=FALSE], i=i)
##		}
##	}
##	x
##})
setMethod("close", "CNSet", function(con, ...){
	callNextMethod(con,...)
	close(linearModelParam(con))
})
setMethod("open", "CNSet", function(con, ...){
	callNextMethod(con,...)
	open(linearModelParam(con))
})
setMethod("lM", "CNSet", function(object) object@lM)
##setMethod("linearModelParam", "AffymetrixCNSet", function(object) object@linearModelParam)
setReplaceMethod("lM", c("CNSet", "list_or_ffdf"), function(object, value){
	object@lM <- value
	object
})
##setReplaceMethod("linearModelParam", c("IlluminaCNSet", "list_or_ffdf"), function(object, value){
##	object@linearModelParam <- value
##	object
##})
##getLinearModelParamFromOptions <- function(object, ...){
##	batchnames <- unique(object$batch)
##	getLinearModelParam(crlmmOptions(object), list(featureNames(object), batchnames))
##}


##setMethod("getLinearModelParam", "IlluminaCNSet", getLinearModelParamFromOptions)

##setMethod("computeCopynumber", "AffymetrixCNSet", function(object){
setMethod("computeCopynumber", "CNSet", function(object){
	ops <- crlmmOptions(object)$cnOpts
	bias.adj <- ops$bias.adj
	## to do the bias adjustment, initial estimates of the parameters are needed
	##  The initial estimates are gotten by running computeCopynumber with cnOptions[["bias.adj"]]=FALSE
	##bias.adj <- cnOptions[["bias.adj"]]
	if(bias.adj & all(is.na(CA(object[, 1])))){
		ops$bias.adj <- FALSE
	}
	crlmmOptions(object)$cnOpts <- ops
	object <- computeCopynumber.CNSet(object)
	object
})

computeCopynumber.CNSet <- function(object){
	batch <- object$batch
	if(is.null(batch) & ncol(object) < 150){
		message("You did not specify the batch variable. Assuming there is no batch effect.")
		batch <- object$batch <- rep("samebatch", ncol(object))
	}
	if(is.null(batch) & ncol(object) >= 150){
		stop("You did not specify the batch variable. ")
	}
	PLATE <- unique(batch)
	if(length(batch) != ncol(object)) stop("Batch must be the same length as the number of samples")
	ops <- crlmmOptions(object)
	cnOptions <- ops$cnOpts
	MIN.SAMPLES <- cnOptions$MIN.SAMPLES
	verbose <- ops$verbose
	use.poe <- cnOptions$use.poe
	bias.adj <- cnOptions$bias.adj
	
	tmp.objects <- instantiateObjects(object)
	if(bias.adj & ncol(object) <= 15){
		warning(paste("bias.adj is TRUE, but too few samples to perform this step"))
		cnOptions$bias.adj <- bias.adj <- FALSE
	}
	if(bias.adj){
		if(verbose) message("Dropping samples with low posterior prob. of normal copy number (samples dropped is locus-specific)")
		if(!use.poe)
			tmp.objects <- biasAdjNP(object, cnOptions, tmp.objects)
		tmp.objects <- biasAdj(object, cnOptions, tmp.objects)
		if(verbose) message("Recomputing location and scale parameters")
	}
	##update tmp.objects
	tmp.objects <- withinGenotypeMoments(object,
					     cnOptions=cnOptions,
					     tmp.objects=tmp.objects)
	object <- locationAndScale(object, cnOptions, tmp.objects)
	tmp.objects <- oneBatch(object,
				cnOptions=cnOptions,
				tmp.objects=tmp.objects)
	##coefs calls nuphiAllele.
	object <- coefs(object, cnOptions, tmp.objects)
	##nuA=getParam(object, "nuA", PLATE)
	THR.NU.PHI <- cnOptions$THR.NU.PHI
	if(THR.NU.PHI){
		if(verbose) message("Thresholding nu and phi")
		object <- thresholdModelParams(object, cnOptions)
	}		
	if(verbose) message("\nAllele specific copy number")	
	object <- polymorphic(object, cnOptions, tmp.objects)
	if(any(!isSnp(object))){ ## there are nonpolymorphic probes
		if(verbose) message("\nCopy number for nonpolymorphic probes...")
		if(!use.poe){
			object <- nonpolymorphic(object, cnOptions, tmp.objects)
		} else {
			object <- nonpolymorphic.poe(object, cnOptions, tmp.objects)
		}
	}
	##---------------------------------------------------------------------------
	##Note: the replacement method multiples by 100
	##CA(object)[, batch==PLATE] <- CA(object)
	##CB(object)[, batch==PLATE] <- CB(object)
	##---------------------------------------------------------------------------
	##update-the plate-specific parameters for copy number
	object <- pr(object, "nuA", PLATE, getParam(object, "nuA", PLATE))
	object <- pr(object, "nuA.se", PLATE, getParam(object, "nuA.se", PLATE))
	object <- pr(object, "nuB", PLATE, getParam(object, "nuB", PLATE))
	object <- pr(object, "nuB.se", PLATE, getParam(object, "nuB.se", PLATE))
	object <- pr(object, "phiA", PLATE, getParam(object, "phiA", PLATE))
	object <- pr(object, "phiA.se", PLATE, getParam(object, "phiA.se", PLATE))
	object <- pr(object, "phiB", PLATE, getParam(object, "phiB", PLATE))
	object <- pr(object, "phiB.se", PLATE, getParam(object, "phiB.se", PLATE))
	object <- pr(object, "tau2A", PLATE, getParam(object, "tau2A", PLATE))
	object <- pr(object, "tau2B", PLATE, getParam(object, "tau2B", PLATE))				
	object <- pr(object, "sig2A", PLATE, getParam(object, "sig2A", PLATE))
	object <- pr(object, "sig2B", PLATE, getParam(object, "sig2B", PLATE))		
	object <- pr(object, "phiAX", PLATE, as.numeric(getParam(object, "phiAX", PLATE)))
	object <- pr(object, "phiBX", PLATE, as.numeric(getParam(object, "phiBX", PLATE)))
	object <- pr(object, "corr", PLATE, getParam(object, "corr", PLATE))
	object <- pr(object, "corrA.BB", PLATE, getParam(object, "corrA.BB", PLATE))
	object <- pr(object, "corrB.AA", PLATE, getParam(object, "corrB.AA", PLATE))		
	rm(object, tmp.objects); gc();
	##object <- object[order(chromosome(object), position(object)), ]
	if(cnOptions[["thresholdCopynumber"]])	object <- thresholdCopynumber(object)
	return(object)
}






