setAs("CallSet", "CNSet", function(from, to){
	if(is.null(from$batch)) stop("Before coercing to CNSet, specify 'batch' in the phenoData")
	lmp <- getLinearModelParam(from)
	stopifnot(class(lmp) != "list" || class(lmp) != "ffdf")
	CA <- initializeBigMatrix(list(featureNames(from), sampleNames(from)))
	CB <- initializeBigMatrix(list(featureNames(from), sampleNames(from)))
	dimnames(CA) <- dimnames(A(from))
	dimnames(CB) <- dimnames(A(from))
	object <- new("CNSet",
		      alleleA=A(from),
		      alleleB=B(from),		  
		      call=calls(from),
		      callProbability=confs(from),
		      CA=CA,
		      CB=CB,
		      phenoData=phenoData(from),
		      experimentData=experimentData(from),
		      annotation=annotation(from),
		      protocolData=protocolData(from),
		      featureData=featureData(from),
		      genomeAnnotation=genomeAnnotation(from),
		      options=crlmmOptions(from),
		      lm=lmp)
	return(object)
})

setMethod("getLinearModelParam", "CallSet", function(object, ...){
	if(is.null(object$batch)) stop("batch must be specified in phenoData")
	dns <- list(featureNames(object), unique(object$batch))
        nr <- length(dns[[1]])
        nc <- length(dns[[2]])             
        ll <- vector("list", 17)
	if(isPackageLoaded("ff")){
		for(i in 1:17) ll[[i]] <- ff(dim=c(nr,nc), vmode="double", finalizer="close", dimnames=dns, overwrite=TRUE)
		names(ll) <- paramNames()
		ll <- do.call(ffdf, ll)
	} else {
		for(i in 1:17) ll[[i]] <- matrix(NA, nr, nc, dimnames=dns)
	}
        return(ll)
})

.getConfs <- function(object) assayDataElement(object, "callProbability")
	##if(class(value)[1] == "matrix") value <- 1-exp(-value/1000)
	##return(value)
##}
setMethod("confs", "CallSet", .getConfs)
##setMethod("confs", "AffymetrixCallSet", .getConfs)
setMethod("calls", "CallSet", function(object) assayDataElement(object, "call"))
##setMethod("calls", "IlluminaCallSet", function(object) assayDataElement(object, "call"))
setMethod("snpCall", "CallSet", calls)
##setMethod("snpCall", "IlluminaCallSet", calls)
setReplaceMethod("calls", signature(object="CallSet", value="matrix_or_ff"),
		 function(object, value) assayDataElementReplace(object, "call", value))
##setReplaceMethod("calls", signature(object="IlluminaCallSet", value="matrix_or_ff"),
##		 function(object, value) assayDataElementReplace(object, "call", value))
setReplaceMethod("confs", "CallSet", value="matrix_or_ff",
		 function(object, value) assayDataElementReplace(object, "callProbability", value))
##setReplaceMethod("confs", "IlluminaCallSet", value="matrix_or_ff",
##		 function(object, value) assayDataElementReplace(object, "callProbability", value))
##setReplaceMethod("confs", signature(object="CrlmmContainer", value="matrix"),
##		 function(object, value){
##			 ##convert probability to integer
##			 P <- value
##			 dns <- dimnames(P)
##			 X <- -1000*log(1-P)
##			 X <- matrix(as.integer(X), nrow(X), ncol(X))
##			 dimnames(X) <- dns
##			 assayDataElementReplace(object, "callProbability", X)
##		 })
crlmmCopynumber <- function(object){
	object <- as(object, "CNSet")
	ops <- crlmmOptions(object)
	verbose <- ops$verbose
	fns <- featureNames(object)
	SNRmin <- ops$SNRMin
	batch <- object$batch
	whichBatch <- ops$cnOpts$whichBatch
	if(is.null(whichBatch)) whichBatch <- seq(along=unique(batch))
	chromosome <- ops$cnOpts$chromosome
	MIN.SAMPLES <- ops$cnOpts$MIN.SAMPLES
	for(CHR in chromosome){
		cat("Chromosome ", CHR, "\n")
		if(CHR > 23) next()
		for(i in whichBatch){
			PLATE <- unique(batch)[i]
			message("Plate: ", PLATE)
			col.index <- which(batch==PLATE)
			if(length(col.index) < MIN.SAMPLES) {
				warning("Plate ", PLATE, " has fewer than 10 samples.  Skipping this plate.")
				next()
			}
			row.index <- which(chromosome(object) == CHR)
			object <- object[row.index, col.index]
			tmp <- computeCopynumber(object[row.index, col.index])
			CA(object)[row.index, col.index] <- CA(tmp)
			CB(object)[row.index, col.index] <- CB(tmp)
		} ## end of batch loop
	} ## end of chromosome loop
	new("CNSet",
	    alleleA=A(cnSet),
	    alleleB=B(cnSet),
	    call=calls(cnSet),
	    callProbability=confs(cnSet),
	    CA=CA(cnSet),
	    CB=CB(cnSet),
	    experimentData=experimentData(cnSet),
	    phenoData=phenoData(cnSet),
	    featureData=featureData(cnSet),
	    protocolData=protocolData(cnSet),
	    annotation=annotation(cnSet))
}
