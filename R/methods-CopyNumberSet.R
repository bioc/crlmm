setValidity("CopyNumberSet", function(object) {
	msg <- validMsg(NULL, assayDataValidMembers(assayData(object), c("CA", "CB")))
	if (is.null(msg)) TRUE else msg
})
##may want to allow thresholding here (... arg)
setMethod("CA", "CopyNumberSet", function(object, ...) assayData(object)[["CA"]]/100)
setMethod("CB", "CopyNumberSet", function(object, ...) assayData(object)[["CB"]]/100)

setReplaceMethod("CA", signature(object="CopyNumberSet", value="matrix"), function(object, value){
	dns <- dimnames(value)
	value <- matrix(as.integer(value*100), nrow(value), ncol(value))
	dimnames(value) <- dns
	assayDataElementReplace(object, "CA", value)
})

setReplaceMethod("CB", signature(object="CopyNumberSet", value="matrix"), function(object, value){
	dns <- dimnames(value)	
	value <- matrix(as.integer(value*100), nrow(value), ncol(value))
	dimnames(value) <- dns	
	assayDataElementReplace(object, "CB", value)
})




setMethod("batch", "CopyNumberSet", function(object){
	if("batch" %in% varLabels(object)){
		result <- object$batch
	} else {
		stop("batch not in varLabels of CopyNumberSet")
	}
	return(result)
})

setMethod("copyNumber", "CopyNumberSet", function(object){
	require(paste(annotation(object), "Crlmm", sep=""), character.only=TRUE) || stop(paste("Annotation package ", annotation(object), "Crlmm not available", sep=""))
	##ensure that 2 + NA = 2 by replacing NA's with zero
	##the above results in copy number 0, 1, or 2 depending on the genotype....safer just to drop
	CA <- CA(object)
	CB <- CB(object)
	##nas <- is.na(CA) & is.na(CB)
	##CA[is.na(CA)] <- 0
	##CB[is.na(CB)] <- 0
	CN <- CA + CB
	##For nonpolymorphic probes, CA is the total copy number
	CN[cnIndex(object, annotation(object)), ] <- CA(object)[cnIndex(object, annotation(object)), ]
	##if both CA and CB are NA, report NA
	##CN[nas] <- NA
	CN
})



##setMethod("ellipse", "CopyNumberSet", function(x, copynumber, ...){
ellipse.CopyNumberSet <- function(x, copynumber, ...){
##setMethod("ellipse", "CopyNumberSet", function(x, copynumber, ...){
	##fittedOrder <- unique(sapply(basename(celFiles), function(x) strsplit(x, "_")[[1]][2]))
	##index <- match(plates, fittedOrder)
	if(nrow(x) > 1) stop("only 1 snp at a time")
	##batch <- unique(x$batch)
	args <- list(...)
	if(!"batch" %in% names(args)){
		jj <- match("batch", varLabels(x))
		if(length(jj) < 1) stop("batch not in varLabels")
		batch <- unique(pData(x)[, jj])
	} else{
		batch <- unique(args$batch)
	}
	if(length(batch) > 1) stop("batch variable not unique")
	nuA <- as.numeric(fData(x)[, match(paste("nuA", batch, sep="_"), fvarLabels(x))])
	nuB <- as.numeric(fData(x)[, match(paste("nuB", batch, sep="_"), fvarLabels(x))])	
	phiA <- as.numeric(fData(x)[, match(paste("phiA", batch, sep="_"), fvarLabels(x))])
	phiB <- as.numeric(fData(x)[, match(paste("phiB", batch, sep="_"), fvarLabels(x))])	
	tau2A <- as.numeric(fData(x)[, match(paste("tau2A", batch, sep="_"), fvarLabels(x))])
	tau2B <- as.numeric(fData(x)[, match(paste("tau2B", batch, sep="_"), fvarLabels(x))])
	sig2A <- as.numeric(fData(x)[, match(paste("sig2A", batch, sep="_"), fvarLabels(x))])
	sig2B <- as.numeric(fData(x)[, match(paste("sig2B", batch, sep="_"), fvarLabels(x))])
	corrA.BB <- as.numeric(fData(x)[, match(paste("corrA.BB", batch, sep="_"), fvarLabels(x))])
	corrB.AA <- as.numeric(fData(x)[, match(paste("corrB.AA", batch, sep="_"), fvarLabels(x))])
	corr <- as.numeric(fData(x)[, match(paste("corr", batch, sep="_"), fvarLabels(x))])
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




