linearParamElementReplace <- function(obj, elt, value) {
    storage.mode <- storageMode(lM(obj))
    switch(storage.mode,
           "lockedEnvironment" = {
               aData <- copyEnv(lM(obj))
               if (is.null(value)) rm(list=elt, envir=aData)
               else aData[[elt]] <- value
               Biobase:::assayDataEnvLock(aData)
               lM(obj) <- aData
           },
           "environment" = {
               if (is.null(value)) rm(list=elt, envir=lM(obj))
               else lM(obj)[[elt]] <- value
           },
           list = lM(obj)[[elt]] <- value)
    obj
}


setMethod("nu", c("CNSet", "character"), function(object, allele) nu(lM(object), allele))
setMethod("phi", c("CNSet", "character"), function(object, allele) phi(lM(object), allele))
setMethod("sigma2", c("CNSet", "character"), function(object, allele) phi(lM(object), allele))
setMethod("tau2", c("CNSet", "character"), function(object, allele) phi(lM(object), allele))
setMethod("corr", c("CNSet", "character"), function(object, allele) phi(lM(object), allele))

setMethod("nuA", signature=signature(object="CNSet"), function(object) nu(object, "A"))
setMethod("nuB", signature=signature(object="CNSet"), function(object) nu(object, "B"))
setMethod("phiA", signature=signature(object="CNSet"), function(object) phi(object, "A"))
setMethod("phiB", signature=signature(object="CNSet"), function(object) phi(object, "B"))
setMethod("sigma2A", signature=signature(object="CNSet"), function(object) sigma2(object, "A"))
setMethod("sigma2B", signature=signature(object="CNSet"), function(object) sigma2(object, "B"))
setMethod("tau2A", signature=signature(object="CNSet"), function(object) tau2(object, "A"))
setMethod("tau2B", signature=signature(object="CNSet"), function(object) tau2(object, "B"))
setMethod("corrAA", signature=signature(object="CNSet"), function(object) corr(object, "AA"))
setMethod("corrAB", signature=signature(object="CNSet"), function(object) corr(object, "AB"))
setMethod("corrBB", signature=signature(object="CNSet"), function(object) corr(object, "BB"))

setReplaceMethod("nuA", signature=signature(object="CNSet", value="matrix"), 
	  function(object, value){
		  linearParamElementReplace(object, "nuA", value)
	  })

setReplaceMethod("nuB", signature=signature(object="CNSet", value="matrix"), 
	  function(object, value){
		  linearParamElementReplace(object, "nuB", value)		  
})

setReplaceMethod("phiA", signature=signature(object="CNSet", value="matrix"), 
	  function(object, value){
		  linearParamElementReplace(object, "phiA", value)		  
})

setReplaceMethod("phiB", signature=signature(object="CNSet", value="matrix"), 
	  function(object, value){
		  linearParamElementReplace(object, "phiB", value)		  
})

setReplaceMethod("sigma2A", signature=signature(object="CNSet", value="matrix"), 
	  function(object, value){
		  linearParamElementReplace(object, "sig2A", value)		  
})

setReplaceMethod("sigma2B", signature=signature(object="CNSet", value="matrix"), 
	  function(object, value){
		  linearParamElementReplace(object, "sig2B", value)		  
})

setReplaceMethod("tau2A", signature=signature(object="CNSet", value="matrix"), 
	  function(object, value){
		  linearParamElementReplace(object, "tau2A", value)		  
})

setReplaceMethod("tau2B", signature=signature(object="CNSet", value="matrix"), 
	  function(object, value){
		  linearParamElementReplace(object, "tau2B", value)		  
})

setReplaceMethod("corrAA", signature=signature(object="CNSet", value="matrix"), 
	  function(object, value){
		  linearParamElementReplace(object, "corrAA", value)		  
})

setReplaceMethod("corrAB", signature=signature(object="CNSet", value="matrix"), 
	  function(object, value){
		  linearParamElementReplace(object, "corrAB", value)		  
})

setReplaceMethod("corrBB", signature=signature(object="CNSet", value="matrix"), 
	  function(object, value){
		  linearParamElementReplace(object, "corrBB", value)		  
})


##setValidity("CNSet",
##	    function(object){
##		    if(!"batch" %in% varLabels(protocolData(object)))
##			    return("'batch' not defined in protocolData")
##		    if(!"chromosome" %in% fvarLabels(object))
##			    return("'chromosome' not defined in featureData")
##		    if(!"position" %in% fvarLabels(object))
##			    return("'position' not defined in featureData")
##		    if(!"isSnp" %in% fvarLabels(object))
##			    return("'isSnp' not defined in featureData")
##		    return(TRUE)
##	    })

setMethod("totalCopyNumber", "CNSet", function(object, i, j){
	if(missing(i) & missing(j)){
		if(inherits(CA(object), "ff") | inherits(CA(object), "ffdf")) stop("Must specify i and/or j for ff objects")
	}
	if(missing(i) & !missing(j)){
		snp.index <- which(isSnp(object))	
		cn.total <- as.matrix(CA(cnSet)[, j])
		cb <- as.matrix(CB(cnSet)[snp.index, j]	)
		cn.total[snp.index, ] <- cn.total[snp.index, ] + cb		
	}
	if(!missing(i) & missing(j)){
		snp.index <- intersect(which(isSnp(object)), i)
		cn.total <- as.matrix(CA(cnSet)[i, ])
		cb <- as.matrix(CB(cnSet)[snp.index, ])	
		cn.total[snp.index, ] <- cn.total[snp.index, ] + cb				
	}
	if(!missing(i) & !missing(j)){
		snp.index <- intersect(which(isSnp(object)), i)		
		cn.total <- as.matrix(CA(cnSet)[i, j])	
		cb <- as.matrix(CB(cnSet)[snp.index, j])
		cn.total[snp.index, ] <- cn.total[snp.index, ] + cb
	}
	cn.total <- cn.total/100
	dimnames(cn.total) <- NULL
	return(cn.total)
})

##setMethod("ellipse", "CNSet", function(x, copynumber, batch, ...){
##	ellipse.CNSet(x, copynumber, batch, ...)
##})



ACN <- function(object, allele, i , j){
	if(missing(i) & missing(j)){
		if(inherits(A(object), "ff") | inherits(A(object), "ffdf")) stop("Must specify i and/or j for ff objects")
	}
	if(missing(i) & !missing(j)){
		## calculate ca only for batches indexed by j
		ubatch <- unique(batch(object))
		batches <- unique(batch(object)[j])
		for(k in seq_along(batches)){
			l <- match(batches[k], ubatch)
			bg <- nu(object, allele)[, l]
			sl <- phi(object, allele)[, l]
			I <- allele(object, allele)[, j]
			acn <- 1/sl*(I - bg)				  
		}
	}
	if(!missing(i) & missing(j)){
		## calculate ca, cb for all batches
		batches <- unique(batch(object))
		for(k in seq_along(batches)){
			##bb <- batches[k]
			bg <- nu(object, allele)[i, k]
			sl <- phi(object, allele)[i, k]
			I <- allele(object, allele)[i, j]
			acn <- 1/sl*(I - bg)
		}
	}
	if(!missing(i) & !missing(j)){
		ubatch <- unique(batch(object))
		batches <- unique(batch(object)[j])
		for(k in seq_along(batches)){
			l <- match(batches[k], ubatch)
			bg <- nu(object, allele)[i, l]
			sl <- phi(object, allele)[i, l]
			I <- allele(object, allele)[i, j]
			acn <- 1/bg*(I - sl)				  
		}			  
	}
	return(acn)
}

setMethod("CA",
	  signature=signature(object="CNSet", i="integerOrMissing", j="integerOrMissing"),
	  function(object, i, j) {
		  ca <- ACN(object, allele="A", i, j)
		  return(ca)
	  })
setMethod("CB",
	  signature=signature(object="CNSet", i="integerOrMissing", j="integerOrMissing"),
	  function(object, i, j) {
		  cb <- ACN(object, allele="B", i, j)
		  return(cb)
	  })

setMethod("totalCopyNumber",
	  signature=signature(object="CNSet", i="integerOrMissing", j="integerOrMissing"),
	  function(object, i, j, ...){
	if(missing(i) & missing(j)){
		if(inherits(A(object), "ff") | inherits(A(object), "ffdf")) stop("Must specify i and/or j for ff objects")
	}
	if(missing(i) & !missing(j)){
		snp.index <- which(isSnp(object))	
		cn.total <- as.matrix(CA(object)[, j])
		if(length(snp.index) > 0){
			cb <- as.matrix(CB(object)[snp.index, j])
			snps <- (1:nrow(cn.total))[i %in% snp.index]
			cn.total[snps, ] <- cn.total[snps, j] + cb				
		}
	}
	if(!missing(i) & missing(j)){
		snp.index <- intersect(which(isSnp(object)), i)
		cn.total <- as.matrix(CA(object)[i, ])
		if(length(snp.index) > 0){
			cb <- as.matrix(CB(object)[snp.index, ])
			snps <- (1:nrow(cn.total))[i %in% snp.index]
			cn.total[snps, ] <- cn.total[snps, ] + cb				
		}
	}
	if(!missing(i) & !missing(j)){
		snp.index <- intersect(which(isSnp(object)), i)		
		cn.total <- as.matrix(CA(object)[i, j])
		if(length(snp.index) > 0){
			cb <- as.matrix(CB(object)[snp.index, j])
			snps <- (1:nrow(cn.total))[i %in% snp.index]
			cn.total[snps, ] <- cn.total[snps, ] + cb
		}
	}
	##cn.total <- cn.total/100
	dimnames(cn.total) <- NULL
	return(cn.total)
})

setReplaceMethod("snpCall", c("CNSet", "ff_or_matrix"),
                 function(object, ..., value){
			 assayDataElementReplace(object, "call", value)
		 })
setReplaceMethod("snpCallProbability", c("CNSet", "ff_or_matrix"),
                 function(object, ..., value){
			 assayDataElementReplace(object, "callProbability", value)
		 })
