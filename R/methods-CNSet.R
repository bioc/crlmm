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


##setMethod("ellipse", "CNSet", function(x, copynumber, batch, ...){
##	ellipse.CNSet(x, copynumber, batch, ...)
##})


ACN <- function(object, allele, i , j){
	bns <- batchNames(object)
	acn <- list()
	if(missing(i) & !missing(j)){
		## calculate ca only for batches indexed by j
		batches <- unique(batch(object)[j])
		for(k in seq_along(batches)){
			l <- match(batches[k], bns)
			bg <- nu(object, allele)[, l]
			slope <- phi(object, allele)[, l]
			I <- allele(object, allele)[, j]
			acn[[k]] <- 1/slope*(I - bg)
		}
	}
	if(!missing(i) & missing(j)){
		## calculate ca, cb for all batches
		batches <- unique(batch(object))
		for(k in seq_along(batches)){
			##bb <- batches[k]
			bg <- nu(object, allele)[i, k]
			slope <- phi(object, allele)[i, k]
			I <- allele(object, allele)[i, j]
			acn[[k]] <- 1/slope*(I - bg)
		}
	}
	if(!missing(i) & !missing(j)){
		ubatch <- unique(batch(object))
		batches <- unique(batch(object)[j])
		acn <- list()
		for(k in seq_along(batches)){
			l <- match(batches[k], ubatch)
			bg <- nu(object, allele)[i, l]
			slope <- phi(object, allele)[i, l]
			I <- allele(object, allele)[i, j]
			acn[[k]] <- 1/slope*(I - bg)
		}
	}
	if(length(acn) > 1) acn <- do.call("cbind", acn)
	if(length(acn) == 1) acn <- acn[[1]]
	return(acn)
}

setMethod("CA", signature=signature(object="CNSet"),
	  function(object, ...){
		  ca <- ACN(object, allele="A", ...)
		  return(ca)
	  })
setMethod("CB", signature=signature(object="CNSet"),
	  function(object, ...) {
		  cb <- ACN(object, allele="B", ...)
		  return(cb)
	  })

##setMethod("totalCopyNumber", signature=signature(object="CNSet"),
totalCopyNumber <- function(object, ..., verbose=TRUE, dimnames=FALSE){
	ca <- CA(object, ...)
	cb <- CB(object, ...)
	is.snp <- isSnp(object)
	dotArgs <- list(...)
	if("i" %in% names(dotArgs)){
		i <- dotArgs[["i"]]
		np.index <- which(!is.snp[i])
		if(length(np.index) > 0) cb[np.index, ] <- 0
	} else {
		np.index <- which(!is.snp)
		if(length(np.index) > 0) cb[np.index, ] <- 0
	}
	return(ca+cb)
}

##
##
##
##
##		  if(missing(i) & !missing(j)){
##			  ca <- CA(object, i, j)
##			  snp.index <- which(isSnp(object))
##			  cn.total <- as.matrix(CA(object)[, j])
##			  if(length(snp.index) > 0){
##				  cb <- as.matrix(CB(object)[snp.index, j])
##				  snps <- (1:nrow(cn.total))[i %in% snp.index]
##				  cn.total[snps, ] <- cn.total[snps, j] + cb
##			  }
##		  }
##		  if(!missing(i) & missing(j)){
##			  snp.index <- intersect(which(isSnp(object)), i)
##			  cn.total <- as.matrix(CA(object)[i, ])
##			  if(length(snp.index) > 0){
##				  cb <- as.matrix(CB(object)[snp.index, ])
##				  snps <- (1:nrow(cn.total))[i %in% snp.index]
##				  cn.total[snps, ] <- cn.total[snps, ] + cb
##			  }
##		  }
##		  if(!missing(i) & !missing(j)){
##			  snp.index <- intersect(which(isSnp(object)), i)
##			  cn.total <- as.matrix(CA(object)[i, j])
##			  if(length(snp.index) > 0){
##				  cb <- as.matrix(CB(object)[snp.index, j])
##				  snps <- (1:nrow(cn.total))[i %in% snp.index]
##				  cn.total[snps, ] <- cn.total[snps, ] + cb
##			  }
##		  }
##		  ##cn.total <- cn.total/100
##		  dimnames(cn.total) <- NULL
##		  return(cn.total)
##	  })


setReplaceMethod("snpCall", c("CNSet", "ff_or_matrix"),
                 function(object, ..., value){
			 assayDataElementReplace(object, "call", value)
		 })
setReplaceMethod("snpCallProbability", c("CNSet", "ff_or_matrix"),
                 function(object, ..., value){
			 assayDataElementReplace(object, "callProbability", value)
		 })
