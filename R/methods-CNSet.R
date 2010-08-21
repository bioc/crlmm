linearParamElementReplace <- function(obj, elt, value) {
    storage.mode <- storageMode(batchStatistics(obj))
    switch(storage.mode,
           "lockedEnvironment" = {
               aData <- copyEnv(batchStatistics(obj))
               if (is.null(value)) rm(list=elt, envir=aData)
               else aData[[elt]] <- value
               Biobase:::assayDataEnvLock(aData)
               batchStatistics(obj) <- aData
           },
           "environment" = {
               if (is.null(value)) rm(list=elt, envir=batchStatistics(obj))
               else batchStatistics(obj)[[elt]] <- value
           },
           list = batchStatistics(obj)[[elt]] <- value)
    obj
}

## parameters
setMethod("nuA", signature=signature(object="CNSet"), function(object) nu(object, "A"))
setMethod("nuB", signature=signature(object="CNSet"), function(object) nu(object, "B"))
setMethod("phiA", signature=signature(object="CNSet"), function(object) phi(object, "A"))
setMethod("phiB", signature=signature(object="CNSet"), function(object) phi(object, "B"))
setMethod("phiPrimeA", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "phiPrimeA")
})
setMethod("phiPrimeB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "phiPrimeB")
})
setMethod("sigma2A", signature=signature(object="CNSet"), function(object) sigma2(object, "A"))
setMethod("sigma2B", signature=signature(object="CNSet"), function(object) sigma2(object, "B"))
setMethod("tau2A", signature=signature(object="CNSet"), function(object) tau2(object, "A"))
setMethod("tau2B", signature=signature(object="CNSet"), function(object) tau2(object, "B"))


setMethod("Ns", signature=signature(object="CNSet"),
	   function(object, ...){
		   Ns(batchStatistics(object), ...)
	   })
setMethod("corr", signature=signature(object="CNSet"),
	   function(object, ...){
		   corr(batchStatistics(object), ...)
	   })
setMethod("medians", signature=signature(object="CNSet"),
	   function(object, ...){
		   medians(batchStatistics(object), ...)
	   })
setMethod("mads", signature=signature(object="CNSet"),
	   function(object, ...){
		   mads(batchStatistics(object), ...)
	   })
setMethod("tau2", signature=signature(object="CNSet"),
	   function(object, ...){
		   tau2(batchStatistics(object), ...)
	   })
##---------------------------------------------------------------------------
## Number of samples with Genotype AA, AB, or BB by batch
setMethod("N.AA", signature=signature(object="CNSet"), function(object){
	assayDataElement(batchStatistics(object), "N.AA")
})
setMethod("N.AB", signature=signature(object="CNSet"), function(object){
	assayDataElement(batchStatistics(object), "N.AB")
})
setMethod("N.BB", signature=signature(object="CNSet"), function(object){
	assayDataElement(batchStatistics(object), "N.BB")
})
setReplaceMethod("N.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "N.AA", value)
	  })
setReplaceMethod("N.AB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "N.AB", value)
	  })
setReplaceMethod("N.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "N.BB", value)
	  })


##---------------------------------------------------------------------------
##  median intensity by genotype cluster for each allele
##should we update the entire matrix rather than one column...
setMethod("medianA.AA", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "medianA.AA")
})
setMethod("medianA.AB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "medianA.AB")
})
setMethod("medianA.BB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "medianA.BB")
})
setMethod("medianB.AA", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "medianB.AA")
})
setMethod("medianB.AB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "medianB.AB")
})
setMethod("medianB.BB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "medianB.BB")
})
setReplaceMethod("medianA.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "medianA.AA", value)
	  })
setReplaceMethod("medianA.AB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "medianA.AB", value)
	  })
setReplaceMethod("medianA.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "medianA.BB", value)
	  })
setReplaceMethod("medianB.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "medianB.AA", value)
	  })
setReplaceMethod("medianB.AB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "medianB.AB", value)
	  })
setReplaceMethod("medianB.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "medianB.BB", value)
	  })

##---------------------------------------------------------------------------
##  mad intensity by genotype cluster for each allele
setMethod("madA.AA", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "madA.AA")
})
setMethod("madA.AB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "madA.AB")
})
setMethod("madA.BB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "madA.BB")
})
setMethod("madB.AA", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "madB.AA")
})
setMethod("madB.AB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "madB.AB")
})
setMethod("madB.BB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "madB.BB")
})
setReplaceMethod("madA.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "madA.AA", value)
	  })
setReplaceMethod("madA.AB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "madA.AB", value)
	  })
setReplaceMethod("madA.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "madA.BB", value)
	  })
setReplaceMethod("madB.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "madB.AA", value)
	  })
setReplaceMethod("madB.AB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "madB.AB", value)
	  })
setReplaceMethod("madB.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "madB.BB", value)
	  })

##---------------------------------------------------------------------------
##  mad of log(intensities) by genotype cluster for each allele
setMethod("tau2A.AA", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "tau2A.AA")
})
setMethod("tau2A.BB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "tau2A.BB")
})
setMethod("tau2B.AA", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "tau2B.AA")
})
setMethod("tau2B.BB", signature=signature(object="CNSet"), function(object) {
	assayDataElement(batchStatistics(object), "tau2B.BB")
})
setReplaceMethod("tau2A.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "tau2A.AA", value)
	  })
setReplaceMethod("tau2A.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "tau2A.BB", value)
	  })
setReplaceMethod("tau2B.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "tau2B.AA", value)
	  })
setReplaceMethod("tau2B.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "tau2B.BB", value)
	  })

## correlation of log2A and log2B within each genotype cluster
setMethod("corrAA", signature=signature(object="CNSet"), function(object){
	assayDataElement(batchStatistics(object), "corrAA")
  })
setMethod("corrAB", signature=signature(object="CNSet"), function(object){
	assayDataElement(batchStatistics(object), "corrAB")
  })
setMethod("corrBB", signature=signature(object="CNSet"), function(object){
	assayDataElement(batchStatistics(object), "corrBB")
  })
setReplaceMethod("corrAA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "corrAA", value)
})

setReplaceMethod("corrAB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "corrAB", value)
})

setReplaceMethod("corrBB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "corrBB", value)
})

setReplaceMethod("phiPrimeA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "phiPrimeA", value)
})

setReplaceMethod("phiPrimeB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "phiPrimeB", value)
})

##setMethod("mad.AA", signature=signature(object="CNSet"), function(object) mads(object, "AA"))
##setMethod("mad.AB", signature=signature(object="CNSet"), function(object) mads(object, "AB"))
##setMethod("mad.BB", signature=signature(object="CNSet"), function(object) mads(object, "BB"))
##

##setReplaceMethod("median.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
##	  function(object, value){
##		  linearParamElementReplace(object, "median.AA", value)
##	  })
##setReplaceMethod("median.AB", signature=signature(object="CNSet", value="ff_or_matrix"),
##	  function(object, value){
##		  linearParamElementReplace(object, "median.AB", value)
##	  })
##setReplaceMethod("median.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
##	  function(object, value){
##		  linearParamElementReplace(object, "median.BB", value)
##	  })
##setReplaceMethod("mad.AA", signature=signature(object="CNSet", value="ff_or_matrix"),
##	  function(object, value){
##		  linearParamElementReplace(object, "mad.AA", value)
##	  })
##setReplaceMethod("mad.AB", signature=signature(object="CNSet", value="ff_or_matrix"),
##	  function(object, value){
##		  linearParamElementReplace(object, "mad.AB", value)
##	  })
##setReplaceMethod("mad.BB", signature=signature(object="CNSet", value="ff_or_matrix"),
##	  function(object, value){
##		  linearParamElementReplace(object, "mad.BB", value)
##	  })

setReplaceMethod("nuA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "nuA", value)
	  })

setReplaceMethod("nuB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "nuB", value)
})

setReplaceMethod("phiA", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "phiA", value)
})

setReplaceMethod("phiB", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "phiB", value)
})

setReplaceMethod("sigma2A", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "sig2A", value)
})

setReplaceMethod("sigma2B", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "sig2B", value)
})

setReplaceMethod("tau2A", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "tau2A", value)
})

setReplaceMethod("tau2B", signature=signature(object="CNSet", value="ff_or_matrix"),
	  function(object, value){
		  linearParamElementReplace(object, "tau2B", value)
})



setReplaceMethod("flags", signature=signature(object="CNSet", value="ff_or_matrix"),
		 function(object, value){
			 linearParamElementReplace(object, "flags", value)
})
setReplaceMethod("flags", signature=signature(object="CNSet", value="ff_matrix"),
		 function(object, value){
			 linearParamElementReplace(object, "flags", value)
})


##setMethod("ellipse", "CNSet", function(x, copynumber, batch, ...){
##	ellipse.CNSet(x, copynumber, batch, ...)
##})

ACN.X <- function(object, allele, i, j){
	acn <- list()
	batches <- unique(batch(object)[j])
	for(k in seq_along(batches)){
		l <- match(batches[k], bns)
		nuA <- nuA(object)[ii, l]
		nuB <- nuB(object)[ii, l]
		phiA <- phiA(object)[ii, l]
		phiB <- phiB(object)[ii, l]
		phiPrimeA <- phiPrimeA(object)[ii, l]
		phiPrimeB <- phiPrimeB(object)[ii, l]
		if(all(is.na(phiPrimeA))){
			I <- allele(object, allele)[ii, j]
			if(allele=="A") acn[[k]] <- 1/phiA*(I - nuA)
			if(allele=="B") acn[[k]] <- 1/phiB*(I - nuB)
		} else {
			A <- A(object)[ii, j]
			B <- B(object)[ii, j]
			phistar <- phiPrimeB/phiA
			tmp <- (B-nuB - phistar*A + phistar*nuA)/phiB
			cb <- tmp/(1-phistar*phiPrime/phiB)
			ca <- A-nuA-phiPrimeA*cb/phiA
			if(allele=="A") acn[[k]] <- ca
			if(allele=="B") acn[[k]] <- cb
		}
	}
	return(acn)
}


ACN <- function(object, allele, i , j){
	if(missing(i) & missing(j)) stop("must specify rows (i) or columns (j)")
	bns <- batchNames(object)
	acn <- list()
	if(missing(i) & !missing(j)){
		ii <- which(chromosome(object) < 23)
		if(length(ii) > 0){
			## calculate ca only for batches indexed by j
			batches <- unique(as.character(batch(object))[j])
			for(k in seq_along(batches)){
				this.batch <- batches[k]
				jj <- j[batch(object)[j] %in% this.batch]
				l <- match(this.batch, bns)
				bg <- nu(object, allele)[ii, l]
				slope <- phi(object, allele)[ii, l]
				I <- allele(object, allele)[ii, jj]
				acn[[k]] <- 1/slope*(I - bg)
			}
		}
		ii <- which(chromosome(object) == 23)
		if(length(ii) > 0){
			acn[[k]] <- ACN.X(object, allele, i, j)
		}
	}
	if(!missing(i) & missing(j)){
		## calculate ca, cb for all batches
		batches <- batchNames(object)
		for(k in seq_along(batches)){
			this.batch <- batches[k]
			l <- match(this.batch, bns)
			##bb <- batches[k]
			bg <- nu(object, allele)[i, l]
			slope <- phi(object, allele)[i, l]
			I <- allele(object, allele)[i, batch(object) == this.batch]
			acn[[k]] <- 1/slope*(I - bg)
		}
	}
	if(!missing(i) & !missing(j)){
		batches <- unique(as.character(batch(object)[j]))
		acn <- list()
		for(k in seq_along(batches)){
			this.batch <- batches[k]
			jj <- j[batch(object)[j] %in% this.batch]
			l <- match(this.batch, bns)
			bg <- nu(object, allele)[i, l]
			slope <- phi(object, allele)[i, l]
			I <- allele(object, allele)[i, jj]
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
		  ca[ca < 0.05] <- 0.05
		  ca[ca > 5] <- 5
		  return(ca)
	  })
setMethod("CB", signature=signature(object="CNSet"),
	  function(object, ...) {
		  cb <- ACN(object, allele="B", ...)
		  cb[cb < 0.05] <- 0.05
		  cb[cb > 5] <- 5
		  return(cb)
	  })

setMethod("totalCopynumber", signature=signature(object="CNSet"),
	  function(object, ...){
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
	  })

setReplaceMethod("snpCall", c("CNSet", "ff_or_matrix"),
                 function(object, ..., value){
			 assayDataElementReplace(object, "call", value)
		 })
setReplaceMethod("snpCallProbability", c("CNSet", "ff_or_matrix"),
                 function(object, ..., value){
			 assayDataElementReplace(object, "callProbability", value)
		 })
