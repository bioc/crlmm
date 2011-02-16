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


## allele A
##   autosome SNPs
##   autosome NPs
##   chromosome X NPs for women
C1 <- function(object, marker.index, batch.index, sample.index){
##	acn <- vector("list", length(batch.index))
	acn <- matrix(NA, nrow=length(marker.index), ncol=length(sample.index))
	for(k in seq_along(batch.index)){
		l <- batch.index[k]
		jj <- sample.index[as.character(batch(object))[sample.index] == batchNames(object)[l]]
		bg <- nuA(object)[marker.index, l]
		slope <- phiA(object)[marker.index, l]
		I <- A(object)[marker.index, jj]
		acn[, match(jj, sample.index)] <- 1/slope * (I - bg)
	}
##	if(length(acn) > 1){
##		acn <- do.call("cbind", acn)
##	} else acn <- acn[[1]]
	return(as.matrix(acn))
}

## allele B  (treated allele 'A' for chromosome X NPs)
##   autosome SNPs
##   chromosome X for male nonpolymorphic markers
C2 <- function(object, marker.index, batch.index, sample.index, NP.X=FALSE){
	acn <- matrix(NA, nrow=length(marker.index), ncol=length(sample.index))
	for(k in seq_along(batch.index)){
		l <- batch.index[k]
		jj <- sample.index[as.character(batch(object))[sample.index] == batchNames(object)[l]]
		bg <- nuB(object)[marker.index, l]
		slope <- phiB(object)[marker.index, l]
		if(!NP.X){
			I <- B(object)[marker.index, jj]
		} else I <- A(object)[marker.index, jj]
		acn[, match(jj, sample.index)] <- 1/slope * (I - bg)
	}
##	if(length(acn) > 1){
##		acn <- do.call("cbind", acn)
##	} else acn <- acn[[1]]
	return(as.matrix(acn))
}

## Chromosome X SNPs
C3 <- function(object, allele, marker.index, batch.index, sample.index){
##	acn <- vector("list", length(batch.index))
	acn <- matrix(NA, nrow=length(marker.index), ncol=length(sample.index))
	for(k in seq_along(batch.index)){
		l <- batch.index[k]
		##j <- which(as.character(batch(object))[sample.index] == batchNames(object)[l])
		jj <- sample.index[as.character(batch(object))[sample.index] == batchNames(object)[l]]
		phiA2 <- phiPrimeA(object)[marker.index, l]
		phiB2 <- phiPrimeB(object)[marker.index, l]
		phiA <- phiA(object)[marker.index, l]
		phiB <- phiB(object)[marker.index, l]
		nuA <- nuA(object)[marker.index, l]
		nuB <- nuB(object)[marker.index, l]
		IA <- A(object)[marker.index, jj]
		IB <- B(object)[marker.index, jj]
		phistar <- phiB2/phiA
		tmp <- (IB - nuB - phistar*IA + phistar*nuA)/phiB
		CB <- tmp/(1-phistar*phiA2/phiB)
		##CB <- 1/(1-phiA2*phiB2/(phiA*phiB)) * 1/phiB * (IA-nuB-phiB2/phiA*(IA-nuA))
		CA <- (IA-nuA-phiA2*CB)/phiA
		if(allele == "B"){
			acn[, match(jj, sample.index)] <- CB
			##acn[[k]] <- CB
		}
		if(allele == "A"){
			acn[, match(jj, sample.index)] <- (IA-nuA-phiA2*CB)/phiA
		}
		if(allele == "AandB"){
			CA <- tmp/(1-phistar*phiA2/phiB)
			CB <- (IA-nuA-phiA2*CB)/phiA
			acn[, match(jj, sample.index)] <- (IA-nuA-phiA2*CB)/phiA
		}
##		if(allele=="AandB")
##			CA <- tmp/(1-phistar*phiA2/phiB)
##			CB <- (IA-nuA-phiA2*CB)/phiA
##			acn[[k]] <- CA+CB
##		}
	}
##	if(length(acn) > 1){
##		acn <- do.call("cbind", acn)
##	} else acn <- acn[[1]]
	return(as.matrix(acn))
}




ACN <- function(object, allele, i , j){
	if(missing(i) & missing(j)) stop("must specify rows (i) or columns (j)")
	is.ff <- is(calls(object), "ff") | is(calls(object), "ffdf")
	missing.i <- missing(i)
	missing.j <- missing(j)
	if(!missing.i){
		is.ann <- !is.na(chromosome(object)[i])
		is.X <- chromosome(object)[i]==23 & is.ann
		is.auto <- chromosome(object)[i] < 23 & is.ann
		is.snp <- isSnp(object)[i] & is.ann
	} else{
		is.ann <- !is.na(chromosome(object))
		is.X <- chromosome(object)==23 & is.ann
		is.auto <- chromosome(object) < 23 & is.ann
		is.snp <- isSnp(object) & is.ann
		i <- 1:nrow(object)
	}
	## Define batch.index and sample.index
	if(!missing.j) {
		batches <- unique(as.character(batch(object))[j])
		##batches <- as.character(batch(object)[j])
		batch.index <- match(batches, batchNames(object))
	} else {
		batch.index <- seq_along(batchNames(object))
		j <- 1:ncol(object)
	}
	nr <- length(i)
	nc <- length(j)
	acn <- matrix(NA, nr, nc)
	dimnames(acn) <- list(featureNames(object)[i],
			      sampleNames(object)[j])
	if(allele == "A"){
		if(is.ff){
			open(nuA(object))
			open(phiA(object))
			open(A(object))
		}
		## --
		## 4 types of markers for allele A
		##--
		## 1. autosomal SNPs or autosomal NPs
		if(any(is.auto)){
			auto.index <- which(is.auto)
			marker.index <- i[is.auto]
			acn[auto.index, ] <- C1(object, marker.index, batch.index, j)
		}
		if(any(is.X)){
			##2. CHR X SNPs (men and women)
			if(any(is.snp)){
				if(is.ff) {
					open(phiPrimeA(object))
					open(phiPrimeB(object))
					open(phiB(object))
					open(nuB(object))
					open(B(object))
				}
				marker.index <- i[is.X & is.snp]
				acn.index <- which(is.X & is.snp)
				acn[acn.index, ] <- C3(object, allele="A", marker.index, batch.index, j)
				if(is.ff) {
					close(phiPrimeA(object))
					close(phiPrimeB(object))
					close(phiB(object))
					close(nuB(object))
					close(B(object))
				}
			}
			if(any(!is.snp)){
				marker.index <- i[is.X & !is.snp]
				acn.index <- which(is.X & !is.snp)
				acn[acn.index, ] <- NA
				female.index <- j[object$gender[j] == 2]
				## 3. CHR X NPs: women
				if(length(female.index) > 0){
					female.batch.index <- match(unique(as.character(batch(object))[female.index]), batchNames(object))
					jj <- which(object$gender[j] == 2)
					acn[acn.index, jj] <- C1(object, marker.index, female.batch.index, female.index)
				}
				male.index <- j[object$gender[j] == 1]
				if(length(male.index) > 0){
					if(is.ff){
						open(nuB(object))
						open(phiB(object))
					}
					male.batch.index <- match(unique(as.character(batch(object))[male.index]), batchNames(object))
					jj <- which(object$gender[j] == 1)
					acn[acn.index, jj] <- C2(object, marker.index, male.batch.index, male.index, NP.X=TRUE)
					if(is.ff){
						close(nuB(object))
						close(phiB(object))
					}
				}
			}
		}
		if(is.ff){
			close(nuA(object))
			close(phiA(object))
			close(A(object))
		}
	}
	if(allele == "B"){
		if(is.ff){
			open(nuB(object))
			open(phiB(object))
			open(B(object))
		}
		if(any(!is.snp)){
			acn.index <- which(!is.snp)
			acn[acn.index, ] <- 0
		}
		if(any(is.auto)){
			auto.index <- which(is.auto & is.snp)
			if(length(auto.index) > 0){
				marker.index <- i[is.auto]
				acn[auto.index, ] <- C2(object, marker.index, batch.index, j)
			}
		}
		if(any(is.X)){
			if(is.ff){
				open(phiPrimeA(object))
				open(phiPrimeB(object))
				open(phiA(object))
				open(nuA(object))
				open(A(object))
			}
			marker.index <- i[is.X & is.snp]
			acn.index <- which(is.X & is.snp)
			acn[acn.index, ] <- C3(object, allele="B", marker.index, batch.index, j)
			if(is.ff){
				close(phiPrimeA(object))
				close(phiPrimeB(object))
				close(phiA(object))
				close(nuA(object))
				close(A(object))
			}
			if(any(!is.snp)){
				acn.index <- which(!is.snp)
				marker.index <- i[!is.snp]
				acn[acn.index, ] <- 0
			}
		}
	}
	if(is.ff){
		close(nuB(object))
		close(phiB(object))
		close(B(object))
	}
	return(acn)
}

setMethod("CA", signature=signature(object="CNSet"),
	  function(object, ...){
		  ca <- ACN(object, allele="A", ...)
		  ca[ca < 0] <- 0
		  ca[ca > 5] <- 5
		  return(ca)
	  })
setMethod("CB", signature=signature(object="CNSet"),
	  function(object, ...) {
		  cb <- ACN(object, allele="B", ...)
		  cb[cb < 0] <- 0
		  cb[cb > 5] <- 5
		  return(cb)
	  })

setMethod("totalCopynumber", signature=signature(object="CNSet"),
	  function(object, ...){
		  ca <- CA(object, ...)
		  cb <- CB(object, ...)
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
