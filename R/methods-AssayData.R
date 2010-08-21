setMethod("Ns", signature(object="AssayData"),
	  function(object, i, j, ...){
		  if(!missing(j)){
			  batchnames <- unique(as.character(batch(object)[j]))
		  } else batchnames <- batchNames(object)
		  nc <- length(batchnames)
		  if(!missing(i)) nr <- length(i) else nr <- nrow(object)
		  res <- array(NA, dim=c(nr, 3, nc))
		  dimnames(res)[[2]] <- c("AA", "AB", "BB")
		  dimnames(res)[[3]] <- batchnames
		  if(missing(i) & missing(j)){
			  N.AA <- as.matrix(assayDataElement(object, "N.AA"))
			  N.AB <- as.matrix(assayDataElement(object, "N.AB"))
			  N.BB <- as.matrix(assayDataElement(object, "N.BB"))
		  }
		  if(missing(i) & !missing(j)){
			  J <- match(batchnames, batchNames(object))
			  stopifnot(length(J) > 0 & !all(is.na(J)))
			  N.AA <- as.matrix(assayDataElement(object, "N.AA"))[, J, ...]
			  N.AB <- as.matrix(assayDataElement(object, "N.AB"))[, J, ...]
			  N.BB <- as.matrix(assayDataElement(object, "N.BB"))[, J, ...]
		  }
		  if(!missing(i) & !missing(j)){
			  J <- match(batchnames, batchNames(object))
			  stopifnot(length(J) > 0 & !all(is.na(J)))
			  N.AA <- as.matrix(assayDataElement(object, "N.AA"))[i, J, ...]
			  N.AB <- as.matrix(assayDataElement(object, "N.AB"))[i, J, ...]
			  N.BB <- as.matrix(assayDataElement(object, "N.BB"))[i, J, ...]
		  }
		  if(!missing(i) & missing(j)){
			  N.AA <- as.matrix(assayDataElement(object, "N.AA"))[i, , ...]
			  N.AB <- as.matrix(assayDataElement(object, "N.AB"))[i, , ...]
			  N.BB <- as.matrix(assayDataElement(object, "N.BB"))[i, , ...]
		  }
		  res[, "AA", ] <- N.AA
		  res[, "AB", ] <- N.AB
		  res[, "BB", ] <- N.BB
		  return(res)
	  })
setMethod("corr", signature(object="AssayData"),
	  function(object, i, j, ...){
		  if(!missing(j)){
			  batchnames <- unique(as.character(batch(object)[j]))
		  } else batchnames <- batchNames(object)
		  nc <- length(batchnames)
		  if(!missing(i)) nr <- length(i) else nr <- nrow(object)
		  res <- array(NA, dim=c(nr, 3, nc))
		  dimnames(res)[[2]] <- c("AA", "AB", "BB")
		  dimnames(res)[[3]] <- batchnames
		  if(missing(i) & missing(j)){
			  corrAA <- as.matrix(assayDataElement(object, "corrAA"))
			  corrAB <- as.matrix(assayDataElement(object, "corrAB"))
			  corrBB <- as.matrix(assayDataElement(object, "corrBB"))
		  }
		  if(missing(i) & !missing(j)){
			  J <- match(batchnames, batchNames(object))
			  stopifnot(length(J) > 0 & !all(is.na(J)))
			  corrAA <- as.matrix(assayDataElement(object, "corrAA"))[, J, ...]
			  corrAB <- as.matrix(assayDataElement(object, "corrAB"))[, J, ...]
			  corrBB <- as.matrix(assayDataElement(object, "corrBB"))[, J, ...]
		  }
		  if(!missing(i) & !missing(j)){
			  J <- match(batchnames, batchNames(object))
			  stopifnot(length(J) > 0 & !all(is.na(J)))
			  corrAA <- as.matrix(assayDataElement(object, "corrAA"))[i, J, ...]
			  corrAB <- as.matrix(assayDataElement(object, "corrAB"))[i, J, ...]
			  corrBB <- as.matrix(assayDataElement(object, "corrBB"))[i, J, ...]
		  }
		  if(!missing(i) & missing(j)){
			  corrAA <- as.matrix(assayDataElement(object, "corrAA"))[i, , ...]
			  corrAB <- as.matrix(assayDataElement(object, "corrAB"))[i, , ...]
			  corrBB <- as.matrix(assayDataElement(object, "corrBB"))[i, , ...]
		  }
		  res[, "AA", ] <- corrAA
		  res[, "AB", ] <- corrAB
		  res[, "BB", ] <- corrBB
		  return(res)
	  })

setMethod("medians", signature(object="AssayData"),
	  function(object, i, j, ...){
		  if(!missing(j)){
			  batchnames <- unique(as.character(batch(object)[j]))
		  } else batchnames <- batchNames(object)
		  nc <- length(batchnames)
		  if(!missing(i)) nr <- length(i) else nr <- nrow(object)
		  res <- array(NA, dim=c(nr, 2, 3, nc))
		  dimnames(res)[[2]] <- c("A", "B")
		  dimnames(res)[[3]] <- c("AA", "AB", "BB")
		  dimnames(res)[[4]] <- batchnames
		  if(missing(i) & missing(j)){
			  medianA.AA <- as.matrix(assayDataElement(object, "medianA.AA"))
			  medianA.AB <- as.matrix(assayDataElement(object, "medianA.AB"))
			  medianA.BB <- as.matrix(assayDataElement(object, "medianA.BB"))
			  medianB.AA <- as.matrix(assayDataElement(object, "medianB.AA"))
			  medianB.AB <- as.matrix(assayDataElement(object, "medianB.AB"))
			  medianB.BB <- as.matrix(assayDataElement(object, "medianB.BB"))
		  }
		  if(missing(i) & !missing(j)){
			  J <- match(batchnames, batchNames(object))
			  stopifnot(length(J) > 0 & !all(is.na(J)))
			  medianA.AA <- as.matrix(assayDataElement(object, "medianA.AA"))[, J, ...]
			  medianA.AB <- as.matrix(assayDataElement(object, "medianA.AB"))[, J, ...]
			  medianA.BB <- as.matrix(assayDataElement(object, "medianA.BB"))[, J, ...]
			  medianB.AA <- as.matrix(assayDataElement(object, "medianB.AA"))[, J, ...]
			  medianB.AB <- as.matrix(assayDataElement(object, "medianB.AB"))[, J, ...]
			  medianB.BB <- as.matrix(assayDataElement(object, "medianB.BB"))[, J, ...]

		  }
		  if(!missing(i) & !missing(j)){
			  J <- match(batchnames, batchNames(object))
			  stopifnot(length(J) > 0 & !all(is.na(J)))
			  medianA.AA <- as.matrix(assayDataElement(object, "medianA.AA"))[i, J, ...]
			  medianA.AB <- as.matrix(assayDataElement(object, "medianA.AB"))[i, J, ...]
			  medianA.BB <- as.matrix(assayDataElement(object, "medianA.BB"))[i, J, ...]
			  medianB.AA <- as.matrix(assayDataElement(object, "medianB.AA"))[i, J, ...]
			  medianB.AB <- as.matrix(assayDataElement(object, "medianB.AB"))[i, J, ...]
			  medianB.BB <- as.matrix(assayDataElement(object, "medianB.BB"))[i, J, ...]
		  }
		  if(!missing(i) & missing(j)){
			  medianA.AA <- as.matrix(assayDataElement(object, "medianA.AA"))[i, ...]
			  medianA.AB <- as.matrix(assayDataElement(object, "medianA.AB"))[i, ...]
			  medianA.BB <- as.matrix(assayDataElement(object, "medianA.BB"))[i, ...]
			  medianB.AA <- as.matrix(assayDataElement(object, "medianB.AA"))[i, ...]
			  medianB.AB <- as.matrix(assayDataElement(object, "medianB.AB"))[i, ...]
			  medianB.BB <- as.matrix(assayDataElement(object, "medianB.BB"))[i, ...]
		  }
		  res[, "A", "AA", ] <- medianA.AA
		  res[, "A", "AB", ] <- medianA.AB
		  res[, "A", "BB", ] <- medianA.BB
		  res[, "B", "AA", ] <- medianB.AA
		  res[, "B", "AB", ] <- medianB.AB
		  res[, "B", "BB", ] <- medianB.BB
		  return(res)
})

setMethod("medians", signature(object="AssayData"),
	  function(object, i, j, ...){
		  if(!missing(j)){
			  batchnames <- unique(as.character(batch(object)[j]))
		  } else batchnames <- batchNames(object)
		  nc <- length(batchnames)
		  if(!missing(i)) nr <- length(i) else nr <- nrow(object)
		  res <- array(NA, dim=c(nr, 2, 3, nc))
		  dimnames(res)[[2]] <- c("A", "B")
		  dimnames(res)[[3]] <- c("AA", "AB", "BB")
		  dimnames(res)[[4]] <- batchnames
		  if(missing(i) & missing(j)){
			  madA.AA <- as.matrix(assayDataElement(object, "madA.AA"))
			  madA.AB <- as.matrix(assayDataElement(object, "madA.AB"))
			  madA.BB <- as.matrix(assayDataElement(object, "madA.BB"))
			  madB.AA <- as.matrix(assayDataElement(object, "madB.AA"))
			  madB.AB <- as.matrix(assayDataElement(object, "madB.AB"))
			  madB.BB <- as.matrix(assayDataElement(object, "madB.BB"))
		  }
		  if(missing(i) & !missing(j)){
			  J <- match(batchnames, batchNames(object))
			  stopifnot(length(J) > 0 & !all(is.na(J)))
			  madA.AA <- as.matrix(assayDataElement(object, "madA.AA"))[, J, ...]
			  madA.AB <- as.matrix(assayDataElement(object, "madA.AB"))[, J, ...]
			  madA.BB <- as.matrix(assayDataElement(object, "madA.BB"))[, J, ...]
			  madB.AA <- as.matrix(assayDataElement(object, "madB.AA"))[, J, ...]
			  madB.AB <- as.matrix(assayDataElement(object, "madB.AB"))[, J, ...]
			  madB.BB <- as.matrix(assayDataElement(object, "madB.BB"))[, J, ...]

		  }
		  if(!missing(i) & !missing(j)){
			  J <- match(batchnames, batchNames(object))
			  stopifnot(length(J) > 0 & !all(is.na(J)))
			  madA.AA <- as.matrix(assayDataElement(object, "madA.AA"))[i, J, ...]
			  madA.AB <- as.matrix(assayDataElement(object, "madA.AB"))[i, J, ...]
			  madA.BB <- as.matrix(assayDataElement(object, "madA.BB"))[i, J, ...]
			  madB.AA <- as.matrix(assayDataElement(object, "madB.AA"))[i, J, ...]
			  madB.AB <- as.matrix(assayDataElement(object, "madB.AB"))[i, J, ...]
			  madB.BB <- as.matrix(assayDataElement(object, "madB.BB"))[i, J, ...]
		  }
		  if(!missing(i) & missing(j)){
			  madA.AA <- as.matrix(assayDataElement(object, "madA.AA"))[i, ...]
			  madA.AB <- as.matrix(assayDataElement(object, "madA.AB"))[i, ...]
			  madA.BB <- as.matrix(assayDataElement(object, "madA.BB"))[i, ...]
			  madB.AA <- as.matrix(assayDataElement(object, "madB.AA"))[i, ...]
			  madB.AB <- as.matrix(assayDataElement(object, "madB.AB"))[i, ...]
			  madB.BB <- as.matrix(assayDataElement(object, "madB.BB"))[i, ...]
		  }
		  res[, "A", "AA", ] <- madA.AA
		  res[, "A", "AB", ] <- madA.AB
		  res[, "A", "BB", ] <- madA.BB
		  res[, "B", "AA", ] <- madB.AA
		  res[, "B", "AB", ] <- madB.AB
		  res[, "B", "BB", ] <- madB.BB
		  return(res)
})



setMethod("tau2", signature(object="AssayData"),
	  function(object, allele, batchname){
		  stopifnot(!missing(allele))
		  if(missing(batchname) & length(sampleNames(object)) > 1)
			  stop("must supply batchname")
		  if(!missing(batchname)){
			  stopifnot(batchname %in% batchNames(object))
			  j <- match(batchname, sampleNames(object))
		  } else j <- 1
		  getTaus <- function(allele){
			  switch(allele,
				 A=cbind(assayDataElement(object, "tau2A.AA")[, j],
				         assayDataElement(object, "tau2A.AB")[, j],
         			         assayDataElement(object, "tau2A.BB")[, j]),
				 B=cbind(assayDataElement(object, "tau2B.AA")[, j],
				         assayDataElement(object, "tau2B.AB")[, j],
				         assayDataElement(object, "tau2B.BB")[, j]),
				 stop("allele must be 'A' or 'B'")
				 )
		  }
		  getTaus(allele)
})


