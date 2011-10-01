setMethod("Ns", signature(object="AssayData"),
	  function(object, ...){
		  dotArgs <- names(list(...))
		  missing.i <- !("i" %in% dotArgs)
		  missing.j <- !("j" %in% dotArgs)
		  if(!missing.j) j <- list(...)[["j"]]
		  batchnames <- sampleNames(object)
		  if(!missing.j) batchnames <- batchnames[j]
		  if(missing.i & missing.j) stop("Must specify either the rows i or batches j")
		  is.ff <- is(assayDataElement(object, "N.AA"), "ff") | is(assayDataElement(object, "N.AA"), "ffdf")
		  if(is.ff){
			  open(assayDataElement(object, "N.AA"))
			  open(assayDataElement(object, "N.AB"))
			  open(assayDataElement(object, "N.BB"))
		  }
		  N.AA <- as.matrix(assayDataElement(object, "N.AA")[...])
		  N.AB <- as.matrix(assayDataElement(object, "N.AB")[...])
		  N.BB <- as.matrix(assayDataElement(object, "N.BB")[...])
		  if(is.ff){
			  close(assayDataElement(object, "N.AA"))
			  close(assayDataElement(object, "N.AB"))
			  close(assayDataElement(object, "N.BB"))
		  }
		  res <- array(NA, dim=c(nrow(N.AA), 3, ncol(N.AA)))
		  dimnames(res)[[2]] <- c("AA", "AB", "BB")
		  dimnames(res)[[3]] <- batchnames
		  res[, "AA", ] <- N.AA
		  res[, "AB", ] <- N.AB
		  res[, "BB", ] <- N.BB
		  return(res)
	  })
setMethod("corr", signature(object="AssayData"),
	  function(object, ...){
		  dotArgs <- names(list(...))
		  missing.i <- !("i" %in% dotArgs)
		  missing.j <- !("j" %in% dotArgs)
		  if(!missing.j) j <- list(...)[["j"]]
		  batchnames <- sampleNames(object)
		  if(!missing.j) batchnames <- batchnames[j]

		  ##if(missing.i & missing.j) stop("Must specify either the rows i or batches j")
		  is.ff <- is(assayDataElement(object, "corrAA"), "ff") | is(assayDataElement(object, "corrAA"), "ffdf")
		  if(is.ff){
			  open(assayDataElement(object, "corrAA"))
			  open(assayDataElement(object, "corrAB"))
			  open(assayDataElement(object, "corrBB"))
		  }
		  corrAA <- as.matrix(assayDataElement(object, "corrAA")[...])
		  corrAB <- as.matrix(assayDataElement(object, "corrAB")[...])
		  corrBB <- as.matrix(assayDataElement(object, "corrBB")[...])
		  if(is.ff){
			  close(assayDataElement(object, "corrAA"))
			  close(assayDataElement(object, "corrAB"))
			  close(assayDataElement(object, "corrBB"))
		  }
		  res <- array(NA, dim=c(nrow(corrAA), 3, ncol(corrAA)))
		  dimnames(res)[[2]] <- c("AA", "AB", "BB")
		  dimnames(res)[[3]] <- batchnames
		  res[, "AA", ] <- corrAA
		  res[, "AB", ] <- corrAB
		  res[, "BB", ] <- corrBB
		  return(res)
	  })

setMethod("medians", signature(object="AssayData"),
	  function(object, ...){
		  dotArgs <- names(list(...))
		  missing.i <- !("i" %in% dotArgs)
		  missing.j <- !("j" %in% dotArgs)
		  if(!missing.j) j <- list(...)[["j"]]
		  batchnames <- sampleNames(object)
		  if(!missing.j) batchnames <- batchnames[j]
		  if(missing.i & missing.j) stop("Must specify either the rows i or batches j")
		  medianA.AA <- as.matrix(assayDataElement(object, "medianA.AA")[...])
		  medianA.AB <- as.matrix(assayDataElement(object, "medianA.AB")[...])
		  medianA.BB <- as.matrix(assayDataElement(object, "medianA.BB")[...])
		  medianB.AA <- as.matrix(assayDataElement(object, "medianB.AA")[...])
		  medianB.AB <- as.matrix(assayDataElement(object, "medianB.AB")[...])
		  medianB.BB <- as.matrix(assayDataElement(object, "medianB.BB")[...])
		  res <- array(NA, dim=c(nrow(medianA.AA), 2, 3, ncol(medianA.AA)))
		  dimnames(res)[[2]] <- c("A", "B")
		  dimnames(res)[[3]] <- c("AA", "AB", "BB")
		  dimnames(res)[[4]] <- batchnames
		  res[, "A", "AA", ] <- medianA.AA
		  res[, "A", "AB", ] <- medianA.AB
		  res[, "A", "BB", ] <- medianA.BB
		  res[, "B", "AA", ] <- medianB.AA
		  res[, "B", "AB", ] <- medianB.AB
		  res[, "B", "BB", ] <- medianB.BB
		  return(res)
})

setMethod("mads", signature(object="AssayData"),
	  function(object, ...){
		  dotArgs <- names(list(...))
		  missing.i <- !("i" %in% dotArgs)
		  missing.j <- !("j" %in% dotArgs)
		  batchnames <- sampleNames(object)
		  if(!missing.j) j <- list(...)[["j"]]
		  if(!missing.j) batchnames <- batchnames[j]
		  if(missing.i & missing.j) stop("Must specify either the rows i or batches j")
		  madA.AA <- as.matrix(assayDataElement(object, "madA.AA")[...])
		  madA.AB <- as.matrix(assayDataElement(object, "madA.AB")[...])
		  madA.BB <- as.matrix(assayDataElement(object, "madA.BB")[...])
		  madB.AA <- as.matrix(assayDataElement(object, "madB.AA")[...])
		  madB.AB <- as.matrix(assayDataElement(object, "madB.AB")[...])
		  madB.BB <- as.matrix(assayDataElement(object, "madB.BB")[...])
		  res <- array(NA, dim=c(nrow(madA.AA), 2, 3, ncol(madA.AA)))
		  dimnames(res)[[2]] <- c("A", "B")
		  dimnames(res)[[3]] <- c("AA", "AB", "BB")
		  dimnames(res)[[4]] <- batchnames
		  res[, "A", "AA", ] <- madA.AA
		  res[, "A", "AB", ] <- madA.AB
		  res[, "A", "BB", ] <- madA.BB
		  res[, "B", "AA", ] <- madB.AA
		  res[, "B", "AB", ] <- madB.AB
		  res[, "B", "BB", ] <- madB.BB
		  return(res)
})



setMethod("tau2", signature(object="AssayData"),
	  function(object, ...){
		  dotArgs <- names(list(...))
		  missing.i <- !("i" %in% dotArgs)
		  missing.j <- !("j" %in% dotArgs)
		  batchnames <- sampleNames(object)
		  if(!missing.j) j <- list(...)[["j"]]
		  if(!missing.j) batchnames <- batchnames[j]
		  if(missing.i & missing.j) stop("Must specify either the rows i or batches j")
		  tau2A.AA <- as.matrix(assayDataElement(object, "tau2A.AA")[...])
		  tau2A.BB <- as.matrix(assayDataElement(object, "tau2A.BB")[...])
		  tau2B.AA <- as.matrix(assayDataElement(object, "tau2B.AA")[...])
		  tau2B.BB <- as.matrix(assayDataElement(object, "tau2B.BB")[...])
		  res <- array(NA, dim=c(nrow(tau2A.AA), 2, 2, ncol(tau2A.AA)))
		  dimnames(res)[[2]] <- c("A", "B")
		  dimnames(res)[[3]] <- c("AA", "BB")
		  dimnames(res)[[4]] <- batchnames
		  res[, "A", "AA", ] <- tau2A.AA
		  res[, "A", "BB", ] <- tau2A.BB
		  res[, "B", "AA", ] <- tau2B.AA
		  res[, "B", "BB", ] <- tau2B.BB
		  return(res)
})


