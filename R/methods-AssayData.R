setMethod("Ns", signature(object="AssayData"),
	  function(object, genotype, batchname){
		  if(missing(batchname) & length(sampleNames(object)) > 1)
			  stop("must supply batchname")
		  if(!missing(batchname)){
			  stopifnot(batchname %in% batchNames(object))
			  j <- match(batchname, sampleNames(object))
		  } else j <- 1
		  if(missing(genotype)){
			  res <- cbind(assayDataElement(object, "N.AA")[, j],
				       assayDataElement(object, "N.AB")[, j],
				       assayDataElement(object, "N.BB")[, j])
			  return(res)
		  } else{
			  getValue <- function(genotype){
				  switch(genotype,
					 AA="N.AA",
					 AB="N.AB",
					 BB="N.BB",
					 stop("allele must be 'AA', 'AB', or 'BB'"))
			  }
			  val <- getValue(genotype)
			  return(assayDataElement(object, val)[, j])
		  }
	  })
setMethod("corr", signature(object="AssayData"),
	  function(object, genotype, batchname){
		  if(missing(batchname) & length(sampleNames(object)) > 1)
			  stop("must supply batchname")
		  if(!missing(batchname)){
			  stopifnot(batchname %in% batchNames(object))
			  j <- match(batchname, sampleNames(object))
		  } else j <- 1
		  if(missing(genotype)){
			  res <- cbind(assayDataElement(object, "corrAA")[, j],
				       assayDataElement(object, "corrAB")[, j],
				       assayDataElement(object, "corrBB")[, j])
			  return(res)
		  } else{
			  getValue <- function(genotype){
				  switch(genotype,
					 AA="corrAA",
					 AB="corrAB",
					 BB="corrBB",
					 stop("allele must be 'AA', 'AB', or 'BB'"))
			  }
			  val <- getValue(genotype)
			  return(assayDataElement(object, val)[, j])
		  }
	  })

setMethod("medians", signature(object="AssayData"),
	  function(object, allele, batchname){
		  stopifnot(!missing(allele))
		  if(missing(batchname) & length(sampleNames(object)) > 1)
			  stop("must supply batchname")
		  if(!missing(batchname)){
			  stopifnot(batchname %in% batchNames(object))
			  j <- match(batchname, sampleNames(object))
		  } else j <- 1
		  getMedians <- function(allele){
			  switch(allele,
				 A=cbind(assayDataElement(object, "medianA.AA")[, j],
				         assayDataElement(object, "medianA.AB")[, j],
         			         assayDataElement(object, "medianA.BB")[, j]),
				 B=cbind(assayDataElement(object, "medianB.AA")[, j],
				         assayDataElement(object, "medianB.AB")[, j],
				         assayDataElement(object, "medianB.BB")[, j]),
				 stop("allele must be 'A' or 'B'")
				 )
		  }
		  getMedians(allele)
})
setMethod("mads", signature(object="AssayData"),
	  function(object, allele, batchname){
		  stopifnot(!missing(allele))
		  if(missing(batchname) & length(sampleNames(object)) > 1)
			  stop("must supply batchname")
		  if(!missing(batchname)){
			  stopifnot(batchname %in% batchNames(object))
			  j <- match(batchname, sampleNames(object))
		  } else j <- 1
		  getMads <- function(allele){
			  switch(allele,
				 A=cbind(assayDataElement(object, "madA.AA")[, j],
				         assayDataElement(object, "madA.AB")[, j],
         			         assayDataElement(object, "madA.BB")[, j]),
				 B=cbind(assayDataElement(object, "madB.AA")[, j],
				         assayDataElement(object, "madB.AB")[, j],
				         assayDataElement(object, "madB.BB")[, j]),
				 stop("allele must be 'A' or 'B'")
				 )
		  }
		  getMads(allele)
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


