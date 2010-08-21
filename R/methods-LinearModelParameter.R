setMethod("nu", c("LinearModelParameter", "character"), 
	  function(object, allele){
		  getValue <- function(allele){
			  switch(allele,
				 A="nuA",
				 B="nuB",
				 stop("allele must be 'A' or 'B'"))
		  }	
		  val <- getValue(allele)
		  assayDataElement(object, val)
})



setMethod("phi", c("LinearModelParameter", "character"),
	  function(object, allele){
		  getValue <- function(allele){
			  switch(allele,
				 A="phiA",
				 B="phiB",
				 stop("allele must be 'A' or 'B'"))
		  }	
		  val <- getValue(allele)
		  assayDataElement(object, val)
	  })

setMethod("sigma2", c("LinearModelParameter", "character"),
	  function(object, allele){
		  getValue <- function(allele){
			  switch(allele,
				 A="sigma2A",
				 B="sigma2B",
				 stop("allele must be 'A' or 'B'"))
		  }
		  val <- getValue(allele)
		  assayDataElement(object, val)
	  })

setMethod("tau2", c("LinearModelParameter", "character"),
	  function(object, allele){
		  getValue <- function(allele){
			  switch(allele,
				 A="tau2A",
				 B="tau2B",
				 stop("allele must be 'A' or 'B'"))
		  }
		  val <- getValue(allele)		  
		  assayDataElement(object, val)
	  })

setMethod("corr", c("LinearModelParameter", "character"),
	  function(object, allele){
		  getValue <- function(allele){
			  switch(allele,
				 AA="corrAA",
				 AB="corrAB",
				 BB="corrBB",
				 stop("allele must be 'AA', 'AB', or 'BB'"))
		  }
		  val <- getValue(allele)		  
		  assayDataElement(object, val)
	  })
