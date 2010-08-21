setMethod("show", "CNSetLM", function(object){
	callNextMethod(object)
	cat("lM: ", length(lM(object)), " elements \n")
	print(names(lM(object)))
})

setMethod("[", "CNSetLM", function(x, i, j, ..., drop=FALSE){
	x <- callNextMethod(x, i, j, ..., drop=drop)
	if(!missing(i)){
		if(class(lM(x)) == "ffdf"){
			lM(x) <- lapply(physical(lM(x)), function(x, i){open(x); x[i, ]}, i=i)
		} else {
			lM(x) <- lapply(lM(x), function(x, i) x[i, , drop=FALSE], i=i)
		}
	}
	x
})

setMethod("[", "CNSetLM", function(x, i, j, ..., drop=FALSE){
	x <- callNextMethod(x, i, j, ..., drop=drop)
##	if(!missing(i)){
##		if(class(lM(x)) == "ffdf"){
##			lM(x) <- lapply(physical(lM(x)), function(x, i){open(x); x[i, ]}, i=i)
##		} else {
##			lM(x) <- lapply(lM(x), function(x, i) x[i, , drop=FALSE], i=i)
##		}
##	}
	x
})


setMethod("lM", "CNSetLM", function(object) object@lM)
setReplaceMethod("lM", c("CNSetLM", "list_or_ffdf"), function(object, value){
	object@lM <- value
	object
})



setMethod("open", "CNSetLM", function(con,...){
	callNextMethod(con,...)
	physical <- get("physical")
	lapply(physical(lM(con)), open)
})

setAs("SnpSuperSet", "CNSetLM", function(from, to){
	stopifnot("batch" %in% varLabels(protocolData(from)))
	cnSet <- new("CNSetLM",
		     alleleA=A(from),
		     alleleB=B(from),
		     call=snpCall(from),
		     callProbability=snpCallProbability(from),
##		     CA=initializeBigMatrix("CA", nrow(from), ncol(from)),
##		     CB=initializeBigMatrix("CB", nrow(from), ncol(from)),
		     annotation=annotation(from),
		     featureData=featureData(from),
		     experimentData=experimentData(from),
		     protocolData=protocolData(from),
		     phenoData=phenoData(from))
	lM(cnSet) <- initializeParamObject(list(featureNames(cnSet), unique(protocolData(from)$batch)))
	return(cnSet)
})

setMethod("computeCopynumber", "CNSet",
	  function(object,
		   MIN.OBS,
		   DF.PRIOR,
		   bias.adj,
		   prior.prob,
		   seed,
		   verbose,
		   GT.CONF.THR,
		   PHI.THR,
		   nHOM.THR,
		   MIN.NU,
		   MIN.PHI,
		   THR.NU.PHI,
		   thresholdCopynumber){
	## to do the bias adjustment, initial estimates of the parameters are needed
	##  The initial estimates are gotten by running computeCopynumber with cnOptions[["bias.adj"]]=FALSE
		  cnOptions <- list(
				    MIN.OBS=MIN.OBS,
				    DF.PRIOR=DF.PRIOR,
				    bias.adj=bias.adj,
				    prior.prob=prior.prob,
				    seed=seed,
				    verbose=verbose,
				    GT.CONF.THR=GT.CONF.THR,
				    PHI.THR=PHI.THR,
				    nHOM.THR=nHOM.THR,
				    MIN.NU=MIN.NU,
				    MIN.PHI=MIN.PHI,
				    THR.NU.PHI=THR.NU.PHI,
				    thresholdCopynumber=thresholdCopynumber)
	bias.adj <- cnOptions[["bias.adj"]]
	if(bias.adj & all(is.na(nu(object, "A")[, 1])){
		cnOptions[["bias.adj"]] <- FALSE
	}
	object <- computeCopynumber.CNSet(object, cnOptions)				
	if(bias.adj & !cnOptions[["bias.adj"]]){
		## Do a second iteration with bias adjustment
		cnOptions[["bias.adj"]] <- TRUE
		object <- computeCopynumber.CNSet(object, cnOptions)
	}
	object
})

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

setMethod("ellipse", "CNSet", function(x, copynumber, batch, ...){
	ellipse.CNSet(x, copynumber, batch, ...)
})

setMethod("nu", c("CNSetLM", "character"), function(object, allele){
	getValue <- function(allele){
		switch(allele,
		       A="nuA",
		       B="nuB",
		       stop("allele must be 'A' or 'B'"))
	}	
	val <- getValue(allele)
	class.lm <- class(lM(object)) 
	if(class.lm == "ffdf"){
		physical <- get("physical")
		res <- physical(lM(object))[[val]]

	} else {
		if(class.lm != "list") stop("lM() must be matrix or ffdf")
		res <- lM(object)[[val]]
	}
	return(res)
})

setMethod("phi", c("CNSetLM", "character"), function(object, allele){
	getValue <- function(allele){
		switch(allele,
		       A="phiA",
		       B="phiB",
		       stop("allele must be 'A' or 'B'"))
	}
	val <- getValue(allele)	
	class.lm <- class(lM(object)) 
	if(class.lm == "ffdf"){
		physical <- get("physical")
		res <- physical(lM(object))[[val]]

	} else {
		if(class.lm != "list") stop("lM() must be matrix or ffdf")
		res <- lM(object)[[val]]
	}
	return(res)
})

setMethod("sigma2", c("CNSetLM", "character"), function(object, allele){
	getValue <- function(allele){
		switch(allele,
		       A="sig2A",
		       B="sig2B",
		       stop("allele must be 'A' or 'B'"))
	}
	val <- getValue(allele)	
	class.lm <- class(lM(object))
	if(class.lm == "ffdf"){
		physical <- get("physical")
		res <- physical(lM(object))[[val]]

	} else {
		if(class.lm != "list") stop("lM() must be matrix or ffdf")
		res <- lM(object)[[val]]
	}
	return(res)
})

setMethod("tau2", c("CNSetLM", "character"), function(object, allele){
	getValue <- function(allele){
		switch(allele,
		       A="tau2A",
		       B="tau2B",
		       stop("allele must be 'A' or 'B'"))
	}
	val <- getValue(allele)
	class.lm <- class(lM(object))
	if(class.lm == "ffdf"){
		physical <- get("physical")
		res <- physical(lM(object))[[val]]

	} else {
		if(class.lm != "list") stop("lM() must be matrix or ffdf")
		res <- lM(object)[[val]]
	}
	return(res)
})

setMethod("corr", c("CNSetLM", "character"), function(object, allele){
	getValue <- function(allele){
		switch(allele,
		       AA="corrAA",
		       AB="corrAB",
		       BB="corrBB",
		       stop("must be AA, AB, or BB"))
	}
	val <- getValue(allele)
	class.lm <- class(lM(object))
	if(class.lm == "ffdf"){
		physical <- get("physical")
		res <- physical(lM(object))[[val]]

	} else {
		if(class.lm != "list") stop("lM() must be matrix or ffdf")
		res <- lM(object)[[val]]
	}
	return(res)
})

