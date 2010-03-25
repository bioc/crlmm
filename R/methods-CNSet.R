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
setMethod("lM", "CNSetLM", function(object) object@lM)
setReplaceMethod("lM", c("CNSetLM", "list_or_ffdf"), function(object, value){
	object@lM <- value
	object
})

setAs("SnpSuperSet", "CNSetLM", function(from, to){
	stopifnot("batch" %in% varLabels(from))
	cnSet <- new("CNSetLM",
		     alleleA=A(from),
		     alleleB=B(from),
		     call=snpCall(from),
		     callProbability=snpCallProbability(from),
		     CA=initializeBigMatrix("CA", nrow(from), ncol(from)),
		     CB=initializeBigMatrix("CB", nrow(from), ncol(from)),
		     annotation=annotation(from),
		     featureData=featureData(from),
		     experimentData=experimentData(from),
		     phenoData=phenoData(from))
	lM(cnSet) <- initializeParamObject(list(featureNames(cnSet), unique(from$batch)))
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
	if(bias.adj & all(is.na(CA(object)))){
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


setMethod("copyNumber", "CNSet", function(object){
	I <- isSnp(object)
	CA <- CA(object)
	CB <- CB(object)
	CN <- CA + CB
	##For nonpolymorphic probes, CA is the total copy number
	CN[!I, ] <- CA(object)[!I, ]
	CN
})


setMethod("ellipse", "CNSet", function(x, copynumber, batch, ...){
	ellipse.CNSet(x, copynumber, batch, ...)
})

##setMethod("ellipse", "CNSet", function(x, copynumber, ...){
ellipse.CNSet <- function(x, copynumber, batch, ...){
	if(nrow(x) > 1) stop("only 1 snp at a time")
	##batch <- unique(x$batch)
	if(missing(batch)){
		stop("must specify batch")
	}
	if(length(batch) > 1) stop("batch variable not unique")
	nuA <- getParam(x, "nuA", batch)
	nuB <- getParam(x, "nuB", batch)
	phiA <- getParam(x, "phiA", batch)
	phiB <- getParam(x, "phiB", batch)
	tau2A <- getParam(x, "tau2A", batch)
	tau2B <- getParam(x, "tau2B", batch)
	sig2A <- getParam(x, "sig2A", batch)
	sig2B <- getParam(x, "sig2B", batch)
	corrA.BB <- getParam(x, "corrA.BB", batch)
	corrB.AA <- getParam(x, "corrB.AA", batch)
	corr <- getParam(x, "corr", batch)
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


