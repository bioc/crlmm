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
setReplaceMethod("lM", c("CNSetLM", "ffdf"), function(object, value){
	object@lM <- value
	object
})
setReplaceMethod("lM", c("CNSetLM", "list"), function(object, value){
	object@lM <- value
	object
})



setMethod("open", "CNSetLM", function(con,...){
	callNextMethod(con,...)
	lapply(physical(lM(con)), open)
})

setAs("SnpSuperSet", "CNSetLM", function(from, to){
	stopifnot("batch" %in% varLabels(protocolData(from)))
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
	message("This accessor will be deprecated.")
	I <- isSnp(object)
	ffIsLoaded <- inherits(CA(object), "ff")
	CA <- CA(object)
	CB <- CB(object)
	if(ffIsLoaded){
		open(CA)
		open(CB)
		CA <- as.matrix(CA[,])
		CB <- as.matrix(CB[,])
	}
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


setMethod("totalCopyNumber",
	  signature=signature(object="CNSet"),
	  function(object, ...){
		  is.ff <- is(CA(object), "ff") | is(CA(object), "ffdf")
		  dotArgs <- list(...)
		  missing.i <- !"i" %in% names(dotArgs)
		  missing.j <- !"j" %in% names(dotArgs)
		  if(missing.i & missing.j){
			  if(is.ff) stop("Must specify i and/or j for ff objects")
		  }
		  if(!missing.i) {
			  i <- dotArgs[["i"]]
			  snp.index <- intersect(i, which(isSnp(object)))
			  ##which rows in the return matrix are snps
			  snp.index2 <- which(isSnp(object)[i])
		  } else {
			  i <- 1:nrow(object)
			  snp.index <- which(isSnp(object))
			  snp.index2 <- snp.index
		  }
		  if(!missing.j){
			  j <- dotArgs[["j"]]
		  } else j <- 1:ncol(object)
		  cn.total <- as.matrix(CA(object)[i, j])
		  if(length(snp.index) > 0){
			  cb <- as.matrix(CB(object)[snp.index, j])
			  cn.total[snp.index2, ] <- cn.total[snp.index2, ] + cb
		  }
		  cn.total <- cn.total/100
		  return(cn.total)
##		  if(missing.i){
####		  if(missing.i & !missing.j){
##			 snp.index <- which(isSnp(object))
##			 cn.total <- as.matrix(CA(object)[, j])
##			 if(length(snp.index) > 0){
##				 cb <- as.matrix(CB(object)[snp.index, j])
####				 snps <- (1:nrow(cn.total))[i %in% snp.index]
##				 cn.total[snp.index, ] <- cn.total[snp.index, ] + cb
##			 }
##		 } else{
##			 snp.index <- intersect(which(isSnp(object)), i)
##			 cn.total <- as.matrix(CA(object)[i, ])
##			 if(length(snp.index) > 0){
##				 cb <- as.matrix(CB(object)[snp.index, ])
##				 snps <- (1:nrow(cn.total))[i %in% snp.index]
##				 cn.total[snps, ] <- cn.total[snps, ] + cb
##			 }
##		 }
##		 if(!missing(i) & missing(j)){
##			 snp.index <- intersect(which(isSnp(object)), i)
##			 cn.total <- as.matrix(CA(object)[i, ])
##			 if(length(snp.index) > 0){
##				 cb <- as.matrix(CB(object)[snp.index, ])
##				 snps <- (1:nrow(cn.total))[i %in% snp.index]
##				 cn.total[snps, ] <- cn.total[snps, ] + cb
##			 }
##		 }
##		 if(!missing(i) & !missing(j)){
##			 snp.index <- intersect(which(isSnp(object)), i)
##			 cn.total <- as.matrix(CA(object)[i, j])
##			 if(length(snp.index) > 0){
##				 cb <- as.matrix(CB(object)[snp.index, j])
##				 snps <- (1:nrow(cn.total))[i %in% snp.index]
##				 cn.total[snps, ] <- cn.total[snps, ] + cb
##			 }
##		 }
##		 cn.total <- cn.total/100
##		 dimnames(cn.total) <- NULL
##		 return(cn.total)
	 })
