setMethod("getOptions", "AffymetrixAlleleSet", function(object){
	copynumber <- TRUE
	nonpolymorphic <- TRUE
	verbose <- TRUE
	seed <- 1
	snprmaOpts <- snprmaOptions.AffymetrixAlleleSet()
	crlmmOpts <- crlmmOptions.AffymetrixAlleleSet()
	cnOpts <- cnOptions()
	list(verbose=verbose,
	     copynumber=copynumber,
	     nonpolymorphic=nonpolymorphic,
	     seed=seed,
	     snprmaOpts=snprmaOpts,
	     crlmmOpts=crlmmOpts)
})

construct.AffymetrixAlleleSet <- function(object, filenames){
	protocolData <- getProtocolData(object, filenames)
	M <- .featureInfo(annotation(object))
	if(isPackageLoaded("ff")){
		M <- ff(M, dim=c(nrow(M), ncol(M)),
			vmode="integer", finalizer="close",
			overwrite=TRUE,
			dimnames=list(rownames(M), c("chromosome", "position", "isSnp")))
	} 
	dns <- list(rownames(M), basename(filenames))
	nr <- nrow(M)
	alleleSet <- new("AffymetrixAlleleSet", 
			 alleleA=initializeBigMatrix(dns),
			 alleleB=initializeBigMatrix(dns),
			 genomeAnnotation=M,
			 options=crlmmOptions(object),
			 annotation=annotation(object))
	protocolData(alleleSet) <- protocolData
	sampleNames(alleleSet) <- basename(filenames)
	featureNames(alleleSet) <- dns[[1]]
	return(alleleSet)
}
setMethod("construct", "AffymetrixAlleleSet", function(object, filenames){
	construct.AffymetrixAlleleSet(object, filenames)
})

snprmaOptions.AffymetrixAlleleSet <- function(mixtureSampleSize=10^5,
					      fitMixture=TRUE,
					      eps=0.1,
					      seed=1){
	list(mixtureSampleSize=mixtureSampleSize,
	     fitMixture=fitMixture,
	     eps=eps,
	     seed=seed)
}

crlmmOptions.AffymetrixAlleleSet <- function(batchSize=1000,
					     probs=rep(1/3,3),
					     DF=6,
					     SNRMin=5,
					     gender=NULL,
					     recallMin=10,
					     recallRegMin=1000,
					     mixtureParams=NA,
					     badSNP=0.7){
	list(batchSize=batchSize,
	     probs=probs,
	     DF=DF,
	     SNRMin=SNRMin,
	     gender=gender,
	     recallMin=recallMin,
	     recallRegMin=recallRegMin,
	     mixtureParams=mixtureParams,
	     badSNP=badSNP)
}

cnOptions <- function(MIN.OBS=3,
		      MIN.SAMPLES=10,
		      DF.PRIOR=50,
		      bias.adj=FALSE,
		      prior.prob=rep(1/4, 4),
		      SNRMin=4,
		      chromosome=1:23,
		      GT.CONF.THR=0.99,
		      PHI.THR=2^6,##used in nonpolymorphic fxn for training
		      nHOM.THR=5, ##used in nonpolymorphic fxn for training
		      MIN.NU=2^3,
		      MIN.PHI=2^3,
		      THR.NU.PHI=TRUE,
		      thresholdCopynumber=TRUE,
		      use.poe=FALSE,
		      whichBatch=NULL){
	list(MIN.OBS=MIN.OBS,
	     MIN.SAMPLES=MIN.SAMPLES,
	     DF.PRIOR=DF.PRIOR,
	     bias.adj=bias.adj,
	     prior.prob=prior.prob,
	     SNRMin=SNRMin,
	     chromosome=chromosome,
	     GT.CONF.THR=GT.CONF.THR,
	     PHI.THR=PHI.THR,##used in nonpolymorphic fxn for training
	     nHOM.THR=nHOM.THR, ##used in nonpolymorphic fxn for training
	     MIN.NU=MIN.NU,
	     MIN.PHI=MIN.PHI,
	     THR.NU.PHI=THR.NU.PHI,
	     thresholdCopynumber=thresholdCopynumber,
	     use.poe=use.poe,
	     whichBatch=whichBatch)
}

setMethod("getProtocolData", "AffymetrixAlleleSet", function(object, filenames){
	getProtocolData.AffymetrixAlleleSet(filenames)
})

getProtocolData.AffymetrixAlleleSet <- function(filenames){
	scanDates <- data.frame(ScanDate=sapply(filenames, celfileDate))
	rownames(scanDates) <- basename(rownames(scanDates))
	protocoldata <- new("AnnotatedDataFrame",
			    data=scanDates,
			    varMetadata=data.frame(labelDescription=colnames(scanDates),
			                           row.names=colnames(scanDates)))
	return(protocoldata)
}

##setAs("AffymetrixAlleleSet", "AffymetrixCallSet", function(from, to){
##	calls <- getCalls(crlmmOptions(from))
##	dimnames(calls) <- list(featureNames(from), sampleNames(from))
##	confs <- getConfs(crlmmOptions(from))
##	dimnames(confs) <- list(featureNames(from), sampleNames(from))
##	new("AffymetrixCallSet",
##	    alleleA=A(from),
##	    alleleB=B(from),		  
##	    call=calls,
##	    callProbability=confs,
##	    phenoData=phenoData(from),
##	    experimentData=experimentData(from),
##	    annotation=annotation(from),
##	    protocolData=protocolData(from),
##	    featureData=featureData(from),
##	    genomeAnnotation=genomeAnnotation(from),
##	    options=crlmmOptions(from))
##})



##
##setMethod("nFeatures", "AffymetrixAlleleSet", function(object){
##	cdfname <- annotation(object)
##	getdimensions <- function(cdfname){
##		switch(cdfname, 
##		       genomewidesnp6=list(snpRange=c(1, 906600), npRange=c(906601, 1852406), nfeatures=1852406),
##		       genomewidesnp5=list(snpRange=c(1, 500568), npRange=c(500569, 917837), nfeatures=917837))
##	}
##	info <- getdimensions(cdfname)
##})
##
##setMethod("getGenomeAnnotation", "AffymetrixAlleleSet", function(object){
##	genomeAnnotation(object) <- .featureInfo(annotation(object))
##	featureNames(object) <- rownames(genomeAnnotation(object))
##	return(object)
##
##})
##setMethod("getGenomeAnnotation", "AffymetrixBigData", function(object){
##	M <- .featureInfo(annotation(object))
##	for(j in 1:3){
##		genomeAnnotation(object)[, j] <- M[, j]
##	}
##	rownames(genomeAnnotation(object)) <- rownames(M)
##	featureNames(object) <- rownames(M)
##	return(object)
##})


.featureInfo <- function(cdfName){
	pkgname <- getCrlmmAnnotationName(cdfName)
	if(!require(pkgname, character.only=TRUE)){
		suggCall <- paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep="")
		msg <- paste("If", pkgname, "is installed on an alternative location, please load it manually by using", suggCall)
		message(strwrap(msg))
		stop("Package ", pkgname, " could not be found.")
		rm(suggCall, msg)
	}
	loader("preprocStuff.rda", .crlmmPkgEnv, pkgname)
	loader("genotypeStuff.rda", .crlmmPkgEnv, pkgname)
	loader("mixtureStuff.rda", .crlmmPkgEnv, pkgname)
	gns <- getVarInEnv("gns")
	path <- system.file("extdata", package=paste(cdfName, "Crlmm", sep=""))
	load(file.path(path, "snpProbes.rda"))
	load(file.path(path, "cnProbes.rda"))
	cnProbes <- get("cnProbes")
	snpIndex <- seq(along=gns)
	npIndex <- seq(along=rownames(cnProbes)) + max(snpIndex) 
	featurenames <- c(gns, rownames(cnProbes))

	fvarlabels=c("chromosome", "position", "isSnp")
	M <- matrix(NA, length(featurenames), 3, dimnames=list(featurenames, fvarlabels))
	index <- match(rownames(snpProbes), rownames(M)) #only snp probes in M get assigned position
	M[index, "position"] <- snpProbes[, grep("pos", colnames(snpProbes))]
	M[index, "chromosome"] <- snpProbes[, grep("chr", colnames(snpProbes))]
	M[index, "isSnp"] <- 1L
	index <- which(is.na(M[, "isSnp"]))
	M[index, "isSnp"] <- 1L
	
	index <- match(rownames(cnProbes), rownames(M)) #only snp probes in M get assigned position
	M[index, "position"] <- cnProbes[, grep("pos", colnames(cnProbes))]
	M[index, "chromosome"] <- cnProbes[, grep("chr", colnames(cnProbes))]
	M[index, "isSnp"] <- 0L
	return(M)
	##list(snpIndex, npIndex, fns)
	##crlmmOpts$snpRange <- range(snpIndex)
	##crlmmOpts$npRange <- range(npIndex)
}
