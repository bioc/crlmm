rowCovs <- function(x, y, ...){
	notna <- !is.na(x)
	N <- rowSums(notna)
	x <- suppressWarnings(log2(x))
	if(!missing(y)) y <- suppressWarnings(log2(y))
	return(rowSums((x - rowMeans(x, ...)) * (y - rowMeans(y, ...)), ...)/(N-1))
}

rowMAD <- function(x, y, ...){
	notna <- !is.na(x)
	mad <- 1.4826*rowMedians(abs(x-rowMedians(x, ...)), ...)
	return(mad)
}

rowCors <- function(x, y, ...){
	N <- rowSums(!is.na(x))
	x <- suppressWarnings(log2(x))
	y <- suppressWarnings(log2(y))
	sd.x <- rowSds(x, ...)
	sd.y <- rowSds(y, ...)
	covar <- rowSums((x - rowMeans(x, ...)) * (y - rowMeans(y, ...)), ...)/(N-1)
	return(covar/(sd.x*sd.y))
}

generateX <- function(w, X) as.numeric(diag(w) %*% X)
generateIXTX <- function(x, nrow=3) {
	X <- matrix(x, nrow=nrow)
	XTX <- crossprod(X)
	solve(XTX)
}


nuphiAllele <- function(object, allele, Ystar, W, tmp.objects, cnOptions){
	I <- isSnp(object)
	Ystar <- Ystar[I, ]
	rownames(Ystar) <- featureNames(object)[isSnp(object)]
	complete <- rowSums(is.na(W)) == 0 & I
	W <- W[I, ]
	if(any(!is.finite(W))){## | any(!is.finite(V))){
		i <- which(rowSums(!is.finite(W)) > 0)
		##browser()
		stop("Possible zeros in the within-genotype estimates of the spread (vA, vB). ")
	}	
	Ns <- tmp.objects[["Ns"]]	
	Ns <- Ns[I, ]
	
	CHR <- unique(chromosome(object))
	batch <- unique(object$batch)
	nuA <- getParam(object, "nuA", batch)
	nuB <- getParam(object, "nuB", batch)
	nuA.se <- getParam(object, "nuA.se", batch)
	nuB.se <- getParam(object, "nuB.se", batch)
	phiA <- getParam(object, "phiA", batch)
	phiB <- getParam(object, "phiB", batch)
	phiA.se <- getParam(object, "phiA.se", batch)
	phiB.se <- getParam(object, "phiB.se", batch)	
	if(CHR==23){
		phiAX <- getParam(object, "phiAX", batch)
		phiBX <- getParam(object, "phiBX", batch)		
	}
	NOHET <- mean(Ns[, "AB"], na.rm=TRUE) < 0.05
	if(missing(allele)) stop("must specify allele")
	if(CHR == 23){
		##Design matrix for X chromosome depends on whether there was a sufficient number of AB genotypes
		if(length(grep("AB", colnames(W))) > 0){
			if(allele == "A") X <- cbind(1, c(1, 0, 2, 1, 0), c(0, 1, 0, 1, 2))
			if(allele == "B") X <- cbind(1, c(0, 1, 0, 1, 2), c(1, 0, 2, 1, 0))
		} else{
			if(allele == "A") X <- cbind(1, c(1, 0, 2, 0), c(0, 1, 0, 2))
			if(allele == "B") X <- cbind(1, c(0, 1, 0, 2), c(1, 0, 2, 0))
		}
	} else {##autosome
		if(allele == "A") X <- cbind(1, 2:0) else X <- cbind(1, 0:2)
		if(NOHET) X <- X[-2, ] ##more than 1 X chromosome, but all homozygous		
	}
	##How to quickly generate Xstar, Xstar = diag(W) %*% X
	Xstar <- apply(W, 1, generateX, X)
	IXTX <- apply(Xstar, 2, generateIXTX, nrow=nrow(X))
	##as.numeric(diag(W[1, ]) %*% X)
	if(CHR == 23){
		betahat <- matrix(NA, 3, nrow(Ystar))
		ses <- matrix(NA, 3, nrow(Ystar))		
	} else{
		betahat <- matrix(NA, 2, nrow(Ystar))
		ses <- matrix(NA, 2, nrow(Ystar))
	}
	for(i in 1:nrow(Ystar)){
		betahat[, i] <- crossprod(matrix(IXTX[, i], ncol(X), ncol(X)), crossprod(matrix(Xstar[, i], nrow=nrow(X)), Ystar[i, ]))
		ssr <- sum((Ystar[i, ] - matrix(Xstar[, i], nrow(X), ncol(X)) %*% matrix(betahat[, i], ncol(X), 1))^2)
		ses[, i] <- sqrt(diag(matrix(IXTX[, i], ncol(X), ncol(X)) * ssr))
	}
	if(allele == "A"){
		nuA[complete] <- betahat[1, ]
		phiA[complete] <- betahat[2, ]
		nuA.se[complete] <- ses[1, ]
		phiA.se[complete] <- ses[2, ]
		object <- pr(object, "nuA", batch, nuA)
		object <- pr(object, "nuA.se", batch, nuA.se)
		object <- pr(object, "phiA", batch, phiA)
		object <- pr(object, "phiA.se", batch, phiA.se)
		if(CHR == 23){
			phiAX[complete] <- betahat[3, ]
			object <- pr(object, "phiAX", batch, phiAX)			
		}
	}
	if(allele=="B"){
		nuB[complete] <- betahat[1, ]
		phiB[complete] <- betahat[2, ]
		nuB.se[complete] <- ses[1, ]
		phiB.se[complete] <- ses[2, ]
		object <- pr(object, "nuB", batch, nuB)
		object <- pr(object, "nuB.se", batch, nuB.se)
		object <- pr(object, "phiB", batch, phiB)		
		object <- pr(object, "phiB.se", batch, phiB.se)
		if(CHR == 23){
			phiBX[complete] <- betahat[3, ]
			object <- pr(object, "phiBX", batch, phiBX)
		}
	}
##	THR.NU.PHI <- cnOptions$THR.NU.PHI
##	if(THR.NU.PHI){
##		verbose <- cnOptions$verbose
##		if(verbose) message("Thresholding nu and phi")
##		object <- thresholdModelParams(object, cnOptions)
##	}
	return(object)
}



predictGender <- function(res, cdfName="genomewidesnp6", SNRMin=5){
	pkgname <- paste(cdfName, "Crlmm", sep="")
	path <- system.file("extdata", package=pkgname)
        loader("23index.rda", .crlmmPkgEnv, pkgname)
	index <- getVarInEnv("index")
	##load(file.path(path, "23index.rda"))
	XIndex <- index[[1]]	
	if(class(res) != "ABset"){
		A <- res$A
		B <- res$B
	} else{
		A <- res@assayData[["A"]]
		B <- res@assayData[["B"]]
	}
	tmp <- which(rowSums(is.na(A)) > 0)
	XMedian <- apply(log2(A[XIndex,, drop=FALSE])+log2(B[XIndex,, drop=FALSE]), 2, median, na.rm=TRUE)/2
	SNR <- res$SNR
	if(sum(SNR>SNRMin)==1){
		gender <- which.min(c(abs(XMedian-8.9), abs(XMedian-9.5)))
	}else{
		gender <- kmeans(XMedian, c(min(XMedian[SNR>SNRMin]), max(XMedian[SNR>SNRMin])))[["cluster"]]
	}
	return(gender)
}


##initializeFFObjects <- function(filenames, cnOptions){
##	outdir <- cnOptions[["outdir"]]
##	cdfName <- cnOptions[["cdfName"]]
##	AFile <- cnOptions[["AFile"]]
##	BFile <- cnOptions[["BFile"]]
##	callsFile <- cnOptions[["callsFile"]]
##	confsFile <- cnOptions[["confsFile"]]
##	snprmaFile <- cnOptions[["snprmaFile"]]
##	cnrmaFile <- cnOptions[["cnrmaFile"]]
##	CAFile <- cnOptions[["CAFile"]]
##	CBFile <- cnOptions[["CBFile"]]
##	load.it <- cnOptions[["load.it"]]
##	fileExists <- list(A=file.exists(AFile),
##			   B=file.exists(BFile),
##			   calls=file.exists(callsFile),
##			   confs=file.exists(confsFile),
##			   CA=file.exists(CAFile),
##			   CB=file.exists(CBFile))
##	allExists <- all(unlist(fileExists))
##	##if files already exist, check that the files have the appropriate dimension
##	if(allExists){
##		load(AFile)
##		open(A)
##		sns <- dimnames(A)[[2]]
##		if(!identical(sns, basename(filenames)) | !load.it){
##		## if not of the same dimension, clean up
##			message("Sample names in previously saved objects differ from the filenames. Removing previously saved objects.")
##			delete(A); gc()
##			unlink(AFile)
##			load(BFile); delete(B); unlink(BFile)
##			unlink(snprmaFile)
##			unlink(cnrmaFile)
##			if(file.exists(file.path(outdir, "cnParams.rda"))){
##				load(file.path(outdir, "cnParams.rda"))
##				delete(cnParams); gc()
##				unlink(file.path(outdir, "cnParams.rda"))
##			}
##			load(callsFile); delete(calls); unlink(callsFile)
##			load(confsFile); delete(confs); unlink(confsFile)
##			load(CAFile); delete(CA); unlink(CAFile)
##			load(CBFile); delete(CB); unlink(CBFile)
##			allExists <- FALSE
##		} 
##	}
##	if(!allExists) 	{
##		message("Initializing ff objects for A, B, confs, calls, CA, and CB.")
##		dns <- .dimnames(filenames, cnOptions[["cdfName"]], cnOptions[["verbose"]])
##		fns <- dns[[1]]
##	}
##	if(!file.exists(AFile)) {A <- initializeBigMatrix(dns); save(A, file=AFile); close(A)}
##	if(!file.exists(BFile)) {B <- initializeBigMatrix(dns); save(B, file=BFile); close(B)}
##	if(!file.exists(confsFile)) {confs <- initializeBigMatrix(dns); save(confs, file=confsFile); close(confs)}
##	if(!file.exists(callsFile)) {calls <- initializeBigMatrix(dns); save(calls, file=callsFile); close(calls)}
##	if(!file.exists(CAFile)) {CA <- initializeBigMatrix(dns); save(CA, file=CAFile); close(CA)}
##	if(!file.exists(CBFile)) {CB <- initializeBigMatrix(dns); save(CB, file=CBFile); close(CB)}
##	featureDataFile <- file.path(outdir, "featureDataFF.rda")
##	if(!file.exists(featureDataFile)){
##		path <- system.file("extdata", package=paste(cnOptions[["cdfName"]], "Crlmm", sep=""))
##		load(file.path(path, "snpProbes.rda"))
##		snpProbes <- get("snpProbes")
##		load(file.path(path, "cnProbes.rda"))
##		cnProbes <- get("cnProbes")
##		message("Initializing featureDataFF.")
##		fvarlabels <- c("chromosome", "position",   "isSnp")
##		M <- matrix(NA, length(fns), 3, dimnames=list(fns, fvarlabels))
##		index <- match(rownames(snpProbes), rownames(M)) #only snp probes in M get assigned position
##		M[index, "position"] <- snpProbes[, grep("pos", colnames(snpProbes))]
##		M[index, "chromosome"] <- snpProbes[, grep("chr", colnames(snpProbes))]
##		M[index, "isSnp"] <- 1L		
##		index <- match(rownames(cnProbes), rownames(M)) #only snp probes in M get assigned position
##		M[index, "position"] <- cnProbes[, grep("pos", colnames(cnProbes))]
##		M[index, "chromosome"] <- cnProbes[, grep("chr", colnames(cnProbes))]
##		M[index, "isSnp"] <- 0L
##		featureDataFF <- ff(M, dim=c(nrow(M), ncol(M)),
##				    vmode="integer", finalizer="close",
##				    overwrite=TRUE,
##				    dimnames=list(fns, fvarlabels))
##		save(featureDataFF, file=file.path(outdir, "featureDataFF.rda"))
##		close(featureDataFF)
##		rm(M, cnProbes, snpProbes, featureDataFF); gc()
##	}
##	## parameters file
##	parameterFile <- file.path(outdir, "cnParams.rda")
##	if(!file.exists(parameterFile)) {
##		message("Initializing parameters file")
##		batch <- cnOptions[["batch"]]
##		dns.batch <- list(fns, unique(batch)) 
##		cnParams <- initializeParamObject(dns.batch)
##		save(cnParams, file=file.path(outdir, "cnParams.rda"))
##		close(cnParams)
##	}
##}

##preprocessAndGenotype <- function(filenames, cnOptions, ...){
##	set.seed(cnOptions[["seed"]])  ##for reproducibility
##	protocolFile <- cnOptions[["protocolFile"]]
##	cdfName <- cnOptions[["cdfName"]]
##	verbose <- cnOptions[["verbose"]]
##	if(file.exists(protocolFile)){
##		## check that file is the same dimension
##		load(protocolFile)
##		if(!identical(sampleNames(protocoldata), basename(filenames)))
##			unlink(protocolFile)
##	}
##	if(!file.exists(protocolFile)){
##		platform <- whichPlatform(paste(cdfName, "Crlmm", sep=""))
##		if(platform=="affymetrix"){
##			if(verbose) message("Creating protocol file with scan dates for the affy arrays")
##			scanDates <- data.frame(ScanDate=sapply(filenames, celfileDate))
##			rownames(scanDates) <- basename(rownames(scanDates))
##			protocoldata <- new("AnnotatedDataFrame",
##					    data=scanDates,
##					    varMetadata=data.frame(labelDescription=colnames(scanDates),
##					    row.names=colnames(scanDates)))
##			save(protocoldata, file=protocolFile)
##		}
##		## protocol file for illumina saved during the readIdatFile step
##	} 
##	if(isPackageLoaded("ff")) initializeFFObjects(filenames, cnOptions)
##	crlmmWrapper(filenames, cnOptions, ...)
##	message("Checking for required files...")
##	message(cnOptions[["AFile"]], ": ", file.exists(cnOptions[["AFile"]]))
##	message(cnOptions[["BFile"]], ": ", file.exists(cnOptions[["BFile"]]))
##	message(cnOptions[["callsFile"]], ": ", file.exists(cnOptions[["callsFile"]]))
##	message(cnOptions[["confsFile"]], ": ", file.exists(cnOptions[["confsFile"]]))
##	message(cnOptions[["snprmaFile"]], ": ", file.exists(cnOptions[["snprmaFile"]]))
##	message(cnOptions[["protocolFile"]], ": ", file.exists(cnOptions[["protocolFile"]]))
##}
	
##crlmmCopynumber <- function(cnOptions, ...){
##crlmmCopynumber <- function(object){
##	ops <- crlmmOptions(object)
##	verbose <- ops$verbose
##	calls <- snpCall(object)
##	confs <- confs(object)
##	fns <- featureNames(object)
##	SNRmin <- ops$SNRMin
##	batch <- object$batch
##	whichBatch <- ops$cnOpts$whichBatch
##	chromosome <- ops$cnOpts$chromosome
##	MIN.SAMPLES <- ops$cnOpts$MIN.SAMPLES
##	##k <- grep("chr", colnames(snpProbes))
##	for(CHR in chromosome){
##		##annotated snp and cn probes
##		##snps <- rownames(snpProbes)[snpProbes[, k] == CHR]
##		##cns <- rownames(cnProbes)[cnProbes[, k] == CHR]
##		##where are the annotated snps in the fns file
##		##index.snps <- match(snps, fns)
##		##index.cn <- match(cns, fns)
##		##row.index <- c(index.snps, index.cn)
##		cat("Chromosome ", CHR, "\n")
##		for(i in whichBatch){
##			PLATE <- unique(batch)[i]
##			message("Plate: ", PLATE)
##			sample.index <- which(batch==PLATE)
##			if(length(sample.index) < MIN.SAMPLES) {
##				warning("Plate ", PLATE, " has fewer than 10 samples.  Skipping this plate.")
##				next()
##			}
##			##cnOptions[["batch"]] <- cnOptions[["batch"]][snpI[["SNR"]]  >= SNRmin]
####			if(isPackageLoaded("ff")){
####				ca <- as.matrix(CA[row.index, sample.index])
####				cb <- as.matrix(CB[row.index, sample.index])
####			} else{
####				dns <- dimnames(A[row.index, sample.index])
####				cb <- ca <- matrix(NA, nr=length(row.index), nc=length(sample.index), dimnames=dns)
####			}
####			cnSet <- new("CNSet",
####				     alleleA=as.matrix(A[row.index, sample.index]),
####				     alleleB=as.matrix(B[row.index, sample.index]),
####				     call=as.matrix(calls[row.index, sample.index]),
####				     callProbability=as.matrix(confs[row.index, sample.index]),
####				     CA=ca,
####				     CB=cb,
####				     featureData=annotatedDataFrameFrom(as.matrix(A[row.index, sample.index]), byrow=TRUE),
####				     phenoData=pD[sample.index, ],
####				     protocolData=protocoldata[sample.index, ])
##			##Verify this is correct
####			annotation(cnSet) <- cnOptions[["cdfName"]]
####			featureNames(cnSet) <- fns[row.index]
##			##add chromosome, position, isSnp
####			cnSet <- annotate(cnSet)
####			if(any(cnSet$SNR > SNRmin)){
####				if(CHR == chromosome[1]) message(paste("Excluding samples with SNR < ", SNRmin))
####				cnSet <- cnSet[, cnSet$SNR >= SNRmin]
####			}
####			featureData(cnSet) <- lm.parameters(cnSet, cnOptions)
##			if(CHR > 23) next()
##			cnSet <- computeCopynumber(object[chromosome(object) == CHR, sample.index])
####			if(!isPackageLoaded("ff") & i == whichBatch[1]) cnParams <- initializeParamObject(list(featureNames(cnSet), unique(cnOptions[["batch"]])[whichBatch]))
####			if(!isPackageLoaded("ff")) {
####				row.index <- 1:nrow(cnSet)
####			} else {
####				##Warning message:
####				##In d[[1]] * d[[2]] : NAs produced by integer overflow
####				CA[row.index, sample.index] <- cnSet@assayData[["CA"]]
####				CB[row.index, sample.index] <- cnSet@assayData[["CB"]]
####			}
####			cnParams <- updateParams(cnParams, cnSet, row.index, batch=unique(batch)[i])
##			## keep only chromosome, position, and isSnp
####			featureData(cnSet) <- featureData(cnSet)[, 1:3]
####			if(!isPackageLoaded("ff")){
####				save(cnSet, file=paste(cnFile, "_", PLATE, "_", CHR, ".rda", sep=""))
####				save(cnParams, file=paste(outdir, "cnParams_", PLATE, "_", CHR, ".rda", sep=""))
####			}
##		} ## end of batch loop
##	} ## end of chromosome loop
####	if(isPackageLoaded("ff")) {
####		close(cnParams)
####		close(A); close(B)		
####		close(CA); close(CB)
####		save(CA, file=CAFile)
####		save(CB, file=CBFile)
####		close(calls); close(confs)
####		return()
####	}
##	return(cnSet)
##}







##loadIlluminaRG <- function(rgFile, crlmmFile, load.it, save.it,...){
####	if(missing(rgFile)){
####		##stop("must specify 'rgFile'.")
####		rgFile <- file.path(dirname(crlmmFile), "rgFile.rda")
####		message("rgFile not specified.  Using ", rgFile)
####	}
##	if(!load.it){
##		RG <- readIdatFiles(...)
##		if(save.it) save(RG, file=rgFile)
##	}
##	if(load.it & !file.exists(rgFile)){
##		message("load.it is TRUE, but rgFile not present.  Attempting to read the idatFiles.")
##		RG <- readIdatFiles(...)
##		if(save.it) save(RG, file=rgFile)
##	}
##	if(load.it & file.exists(rgFile)){
##		message("Loading RG file")
##		load(rgFile)
##		RG <- get("RG")
##	}
##	return(RG)
##}
##
##loadIlluminaCallSet <- function(crlmmFile, snprmaFile, RG, load.it, save.it, cdfName){
##	if(!file.exists(crlmmFile) | !load.it){		
##		callSet <- crlmmIllumina(RG=RG,
##					 cdfName=cdfName,
##					 sns=sampleNames(RG),
##					 returnParams=TRUE,
##					 save.it=TRUE,
##					 intensityFile=snprmaFile)
##		if(save.it) save(callSet, file=crlmmFile)
##	} else {
##		message("Loading ", crlmmFile, "...")
##		load(crlmmFile)
##		callSet <- get("callSet")
##	}
##	protocolData(callSet) <- protocolData(RG)
##	return(callSet)
##}


##loadAffyCallSet <- function(filenames, confsFile, callsFile, snprmaFile, load.it, save.it,  cdfName){
##
####		if(save.it){
####			message("Saving callSet to", callsFile)
####			##save(callSet, cnrmaResult, file=callsFile)
####			save(callSet, file=callsFile)
####		}
####	} else {
####		message("Loading ", callsFile, "...")
####		##load(snprmaFile)				
####		load(callsFile)
####		callSet <- get("callSet")
####		##cnrmaResult <- get("cnrmaResult")
####	}
##
##}
##	if(platform=="affymetrix") {
##		protocolData(callSet)[["ScanDate"]] <- as.character(celDates(filenames))
##		sampleNames(protocolData(callSet)) <- sampleNames(callSet)
##	}

##loadAffyCnrma <- function(filenames, cnrmaFile, cdfName, outdir, load.it, save.it, use.bigmemory=FALSE){
##	if(!file.exists(cnrmaFile) | !load.it){	
##		message("Quantile normalizing the copy number probes...")
##		cnrmaResult <- cnrma2(filenames=filenames, cdfName=cdfName, outdir=outdir, cnrmaFile=cnrmaFile)
##	} else cnrmaResult <- attach.big.matrix("NP.desc", outdir)
##	return(cnrmaResult)
##}

##loadIlluminaCnrma <- function(){
##	if(exists("cnAB")){
##		np.A <- as.integer(cnAB$A)
##		np.B <- as.integer(cnAB$B)
##		np <- ifelse(np.A > np.B, np.A, np.B)
##		np <- matrix(np, nrow(cnAB$A), ncol(cnAB$A))
##		rownames(np) <- cnAB$gns
##		colnames(np) <- cnAB$sns
##		cnAB$NP <- np
##		##sampleNames(callSet) <- res$sns
##		sampleNames(callSet) <- cnAB$sns
##		cnrmaResult <- get("cnAB")
##	} else cnrmaResult <- NULL
##	return(cnrmaResult)
##}
##
##crlmmWrapper <- function(filenames, cnOptions, ...){
##	crlmmBatchSize <- cnOptions[["crlmmBatchSize"]]
##	cdfName <- cnOptions[["cdfName"]]
##	load.it <- cnOptions[["load.it"]]
##	save.it <- cnOptions[["save.it"]]
##	callsFile <- cnOptions[["callsFile"]]
##	confsFile <- cnOptions[["confsFile"]]
##	AFile=cnOptions[["AFile"]]
##	BFile=cnOptions[["BFile"]]	
##	snprmaFile=cnOptions[["snprmaFile"]]
##	cnrmaFile=cnOptions[["cnrmaFile"]]
##	rgFile=cnOptions[["rgFile"]]
##	protocolFile <- cnOptions[["protocolFile"]]
##	outdir <- cnOptions[["outdir"]]
##	if(missing(cdfName)) stop("cdfName is missing -- a valid cdfName is required.  See crlmm:::validCdfNames()")
##	platform <- whichPlatform(cdfName)
##	if(!(platform %in% c("affymetrix", "illumina"))){
##		stop("Only 'affymetrix' and 'illumina' platforms are supported at this time.")
##	} else {
##		if(!isValidCdfName(cdfName)){
##			stop(cdfName, " is not a valid entry.  See crlmm:::validCdfNames(platform)")
##		} else  message("Using the annotation package ", cdfName, " for this ", platform, " platform")
##	}
##	if(platform == "illumina") {
##		if(!file.exists(rgFile)){
##			if(load.it) message(rgFile, " does not exist and you chose to load.it.  Re-reading the R and G intensities from the IDAT files")
##			sampleSheet <- cnOptions$sampleSheet
##			ids <- cnOptions$ids
##			arrayInfoColNames <- cnOptions$arrayInfoColNames
##			highDensity <- cnOptions$highDensity
##			##this is either an NChannelSet object, or a list of pointers to the ff objects
##			RG <- readIdatFiles(sampleSheet=sampleSheet,
##					    arrayNames=basename(filenames),
##					    ids=ids,
##					    path=dirname(filenames),
##					    highDensity=highDensity,
##					    fileExt=cnOptions$fileExt[1:2],
##					    sep=cnOptions$fileExt[[3]],
##					    saveDate=FALSE,  ## I do this earlier
##					    verbose=cnOptions[["verbose"]],
##					    protocolFile=protocolFile)
##			if(save.it) save(RG, file=rgFile)
##			##RG <- loadIlluminaRG(rgFile, callsFile, load.it, save.it)
##		} else{
##			if(!isPackageLoaded("ff")) {load(rgFile); RG <- get("RG")}
##		}
##	}
##	if(!(file.exists(dirname(callsFile)))) stop(dirname(callsFile), " does not exist.")
##	if(!(file.exists(dirname(snprmaFile)))) stop(dirname(snprmaFile), " does not exist.")
##	if(platform == "affymetrix"){
##		crlmm(filenames=filenames,
##		      cdfName=cdfName,
##		      save.it=TRUE,
##		      load.it=load.it,
##		      snprmaFile=snprmaFile,
##		      callsFile=callsFile,
##		      confsFile=confsFile,
##		      AFile=AFile,
##		      BFile=BFile,
##		      crlmmBatchSize=crlmmBatchSize,
##		      SNRMin=cnOptions[["SNRMin"]])
##	}
##	gc()
##	if(platform == "illumina") {
##		callSet <- crlmmIllumina(RG=RG,
##					 cdfName=cdfName,
##					 sns=sampleNames(RG),
##					 returnParams=TRUE,
##					 save.it=TRUE,
##					 snprmaFile=snprmaFile,
##					 callsFile=callsFile,
##					 confsFile=confsFile,
##					 AFile=AFile,
##					 BFile=BFile)
##		##callSet <- loadIlluminaCallSet(callsFile, snprmaFile, RG, load.it, save.it, cdfName)
##	}
##	if(platform == "affymetrix"){
##		if(!file.exists(cnrmaFile) | !load.it){	
##			message("Quantile normalizing the copy number probes...")
##			## updates A matrix and saves cnrmaFile
##			cnrma(filenames=filenames, cdfName=cdfName, outdir=outdir, verbose=cnOptions[["verbose"]], cnrmaFile=cnrmaFile, AFile=AFile, snprmaFile=snprmaFile)
##		} 
##	}
####	if(!is.null(cnrmaResult)){
####		for(CHR in chromosome){
####			cat(CHR, " ")
####			cnps <- rownames(cnProbes)[cnProbes[, k] == CHR]
####			index.nps <- match(cnps, rownames(cnrmaResult[["NP"]]))
####			NP <- cnrmaResult$NP[index.nps, ]
####			save(NP, file=file.path(tmpdir, paste("NP_", CHR, ".rda", sep="")))
####			rm(NP); gc()
####		}
####	}
##	if(!save.it){
##		message("Cleaning up")		
##		unlink(snprmaFile); unlink(cnrmaFile)
##	}
##}



# steps: quantile normalize hapmap: create 1m_reference_cn.rda object
##cnrma <- function(filenames, cdfName, sns, seed=1, verbose=FALSE, outdir){
##	if(missing(cdfName)) stop("must specify cdfName")
##	pkgname <- getCrlmmAnnotationName(cdfName)
##	require(pkgname, character.only=TRUE) || stop("Package ", pkgname, " not available")
##	if (missing(sns)) sns <- basename(filenames)
##        loader("npProbesFid.rda", .crlmmPkgEnv, pkgname)
##	fid <- getVarInEnv("npProbesFid")
##	set.seed(seed)
##	idx2 <- sample(length(fid), 10^5) ##for skewness. no need to do everything
##	SKW <- vector("numeric", length(filenames))
##
##	
##	NP <- matrix(NA, length(fid), length(filenames))
##	verbose <- TRUE
##	if(verbose){
##		message("Processing ", length(filenames), " files.")
##		if (getRversion() > '2.7.0') pb <- txtProgressBar(min=0, max=length(filenames), style=3)
##	}
##	if(cdfName=="genomewidesnp6"){
##		loader("1m_reference_cn.rda", .crlmmPkgEnv, pkgname)
##	}
##	if(cdfName=="genomewidesnp5"){
##		loader("5.0_reference_cn.rda", .crlmmPkgEnv, pkgname)
##	}
##	reference <- getVarInEnv("reference")
##	##if(!is.matrix(reference)) stop("target distribution for quantile normalization not available.")
##	for(i in seq(along=filenames)){
##		y <- as.matrix(read.celfile(filenames[i], intensity.means.only=TRUE)[["INTENSITY"]][["MEAN"]][fid])
##		x <- log2(y[idx2])
##		SKW[i] <- mean((x-mean(x))^3)/(sd(x)^3)
##		rm(x)
##		NP[, i] <- as.integer(normalize.quantiles.use.target(y, target=reference))
##		if (verbose)
##			if (getRversion() > '2.7.0') setTxtProgressBar(pb, i)
##			else cat(".")
##	}
##	dimnames(NP) <- list(names(fid), sns)
##	##dimnames(NP) <- list(map[, "man_fsetid"], sns)
##	res3 <- list(NP=NP, SKW=SKW)
##	cat("\n")
##	return(res3)
##}

cnrma <- function(object, filenames){
	ops <- crlmmOptions(object)
	cdfName <- annotation(object)
	seed <- ops$seed
	verbose <- ops$verbose
	##cnrmaFile <- ops$cnrmaFile
	A <- A(object)
	if(missing(cdfName)) stop("must specify cdfName")
	pkgname <- getCrlmmAnnotationName(cdfName)
	require(pkgname, character.only=TRUE) || stop("Package ", pkgname, " not available")
	sns <- basename(filenames)
        loader("npProbesFid.rda", .crlmmPkgEnv, pkgname)
	fid <- getVarInEnv("npProbesFid")
	set.seed(seed)
	idx2 <- sample(length(fid), 10^5) ##for skewness. no need to do everything
	SKW <- vector("numeric", length(filenames))
	index <- match(names(fid), featureNames(object))
	stopifnot(identical(featureNames(object)[index], names(fid)))
	if(length(index) < 1) stop("None of the names for the nonpolymorphic probes in the annotation package match the names stored in the snprmaFile.")
	if(verbose){
		message("Processing ", length(filenames), " files.")
		if (getRversion() > '2.7.0') pb <- txtProgressBar(min=0, max=length(filenames), style=3)
	}
	if(cdfName=="genomewidesnp6"){
		loader("1m_reference_cn.rda", .crlmmPkgEnv, pkgname)
	}
	if(cdfName=="genomewidesnp5"){
		loader("5.0_reference_cn.rda", .crlmmPkgEnv, pkgname)
	}
	reference <- getVarInEnv("reference")
	for(i in seq(along=filenames)){
		y <- as.matrix(read.celfile(filenames[i], intensity.means.only=TRUE)[["INTENSITY"]][["MEAN"]][fid])
		x <- log2(y[idx2])
		SKW[i] <- mean((x-mean(x))^3)/(sd(x)^3)
		rm(x); gc()
		A[index, i] <- as.integer(normalize.quantiles.use.target(y, target=reference))
		if (verbose)
			if (getRversion() > '2.7.0') setTxtProgressBar(pb, i)
			else cat(".")
		rm(y); gc()
	}
	cat("\nDone\n")
	pData(object)$SKW_nonpolymorphic <- SKW
	object@assayData[["alleleA"]] <- A
	return(object)
}

getFlags <- function(object, PHI.THR){
	batch <- unique(object$batch)
	nuA <- getParam(object, "nuA", batch)
	nuB <- getParam(object, "nuB", batch)
	phiA <- getParam(object, "phiA", batch)
	phiB <- getParam(object, "phiB", batch)
	negativeNus <- nuA < 1 | nuB < 1
	negativePhis <- phiA < PHI.THR | phiB < PHI.THR
	negativeCoef <- negativeNus | negativePhis
	notfinitePhi <- !is.finite(phiA) | !is.finite(phiB)
	flags <- negativeCoef | notfinitePhi
	return(flags)
}


instantiateObjects <- function(object){
	Ns <- matrix(NA, nrow(object), 5)
	colnames(Ns) <- c("A", "B", "AA", "AB", "BB")
	vA <- vB <- muB <- muA <- Ns
	normal <- matrix(TRUE, nrow(object), ncol(object))
	dimnames(normal) <- list(featureNames(object), sampleNames(object))
	tmp.objects <- list(vA=vA,
			    vB=vB,
			    muB=muB,
			    muA=muA,
			    Ns=Ns,
			    normal=normal)
        return(tmp.objects)
}

thresholdCopynumber <- function(object){
	ca <- CA(object)
	cb <- CB(object)
	ca[ca < 0.05] <- 0.05
	ca[ca > 5] <- 5
	cb[cb < 0.05] <- 0.05
	cb[cb > 5] <- 5
	CA(object) <- ca
	CB(object) <- cb
	return(object)
}

##linear model parameters
##lm.parameters <- function(object, cnOptions){
##	fD <- fData(object)
##	batch <- object$batch
##	uplate <- unique(batch)
##	parameterNames <- c(paste("tau2A", uplate, sep="_"),
##			    paste("tau2B", uplate, sep="_"),
##			    paste("sig2A", uplate, sep="_"),
##			    paste("sig2B", uplate, sep="_"),
##			    paste("nuA", uplate, sep="_"),
##			    paste("nuA.se", uplate, sep="_"),			    
##			    paste("nuB", uplate, sep="_"),
##			    paste("nuB.se", uplate, sep="_"),			    			    
##			    paste("phiA", uplate, sep="_"),
##			    paste("phiA.se", uplate, sep="_"),			    
##			    paste("phiB", uplate, sep="_"),
##			    paste("phiB.se", uplate, sep="_"),			    
##			    paste("phiAX", uplate, sep="_"),
##			    paste("phiBX", uplate, sep="_"),			    
##			    paste("corr", uplate, sep="_"),
##			    paste("corrA.BB", uplate, sep="_"),
##			    paste("corrB.AA", uplate, sep="_"))
##	pMatrix <- data.frame(matrix(numeric(0),
##				     nrow(object),
##				     length(parameterNames)),
##				     row.names=featureNames(object))
##	colnames(pMatrix) <- parameterNames
##	fD2 <- cbind(fD, pMatrix)
##	new("AnnotatedDataFrame", data=fD2,
##	    varMetadata=data.frame(labelDescription=colnames(fD2),
##	    row.names=colnames(fD2)))
##}

nonpolymorphic <- function(object, cnOptions, tmp.objects){
	chromosome <- cnOptions[["chromosome"]]
	batch <- unique(object$batch)
	CHR <- unique(chromosome(object))
	verbose <- cnOptions[["verbose"]]
	if(CHR != chromosome[1]) verbose <- FALSE
	if(batch != unique(cnOptions[["batch"]])[1]) verbose <- FALSE
	goodSnps <- function(object, PHI.THR, tmp.objects, nAA.THR, nBB.THR){
		Ns <- tmp.objects[["Ns"]]
		##Ns <- get("Ns", envir)
		flags <- getFlags(object, PHI.THR)
		fewAA <- Ns[, "AA"] < nAA.THR
		fewBB <- Ns[, "BB"] < nBB.THR
		flagsA <- flags | fewAA
		flagsB <- flags | fewBB
		flags <- list(A=flagsA, B=flagsB)
		return(flags)
	}
	nAA.THR <- cnOptions$nHOM.THR
	nBB.THR <- cnOptions$nHOM.THR
	PHI.THR <- cnOptions$PHI.THR
	snpflags <- goodSnps(object, PHI.THR, tmp.objects, nAA.THR, nBB.THR)
	flagsA <- snpflags$A
	flagsB <- snpflags$B
##	if(all(flagsA) | all(flagsB)) stop("all snps are flagged")
	nuA <- getParam(object, "nuA", batch)
	nuB <- getParam(object, "nuB", batch)
	phiA <- getParam(object, "phiA", batch)
	phiB <- getParam(object, "phiB", batch)	
	sns <- sampleNames(object)
	muA <- tmp.objects[["muA"]]
	muB <- tmp.objects[["muB"]]
	A <- A(object)
	B <- B(object)
	CA <- CA(object)
	CB <- CB(object)
	if(CHR == 23){
		phiAX <- getParam(object, "phiAX", batch)
		phiBX <- getParam(object, "phiBX", batch)
	}
	##---------------------------------------------------------------------------
	## Train on unflagged SNPs
	##---------------------------------------------------------------------------
	##Might be best to train using the X chromosome, since for the
	##X phi and nu have been adjusted for cross-hybridization
	##plateInd <- plate == uplate[p]
	##muA <- muA[!flagsA, p, c("A", "AA")]
	##muB <- muB[!flagsB, p, c("B", "BB")]
	muA <- muA[!flagsA, "AA"]
	muB <- muB[!flagsB, "BB"]
	X <- cbind(1, log2(c(muA, muB)))
	Y <- log2(c(phiA[!flagsA], phiB[!flagsB]))
	if(nrow(X) > 5000){
		ix <- sample(1:nrow(X), 5000)
	} else {
		ix <- 1:nrow(X)
	}
	betahat <- solve(crossprod(X[ix, ]), crossprod(X[ix, ], Y[ix]))
	normal <- tmp.objects[["normal"]][!isSnp(object), ]	
	if(CHR == 23){
		##normalNP <- envir[["normalNP"]]
		##normalNP <- normalNP[, plate==uplate[p]]
		##nuT <- envir[["nuT"]]
		##phiT <- envir[["phiT"]]
		
		##cnvs <- envir[["cnvs"]]
                ##loader("cnProbes.rda", pkgname=pkgname, envir=.crlmmPkgEnv)
                ##cnProbes <- get("cnProbes", envir=.crlmmPkgEnv)
		##cnProbes <- cnProbes[match(cnvs, rownames(cnProbes)), ]

		##For build Hg18
		##http://genome.ucsc.edu/cgi-bin/hgGateway
		##pseudo-autosomal regions on X
		##chrX:1-2,709,520 and chrX:154584237-154913754, respectively
		##par:pseudo-autosomal regions
		pseudoAR <- position(object) < 2709520 | (position(object) > 154584237 & position(object) < 154913754)
		##pseudoAR <- cnProbes[, "position"] < 2709520 | (cnProbes[, "position"] > 154584237 & cnProbes[, "position"] < 154913754)
		##in case some of the cnProbes are not annotated
		pseudoAR[is.na(pseudoAR)] <- FALSE
		pseudoAR <- pseudoAR[!isSnp(object)]
		##gender <- envir[["gender"]]
		gender <- object$gender
		obj1 <- object[!isSnp(object), ]
		A.male <- A(obj1[, gender==1])
		mu1 <- rowMedians(A.male, na.rm=TRUE)
		##mu1 <- rowMedians(NP[, gender=="male"], na.rm=TRUE)
		##mu2 <- rowMedians(NP[, gender=="female"], na.rm=TRUE)
		A.female <- A(obj1[, gender==2])
		mu2 <- rowMedians(A.female, na.rm=TRUE)
		mus <- log2(cbind(mu1, mu2))
		X.men <- cbind(1, mus[, 1])
		X.fem <- cbind(1, mus[, 2])
		
		Yhat1 <- as.numeric(X.men %*% betahat)
		Yhat2 <- as.numeric(X.fem %*% betahat)
		phi1 <- 2^(Yhat1)
		phi2 <- 2^(Yhat2)
		nu1 <- 2^(mus[, 1]) - phi1
		nu2 <- 2^(mus[, 2]) - 2*phi2

		if(any(pseudoAR)){
			nu1[pseudoAR] <- 2^(mus[pseudoAR, 1]) - 2*phi1[pseudoAR]
		}
		CT1 <- 1/phi1*(A.male-nu1)
		CT2 <- 1/phi2*(A.female-nu2)
		##CT2 <- 1/phi2*(NP[, gender=="female"]-nu2)
		##CT1 <- matrix(as.integer(100*CT1), nrow(CT1), ncol(CT1))
		##CT2 <- matrix(as.integer(100*CT2), nrow(CT2), ncol(CT2))
		##CT <- envir[["CT"]]
		CA <- CA(obj1)
		CA[, gender==1] <- CT1
		CA[, gender==2] <- CT2
		CA(object)[!isSnp(object), ] <- CA
		##CT[, plate==uplate[p] & gender=="male"] <- CT1
		##CT[, plate==uplate[p] & gender=="female"] <- CT2
		##envir[["CT"]] <- CT

		##only using females to compute the variance
		##normalNP[, gender=="male"] <- NA
		normal[, gender==1] <- NA
		sig2A <- getParam(object, "sig2A", batch)
		normal.f <- normal[, object$gender==2]
		sig2A[!isSnp(object)] <- rowMAD(log2(A.female*normal.f), na.rm=TRUE)^2
		sig2A[!isSnp(object) & is.na(sig2A)] <- median(sig2A[!isSnp(object)], na.rm=TRUE)
		##sig2T[, p] <- rowMAD(log2(NP*normalNP), na.rm=TRUE)^2
		object <- pr(object, "sig2A", batch, sig2A)

		nuA[!isSnp(object)] <- nu2
		phiA[!isSnp(object)] <- phi2
		
		THR.NU.PHI <- cnOptions$THR.NU.PHI
		if(THR.NU.PHI){
			##Assign values to object
			object <- pr(object, "nuA", batch, nuA)
			object <- pr(object, "phiA", batch, phiA)
			if(verbose) message("Thresholding nu and phi")
			object <- thresholdModelParams(object, cnOptions)
		} else {
			object <- pr(object, "nuA", batch, nuA)		
			object <- pr(object, "phiA", batch, phiA)
		}
	} else {
		A <- A(object)[!isSnp(object), ]
		mus <- rowMedians(A * normal, na.rm=TRUE)
		crosshyb <- max(median(muA) - median(mus), 0)
		X <- cbind(1, log2(mus+crosshyb))
		logPhiT <- X %*% betahat
		phiA[!isSnp(object)] <- 2^(logPhiT)
		nuA[!isSnp(object)] <- mus-2*phiA[!isSnp(object)]

		THR.NU.PHI <- cnOptions$THR.NU.PHI
		if(THR.NU.PHI){
			##Assign values to object
			object <- pr(object, "nuA", batch, nuA)
			object <- pr(object, "phiA", batch, phiA)			
			if(verbose) message("Thresholding nu and phi")
			object <- thresholdModelParams(object, cnOptions)
			##reassign values (now thresholded at MIN.NU and MIN.PHI
			nuA <- getParam(object, "nuA", batch)
			phiA <- getParam(object, "phiA", batch)
		}
		CA(object)[!isSnp(object), ] <- 1/phiA[!isSnp(object)]*(A - nuA[!isSnp(object)])
		sig2A <- getParam(object, "sig2A", batch)
		sig2A[!isSnp(object)] <- rowMAD(log2(A*normal), na.rm=TRUE)^2
		object <- pr(object, "sig2A", batch, sig2A)
		##added
		object <- pr(object, "nuA", batch, nuA)
		object <- pr(object, "phiA", batch, phiA)
	}
	return(object)
}

nonpolymorphic.poe <- function(object, cnOptions, tmp.object){
	require(metaArray)
	nps <- log2(A(object)[!isSnp(object), ])
	nps <- (nps-rowMedians(nps))/rowMAD(nps)
	runAvg <- apply(nps, 2, myfilter, filter=rep(1/10, 10))
	rownames(runAvg) <- featureNames(object)[!isSnp(object)]
	rm.nas <- rowSums(is.na(runAvg)) == 0
	runAvg <- runAvg[rm.nas, ]
	
	poe.scale <- poe.em(runAvg, cl=rep(0, ncol(nps)))$data
	pinegg <- piposg <- poe.scale
	piposg[piposg < 0] <- 0
	pinegg[pinegg > 0] <- 0
	pinegg <- pinegg*-1
	pm.em <- 1*pinegg + 2*(1-pinegg-piposg) + 3*piposg
	rownames(pm.em) <- rownames(runAvg)
	CA(object)[match(rownames(pm.em), featureNames(object)), ] <- pm.em
	##CA(object)[!isSnp(object), ] <- pm.em
	return(object)
}

##sufficient statistics on the intensity scale
withinGenotypeMoments <- function(object, cnOptions, tmp.objects){
	normal <- tmp.objects[["normal"]]
	## muA, muB: robust estimates of the within-genotype center (intensity scale)
	muA <- tmp.objects[["muA"]]
	muB <- tmp.objects[["muB"]]
	## vA, vB: robust estimates of the within-genotype variance (intensity scale)
	vA <- tmp.objects[["vA"]]
	vB <- tmp.objects[["vB"]]
	Ns <- tmp.objects[["Ns"]]
	G <- calls(object) 
	GT.CONF.THR <- cnOptions$GT.CONF.THR
	CHR <- unique(chromosome(object))

	A <- A(object)
	B <- B(object)
##	highConf <- (1-exp(-confs(object)/1000)) > GT.CONF.THR
	highConf <- confs(object) > GT.CONF.THR
	##highConf <- highConf > GT.CONF.THR
	if(CHR == 23){
		gender <- object$gender
##		gender <- envir[["gender"]]
		IX <- matrix(gender, nrow(G), ncol(G), byrow=TRUE)
##		IX <- IX == "female"
		IX <- IX == 2  ##2=female, 1=male
	} else IX <- matrix(TRUE, nrow(G), ncol(G))
	index <- GT.B <- GT.A <- vector("list", 3)
	names(index) <- names(GT.B) <- names(GT.A) <- c("AA", "AB", "BB")
	##--------------------------------------------------
	##within-genotype sufficient statistics
	##--------------------------------------------------
	##GT.B <- GT.A <- list()
	snpIndicator <- matrix(isSnp(object), nrow(object), ncol(object)) ##RS: added
	for(j in 1:3){
		GT <- G==j & highConf & IX & snpIndicator
		GT <- GT * normal
		Ns[, j+2] <- rowSums(GT, na.rm=TRUE)				
		GT[GT == FALSE] <- NA
		GT.A[[j]] <- GT*A
		GT.B[[j]] <- GT*B
		index[[j]] <- which(Ns[, j+2] > 0 & isSnp(object)) ##RS: added		
		muA[, j+2] <- rowMedians(GT.A[[j]], na.rm=TRUE)
		muB[, j+2] <- rowMedians(GT.B[[j]], na.rm=TRUE)
		vA[, j+2] <- rowMAD(GT.A[[j]], na.rm=TRUE)
		vB[, j+2] <- rowMAD(GT.B[[j]], na.rm=TRUE)

		##Shrink towards the typical variance
		DF <- Ns[, j+2]-1
		DF[DF < 1] <- 1
		v0A <- median(vA[, j+2], na.rm=TRUE)
		v0B <- median(vB[, j+2], na.rm=TRUE)
		if(v0A == 0) v0A <- NA
		if(v0B == 0) v0B <- NA
		DF.PRIOR <- cnOptions$DF.PRIOR
		vA[, j+2] <- (vA[, j+2]*DF + v0A*DF.PRIOR)/(DF.PRIOR+DF)
		vA[is.na(vA[, j+2]), j+2] <- v0A
		vB[, j+2] <- (vB[, j+2]*DF + v0B*DF.PRIOR)/(DF.PRIOR+DF)
		vB[is.na(vB[, j+2]), j+2] <- v0B
	}
	if(CHR == 23){
		k <- 1
		for(j in c(1,3)){
			GT <- G==j & highConf & !IX 
			Ns[, k] <- rowSums(GT)
			GT[GT == FALSE] <- NA
			muA[, k] <- rowMedians(GT*A, na.rm=TRUE)
			muB[, k] <- rowMedians(GT*B, na.rm=TRUE)
			vA[, k] <- rowMAD(GT*A, na.rm=TRUE)
			vB[, k] <- rowMAD(GT*B, na.rm=TRUE)
			
			DF <- Ns[, k]-1
			DF[DF < 1] <- 1
			v0A <- median(vA[, k], na.rm=TRUE)
			v0B <- median(vB[, k], na.rm=TRUE)
			vA[, k] <- (vA[, k]*DF + v0A*DF.PRIOR)/(DF.PRIOR+DF)
			vA[is.na(vA[, k]), k] <- v0A
			vB[, k] <- (vB[, k]*DF + v0B*DF.PRIOR)/(DF.PRIOR+DF)
			vB[is.na(vB[, k]), k] <- v0B			
			k <- k+1
		}
	}
	tmp.objects[["Ns"]] <- Ns
	tmp.objects[["vA"]] <- vA
	tmp.objects[["vB"]] <- vB
	tmp.objects[["muA"]] <- muA
	tmp.objects[["muB"]] <- muB
	tmp.objects$index <- index
	tmp.objects$GT.A <- GT.A
	tmp.objects$GT.B <- GT.B
	return(tmp.objects)
}

oneBatch <- function(object, cnOptions, tmp.objects){
	muA <- tmp.objects[["muA"]]
	muB <- tmp.objects[["muB"]]
	Ns <- tmp.objects[["Ns"]]
	CHR <- unique(chromosome(object))
	PLATE <- unique(object$batch)
	##---------------------------------------------------------------------------
	## Impute sufficient statistics for unobserved genotypes (plate-specific)
	##---------------------------------------------------------------------------
	MIN.OBS <- cnOptions$MIN.OBS
	index.AA <- which(Ns[, "AA"] >= MIN.OBS)
	index.AB <- which(Ns[, "AB"] >= MIN.OBS)
	index.BB <- which(Ns[, "BB"] >= MIN.OBS)
	correct.orderA <- muA[, "AA"] > muA[, "BB"]
	correct.orderB <- muB[, "BB"] > muB[, "AA"]
	##For chr X, this will ignore the males 
	nobs <- rowSums(Ns[, 3:5] >= MIN.OBS, na.rm=TRUE) == 3
	index.complete <- which(correct.orderA & correct.orderB & nobs) ##be selective here
	size <- min(5000, length(index.complete))
	if(size == 5000) index.complete <- sample(index.complete, 5000)
	if(length(index.complete) < 200){
		warning("There are too few samples in plate ", PLATE, " to estimate the copy number for chromosome ", CHR, ".  CA,CB values are NAs")
		return(tmp.objects)
	}
	index <- tmp.objects[["index"]]
	index[[1]] <- which(Ns[, "AA"] == 0 & (Ns[, "AB"] >= MIN.OBS & Ns[, "BB"] >= MIN.OBS))
	index[[2]] <- which(Ns[, "AB"] == 0 & (Ns[, "AA"] >= MIN.OBS & Ns[, "BB"] >= MIN.OBS))
	index[[3]] <- which(Ns[, "BB"] == 0 & (Ns[, "AB"] >= MIN.OBS & Ns[, "AA"] >= MIN.OBS))
	mnA <- muA[, 3:5]
	mnB <- muB[, 3:5]
	for(j in 1:3){
		if(length(index[[j]]) == 0) next()
		X <- cbind(1, mnA[index.complete,  -j, drop=FALSE], mnB[index.complete,  -j, drop=FALSE])
		Y <- cbind(mnA[index.complete, j], mnB[index.complete,  j])
		betahat <- solve(crossprod(X), crossprod(X,Y))
		X <- cbind(1, mnA[index[[j]],  -j, drop=FALSE],  mnB[index[[j]],  -j, drop=FALSE])
		mus <- X %*% betahat
		muA[index[[j]], j+2] <- mus[, 1]
		muB[index[[j]], j+2] <- mus[, 2]
	}
	if(CHR == 23){
		nobsA <- Ns[, "A"] > MIN.OBS
		nobsB <- Ns[, "B"] > MIN.OBS
		notMissing <- !(is.na(muA[, "A"]) | is.na(muA[, "B"]) | is.na(muB[, "A"]) | is.na(muB[, "B"]))
		complete <- list()
		complete[[1]] <- which(correct.orderA & correct.orderB & nobsA & notMissing) ##be selective here
		complete[[2]] <- which(correct.orderA & correct.orderB & nobsB & notMissing) ##be selective here
		if(length(complete[[1]]) < 1 | length(complete[[2]]) < 1) stop("Too few observations to estimate the center for 'A' and 'B'  clusters on chrom X")		
		size <- min(5000, length(complete[[1]]))
		if(size == 5000) complete <- lapply(complete, function(x) sample(x, size))
		index <- list()
		index[[1]] <- which(Ns[, "A"] == 0)
		index[[2]] <- which(Ns[, "B"] == 0)
		cols <- 2:1
		for(j in 1:2){
			if(length(index[[j]]) == 0) next()
			X <- cbind(1, muA[complete[[j]], cols[j]], muB[complete[[j]], cols[j]])
			Y <- cbind(muA[complete[[j]], j], muB[complete[[j]], j])
			betahat <- solve(crossprod(X), crossprod(X,Y))
			X <- cbind(1, muA[index[[j]], cols[j]],  muB[index[[j]], cols[j]])
			mus <- X %*% betahat
			muA[index[[j]], j] <- mus[, 1]
			muB[index[[j]], j] <- mus[, 2]
		}
	}
	##missing two genotypes
	noAA <- Ns[, "AA"] < MIN.OBS
	noAB <- Ns[, "AB"] < MIN.OBS
	noBB <- Ns[, "BB"] < MIN.OBS
	index[[1]] <- noAA & noAB
	index[[2]] <- noBB & noAB
	index[[3]] <- noAA & noBB
##	snpflags <- envir[["snpflags"]]
	snpflags <- index[[1]] | index[[2]] | index[[3]]
##	snpflags[, p] <- index[[1]] | index[[2]] | index[[3]]
	##---------------------------------------------------------------------------
	## Two genotype clusters not observed -- would sequence help? (didn't seem to)
	## 1. extract index of complete data
	## 2. Regress  mu_missing ~ sequence + mu_observed
	## 3. solve for nu assuming the median is 2
	##---------------------------------------------------------------------------
	cols <- c(3, 1, 2)
	for(j in 1:3){
		if(sum(index[[j]]) == 0) next()
		k <- cols[j]
		X <- cbind(1, mnA[index.complete, k], mnB[index.complete, k])
		Y <- cbind(mnA[index.complete,  -k],
			   mnB[index.complete,  -k])
		betahat <- solve(crossprod(X), crossprod(X,Y))
		X <- cbind(1, mnA[index[[j]],  k], mnB[index[[j]],  k])
		mus <- X %*% betahat
		muA[index[[j]], -c(1, 2, k+2)] <- mus[, 1:2]
		muB[index[[j]], -c(1, 2, k+2)] <- mus[, 3:4]
	}
	negA <- rowSums(muA < 0) > 0
	negB <- rowSums(muB < 0) > 0	
	snpflags <- snpflags| negA | negB | rowSums(is.na(muA[, 3:5]), na.rm=TRUE) > 0
	tmp.objects$snpflags <- snpflags
	tmp.objects[["muA"]] <- muA
	tmp.objects[["muB"]] <- muB
	tmp.objects[["snpflags"]] <- snpflags
	return(tmp.objects)
}

##Estimate tau2, sigma2, and correlation (updates the object)
locationAndScale <- function(object, cnOptions, tmp.objects){
	Ns <- tmp.objects[["Ns"]]
	index <- tmp.objects[["index"]]
	
	index.AA <- index[[1]]
	index.AB <- index[[2]]
	index.BB <- index[[3]]
	rm(index); gc()

	GT.A <- tmp.objects[["GT.A"]]
	GT.B <- tmp.objects[["GT.B"]]
	AA.A <- GT.A[[1]]
	AB.A <- GT.A[[2]]
	BB.A <- GT.A[[3]]
	
	AA.B <- GT.B[[1]]
	AB.B <- GT.B[[2]]
	BB.B <- GT.B[[3]]
	
	x <- BB.A[index.BB, ]
	batch <- unique(object$batch)
	DF.PRIOR <- cnOptions$DF.PRIOR

	tau2A <- getParam(object, "tau2A", batch)
	tau2A[index.BB] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2	
	DF <- Ns[, "BB"]-1
	DF[DF < 1] <- 1
	med <- median(tau2A, na.rm=TRUE)
	tau2A <- (tau2A * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
	tau2A[is.na(tau2A) & isSnp(object)] <- med
	object <- pr(object, "tau2A", batch, tau2A)

	sig2B <- getParam(object, "sig2B", batch)	
	x <- BB.B[index.BB, ]
	sig2B[index.BB] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2	
	med <- median(sig2B, na.rm=TRUE)
	sig2B <- (sig2B * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
	sig2B[is.na(sig2B) & isSnp(object)] <- med
	object <- pr(object, "sig2B", batch, sig2B)	

	tau2B <- getParam(object, "tau2B", batch)	
	x <- AA.B[index.AA, ]
	tau2B[index.AA] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2		
	DF <- Ns[, "AA"]-1
	DF[DF < 1] <- 1
	med <- median(tau2B, na.rm=TRUE)
	tau2B <- (tau2B * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
	tau2B[is.na(tau2B) & isSnp(object)] <- med
	object <- pr(object, "tau2B", batch, tau2B)

	sig2A <- getParam(object, "sig2A", batch)	
	x <- AA.A[index.AA, ]
	sig2A[index.AA] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2##var(log(IA)|AA)	
	med <- median(sig2A, na.rm=TRUE)
	sig2A <- (sig2A * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
	sig2A[is.na(sig2A) & isSnp(object)] <- med
	object <- pr(object, "sig2A", batch, sig2A)

	if(length(index.AB) > 0){ ##all homozygous is possible
		corr <- getParam(object, "corr", batch)
		x <- AB.A[index.AB, ]
		y <- AB.B[index.AB, ]
		corr[index.AB] <- rowCors(x, y, na.rm=TRUE)
		corr[corr < 0] <- 0
		DF <- Ns[, "AB"]-1
		DF[DF<1] <- 1
		med <- median(corr, na.rm=TRUE)
		corr <- (corr*DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
		corr[is.na(corr) & isSnp(object)] <- med
		object <- pr(object, "corr", batch, corr)		
	}
	corrB.AA <- getParam(object, "corrB.AA", batch)
	backgroundB <- AA.B[index.AA, ]
	signalA <- AA.A[index.AA, ]
	corrB.AA[index.AA] <- rowCors(backgroundB, signalA, na.rm=TRUE)
	DF <- Ns[, "AA"]-1
	DF[DF < 1] <- 1
	med <- median(corrB.AA, na.rm=TRUE)
	corrB.AA <- (corrB.AA*DF + med*DF.PRIOR)/(DF.PRIOR + DF)
	corrB.AA[is.na(corrB.AA) & isSnp(object)] <- med
	object <- pr(object, "corrB.AA", batch, corrB.AA)

	corrA.BB <- getParam(object, "corrA.BB", batch)	
	backgroundA <- BB.A[index.BB, ]
	signalB <- BB.B[index.BB, ]
	corrA.BB[index.BB] <- rowCors(backgroundA, signalB, na.rm=TRUE)
	DF <- Ns[, "BB"]-1
	DF[DF < 1] <- 1
	med <- median(corrA.BB, na.rm=TRUE)
	corrA.BB <- (corrA.BB*DF + med*DF.PRIOR)/(DF.PRIOR + DF)
	corrA.BB[is.na(corrA.BB) & isSnp(object)] <- med
	object <- pr(object, "corrA.BB", batch, corrA.BB)
	return(object)
}

coefs <- function(object, cnOptions, tmp.objects){
	batch <- unique(object$batch)
	CHR <- unique(chromosome(object))
	muA <- tmp.objects[["muA"]]
	muB <- tmp.objects[["muB"]]
	vA <- tmp.objects[["vA"]]
	vB <- tmp.objects[["vB"]]
	Ns <- tmp.objects[["Ns"]]
	if(CHR != 23){
		IA <- muA[, 3:5]
		IB <- muB[, 3:5]
		vA <- vA[, 3:5]
		vB <- vB[, 3:5]
		Np <- Ns[, 3:5]
	} else {
		NOHET <- is.na(median(vA[, "AB"], na.rm=TRUE))
		if(NOHET){
			IA <- muA[, -4]
			IB <- muB[, -4]
			vA <- vA[, -4]
			vB <- vB[, -4]
			Np <- Ns[, -4] 
		} else{
			IA <- muA
			IB <- muB
			vA <- vA
			vB <- vB
			Np <- Ns
		}
		
	}
	Np[Np < 1] <- 1
	vA2 <- vA^2/Np
	vB2 <- vB^2/Np
	wA <- sqrt(1/vA2)
	wB <- sqrt(1/vB2)
	YA <- IA*wA
	YB <- IB*wB
	##update lm.coefficients stored in object
	object <- nuphiAllele(object, allele="A", Ystar=YA, W=wA, tmp.objects=tmp.objects, cnOptions=cnOptions)
	object <- nuphiAllele(object, allele="B", Ystar=YB, W=wB, tmp.objects=tmp.objects, cnOptions=cnOptions)	
	##---------------------------------------------------------------------------
	##Estimate crosshyb using X chromosome and sequence information
	##---------------------------------------------------------------------------
	##browser()
	####data(sequences, package="genomewidesnp6Crlmm")
	##snpflags <- envir[["snpflags"]]
	##muA <- envir[["muA"]][, p, 3:5]
	##muB <- envir[["muB"]][, p, 3:5]
	##Y <- envir[["phiAx"]]
	##load("sequences.rda")
	##seqA <- sequences[, "A", ][, 1]
	##seqA <- seqA[match(snps, names(seqA))]
	##X <- cbind(1, sequenceDesignMatrix(seqA))
	##X <- cbind(X, nuA[, p], phiA[, p], nuB[, p], phiB[, p])
	##missing <- rowSums(is.na(X)) > 0
	##betahat <- solve(crossprod(X[!missing, ]), crossprod(X[!missing, ], Y[!missing]))
	return(object)
}

polymorphic <- function(object, cnOptions, tmp.objects){
	batch <- unique(object$batch)
	CHR <- unique(chromosome(object))
	vA <- tmp.objects[["vA"]]
	vB <- tmp.objects[["vB"]]
	Ns <- tmp.objects[["Ns"]]
	
	nuA <- getParam(object, "nuA", batch)
	nuB <- getParam(object, "nuB", batch)
	nuA.se <- getParam(object, "nuA.se", batch)
	nuB.se <- getParam(object, "nuB.se", batch)

	phiA <- getParam(object, "phiA", batch)
	phiB <- getParam(object, "phiB", batch)
	phiA.se <- getParam(object, "phiA.se", batch)
	phiB.se <- getParam(object, "phiB.se", batch)
	A <- A(object)
	B <- B(object)

	NOHET <- mean(Ns[, "AB"], na.rm=TRUE) < 0.05
	##---------------------------------------------------------------------------
	## Estimate CA, CB
	##---------------------------------------------------------------------------
	if(CHR == 23){
		phiAX <- getParam(object, "phiAX", batch)  ##nonspecific hybridization coef
		phiBX <- getParam(object, "phiBX", batch)  ##nonspecific hybridization coef			
		phistar <- phiBX/phiA  
		tmp <- (B-nuB - phistar*A + phistar*nuA)/phiB
		copyB <- tmp/(1-phistar*phiAX/phiB)
		copyA <- (A-nuA-phiAX*copyB)/phiA
		CB(object) <- copyB  ## multiplies by 100 and converts to integer 
		CA(object) <- copyA
	} else{
		CA(object) <- matrix((1/phiA*(A-nuA)), nrow(A), ncol(A))
		CB(object) <- matrix((1/phiB*(B-nuB)), nrow(B), ncol(B))

	}
	return(object)
}

posteriorProbability.snps <- function(object, cnOptions, tmp.objects=list()){
	I <- isSnp(object)
	gender <- object$gender
	CHR <- unique(chromosome(object))
	batch <- unique(object$batch)
	if(CHR == 23){
		phiAX <- getParam(object, "phiAX", batch)[I]
		phiBX <- getParam(object, "phiBX", batch)[I]
	}
	A <- A(object)[I, ]
	B <- B(object)[I, ]
	sig2A <- getParam(object, "sig2A", batch)[I]
	sig2B <- getParam(object, "sig2B", batch)[I]
	tau2A <- getParam(object, "tau2A", batch)[I]
	tau2B <- getParam(object, "tau2B", batch)[I]	
	corrA.BB <- getParam(object, "corrA.BB", batch)[I]
	corrB.AA <- getParam(object, "corrB.AA", batch)[I]
	corr <- getParam(object, "corr", batch)[I]		
	nuA <- getParam(object, "nuA", batch)[I]
	nuB <- getParam(object, "nuB", batch)[I]	
	phiA <- getParam(object, "phiA", batch)[I]
	phiB <- getParam(object, "phiB", batch)[I]
	normal <- tmp.objects[["normal"]][I, ]
	prior.prob <- cnOptions$prior.prob
	emit <- array(NA, dim=c(nrow(A), ncol(A), 10))##SNPs x sample x 'truth'	
	lA <- log2(A)
	lB <- log2(B)	
	X <- cbind(lA, lB)
	counter <- 1##state counter								
	for(CT in 0:3){
		for(CA in 0:CT){
			cat(".")
			CB <- CT-CA
			A.scale <- sqrt(tau2A*(CA==0) + sig2A*(CA > 0))
			B.scale <- sqrt(tau2B*(CA==0) + sig2B*(CA > 0))
			if(CA == 0 & CB == 0) rho <- 0
			if(CA == 0 & CB > 0) rho <- corrA.BB
			if(CA > 0 & CB == 0) rho <- corrB.AA
			if(CA > 0 & CB > 0) rho <- corr
			if(CHR == 23){
				##means <- cbind(suppressWarnings(log2(nuA[, p]+CA*phiA[, p] + CB*phiAx[, p])), suppressWarnings(log2(nuB[, p]+CB*phiB[, p] + CA*phiBx[, p])))
				meanA <- suppressWarnings(log2(nuA+CA*phiA + CB*phiAX))
				meanB <- suppressWarnings(log2(nuB+CB*phiB + CA*phiBX))				
			} else{
				##means <- cbind(suppressWarnings(log2(nuA+CA*phiA)), suppressWarnings(log2(nuB+CB*phiB)))
				meanA <- suppressWarnings(log2(nuA+CA*phiA))
				meanB <- suppressWarnings(log2(nuB+CB*phiB))
				covs <- rho*A.scale*B.scale
				A.scale2 <- A.scale^2
				B.scale2 <- B.scale^2
			}
			Q.x.y <- 1/(1-rho^2)*(((lA - meanA)/A.scale)^2 + ((lB - meanB)/B.scale)^2 - 2*rho*((lA - meanA)*(lB - meanB))/(A.scale*B.scale))
			f.x.y <- 1/(2*pi*A.scale*B.scale*sqrt(1-rho^2))*exp(-0.5*Q.x.y)
			emit[, , counter] <- f.x.y
			counter <- counter+1
		}
	}
	priorProb <- cnOptions$prior.prob
	homDel <- priorProb[1]*emit[, , 1]
	hemDel <- priorProb[2]*emit[, , c(2, 3)] # + priorProb[3]*emit[, c(4, 5, 6)] + priorProb[4]*emit[, c(7:10)]
	norm <- priorProb[3]*emit[, , 4:6]
	amp <- priorProb[4]*emit[, , 7:10]
	##sum over the different combinations within each copy number state
	hemDel <- hemDel[, , 1] + hemDel[, , 2]
	norm <- norm[, , 1] + norm[, , 2] + norm[, , 3]
	amp <- amp[, , 1] + amp[, , 2] + amp[ , , 3] + amp[, , 4]
	total <- homDel + hemDel + norm + amp
	homDel <- homDel/total
	hemDel <- hemDel/total
	norm <- norm/total
	amp <- amp/total
	tmp.objects$posteriorProb <- list(hemDel=hemDel, norm=norm, amp=amp)
	##envir[["posteriorProb"]] <- list(hemDel=hemDel, norm=norm, amp=amp)
	posteriorProb <- array(NA, dim=c(nrow(A), ncol(A), 4))
	posteriorProb[, , 1] <- homDel
	posteriorProb[, , 2] <- hemDel
	posteriorProb[, , 3] <- norm
	posteriorProb[, , 4] <- amp
	return(list(tmp.objects, posteriorProb))
}



biasAdj <- function(object, cnOptions, tmp.objects){
	gender <- object$gender
	CHR <- unique(chromosome(object))
	I <- isSnp(object)
	A <- A(object)[I, ]
	normal <- tmp.objects[["normal"]][I, ]
	results <- posteriorProbability.snps(object=object, cnOptions=cnOptions, tmp.objects=tmp.objects)
	posteriorProb <- results[[2]]
	tmp.objects <- results[[1]]
	mostLikelyState <- apply(posteriorProb, c(1, 2), function(x) order(x, decreasing=TRUE)[1])
	if(CHR == 23){
		##so state index 3 is the most likely state for men and women
		mostLikelyState[, gender==1] <- mostLikelyState[, gender==1] + 1
	}
	proportionSamplesAltered <- rowMeans(mostLikelyState != 3)
	ii <- proportionSamplesAltered < 0.8 & proportionSamplesAltered > 0.01
	##  only exclude observations from one tail, depending on
	##  whether more are up or down
	moreup <- rowSums(mostLikelyState > 3) > rowSums(mostLikelyState < 3) ##3 is normal
	## if equal, which points have a high posterior probability of altered copy number.  drop those.
	NORM <- matrix(FALSE, nrow(A), ncol(A))
	##NORM[proportionSamplesAltered > 0.8, ] <- FALSE
	ratioUp <- posteriorProb[, , 4]/posteriorProb[, , 3]
	NORM[ii & moreup, ] <- ratioUp[moreup & ii] < 1  ##normal more likely
	ratioDown <- posteriorProb[, , 2]/posteriorProb[, , 3]
	NORM[ii & !moreup, ] <- ratioDown[!moreup & ii] < 1  ##normal more likely
	normal <- NORM*normal
	tmp <- tmp.objects[["normal"]]
	tmp[I, ] <- normal
	tmp.objects[["normal"]] <- tmp
	return(tmp.objects)
}


##biasAdjNP <- function(plateIndex, envir, priorProb){
biasAdjNP <- function(object, cnOptions, tmp.objects){
	batch <- unique(object$batch)
	normalNP <- tmp.objects[["normal"]][!isSnp(object), ]
	CHR <- unique(chromosome(object))
	A <- A(object)[!isSnp(object), ]
	sig2A <- getParam(object, "sig2A", batch)
	gender <- object$gender
	##Assume that on the log-scale, that the background variance is the same...
	tau2A <- sig2A
	nuA <- getParam(object, "nuA", batch)
	phiA <- getParam(object, "phiA", batch)
	prior.prob <- cnOptions$prior.prob
	emit <- array(NA, dim=c(nrow(A), ncol(A), 4))##SNPs x sample x 'truth'	
	lT <- log2(A)
	I <- isSnp(object)
	counter <- 1 ##state counter
##	for(CT in 0:3){
##		sds <- sqrt(tau2A[I]*(CT==0) + sig2A[I]*(CT > 0))
##		means <- suppressWarnings(log2(nuA[I]+CT*phiA[I]))
##		tmp <- dnorm(lT, mean=means, sd=sds)
##		emit[, , counter] <- tmp
##		counter <- counter+1
##	}
##	mostLikelyState <- apply(emit, c(1, 2), function(x) order(x, decreasing=TRUE)[1])
	counter <- 1
	for(CT in c(0,1,2,2.5)){
		sds <- sqrt(tau2A[I]*(CT==0) + sig2A[I]*(CT > 0))
		means <- suppressWarnings(log2(nuA[I]+CT*phiA[I]))
		tmp <- dnorm(lT, mean=means, sd=sds)
		emit[, , counter] <- tmp
		counter <- counter+1
	}
	mostLikelyState <- apply(emit, c(1, 2), function(x) order(x, decreasing=TRUE)[1])
	
	if(CHR == 23){
		## the state index for male on chromosome 23  is 2
		## add 1 so that the state index is 3 for 'normal' state
		mostLikelyState[, gender=="male"] <- mostLikelyState[, gender==1] + 1
	}
	tmp3 <- mostLikelyState != 3
	##Those near 1 have NaNs for nu and phi.  this occurs by NaNs in the muA[,, "A"] or muA[, , "B"] for X chromosome
	proportionSamplesAltered <- rowMeans(tmp3)##prop normal
	ii <- proportionSamplesAltered < 0.75
	moreup <- rowSums(mostLikelyState > 3) > rowSums(mostLikelyState < 3)
	notUp <-  mostLikelyState[ii & moreup, ] <= 3
	notDown <- mostLikelyState[ii & !moreup, ] >= 3
	NORM <- matrix(TRUE, nrow(A), ncol(A))
	NORM[ii & moreup, ] <- notUp
	NORM[ii & !moreup, ] <- notDown
	normalNP <- normalNP*NORM

	##flagAltered <- which(proportionSamplesAltered > 0.5)
	##envir[["flagAlteredNP"]] <- flagAltered
	normal <- tmp.objects[["normal"]]
	normal[!isSnp(object), ] <- normalNP
	tmp.objects[["normal"]] <- normal
	return(tmp.objects)
}


getParams <- function(object, batch){
	batch <- unique(object$batch)
	if(length(batch) > 1) stop("batch variable not unique")		
	nuA <- as.numeric(fData(object)[, match(paste("nuA", batch, sep="_"), fvarLabels(object))])
	nuB <- as.numeric(fData(object)[, match(paste("nuB", batch, sep="_"), fvarLabels(object))])	
	phiA <- as.numeric(fData(object)[, match(paste("phiA", batch, sep="_"), fvarLabels(object))])
	phiB <- as.numeric(fData(object)[, match(paste("phiB", batch, sep="_"), fvarLabels(object))])	
	tau2A <- as.numeric(fData(object)[, match(paste("tau2A", batch, sep="_"), fvarLabels(object))])
	tau2B <- as.numeric(fData(object)[, match(paste("tau2B", batch, sep="_"), fvarLabels(object))])
	sig2A <- as.numeric(fData(object)[, match(paste("sig2A", batch, sep="_"), fvarLabels(object))])
	sig2B <- as.numeric(fData(object)[, match(paste("sig2B", batch, sep="_"), fvarLabels(object))])
	corrA.BB <- as.numeric(fData(object)[, match(paste("corrA.BB", batch, sep="_"), fvarLabels(object))])
	corrB.AA <- as.numeric(fData(object)[, match(paste("corrB.AA", batch, sep="_"), fvarLabels(object))])
	corr <- as.numeric(fData(object)[, match(paste("corr", batch, sep="_"), fvarLabels(object))])
	params <- list(nuA=nuA,
		       nuB=nuB,
		       phiA=phiA,
		       phiB=phiB,
		       tau2A=tau2A,
		       tau2B=tau2B,
		       sig2A=sig2A,
		       sig2B=sig2B,
		       corrA.BB=corrA.BB,
		       corrB.AA=corrB.AA,
		       corr=corr)
	return(params)
}


## Constrain nu and phi to positive values
thresholdModelParams <- function(object, cnOptions){
	MIN.NU <- cnOptions$MIN.NU
	MIN.PHI <- cnOptions$MIN.PHI
	batch <- unique(object$batch)
	nuA <- getParam(object, "nuA", batch)
	nuA[nuA < MIN.NU] <- MIN.NU
	object <- pr(object, "nuA", batch, nuA)
	nuB <- getParam(object, "nuB", batch)
	if(!all(is.na(nuB))){
		nuB[nuB < MIN.NU] <- MIN.NU
		object <- pr(object, "nuB", batch, nuB)
	}
	phiA <- getParam(object, "phiA", batch)
	phiA[phiA < MIN.PHI] <- MIN.PHI
	object <- pr(object, "phiA", batch, phiA)
	phiB <- getParam(object, "phiB", batch)
	if(!all(is.na(phiB))){
		phiB[phiB < MIN.PHI] <- MIN.PHI
		object <- pr(object, "phiB", batch, phiB)		
	}
	phiAX <- as.numeric(getParam(object, "phiAX", batch))	
	if(!all(is.na(phiAX))){
		phiAX[phiAX < MIN.PHI] <- MIN.PHI
		object <- pr(object, "phiAX", batch, phiAX)		
	}
	phiBX <- as.numeric(getParam(object, "phiBX", batch))
	if(!all(is.na(phiBX))){
		phiBX[phiBX < MIN.PHI] <- MIN.PHI
		object <- pr(object, "phiBX", batch, phiBX)			
	}
	return(object)
}

##computeCopynumber.SnpSuperSet <- function(object, cnOptions){
####	use.ff <- cnOptions[["use.ff"]]
####	if(!use.ff){
####		object <- as(object, "CrlmmSet")
####	} else	object <- as(object, "CrlmmSetFF")
##	bias.adj <- cnOptions[["bias.adj"]]
##	##must be FALSE to initialize parameters
##	cnOptions[["bias.adj"]] <- FALSE
##	## Add linear model parameters to the CrlmmSet object
##	featureData(object) <- lm.parameters(object, cnOptions)
##	if(!isValidCdfName(annotation(object))) stop(annotation(object), " not supported.")
##	object <- as(object, "CNSet")
##	object <- computeCopynumber.CNSet(object, cnOptions)
##	if(bias.adj==TRUE){## run a second time
##		object <- computeCopynumber.CNSet(object, cnOptions)
##	}
##	return(object)
##}



## cite metaArray: Choi/Ghosh
crlmm_poe <- function(mat, cl, threshold=0.00001, every=100,use.mad=TRUE) {
  mat <- as.matrix(mat)
  nc <- ncol(mat)
  nr <- nrow(mat)
  if(all(is.null(cl))) cl <- rep(0,dim(mat)[2])
  cat("Number of Samples:", nc, "\n")
  cat("Number of Genes:", nr, "\n")
  cat("This model assumes that the samples are centered and scaled.\n")
  new.mat <- matrix(0,nr,nc)
  med.expr <- apply(mat,1,median,na.rm=TRUE)
  new.mat <- sweep(mat,1,med.expr)
  for(i in 1:nr) {
    if(sum(is.na(as.numeric(mat[i,]))) > .25 * nc ) stop("More than 25% missing values for gene", i, "\n")
    zvec <- crlmm_fit.em(as.numeric(mat[i,]), cl, threshold=threshold, use.mad=use.mad) 
    new.mat[i,] <- zvec$expr
    if(i%%every==0) cat(i, "genes fitted\n")   
  }
  rownames(new.mat) <- rownames(mat)
  colnames(new.mat) <- colnames(mat)
  return(list(data=new.mat))
}

crlmm_fit.em <- function(x, cl, threshold=1e-6, use.mad=TRUE){
	sup <- sum(cl==0) > 0 && sum(cl==1) > 0
	x <- as.numeric(x)
	len <- length(x)
	z <- rep(0,length(x))
	log.lik <- 1000
	lik.rec <- NULL
	num.iter <- 0
	a <- min(x,na.rm=TRUE); b <- max(x,na.rm=TRUE)
	err <- err.old <- 1
	if(!sup) {
		z <- find.init(x)
		Pi <- mean(z)
		mu <- sum((1-z)*x) / sum(1-z)
		sigmasq <- sum((1-z)*(x-mu)^2) / sum(1-z)
		tt <- len
		while(err > threshold) {
			num.iter <- num.iter + 1
			log.lik.old <- log.lik
			err.old <- err
			## E Step
			est.u <- dunif(x,a,b)
			est.u.p <- Pi * est.u
			est.n <- dnorm(x,mu,sqrt(sigmasq))
			est.n.p <- (1-Pi) * est.n
			z <- est.u.p / (est.n.p + est.u.p)

			if(any(is.na(z))) stop("NA occurred in imputation\n")
			## M Step
			mu <- sum((1-z)*x) / sum(1-z)
			sigmasq <- sum((1-z)*((x-mu)^2)) / sum(1-z)
			Pi <- sum(z) / len    
			sgn.z <- ifelse(x < mu, -1, 1)

			## Likelihood
			est.u <- dunif(x,a,b)
			est.u.p <- Pi * est.u
			est.n <- dnorm(x,mu,sqrt(sigmasq))
			est.n.p <- (1-Pi) * est.n   
			log.lik <- sum(log(est.u.p + est.n.p))
			err <- abs(log.lik.old - log.lik)
			if(num.iter != 1) lik.rec[num.iter-1] <- log.lik
		}    
	}  else {
		tt <- sum(cl==0)
		z[cl==0] <- runif(tt,0,1)
		Pi <- mean(z[cl==0])
		mu <- sum((1-z)*x) / sum(1-z)
		sigmasq <- sum((1-z)*(x-mu)^2) / sum(1-z)

		while(err > threshold) {
			num.iter <- num.iter + 1
			log.lik.old <- log.lik
			err.old <- err

			## E Step
			est.u <- dunif(x,a,b)
			est.u.p <- Pi * est.u
			est.n <- dnorm(x,mu,sqrt(sigmasq))
			est.n.p <- (1-Pi) * est.n
			z <- est.u.p / (est.n.p + est.u.p)
			z[cl==1] <- 0
			if(any(is.na(z))) stop("NA occurred in imputation\n")
			## M Step
			mu <- sum((1-z)*x) / sum(1-z)
			sigmasq <- sum((1-z)*((x-mu)^2)) / sum(1-z)
			Pi <- mean(z[cl==0])    
			sgn.z <- ifelse(x < mu, -1, 1)

			## Likelihood
			est.u <- dunif(x,a,b)
			est.u.p <- Pi * est.u
			est.n <- dnorm(x,mu,sqrt(sigmasq))
			est.n.p <- (1-Pi) * est.n   
			log.lik <- sum(log(est.u.p[cl==0] + est.n.p[cl==0])) + sum(log(est.n[cl==1]))
			err <- abs(log.lik.old - log.lik)
			if(num.iter != 1) lik.rec[num.iter-1] <- log.lik
		}    
	}
	est.u.p <- Pi * dunif(x,a,b) 
	est.n.p <- (1-Pi) * dnorm(x,mu,sqrt(sigmasq))
	est.u.mu <- Pi * dunif(mu,a,b)
	est.n.mu <- (1-Pi) * dnorm(mu,mu,sqrt(sigmasq))
	z0 <- est.u.p / (est.n.p + est.u.p)
	zmu <- est.u.mu / (est.u.mu + est.n.mu)
	sgn.z0 <- ifelse(x < mu, -1, 1) 
					#loc <- (max(lik.rec) != lik.rec[length(lik.rec)])
	expr <- rep(0, len)
	expr <- (z0 - zmu) * sgn.z0
	return(list(expr=expr, a=a, b=b, sigmasq=sigmasq, mu=mu, Pi=Pi, lik.rec=lik.rec))
}



.copyNumber <- function(object){
	I <- isSnp(object)
	CA <- CA(object)
	CB <- CB(object)
	CN <- CA + CB
	##For nonpolymorphic probes, CA is the total copy number
	CN[!I, ] <- CA(object)[!I, ]
	CN
}



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


##initializeCNFiles <- function(cnOpts){
##	load.it <- cnOptions[["load.it"]]
##	outdir <- cnOptions[["outdir"]]
##	snprmaFile <- cnOptions[["snprmaFile"]]
##	protocolFile <- cnOptions[["protocolFile"]]
##	CAFile <- cnOpts[["CAFile"]]
##	CBFile <- cnOpts[["CBFile"]]
##	path <- system.file("extdata", package=paste(cnOptions[["cdfName"]], "Crlmm", sep=""))
##	load(file.path(path, "snpProbes.rda"))
##	snpProbes <- get("snpProbes")
##	load(file.path(path, "cnProbes.rda"))
##	cnProbes <- get("cnProbes")		
##	load(snprmaFile)
##	res <- get("res")
##	load(protocolFile)
##	if(isPackageLoaded("ff")){
##		if(file.exists(CAFile)){
##			load(CAFile)
##			if(is.null(dim(CA)) | !all(dim(CA) == dim(A))) {
##				unlink(CAFile)
##				CA <- initializeBigMatrix(nr, nc)
##				save(CA, file=CAFile)
##				CA[,] <- NA
##
##			}
##		} else {
##			CA <- initializeBigMatrix(nr, nc)
##			save(CA, file=CAFile)
##			CA[,] <- NA
##		}
##		if(file.exists(CBFile)){
##			load(CBFile)
##			if(is.null(dim(CB)) | !all(dim(CB) == dim(A))) {			     
##				CB <- initializeBigMatrix(nr, nc)
##				save(CB, file=CBFile)					
##				CB[,] <- NA
##			}
##		} else{
##			CB <- initializeBigMatrix(nr, nc)
##			CB[,] <- NA
##			save(CB, file=CBFile)
##		}
##		open(CA)
##		open(CB)
##		colnames(CA) <- colnames(CB) <- sampleNames(protocoldata)
##	}
##	close(CA); close(CB)
##
##
##}
##collect <- function(cnOptions, CHR, PLATE){
##	snprmaFile <- cnOptions[["snprmaFile"]]
##	cnFile <- cnOptions[["cnFile"]]
##	cnrmaFile <- cnOptions[["cnrmaFile"]]
##	snprmaFile <- cnOptions[["snprmaFile"]]
##	callsFile <- cnOptions[["callsFile"]]
##	confsFile <- cnOptions[["confsFile"]]
##	protocolFile <- cnOptions[["protocolFile"]]
##	AFile <- cnOptions[["AFile"]]
##	BFile <- cnOptions[["BFile"]]
##	CAFile <- cnOptions[["CAFile"]]
##	CBFile <- cnOptions[["CBFile"]]
##	if(isPackageLoaded("ff")){
##		message("collecting ff objects to create an FFSet instance")
##		load(protocolFile)
##		load(snprmaFile)
##		load(file.path(cnOptions[["outdir"]], "cnParams.rda"))
##		res <- get("res")
##		sampleStats <- data.frame(SKW=res$SKW,
##					  SNR=res$SNR,
##					  ##mixtureParams=res$mixtureParams,
##					  gender=res$gender,
##					  batch=cnOptions[["batch"]])##)
##		pD <- new("AnnotatedDataFrame",
##			  data=sampleStats,
##			  varMetadata=data.frame(labelDescription=colnames(sampleStats)))	
##		load(file.path(cnOptions[["outdir"]], "featureDataFF.rda"))
##		open(featureDataFF)
##		load(AFile); open(A)
##		load(BFile); open(B)
##		load(CAFile); open(CA)
##		load(CBFile); open(CB)
##		load(callsFile); open(calls)
##		load(confsFile); open(confs)
##		sampleNames(pD) <- sampleNames(protocoldata)
##
##		if(any(colnames(CA) != colnames(A)) | any(colnames(CA) != sampleNames(protocoldata)))
##			colnames(calls) <- colnames(confs) <- colnames(CA) <- colnames(CB) <- colnames(B) <- colnames(A) <- sampleNames(protocoldata)
####		if(any(rownames(CA) != fns))
####			rownames(CA) <- rownames(CB) <- rownames(B) <- rownames(calls) <- rownames(confs) <- rownames(A)
##		ffSet <- new("FFSet",
##			     alleleA=A,
##			     alleleB=B,
##			     call=calls,
##			     callProbability=confs,
##			     CA=CA,
##			     CB=CB,
##			     phenoData=pD,
##			     protocolData=protocoldata,
##			     ##annotation=cnOptions[["cdfName"]],
##			     genomeAnnotation=featureDataFF,
##			     linearModelParam=cnParams)
##		##featureNames(ffSet) <- rownames(featureDataFF)
##		##sampleNames(ffSet) <- sampleNames(protocoldata)
##		return(ffSet)
##	} else{
##		if(missing(PLATE)) stop("must specify PLATE")
##		load(paste(cnOptions[["cnFile"]], "_", PLATE, "_", CHR, ".rda", sep=""))
##		return(cnSet)
##	}
##}


lmPlot <- function(object, plate, fill,...){
	stopifnot(nrow(object) == 1)
	Ai <- split(A(object)[object$batch==plate], snpCall(object)[object$batch==plate])
	Bi <- split(B(object)[object$batch==plate], snpCall(object)[object$batch==plate])
	##xlim <- ylim <- c(9,12)
	boxplot(Ai, xaxt="n", ylab=expression(I[A]), main=featureNames(object), ...)
	legend("topleft", bty="n", legend=c("AA", "AB", "BB"), fill=fill)
	nuA <- getParam(object, "nuA", plate)
	phiA <- getParam(object, "phiA", plate)
	segments(3, nuA, 1, nuA+2*phiA, lwd=2, col="blue")
	boxplot(Bi, xaxt="n", ylab=expression(I[B]), ...)
	nuB <- getParam(object, "nuB", plate)
	phiB <- getParam(object, "phiB", plate)	
	segments(1, nuB,  3.0, nuB+2*phiB, lwd=2, col="blue")
	axis(1, at=1:3, labels=c("AA", "AB", "BB"))
}
