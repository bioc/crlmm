##---------------------------------------------------------------------------
##---------------------------------------------------------------------------
getProtocolData.Affy <- function(filenames){
	scanDates <- data.frame(ScanDate=sapply(filenames, celfileDate))
	rownames(scanDates) <- basename(rownames(scanDates))
	protocoldata <- new("AnnotatedDataFrame",
			    data=scanDates,
			    varMetadata=data.frame(labelDescription=colnames(scanDates),
			                           row.names=colnames(scanDates)))
	return(protocoldata)
}
getFeatureData <- function(cdfName, copynumber=FALSE){
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
	snpProbes <- get("snpProbes")
	if(copynumber){
		load(file.path(path, "cnProbes.rda"))
		cnProbes <- get("cnProbes")
		snpIndex <- seq(along=gns)
		npIndex <- seq(along=rownames(cnProbes)) + max(snpIndex)
		featurenames <- c(gns, rownames(cnProbes))
	} else featurenames <- gns
	fvarlabels=c("chromosome", "position", "isSnp")
	M <- matrix(NA, length(featurenames), 3, dimnames=list(featurenames, fvarlabels))
	index <- match(rownames(snpProbes), rownames(M)) #only snp probes in M get assigned position
	M[index, "position"] <- snpProbes[, grep("pos", colnames(snpProbes))]
	M[index, "chromosome"] <- chromosome2integer(snpProbes[, grep("chr", colnames(snpProbes))])
	M[index, "isSnp"] <- 1L
	##index <- which(is.na(M[, "isSnp"]))
	##M[index, "isSnp"] <- 1L
	if(copynumber){
		index <- match(rownames(cnProbes), rownames(M)) #only snp probes in M get assigned position
		M[index, "position"] <- cnProbes[, grep("pos", colnames(cnProbes))]
		M[index, "chromosome"] <- chromosome2integer(cnProbes[, grep("chr", colnames(cnProbes))])
		M[index, "isSnp"] <- 0L
	}
	##A few of the snpProbes do not match -- I think it is chromosome Y.
	M[is.na(M[, "isSnp"]), "isSnp"] <- 1L
	M <- data.frame(M)
	return(new("AnnotatedDataFrame", data=M))
}
getFeatureData.Affy <- getFeatureData

construct <- function(filenames,
		      cdfName,
		      copynumber=TRUE,
		      sns, verbose=TRUE, batch, fns){
	if(!missing(batch)){
		stopifnot(length(batch) == length(sns))
	}
	if(missing(sns) & missing(filenames)) stop("one of filenames or samplenames (sns) must be provided")
	if(verbose) message("Initializing container for copy number estimation")
	featureData <- getFeatureData(cdfName, copynumber=copynumber)
	if(!missing(fns)){
		warning("subsetting the featureData can cause the snprma object and CNSet object to be misaligned. This problem is fixed in devel as we match the names prior to assigning values from snprma")
		index <- match(fns, featureNames(featureData))
		if(all(is.na(index))) stop("fns not in featureNames")
		featureData <- featureData[index, ]
	}
	nr <- nrow(featureData); nc <- length(sns)
	cnSet <- new("CNSet",
		     alleleA=initializeBigMatrix(name="A", nr, nc),
		     alleleB=initializeBigMatrix(name="B", nr, nc),
		     call=initializeBigMatrix(name="call", nr, nc),
		     callProbability=initializeBigMatrix(name="callPr", nr,nc),
		     annotation=cdfName,
		     batch=batch)
	sampleNames(cnSet) <- sns
	if(!missing(filenames)){
		if(missing(sns)) sns <- basename(filenames)
		protocolData <- getProtocolData.Affy(filenames)
	} else{
		protocolData <- annotatedDataFrameFrom(A(cnSet), byrow=FALSE)
	}
	rownames(pData(protocolData)) <- sns
	protocolData(cnSet) <- protocolData
	featureData(cnSet) <- featureData
	featureNames(cnSet) <- featureNames(featureData)
	pd <- data.frame(matrix(NA, nc, 3), row.names=sns)
	colnames(pd)=c("SKW", "SNR", "gender")
	phenoData(cnSet) <- new("AnnotatedDataFrame", data=pd)
	gns <- getVarInEnv("gns")
	stopifnot(identical(gns, featureNames(cnSet)[1:length(gns)]))
	return(cnSet)
}

genotype <- function(filenames,
		       cdfName,
		       batch,
		       mixtureSampleSize=10^5,
		       eps=0.1,
		       verbose=TRUE,
		       seed=1,
		       sns,
		       probs=rep(1/3, 3),
		       DF=6,
		       SNRMin=5,
		       recallMin=10,
		       recallRegMin=1000,
		       gender=NULL,
		       returnParams=TRUE,
		       badSNP=0.7){
	is.lds <- ifelse(isPackageLoaded("ff"), TRUE, FALSE)
	if(missing(cdfName)) stop("must specify cdfName")
	if(!isValidCdfName(cdfName)) stop("cdfName not valid.  see validCdfNames")
	if(missing(sns)) sns <- basename(filenames)
	callSet <- construct(filenames=filenames,
			     cdfName=cdfName,
			     copynumber=TRUE,
			     sns=sns,
			     verbose=verbose,
			     batch=batch)
	FUN <- ifelse(is.lds, "snprma2", "snprma")
	snprmaFxn <- function(FUN,...){
		switch(FUN,
		       snprma=snprma(...),
		       snprma2=snprma2(...))
	}
	snprmaRes <- snprmaFxn(FUN,
			       filenames=filenames,
			       mixtureSampleSize=mixtureSampleSize,
			       fitMixture=TRUE,
			       eps=eps,
			       verbose=verbose,
			       seed=seed,
			       cdfName=cdfName,
			       sns=sns)
	gns <- snprmaRes[["gns"]]
	snp.I <- isSnp(callSet)
	is.snp <- which(snp.I)
	snp.index <- match(featureNames(callSet)[is.snp], gns)
	stopifnot(identical(featureNames(callSet)[is.snp], gns[snp.index]))
##	is.snp <- isSnp(callSet)
##	snp.index <- which(is.snp)
	if(is.lds) open(callSet)
	mixtureParams <- matrix(NA, 4, length(filenames))
	##message("Saving snprmaRes file")
	##save(snprmaRes, file=file.path(outdir, "snprmaRes.rda"))
	if(verbose) message("Finished preprocessing.")
	gns <- snprmaRes[["gns"]]
	snp.I <- isSnp(callSet)
	is.snp <- which(snp.I)
	snp.index <- match(featureNames(callSet)[is.snp], gns)
	stopifnot(identical(featureNames(callSet)[is.snp], gns[snp.index]))
	if(is.lds) open(callSet)
	mixtureParams <- matrix(NA, 4, length(filenames))
	if(is.lds){
		open(snprmaRes[["A"]])
		open(snprmaRes[["B"]])
		open(snprmaRes[["SNR"]])
		open(snprmaRes[["mixtureParams"]])
		bb <- getOption("ffbatchbytes")
		ffcolapply(A(callSet)[is.snp, i1:i2] <- snprmaRes[["A"]][snp.index, i1:i2], X=snprmaRes[["A"]],
			   BATCHBYTES=bb)
		ffcolapply(B(callSet)[is.snp, i1:i2] <- snprmaRes[["B"]][snp.index, i1:i2], X=snprmaRes[["B"]],
			   BATCHBYTES=bb)
	} else{
		A(callSet)[is.snp, ] <- snprmaRes[["A"]][snp.index, ]
		B(callSet)[is.snp, ] <- snprmaRes[["B"]][snp.index, ]
	}
	pData(callSet)$SKW <- snprmaRes[["SKW"]]
	pData(callSet)$SNR <- snprmaRes[["SNR"]]
	mixtureParams <- snprmaRes$mixtureParams
	np.index <- which(!snp.I)
	if(verbose) message("Normalizing nonpolymorphic markers")
	FUN <- ifelse(is.lds, "cnrma2", "cnrma")
	## main purpose is to update 'alleleA'
	cnrmaFxn <- function(FUN,...){
		switch(FUN,
		       cnrma=cnrma(...),
		       cnrma2=cnrma2(...))
	}
	## consider passing only A for NPs.
	if(verbose) message("Quantile normalizing nonpolymorphic markers")
	AA <- cnrmaFxn(FUN, A=A(callSet),
		       filenames=filenames,
		       row.names=featureNames(callSet)[np.index],
		       cdfName=cdfName,
		       sns=sns,
		       seed=seed,
		       verbose=verbose)
	if(!is.lds) A(callSet) <- AA
	rm(AA)
	FUN <- ifelse(is.lds, "crlmmGT2", "crlmmGT")
	## genotyping
	crlmmGTfxn <- function(FUN,...){
		switch(FUN,
		       crlmmGT2=crlmmGT2(...),
		       crlmmGT=crlmmGT(...))
	}
	if(verbose) message("Running crlmmGT2")
	tmp <- crlmmGTfxn(FUN,
			  A=snprmaRes[["A"]],
			  B=snprmaRes[["B"]],
			  SNR=snprmaRes[["SNR"]],
			  mixtureParams=snprmaRes[["mixtureParams"]],
			  cdfName=cdfName,
			  row.names=NULL,
			  col.names=sampleNames(callSet),
			  SNRMin=SNRMin,
			  recallMin=recallMin,
			  recallRegMin=recallRegMin,
			  gender=gender,
			  verbose=verbose,
			  returnParams=returnParams,
			  badSNP=badSNP)
	if(verbose) message("Genotyping finished.  Updating container with genotype calls and confidence scores.")
	if(is.lds){
		open(tmp[["calls"]])
		open(tmp[["confs"]])
		ffcolapply(snpCall(callSet)[is.snp, i1:i2] <- tmp[["calls"]][snp.index, i1:i2], X=tmp[["calls"]],
			   BATCHBYTES=bb)
		ffcolapply(snpCallProbability(callSet)[is.snp, i1:i2] <- tmp[["confs"]][snp.index, i1:i2], X=tmp[["confs"]],
			   BATCHBYTES=bb)
		close(tmp[["calls"]])
		close(tmp[["confs"]])
	} else {
		calls(callSet)[is.snp, ] <- tmp[["calls"]][snp.index, ]
		snpCallProbability(callSet)[is.snp, ] <- tmp[["confs"]][snp.index, ]
	}
	message("Finished updating.  Cleaning up.")
	callSet$gender <- tmp$gender
	if(is.lds){
		delete(snprmaRes[["A"]])
		delete(snprmaRes[["B"]])
		delete(snprmaRes[["mixtureParams"]])
		rm(snprmaRes)
	}
	close(callSet)
	return(callSet)
}

genotype2 <- function(){
	.Defunct(msg="The genotype2 function has been deprecated. The function genotype should be used instead.  genotype will support large data using ff provided that the ff package is loaded.")
}
genotypeLD <- genotype2

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

shrink <- function(x, Ns, DF.PRIOR){
	DF <- Ns-1
	DF[DF < 1] <- 1
	x.0 <- apply(x, 2, median, na.rm=TRUE)
	x <- (x*DF + x.0*DF.PRIOR)/(DF.PRIOR + DF)
	for(j in 1:ncol(x)) x[is.na(x[, j]), j] <- x.0[j]
	return(x)
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

dqrlsWrapper <- function(x, y, wts, tol=1e-7){
	n <- NROW(y)
	p <- ncol(x)
	ny <- NCOL(y)
	.Fortran("dqrls", qr=x*wts, n=n, p=p, y=y * wts, ny=ny,
		 tol=as.double(tol), coefficients=mat.or.vec(p, ny),
		 residuals=y, effects=mat.or.vec(n, ny),
		 rank=integer(1L), pivot=1L:p, qraux=double(p),
		 work=double(2 * p), PACKAGE="base")[["coefficients"]]
}


fit.wls <- function(NN, sigma, allele, Y, autosome, X){
	Np <- NN
	Np[Np < 1] <- 1
	W <- (sigma/sqrt(Np))^-1
	Ystar <- Y*W
	complete <- which(rowSums(is.na(W)) == 0 & rowSums(is.na(Ystar)) == 0)
	if(missing(allele)) stop("must specify allele")
	if(autosome & missing(X)){
		if(allele == "A") X <- cbind(1, 2:0)
		if(allele == "B") X <- cbind(1, 0:2)
	}
	if(!autosome & missing(X)){
		if(allele == "A") X <- cbind(1, c(1, 0, 2, 1, 0), c(0, 1, 0, 1, 2))
		if(allele == "B") X <- cbind(1, c(0, 1, 0, 1, 2), c(1, 0, 2, 1, 0))
	}
	betahat <- matrix(NA, ncol(X), nrow(Ystar))
	ww <- rep(1, ncol(Ystar))
	for(i in complete){
		betahat[, i] <- dqrlsWrapper(W[i, ] * X, Ystar[i, ], ww)
		##ssr <- sum((Ystar[i, ] - matrix(Xstar[, i], nrow(X), ncol(X)) %*% matrix(betahat[, i], ncol(X), 1))^2)
	}
	return(betahat)
}

shrinkGenotypeSummaries <- function(strata, index.list, object, MIN.OBS, MIN.SAMPLES, DF.PRIOR,
				    verbose, is.lds){
	if(is.lds) {physical <- get("physical"); open(object)}
	if(verbose) message("      Probe stratum ", strata, " of ", length(index.list))
	marker.index <- index.list[[strata]]
	batches <- split(seq_along(batch(object)), as.character(batch(object)))
	batches <- batches[sapply(batches, length) >= MIN.SAMPLES]
	batchnames <- batchNames(object)
	N.AA <- as.matrix(N.AA(object)[marker.index, ])
	N.AB <- as.matrix(N.AB(object)[marker.index, ])
	N.BB <- as.matrix(N.BB(object)[marker.index, ])
	medianA.AA <- as.matrix(medianA.AA(object)[marker.index,])
	medianA.AB <- as.matrix(medianA.AB(object)[marker.index,])
	medianA.BB <- as.matrix(medianA.BB(object)[marker.index,])
	medianB.AA <- as.matrix(medianB.AA(object)[marker.index,])
	medianB.AB <- as.matrix(medianB.AB(object)[marker.index,])
	medianB.BB <- as.matrix(medianB.BB(object)[marker.index,])
	madA.AA <- as.matrix(madA.AA(object)[marker.index,])
	madA.AB <- as.matrix(madA.AB(object)[marker.index,])
	madA.BB <- as.matrix(madA.BB(object)[marker.index,])
	madB.AA <- as.matrix(madB.AA(object)[marker.index,])
	madB.AB <- as.matrix(madB.AB(object)[marker.index,])
	madB.BB <- as.matrix(madB.BB(object)[marker.index,])
	medianA <- medianB <- shrink.madB <- shrink.madA <- vector("list", length(batchnames))
	shrink.tau2A.AA <- tau2A.AA <- as.matrix(tau2A.AA(object)[marker.index,])
	shrink.tau2B.BB <- tau2B.BB <- as.matrix(tau2B.BB(object)[marker.index,])
	shrink.tau2A.BB <- tau2A.BB <- as.matrix(tau2A.BB(object)[marker.index,])
	shrink.tau2B.AA <- tau2B.AA <- as.matrix(tau2B.AA(object)[marker.index,])
	shrink.corrAA <- corrAA <- as.matrix(corrAA(object)[marker.index, ])
	shrink.corrAB <- corrAB <- as.matrix(corrAB(object)[marker.index, ])
	shrink.corrBB <- corrBB <- as.matrix(corrBB(object)[marker.index, ])
	flags <- as.matrix(flags(object)[marker.index, ])
	for(k in seq(along=batches)){
		sample.index <- batches[[k]]
		this.batch <- unique(as.character(batch(object)[sample.index]))
		medianA[[k]] <- cbind(medianA.AA[, k], medianA.AB[, k], medianA.BB[, k])
		medianB[[k]] <- cbind(medianB.AA[, k], medianB.AB[, k], medianB.BB[, k])
		madA <- cbind(madA.AA[, k], madA.AB[, k], madA.BB[, k])
		madB <- cbind(madB.AA[, k], madB.AB[, k], madB.BB[, k])
		NN <- cbind(N.AA[, k], N.AB[, k], N.BB[, k])
		##RS: estimate DF.PRIOR
		shrink.madA[[k]] <- shrink(madA, NN, DF.PRIOR)
		shrink.madB[[k]] <- shrink(madB, NN, DF.PRIOR)
		## an estimate of the background variance is the MAD
		## of the log2(allele A) intensities among subjects with
		## genotypes BB
		shrink.tau2A.BB[, k] <- shrink(tau2A.BB[, k, drop=FALSE], NN[, 3], DF.PRIOR)[, drop=FALSE]
		shrink.tau2B.AA[, k] <- shrink(tau2B.AA[, k, drop=FALSE], NN[, 1], DF.PRIOR)[, drop=FALSE]
		## an estimate of the signal variance is the MAD
		## of the log2(allele A) intensities among subjects with
		## genotypes AA
		shrink.tau2A.AA[, k] <- shrink(tau2A.AA[, k, drop=FALSE], NN[, 1], DF.PRIOR)[, drop=FALSE]
		shrink.tau2B.BB[, k] <- shrink(tau2B.BB[, k, drop=FALSE], NN[, 3], DF.PRIOR)[, drop=FALSE]
		cor.AA <- corrAA[, k, drop=FALSE]
		cor.AB <- corrAB[, k, drop=FALSE]
		cor.BB <- corrBB[, k, drop=FALSE]
		shrink.corrAA[, k] <- shrink(cor.AA, NN[, 1], DF.PRIOR)
		shrink.corrAB[, k] <- shrink(cor.AB, NN[, 2], DF.PRIOR)
		shrink.corrBB[, k] <- shrink(cor.BB, NN[, 3], DF.PRIOR)
		##
		##---------------------------------------------------------------------------
		## SNPs that we'll use for imputing location/scale of unobserved genotypes
		##---------------------------------------------------------------------------
		index.complete <- indexComplete(NN, medianA[[k]], medianB[[k]], MIN.OBS)
		if(length(index.complete) == 1){
			if(index.complete == FALSE) return()
		}
		##
		##---------------------------------------------------------------------------
		## Impute sufficient statistics for unobserved genotypes (plate-specific)
		##---------------------------------------------------------------------------
		unobservedAA <- NN[, 1] < MIN.OBS
		unobservedAB <- NN[, 2] < MIN.OBS
		unobservedBB <- NN[, 3] < MIN.OBS
		unobserved.index <- vector("list", 3)
		## AB, BB observed
		unobserved.index[[1]] <- which(unobservedAA & (NN[, 2] >= MIN.OBS & NN[, 3] >= MIN.OBS))
		## AA and BB observed (strange)
		unobserved.index[[2]] <- which(unobservedAB & (NN[, 1] >= MIN.OBS & NN[, 3] >= MIN.OBS))
		## AB and AA observed
		unobserved.index[[3]] <- which(unobservedBB & (NN[, 2] >= MIN.OBS & NN[, 1] >= MIN.OBS))
		res <- imputeCenter(medianA[[k]], medianB[[k]], index.complete, unobserved.index)
		medianA[[k]] <- res[[1]]
		medianB[[k]] <- res[[2]]
		rm(res)
		##the NA's in 'medianA' and 'medianB' are monomorphic if MIN.OBS = 1
		##
		## RS: For Monomorphic SNPs a mixture model may be better
		## RS: Further, we can improve estimation by borrowing strength across batch
		unobserved.index[[1]] <- which(unobservedAA & unobservedAB)
		unobserved.index[[2]] <- which(unobservedBB & unobservedAB)
		unobserved.index[[3]] <- which(unobservedAA & unobservedBB) ## strange
		res <- imputeCentersForMonomorphicSnps(medianA[[k]], medianB[[k]],
						       index.complete,
						       unobserved.index)
		medianA[[k]] <- res[[1]]; medianB[[k]] <- res[[2]]
		rm(res)
		negA <- rowSums(medianA[[k]] < 0) > 0
		negB <- rowSums(medianB[[k]] < 0) > 0
		flags[, k] <- as.integer(rowSums(NN == 0) > 0 | negA | negB)
	}
	flags(object)[marker.index, ] <- flags
	medianA.AA(object)[marker.index, ] <- do.call("cbind", lapply(medianA, function(x) x[, 1]))
	medianA.AB(object)[marker.index, ] <- do.call("cbind", lapply(medianA, function(x) x[, 2]))
	medianA.BB(object)[marker.index, ] <- do.call("cbind", lapply(medianA, function(x) x[, 3]))
	medianB.AA(object)[marker.index, ] <- do.call("cbind", lapply(medianB, function(x) x[, 1]))
	medianB.AB(object)[marker.index, ] <- do.call("cbind", lapply(medianB, function(x) x[, 2]))
	medianB.BB(object)[marker.index, ] <- do.call("cbind", lapply(medianB, function(x) x[, 3]))
	##
	madA.AA(object)[marker.index, ] <- do.call("cbind", lapply(shrink.madA, function(x) x[, 1]))
	madA.AB(object)[marker.index, ] <- do.call("cbind", lapply(shrink.madA, function(x) x[, 2]))
	madA.BB(object)[marker.index, ] <- do.call("cbind", lapply(shrink.madA, function(x) x[, 3]))
	madB.AA(object)[marker.index, ] <- do.call("cbind", lapply(shrink.madB, function(x) x[, 1]))
	madB.AB(object)[marker.index, ] <- do.call("cbind", lapply(shrink.madB, function(x) x[, 2]))
	madB.BB(object)[marker.index, ] <- do.call("cbind", lapply(shrink.madB, function(x) x[, 3]))
	##
	corrAA(object)[marker.index, ] <- shrink.corrAA
	corrAB(object)[marker.index, ] <- shrink.corrAB
	corrBB(object)[marker.index, ] <- shrink.corrBB
	tau2A.AA(object)[marker.index,] <- shrink.tau2A.AA
	tau2A.BB(object)[marker.index,] <- shrink.tau2A.BB
	tau2B.AA(object)[marker.index,] <- shrink.tau2B.AA
	tau2B.BB(object)[marker.index,] <- shrink.tau2B.BB
	if(is.lds) return(TRUE) else return(object)
}



fit.lm1 <- function(strata,
		    index.list,
		    object,
		    MIN.SAMPLES,
		    THR.NU.PHI,
		    MIN.NU,
		    MIN.PHI,
		    verbose, is.lds,
		    CHR.X, ...){
	if(is.lds) {physical <- get("physical"); open(object)}
	if(verbose) message("      Probe stratum ", strata, " of ", length(index.list))
	snps <- index.list[[strata]]
	batches <- split(seq_along(batch(object)), as.character(batch(object)))
	batches <- batches[sapply(batches, length) >= MIN.SAMPLES]
	batchnames <- batchNames(object)
	N.AA <- as.matrix(N.AA(object)[snps, ])
	N.AB <- as.matrix(N.AB(object)[snps, ])
	N.BB <- as.matrix(N.BB(object)[snps, ])
	medianA.AA <- as.matrix(medianA.AA(object)[snps,])
	medianA.AB <- as.matrix(medianA.AB(object)[snps,])
	medianA.BB <- as.matrix(medianA.BB(object)[snps,])
	medianB.AA <- as.matrix(medianB.AA(object)[snps,])
	medianB.AB <- as.matrix(medianB.AB(object)[snps,])
	medianB.BB <- as.matrix(medianB.BB(object)[snps,])
	madA.AA <- as.matrix(madA.AA(object)[snps,])
	madA.AB <- as.matrix(madA.AB(object)[snps,])
	madA.BB <- as.matrix(madA.BB(object)[snps,])
	madB.AA <- as.matrix(madB.AA(object)[snps,])
	madB.AB <- as.matrix(madB.AB(object)[snps,])
	madB.BB <- as.matrix(madB.BB(object)[snps,])
	tau2A.AA <- as.matrix(tau2A.AA(object)[snps,])
	tau2B.BB <- as.matrix(tau2B.BB(object)[snps,])
	tau2A.BB <- as.matrix(tau2A.BB(object)[snps,])
	tau2B.AA <- as.matrix(tau2B.AA(object)[snps,])
	corrAA <- as.matrix(corrAA(object)[snps, ])
	corrAB <- as.matrix(corrAB(object)[snps, ])
	corrBB <- as.matrix(corrBB(object)[snps, ])
	nuA <- as.matrix(nuA(object)[snps, ])
	phiA <- as.matrix(phiA(object)[snps, ])
	nuB <- as.matrix(nuB(object)[snps, ])
	phiB <- as.matrix(phiB(object)[snps, ])
	flags <- as.matrix(flags(object)[snps, ])
	for(k in seq(along=batches)){
		B <- batches[[k]]
		if(length(B) < MIN.SAMPLES) next()
		this.batch <- unique(as.character(batch(object)[B]))
		medianA <- cbind(medianA.AA[, k], medianA.AB[, k], medianA.BB[, k])
		medianB <- cbind(medianB.AA[, k], medianB.AB[, k], medianB.BB[, k])
		madA <- cbind(madA.AA[, k], madA.AB[, k], madA.BB[, k])
		madB <- cbind(madB.AA[, k], madB.AB[, k], madB.BB[, k])
		NN <- cbind(N.AA[, k], N.AB[, k], N.BB[, k])
		## we're regressing on the medians using the standard errors (hence the division by N) as weights
		res <- fit.wls(NN=NN, sigma=madA, allele="A", Y=medianA, autosome=!CHR.X)
		nuA[, k] <- res[1, ]
		phiA[, k] <- res[2, ]
		rm(res)
		res <- fit.wls(NN=NN, sigma=madB, allele="B", Y=medianB, autosome=!CHR.X)##allele="B", Ystar=YB, W=wB, Ns=Ns)
		nuB[, k] <- res[1, ]
		phiB[, k] <- res[2, ]
##		cA[, k] <- matrix((1/phiA[, J]*(A-nuA[, J])), nrow(A), ncol(A))
##		cB[, k] <- matrix((1/phiB[, J]*(B-nuB[, J])), nrow(B), ncol(B))
	}
	if(THR.NU.PHI){
		nuA[nuA < MIN.NU] <- MIN.NU
		nuB[nuB < MIN.NU] <- MIN.NU
		phiA[phiA < MIN.PHI] <- MIN.PHI
		phiB[phiB < MIN.PHI] <- MIN.PHI
	}
	nuA(object)[snps, ] <- nuA
	nuB(object)[snps, ] <- nuB
	phiA(object)[snps, ] <- phiA
	phiB(object)[snps, ] <- phiB
	if(is.lds){
		close(object)
		return(TRUE)
	} else{
		return(object)
	}
}

## nonpolymorphic markers
fit.lm2 <- function(strata,
		    index.list,
		    object,
		    MIN.SAMPLES,
		    THR.NU.PHI,
		    MIN.NU,
		    MIN.PHI,
		    verbose, is.lds, CHR.X, ...){
	if(is.lds) {physical <- get("physical"); open(object)}
	if(verbose) message("      Probe stratum ", strata, " of ", length(index.list))
	marker.index <- index.list[[strata]]
	batches <- split(seq_along(batch(object)), as.character(batch(object)))
	batches <- batches[sapply(batches, length) >= MIN.SAMPLES]

	ii <- isSnp(object) & chromosome(object) < 23 & !is.na(chromosome(object))
	flags <- as.matrix(flags(object)[ii, ])
	fns <- featureNames(object)[ii]
	fns.noflags <- fns[rowSums(flags, na.rm=T) == 0]
	snp.index <- sample(match(fns.noflags, featureNames(object)), 5000)

	nuA.np <- as.matrix(nuA(object)[marker.index, ])
	phiA.np <- as.matrix(phiA(object)[marker.index, ])
	tau2A.AA <- as.matrix(tau2A.AA(object)[marker.index, ])

	nuA.snp <- as.matrix(nuA(object)[snp.index, ])
	nuB.snp <- as.matrix(nuB(object)[snp.index, ])
	phiA.snp <- as.matrix(phiA(object)[snp.index, ])
	phiB.snp <- as.matrix(phiB(object)[snp.index, ])
	medianA.AA <- as.matrix(medianA.AA(object)[snp.index,])
	medianB.BB <- as.matrix(medianB.BB(object)[snp.index,])

	medianA.AA.np <- as.matrix(medianA.AA(object)[marker.index,])
	for(k in seq_along(batches)){
		B <- batches[[k]]
		this.batch <- unique(as.character(batch(object)[B]))
		X <- cbind(1, log2(c(medianA.AA[, k], medianB.BB[, k])))
		Y <- log2(c(phiA.snp[, k], phiB.snp[, k]))
		betahat <- solve(crossprod(X), crossprod(X, Y))
		##crosshyb <- max(median(medianA.AA[, k]) - median(medianA.AA.np[, k]), 0)
		crosshyb <- 0
		X <- cbind(1, log2(medianA.AA.np[, k] + crosshyb))
		logPhiT <- X %*% betahat
		phiA.np[, k] <- 2^(logPhiT)
		nuA.np[, k] <- medianA.AA.np[,k]-2*phiA.np[, k]
##		cA[, k] <- 1/phiA.np[, J] * (A.np - nuA.np[, J])
	}
	if(THR.NU.PHI){
		nuA.np[nuA.np < MIN.NU] <- MIN.NU
		phiA.np[phiA.np < MIN.PHI] <- MIN.PHI
	}
	nuA(object)[marker.index, ] <- nuA.np
	phiA(object)[marker.index, ] <- phiA.np
	if(is.lds) { close(object); return(TRUE)}
	return(object)
}

summarizeMaleXNps <- function(marker.index,
			      batches,
			      object, MIN.SAMPLES){
	nr <- length(marker.index)
	nc <- length(batchNames(object))
	NN.Mlist <- imputed.medianA <- imputed.medianB <- shrink.madA <- shrink.madB <- vector("list", nc)
	gender <- object$gender
	AA <- as.matrix(A(object)[marker.index, gender==1])
	madA.AA <- medianA.AA <- matrix(NA, nr, nc)
	numberMenPerBatch <- rep(NA, nc)
	for(k in seq_along(batches)){
		B <- batches[[k]]
		this.batch <- unique(as.character(batch(object)[B]))
		gender <- object$gender[B]
		if(sum(gender==1) < MIN.SAMPLES) next()
		sns.batch <- sampleNames(object)[B]
		##subset GG apppriately
		sns <- colnames(AA)
		J <- sns%in%sns.batch
		numberMenPerBatch[k] <- length(J)
		medianA.AA[, k] <- rowMedians(AA[, J], na.rm=TRUE)
		madA.AA[, k] <- rowMAD(AA[, J], na.rm=TRUE)
	}
	return(list(medianA.AA=medianA.AA,
		    madA.AA=madA.AA))
}


summarizeXGenotypes <- function(marker.index,
				    batches,
				    object,
				    GT.CONF.THR,
				    MIN.OBS,
				    MIN.SAMPLES,
				    verbose,
				    is.lds,
				    DF.PRIOR,
				    gender="male", ...){
	if(gender == "male"){
		sample.index <- which(object$gender==1)
	} else sample.index <- which(object$gender==2)
	nr <- length(marker.index)
	nc <- length(batchNames(object))
##	NN.Mlist <- imputed.medianA <- imputed.medianB <- shrink.madA <- shrink.madB <- vector("list", nc)
	NN.Mlist <- medianA <- medianB <- shrink.madA <- shrink.madB <- vector("list", nc)
	##gender <- object$gender
	GG <- as.matrix(calls(object)[marker.index, sample.index])
	CP <- as.matrix(snpCallProbability(object)[marker.index, sample.index])
	AA <- as.matrix(A(object)[marker.index, sample.index])
	BB <- as.matrix(B(object)[marker.index, sample.index])
	for(k in seq_along(batches)){
		sample.index <- batches[[k]]
		this.batch <- unique(as.character(batch(object)[sample.index]))
		##gender <- object$gender[sample.index]
		if(length(sample.index) < MIN.SAMPLES) next()
		sns.batch <- sampleNames(object)[sample.index]
		##subset GG apppriately
		sns <- colnames(GG)
		J <- sns%in%sns.batch
		G <- GG[, J]
		xx <- CP[, J]
		highConf <- (1-exp(-xx/1000)) > GT.CONF.THR
		G <- G*highConf
		A <- AA[, J]
		B <- BB[, J]
		G.AA <- G==1
		G.AA[G.AA==FALSE] <- NA
		G.AB <- G==2
		G.AB[G.AB==FALSE] <- NA
		G.BB <- G==3
		G.BB[G.BB==FALSE] <- NA
		N.AA <- rowSums(G.AA, na.rm=TRUE)
		if(gender == "female")
			N.AB <- rowSums(G.AB, na.rm=TRUE)
		N.BB <- rowSums(G.BB, na.rm=TRUE)
		summaryStats <- function(X, INT, FUNS){
			tmp <- matrix(NA, nrow(X), length(FUNS))
			for(j in seq_along(FUNS)){
				FUN <- match.fun(FUNS[j])
				tmp[, j] <- FUN(X*INT, na.rm=TRUE)
			}
			tmp
		}
		statsA.AA <- summaryStats(G.AA, A, FUNS=c("rowMedians", "rowMAD"))
		if(gender == "female")
			statsA.AB <- summaryStats(G.AB, A, FUNS=c("rowMedians", "rowMAD"))
		statsA.BB <- summaryStats(G.BB, A, FUNS=c("rowMedians", "rowMAD"))
		statsB.AA <- summaryStats(G.AA, B, FUNS=c("rowMedians", "rowMAD"))
		if(gender == "female")
			statsB.AB <- summaryStats(G.AB, B, FUNS=c("rowMedians", "rowMAD"))
		statsB.BB <- summaryStats(G.BB, B, FUNS=c("rowMedians", "rowMAD"))
		if(gender=="male"){
			medianA[[k]] <- cbind(statsA.AA[, 1], ##statsA.AB[, 1],
					      statsA.BB[, 1])
			medianB[[k]] <- cbind(statsB.AA[, 1], ##statsB.AB[, 1],
					      statsB.BB[, 1])
			madA <- cbind(statsA.AA[, 2],  ##statsA.AB[, 2],
				      statsA.BB[, 2])
			madB <- cbind(statsB.AA[, 2], ##statsB.AB[, 2],
				      statsB.BB[, 2])
			NN <- cbind(N.AA, N.BB)
			rm(statsA.AA, statsA.BB, statsB.AA, statsB.BB)
		} else {
			medianA[[k]] <- cbind(statsA.AA[, 1], statsA.AB[, 1],
					      statsA.BB[, 1])
			medianB[[k]] <- cbind(statsB.AA[, 1], statsB.AB[, 1],
					      statsB.BB[, 1])
			madA <- cbind(statsA.AA[, 2],  statsA.AB[, 2],
				      statsA.BB[, 2])
			madB <- cbind(statsB.AA[, 2], statsB.AB[, 2],
				      statsB.BB[, 2])
			NN <- cbind(N.AA, N.AB, N.BB)
			rm(statsA.AA, statsA.AB, statsA.BB, statsB.AA, statsB.AB, statsB.BB)
		}
		NN.Mlist[[k]] <- NN
		shrink.madA[[k]] <- shrink(madA, NN, DF.PRIOR)
		shrink.madB[[k]] <- shrink(madB, NN, DF.PRIOR)
		##---------------------------------------------------------------------------
		## SNPs that we'll use for imputing location/scale of unobserved genotypes
		##---------------------------------------------------------------------------
		index.complete <- indexComplete(NN, medianA[[k]], medianB[[k]], MIN.OBS)
		##---------------------------------------------------------------------------
		## Impute sufficient statistics for unobserved genotypes (plate-specific)
		##---------------------------------------------------------------------------
		if(gender=="male"){
			res <- imputeCenterX(medianA[[k]], medianB[[k]], NN, index.complete, MIN.OBS)
		} else {
			unobservedAA <- NN[, 1] < MIN.OBS
			unobservedAB <- NN[, 2] < MIN.OBS
			unobservedBB <- NN[, 3] < MIN.OBS
			unobserved.index <- vector("list", 3)
			## AB, BB observed
			unobserved.index[[1]] <- which(unobservedAA & (NN[, 2] >= MIN.OBS & NN[, 3] >= MIN.OBS))
			## AA and BB observed (strange)
			unobserved.index[[2]] <- which(unobservedAB & (NN[, 1] >= MIN.OBS & NN[, 3] >= MIN.OBS))
			## AB and AA observed
			unobserved.index[[3]] <- which(unobservedBB & (NN[, 2] >= MIN.OBS & NN[, 1] >= MIN.OBS))
			res <- imputeCenter(medianA[[k]], medianB[[k]], index.complete, unobserved.index)
		}
		medianA[[k]] <- res[[1]]
		medianB[[k]] <- res[[2]]
	}
	return(list(madA=shrink.madA,
		    madB=shrink.madB,
		    NN=NN.Mlist,
		    medianA=medianA,
		    medianB=medianB))
}

## X chromosome, SNPs
fit.lm3 <- function(strata,
		    index.list,
		    object,
		    SNRMin,
		    MIN.SAMPLES,
		    MIN.OBS,
		    DF.PRIOR,
		    GT.CONF.THR,
		    THR.NU.PHI,
		    MIN.NU,
		    MIN.PHI,
		    verbose, is.lds, CHR.X, ...){
	if(is.lds) {physical <- get("physical"); open(object)}
	if(verbose) message("      Probe stratum ", strata, " of ", length(index.list))
	gender <- object$gender
	enough.males <- sum(gender==1) > MIN.SAMPLES
	enough.females <- sum(gender==2) > MIN.SAMPLES
	if(!enough.males & !enough.females){
		message(paste("fewer than", MIN.SAMPLES, "men and women.  Copy number not estimated for CHR X"))
		return(object)
	}
	marker.index <- index.list[[strata]]
	batches <- split(seq_along(batch(object)), as.character(batch(object)))
	batches <- batches[sapply(batches, length) >= MIN.SAMPLES]
	nuA <- as.matrix(nuA(object)[marker.index, ])
	nuB <- as.matrix(nuB(object)[marker.index, ])
	phiA <- as.matrix(phiA(object)[marker.index, ])
	phiB <- as.matrix(phiB(object)[marker.index, ])
	phiA2 <- as.matrix(phiPrimeA(object)[marker.index, ])
	phiB2 <- as.matrix(phiPrimeB(object)[marker.index, ])
	if(enough.males){
		res <- summarizeXGenotypes(marker.index=marker.index,
					       batches=batches,
					       object=object,
					       GT.CONF.THR=GT.CONF.THR,
					       MIN.SAMPLES=MIN.SAMPLES,
					       MIN.OBS=MIN.OBS,
					       verbose=verbose,
					       is.lds=is.lds,
					       DF.PRIOR=DF.PRIOR/2,
					   gender="male")
		madA.Mlist <- res[["madA"]]
		madB.Mlist <- res[["madB"]]
		medianA.Mlist <- res[["medianA"]]
		medianB.Mlist <- res[["medianB"]]
		NN.Mlist <- res[["NN"]]
		rm(res)
		## Need N, median, mad
	}
	if(enough.females){
		res <- summarizeXGenotypes(marker.index=marker.index,
						 batches=batches,
						 object=object,
						 GT.CONF.THR=GT.CONF.THR,
						 MIN.SAMPLES=MIN.SAMPLES,
						 MIN.OBS=MIN.OBS,
						 verbose=verbose,
						 is.lds=is.lds,
						 DF.PRIOR=DF.PRIOR/2,
					   gender="female")
		madA.Flist <- res[["madA"]]
		madB.Flist <- res[["madB"]]
		medianA.Flist <- res[["medianA"]]
		medianB.Flist <- res[["medianB"]]
		NN.Flist <- res[["NN"]]
		rm(res)
	}
	for(k in seq_along(batches)){
		sample.index <- batches[[k]]
		this.batch <- unique(as.character(batch(object)[sample.index]))
		gender <- object$gender[sample.index]
		enough.men <- sum(gender==1) >= MIN.SAMPLES
		enough.women <- sum(gender==2) >= MIN.SAMPLES
		if(!enough.men & !enough.women) {
			if(verbose) message(paste("fewer than", MIN.SAMPLES, "men and women in batch", this.batch, ". CHR X copy number not available. "))
			next()
		}
		if(enough.women){
			madA.F <- madA.Flist[[k]]
			madB.F <- madB.Flist[[k]]
			medianA.F <- medianA.Flist[[k]]
			medianB.F <- medianB.Flist[[k]]
			NN.F <- NN.Flist[[k]]
		}
		if(enough.men){
			madA.M <- madA.Mlist[[k]]
			madB.M <- madB.Mlist[[k]]
			medianA.M <- medianA.Mlist[[k]]
			medianB.M <- medianB.Mlist[[k]]
			NN.M <- NN.Mlist[[k]]
		}
		if(enough.men & enough.women){
			betas <- fit.wls(cbind(NN.M, NN.F),
					 sigma=cbind(madA.M, madA.F),
					 allele="A",
					 Y=cbind(medianA.M, medianA.F),
					 autosome=FALSE)
			nuA[, k] <- betas[1, ]
			phiA[, k] <- betas[2, ]
			phiA2[, k] <- betas[3, ]
			betas <- fit.wls(cbind(NN.M, NN.F),
					 sigma=cbind(madB.M, madB.F),
					 allele="B",
					 Y=cbind(medianB.M, medianB.F),
					 autosome=FALSE)
			nuB[, k] <- betas[1, ]
			phiB[, k] <- betas[2, ]
			phiB2[, k] <- betas[3, ]
		}
		if(enough.men & !enough.women){
			betas <- fit.wls(NN.M,##[, c(1,3)],
					 sigma=madA.M,##[, c(1,3)],
					 allele="A",
					 Y=medianA.M,##[, c(1,3)],
					 autosome=FALSE,
					 X=cbind(1, c(1, 0)))
			nuA[, k] <- betas[1, ]
			phiA[, k] <- betas[2, ]
			betas <- fit.wls(NN.M,##[, c(1,3)],
					 sigma=madB.M,#[, c(1,3)],
					 allele="B",
					 Y=medianB.M,#[, c(1,3)],
					 autosome=FALSE,
					 X=cbind(1, c(0, 1)))
			nuB[, k] <- betas[1, ]
			phiB[, k] <- betas[2, ]
		}
		if(!enough.men & enough.women){
			betas <- fit.wls(NN.F,
					 sigma=madA.F,
					 allele="A",
					 Y=medianA.F,
					 autosome=TRUE) ## can just use the usual design matrix for the women-only analysis
			nuA[, k] <- betas[1, ]
			phiA[, k] <- betas[2, ]
			betas <- fit.wls(NN.F,
					 sigma=madB.F,
					 allele="B",
					 Y=medianB.F,
					 autosome=TRUE) ## can just use the usual design matrix for the women-only analysis
			nuB[, k] <- betas[1, ]
			phiB[, k] <- betas[2, ]
		}
	}
	if(THR.NU.PHI){
		nuA[nuA < MIN.NU] <- MIN.NU
		nuB[nuB < MIN.NU] <- MIN.NU
		phiA[phiA < MIN.PHI] <- MIN.PHI
		phiA2[phiA2 < 1] <- 1
		phiB[phiB < MIN.PHI] <- MIN.PHI
		phiB2[phiB2 < 1] <- 1
	}
	nuA(object)[marker.index, ] <- nuA
	nuB(object)[marker.index, ] <- nuB
	phiA(object)[marker.index, ] <- phiA
	phiB(object)[marker.index, ] <- phiB
	phiPrimeA(object)[marker.index, ] <- phiA2
	phiPrimeB(object)[marker.index, ] <- phiB2
	if(is.lds) {close(object); return(TRUE)} else return(object)
}

fit.lm4 <- function(strata,
		    index.list,
		    object,
		    MIN.SAMPLES,
		    THR.NU.PHI,
		    MIN.NU,
		    MIN.PHI,
		    verbose, is.lds, ...){
	if(is.lds) {physical <- get("physical"); open(object)}
	gender <- object$gender
	enough.males <- sum(gender==1) > MIN.SAMPLES
	enough.females <- sum(gender==2) > MIN.SAMPLES
	if(!enough.males & !enough.females){
		message(paste("fewer than", MIN.SAMPLES, "men and women.  Copy number not estimated for CHR X"))
		return(object)
	}
	if(verbose) message("      Probe stratum ", strata, " of ", length(index.list))
	marker.index <- index.list[[strata]]
	batches <- split(seq_along(batch(object)), as.character(batch(object)))
	batches <- batches[sapply(batches, length) >= MIN.SAMPLES]
	nc <- length(batchNames(object))
	if(enough.males){
		res <- summarizeMaleXNps(marker.index=marker.index,
					 batches=batches,
					 object=object, MIN.SAMPLES=MIN.SAMPLES)
		medianA.AA.M <- res[["medianA.AA"]]
		madA.AA.M <- res[["madA.AA"]]

	}
	medianA.AA.F <- as.matrix(medianA.AA(object)[marker.index, ]) ## median for women
	madA.AA.F <- as.matrix(madA.AA(object)[marker.index, ]) ## median for women
	split.gender <- split(gender, as.character(batch(object)))
	N.M <- sapply(split.gender, function(x) sum(x==1))
	N.F <- sapply(split.gender, function(x) sum(x==2))
	nuA <- as.matrix(nuA(object)[marker.index, ])
	nuB <- as.matrix(nuB(object)[marker.index, ])
	phiA <- as.matrix(phiA(object)[marker.index, ])
	phiB <- as.matrix(phiB(object)[marker.index, ])
	ii <- isSnp(object) & chromosome(object) < 23 & !is.na(chromosome(object))
	fns <- featureNames(object)[ii]
	flags <- as.matrix(flags(object)[ii, ])
	fns.noflags <- fns[rowSums(flags, na.rm=T) == 0]
	snp.index <- sample(match(fns.noflags, featureNames(object)), 10000)
	N.AA <- as.matrix(N.AA(object)[snp.index, ])
	N.AB <- as.matrix(N.AA(object)[snp.index, ])
	N.BB <- as.matrix(N.AA(object)[snp.index, ])
	enoughAA <- rowSums(N.AA < 5) == 0
	enoughAB <- rowSums(N.AB < 5) == 0
	enoughBB <- rowSums(N.BB < 5) == 0
	snp.index <- snp.index[enoughAA & enoughAB & enoughBB]
	stopifnot(length(snp.index) > 100)
	nuA.snp.notmissing <- rowSums(is.na(as.matrix(nuA(object)[snp.index, ]))) == 0
	nuA.snp.notnegative <- rowSums(as.matrix(nuA(object)[snp.index, ]) < 20) == 0
	snp.index <- snp.index[nuA.snp.notmissing & nuA.snp.notnegative]
	stopifnot(length(snp.index) > 100)

	medianA.AA.snp <- as.matrix(medianA.AA(object)[snp.index,])
	medianB.BB.snp <- as.matrix(medianB.BB(object)[snp.index,])

	nuA.snp <- as.matrix(nuA(object)[snp.index, ])
	nuB.snp <- as.matrix(nuB(object)[snp.index, ])
	phiA.snp <- as.matrix(phiA(object)[snp.index, ])
	phiB.snp <- as.matrix(phiB(object)[snp.index, ])
##	pseudoAR <- position(object)[snp.index] < 2709520 | (position(object)[snp.index] > 154584237 & position(object)[snp.index] < 154913754)
##	pseudoAR[is.na(pseudoAR)] <- FALSE
	for(k in seq_along(batches)){
		B <- batches[[k]]
		this.batch <- unique(as.character(batch(object)[B]))
		gender <- object$gender[B]
		enough.men <- N.M[k] >= MIN.SAMPLES
		enough.women <- N.F[k] >= MIN.SAMPLES
		if(!enough.men & !enough.women) {
			if(verbose) message(paste("fewer than", MIN.SAMPLES, "men and women in batch", this.batch, ". CHR X copy number not available. "))
			next()
		}
		tmp <- cbind(medianA.AA.snp[, k], medianB.BB.snp[,k], phiA.snp[, k], phiB.snp[, k])
		tmp <- tmp[rowSums(is.na(tmp) == 0) & rowSums(tmp < 20) == 0, ]
		if(nrow(tmp) < 100){
			stop("too few markers for estimating nonpolymorphic CN on chromosome X")
		}
		X <- cbind(1, log2(c(tmp[, 1], tmp[, 2])))
		Y <- log2(c(tmp[, 3], tmp[, 4]))
		betahat <- solve(crossprod(X), crossprod(X, Y))
		if(enough.men){
			X.men <- cbind(1, log2(medianA.AA.M[, k]))
			Yhat1 <- as.numeric(X.men %*% betahat)
			## put intercept and slope for men in nuB and phiB
			phiB[, k] <- 2^(Yhat1)
			nuB[, k] <- medianA.AA.M[, k] - 1*phiB[, k]
		}
		X.fem <- cbind(1, log2(medianA.AA.F[, k]))
		Yhat2 <- as.numeric(X.fem %*% betahat)
		phiA[, k] <- 2^(Yhat2)
		nuA[, k] <- medianA.AA.F[, k] - 2*phiA[, k]
	}
	if(THR.NU.PHI){
		nuA[nuA < MIN.NU] <- MIN.NU
		phiA[phiA < MIN.PHI] <- MIN.PHI
		nuB[nuB < MIN.NU] <- MIN.NU
		phiB[phiB < MIN.PHI] <- MIN.PHI
	}
	nuA(object)[marker.index, ] <- nuA
	phiA(object)[marker.index, ] <- phiA
	nuB(object)[marker.index, ] <- nuB
	phiB(object)[marker.index, ] <- phiB
	if(is.lds) {close(object); return(TRUE)} else return(object)
}

whichPlatform <- function(cdfName){
	index <- grep("genomewidesnp", cdfName)
	if(length(index) > 0){
		platform <- "affymetrix"
	} else{
		index <- grep("human", cdfName)
		platform <- "illumina"
	}
	return(platform)
}

cnrma <- function(A, filenames, row.names, verbose=TRUE, seed=1, cdfName, sns){
	if(missing(cdfName)) stop("must specify cdfName")
	pkgname <- getCrlmmAnnotationName(cdfName)
	require(pkgname, character.only=TRUE) || stop("Package ", pkgname, " not available")
	if (missing(sns)) sns <- basename(filenames)
	if(verbose) message("Loading annotations for nonpolymorphic probes")
        loader("npProbesFid.rda", .crlmmPkgEnv, pkgname)
	fid <- getVarInEnv("npProbesFid")
	if(cdfName=="genomewidesnp6"){
		loader("1m_reference_cn.rda", .crlmmPkgEnv, pkgname)
	}
	if(cdfName=="genomewidesnp5"){
		loader("5.0_reference_cn.rda", .crlmmPkgEnv, pkgname)
	}
	reference <- getVarInEnv("reference")
	fid <- fid[match(row.names, names(fid))]
	np.index <- match(row.names, rownames(A))
	gns <- names(fid)
	set.seed(seed)
	##idx2 <- sample(length(fid), 10^5)
	for (k in seq_along(filenames)){
		y <- as.matrix(read.celfile(filenames[k], intensity.means.only=TRUE)[["INTENSITY"]][["MEAN"]][fid])
		A[np.index, k] <- as.integer(normalize.quantiles.use.target(y, target=reference))
	}
	return(A)
}

cnrma2 <- function(A, filenames, row.names, verbose=TRUE, seed=1, cdfName, sns){
	if(missing(cdfName)) stop("must specify cdfName")
	pkgname <- getCrlmmAnnotationName(cdfName)
	require(pkgname, character.only=TRUE) || stop("Package ", pkgname, " not available")
	if (missing(sns)) sns <- basename(filenames)
	sampleBatches <- splitIndicesByNode(seq(along=filenames))
	if(verbose) message("Processing nonpolymorphic probes for ", length(filenames), " files.")
	## updates A
	ocLapply(sampleBatches,
		 processCEL2,
		 row.names=row.names,
		 filenames=filenames,
		 A=A,
		 seed=seed,
		 pkgname=pkgname,
		 cdfName=cdfName,
		 neededPkgs=c("crlmm", pkgname))
	##list(sns=sns, gns=row.names, SKW=SKW, cdfName=cdfName)
	return(A)
}

processCEL2 <- function(i, filenames, row.names, A, seed, cdfName, pkgname){
	if(cdfName=="genomewidesnp6"){
		loader("1m_reference_cn.rda", .crlmmPkgEnv, pkgname)
	}
	if(cdfName=="genomewidesnp5"){
		loader("5.0_reference_cn.rda", .crlmmPkgEnv, pkgname)
	}
	reference <- getVarInEnv("reference")
        loader("npProbesFid.rda", .crlmmPkgEnv, pkgname)
	fid <- getVarInEnv("npProbesFid")
	row.names <- row.names[row.names %in% names(fid)] ##removes chromosome Y
	fid <- fid[match(row.names, names(fid))]
	##stopifnot(all.equal(names(fid), row.names))
	np.index <- match(row.names, rownames(A))
	gns <- names(fid)
	set.seed(seed)
	##idx2 <- sample(length(fid), 10^5)
	for (k in i){
		y <- as.matrix(read.celfile(filenames[k], intensity.means.only=TRUE)[["INTENSITY"]][["MEAN"]][fid])
		A[np.index, k] <- as.integer(normalize.quantiles.use.target(y, target=reference))
	}
	return(TRUE)
}


imputeCenter <- function(muA, muB, index.complete, index.missing){
	index <- index.missing
	mnA <- muA
	mnB <- muB
	if(length(index.complete) >= 100){
		for(j in 1:3){
			if(length(index[[j]]) == 0) next()
			X <- cbind(1, mnA[index.complete,  -j, drop=FALSE], mnB[index.complete,  -j, drop=FALSE])
			Y <- cbind(mnA[index.complete, j], mnB[index.complete,  j])
			betahat <- solve(crossprod(X), crossprod(X,Y))
			X <- cbind(1, mnA[index[[j]],  -j, drop=FALSE],  mnB[index[[j]],  -j, drop=FALSE])
			mus <- X %*% betahat
			muA[index[[j]], j] <- mus[, 1]
			muB[index[[j]], j] <- mus[, 2]
		}
	}
	results <- list(muA, muB)
	return(results)
}

indexComplete <- function(NN, medianA, medianB, MIN.OBS){
	Nindex <- which(rowSums(NN > MIN.OBS) == ncol(NN))
	correct.order <- which(medianA[, 1] > medianA[, ncol(medianA)] & medianB[, ncol(medianB)] > medianB[, 1])
	index.complete <- intersect(Nindex, correct.order)
	size <- min(5000, length(index.complete))
	if(size == 5000) index.complete <- sample(index.complete, 5000, replace=TRUE)
##	if(length(index.complete) < 100){
##		stop("fewer than 100 snps pass criteria for imputing unobserved genotype location/scale")
##	}
	return(index.complete)
}

imputeCentersForMonomorphicSnps <- function(medianA, medianB, index.complete, unobserved.index){
	cols <- c(3, 1, 2)
	if(length(index.complete) < 100){
		##message("index.complete less than 100.  No imputation")
		results <- list(medianA=medianA, medianB=medianB)
		return(results)
	}
	for(j in 1:3){
		if(length(unobserved.index[[j]]) == 0) next()
		kk <- cols[j]
		X <- cbind(1, medianA[index.complete, kk], medianB[index.complete, kk])
		Y <- cbind(medianA[index.complete,  -kk],
			   medianB[index.complete,  -kk])
		betahat <- solve(crossprod(X), crossprod(X,Y))
		X <- cbind(1, medianA[unobserved.index[[j]],  kk], medianB[unobserved.index[[j]],  kk])
		mus <- X %*% betahat
		medianA[unobserved.index[[j]], -kk] <- mus[, 1:2]
		medianB[unobserved.index[[j]], -kk] <- mus[, 3:4]
	}
	results <- list(medianA=medianA, medianB=medianB)
	return(results)
}


imputeCenterX <- function(muA, muB, Ns, index.complete, MIN.OBS){
	##index1 <- which(Ns[, 1] == 0 & Ns[, 3] > MIN.OBS)
	index1 <- which(Ns[, 1] == 0 & Ns[, 2] > MIN.OBS)
	if(length(index1) > 0){
		##X <- cbind(1, muA[index.complete, 3], muB[index.complete, 3])
		X <- cbind(1, muA[index.complete, 2], muB[index.complete, 2])
		Y <- cbind(1, muA[index.complete, 1], muB[index.complete, 1])
##		X <- cbind(1, muA[index.complete[[1]], 3], muB[index.complete[[1]], 3])
##		Y <- cbind(1, muA[index.complete[[1]], 1], muB[index.complete[[1]], 1])
		betahat <- solve(crossprod(X), crossprod(X,Y))
		##now with the incomplete SNPs
		##X <- cbind(1, muA[index1, 3], muB[index1, 3])
		X <- cbind(1, muA[index1, 2], muB[index1, 2])
		mus <- X %*% betahat
		muA[index1, 1] <- mus[, 2]
		muB[index1, 1] <- mus[, 3]
	}
	index1 <- which(Ns[, 2] == 0)
	if(length(index1) > 0){
		X <- cbind(1, muA[index.complete, 1], muB[index.complete, 1])
		Y <- cbind(1, muA[index.complete, 2], muB[index.complete, 2])
##		X <- cbind(1, muA[index.complete[[2]], 1], muB[index.complete[[2]], 1])
##		Y <- cbind(1, muA[index.complete[[2]], 3], muB[index.complete[[2]], 3])
		betahat <- solve(crossprod(X), crossprod(X,Y))
		##now with the incomplete SNPs
		X <- cbind(1, muA[index1, 1], muB[index1, 1])
		mus <- X %*% betahat
		muA[index1, 2] <- mus[, 2]
		muB[index1, 2] <- mus[, 3]
	}
	list(muA, muB)
}





## constructors for Illumina platform
constructIlluminaFeatureData <- function(gns, cdfName){
	pkgname <- paste(cdfName, "Crlmm", sep="")
	path <- system.file("extdata", package=pkgname)
	load(file.path(path, "cnProbes.rda"))
	load(file.path(path, "snpProbes.rda"))
	cnProbes$chr <- chromosome2integer(cnProbes$chr)
	cnProbes <- as.matrix(cnProbes)
	snpProbes$chr <- chromosome2integer(snpProbes$chr)
	snpProbes <- as.matrix(snpProbes)
	mapping <- rbind(snpProbes, cnProbes, deparse.level=0)
	mapping <- mapping[match(gns, rownames(mapping)), ]
	isSnp <- 1L-as.integer(gns %in% rownames(cnProbes))
	mapping <- cbind(mapping, isSnp, deparse.level=0)
	stopifnot(identical(rownames(mapping), gns))
	colnames(mapping) <- c("chromosome", "position", "isSnp")
	new("AnnotatedDataFrame",
	    data=data.frame(mapping),
	    varMetadata=data.frame(labelDescription=colnames(mapping)))
}
constructIlluminaAssayData <- function(np, snp, object, storage.mode="environment", nr){
	stopifnot(identical(snp$gns, featureNames(object)))
	gns <- c(featureNames(object), np$gns)
	sns <- np$sns
	np <- np[1:2]
	snp <- snp[1:2]
	stripnames <- function(x) {
		dimnames(x) <- NULL
		x
	}
	np <- lapply(np, stripnames)
	snp <- lapply(snp, stripnames)
	is.ff <- is(snp[[1]], "ff") | is(snp[[1]], "ffdf")
	if(is.ff){
		lapply(snp, open)
		open(calls(object))
		open(snpCallProbability(object))
##		lapply(np, open)
	}
	##tmp <- rbind(as.matrix(snp[[1]]), as.matrix(np[[1]]), deparse.level=0)
	A.snp <- snp[[1]]
	B.snp <- snp[[2]]
	##Why is np not a ff object?
	A.np <- np[[1]]
	B.np <- np[[2]]
	nc <- ncol(object)
	if(is.ff){
		NA.vec <- rep(NA, nrow(A.np))
		AA <- initializeBigMatrix("A", nr, nc, vmode="integer")
		BB <- initializeBigMatrix("B", nr, nc, vmode="integer")
		GG <- initializeBigMatrix("calls", nr, nc, vmode="integer")
		PP <- initializeBigMatrix("confs", nr, nc, vmode="integer")
		for(j in 1:ncol(object)){
			AA[, j] <- c(snp[[1]][, j], np[[1]][, j])
			BB[, j] <- c(snp[[2]][, j], np[[2]][, j])
			GG[, j] <- c(calls(object)[, j], NA.vec)
			PP[, j] <- c(snpCallProbability(object)[, j], NA.vec)

		}
	} else {
		AA <- rbind(snp[[1]], np[[1]], deparse.level=0)
		BB <- rbind(snp[[2]], np[[2]], deparse.level=0)
		gt <- stripnames(calls(object))
		emptyMatrix <- matrix(integer(), nrow(np[[1]]), ncol(A))
		GG <- rbind(gt, emptyMatrix, deparse.level=0)
		pr <- stripnames(snpCallProbability(object))
		PP <- rbind(pr, emptyMatrix, deparse.level=0)
	}
	assayDataNew(storage.mode,
		     alleleA=AA,
		     alleleB=BB,
		     call=GG,
		     callProbability=PP)
}
constructIlluminaCNSet <- function(crlmmResult,
				   path,
				   snpFile,
				   cnFile){
	load(file.path(path, "snpFile.rda"))
	res <- get("res")
	load(file.path(path, "cnFile.rda"))
	cnAB <- get("cnAB")
	fD <- constructIlluminaFeatureData(c(res$gns, cnAB$gns), cdfName="human370v1c")
	##new.order <- order(fD$chromosome, fD$position)
	##fD <- fD[new.order, ]
	aD <- constructIlluminaAssayData(cnAB, res, crlmmResult, nr=nrow(fD))
	##protocolData(crlmmResult)$batch <- vector("integer", ncol(crlmmResult))
	new("CNSet",
	    call=aD[["call"]],
	    callProbability=aD[["callProbability"]],
	    alleleA=aD[["alleleA"]],
	    alleleB=aD[["alleleB"]],
	    phenoData=phenoData(crlmmResult),
	    protocolData=protocolData(crlmmResult),
	    featureData=fD,
	    batch=batch,
	    annotation="human370v1c")

}



ellipseCenters <- function(object, index, allele, batch, log.it=TRUE){
	ubatch <- unique(batch(object))[batch]
	Nu <- nu(object, allele)[index, batch]
	Phi <- phi(object, allele)[index, batch]
	centers <- list(Nu, Nu+Phi, Nu+2*Phi)
	if(log.it)
		centers <- lapply(centers, log2)
	myLabels <- function(allele){
		switch(allele,
		       A=c("BB", "AB", "AA"),
		       B=c("AA", "AB", "BB"),
		       stop("allele must be 'A' or 'B'")
		       )
	}
	nms <- myLabels(allele)
	names(centers) <- nms
	fns <- featureNames(object)[index]
	centers$fns <- fns
	return(centers)
}


shrinkSummary <- function(object,
			  type=c("SNP", "X.SNP"), ##"X.snps", "X.nps"),
			  MIN.OBS=1,
			  MIN.SAMPLES=10,
			  DF.PRIOR=50,
			  verbose=TRUE,
			  marker.index,
			  is.lds){
	stopifnot(type[[1]] %in% c("SNP", "X.SNP"))
	if(type[[1]] == "X.SNP"){
		gender <- object$gender
		if(sum(gender == 2) < 3) {
			message("too few females to estimate within genotype summary statistics on CHR X")
			return(object)
		}
		CHR.X <- TRUE
	} else CHR.X <- FALSE
	if(missing(marker.index)){
		batch <- batch(object)
		is.snp <- isSnp(object)
		is.autosome <- chromosome(object) < 23
		is.annotated <- !is.na(chromosome(object))
		is.X <- chromosome(object) == 23
		is.lds <- is(calls(object), "ffdf") | is(calls(object), "ff_matrix")
		if(is.lds) stopifnot(isPackageLoaded("ff"))
		marker.index <- whichMarkers(type[[1]], is.snp, is.autosome, is.annotated, is.X)
	}
	if(is.lds){
		index.list <- splitIndicesByLength(marker.index, ocProbesets())
		ocLapply(seq(along=index.list),
			 shrinkGenotypeSummaries,
			 index.list=index.list,
			 object=object,
			 verbose=verbose,
			 MIN.OBS=MIN.OBS,
			 MIN.SAMPLES=MIN.SAMPLES,
			 DF.PRIOR=DF.PRIOR,
			 is.lds=is.lds,
			 neededPkgs="crlmm")
	} else {
		object <- shrinkGenotypeSummaries(strata=1,
			      index.list=list(marker.index),
			      object=object,
			      verbose=verbose,
			      MIN.OBS=MIN.OBS,
			      MIN.SAMPLES=MIN.SAMPLES,
			      DF.PRIOR=DF.PRIOR,
			      is.lds=is.lds)
	}
	return(object)
}

genotypeSummary <- function(object,
			    GT.CONF.THR=0.95,
			    type=c("SNP", "NP", "X.SNP", "X.NP"), ##"X.snps", "X.nps"),
			    verbose=TRUE,
			    marker.index,
			    is.lds){
	if(type == "X.SNP" | type=="X.NP"){
##		gender <- object$gender
##		## the number of women in each batch could be less than 3...
##		if(sum(gender == 2) < 3) {
##			message("too few females to estimate within genotype summary statistics on CHR X")
##			return(object)
##		}
		CHR.X <- TRUE
	} else CHR.X <- FALSE
	if(missing(marker.index)){
		batch <- batch(object)
		is.snp <- isSnp(object)
		is.autosome <- chromosome(object) < 23
		is.annotated <- !is.na(chromosome(object))
		is.X <- chromosome(object) == 23
		is.lds <- is(calls(object), "ffdf") | is(calls(object), "ff_matrix")
		if(is.lds) stopifnot(isPackageLoaded("ff"))
		marker.index <- whichMarkers(type[[1]], is.snp, is.autosome, is.annotated, is.X)
	}
	summaryFxn <- function(type){
		switch(type,
		       SNP="summarizeSnps",
		       NP="summarizeNps",
##		       X.SNP="summarizeSnps",
		       X.NP="summarizeNps")
	}
	myf <- summaryFxn(type[[1]])
	FUN <- get(myf)
	if(is.lds){
		index.list <- splitIndicesByLength(marker.index, ocProbesets())
		ocLapply(seq(along=index.list),
			 FUN,
			 index.list=index.list,
			 object=object,
			 batchSize=ocProbesets(),
			 GT.CONF.THR=GT.CONF.THR,
			 verbose=verbose,
			 is.lds=is.lds,
			 CHR.X=CHR.X,
			 neededPkgs="crlmm")
	} else{
		object <- FUN(strata=1,
			      index.list=list(marker.index),
			      object=object,
			      batchSize=ocProbesets(),
			      GT.CONF.THR=GT.CONF.THR,
			      verbose=verbose,
			      is.lds=is.lds,
			      CHR.X=CHR.X)
	}
	return(object)
}

whichMarkers <- function(type, is.snp, is.autosome, is.annotated, is.X){
	switch(type,
	       SNP=which(is.snp & is.autosome & is.annotated),
	       NP=which(!is.snp & is.autosome),
	       X.SNP=which(is.snp & is.X),
	       X.NP=which(!is.snp & is.X),
	       stop("'type' must be one of 'SNP', 'NP', 'X.SNP', or 'X.NP'")
	       )
}

summarizeNps <- function(strata, index.list, object, batchSize,
			 GT.CONF.THR, verbose, is.lds, CHR.X, ...){
	if(is.lds) {physical <- get("physical"); open(object)}
	if(verbose) message("      Probe stratum ", strata, " of ", length(index.list))
	index <- index.list[[strata]]
##	if(CHR.X) {
##		sample.index <- which(object$gender==2)
##		batches <- split(sample.index, as.character(batch(object))[sample.index])
##	} else {
##		batches <- split(seq_along(batch(object)), as.character(batch(object)))
##	}
	batches <- split(seq_along(batch(object)), as.character(batch(object)))
	batchnames <- batchNames(object)
	nr <- length(index)
	nc <- length(batchnames)
	N.AA <- medianA.AA <- madA.AA <- tau2A.AA <- matrix(NA, nr, nc)
	is.illumina <- whichPlatform(annotation(object)) == "illumina"
	AA <- as.matrix(A(object)[index, ])
	if(is.illumina){
		BB <- as.matrix(B(object)[index, ])
		AVG <- (AA+BB)/2
		A(object)[index, ] <- AVG
		AA <- AVG
		rm(AVG, BB)
	}
	for(k in seq_along(batches)){
		sample.index <- batches[[k]]
		N.AA[, k] <- length(sample.index)
		if(CHR.X){
			gender <- object$gender[sample.index]
			sample.index <- sample.index[gender == 2]
			if(length(sample.index) == 0) next()
		}
		this.batch <- unique(as.character(batch(object)[sample.index]))
		j <- match(this.batch, batchnames)
		I.A <- AA[, sample.index]
		medianA.AA[, k] <- rowMedians(I.A, na.rm=TRUE)
		madA.AA[, k] <- rowMAD(I.A, na.rm=TRUE)
		## log2 Transform Intensities
		I.A <- log2(I.A)
		tau2A.AA[, k] <- rowMAD(I.A, na.rm=TRUE)^2
	}
	N.AA(object)[index,] <- N.AA
	medianA.AA(object)[index,] <- medianA.AA
	madA.AA(object)[index, ] <- madA.AA
	tau2A.AA(object)[index, ] <- tau2A.AA
	if(is.lds) return(TRUE) else return(object)
}

summarizeSnps <- function(strata,
			  index.list,
			  object,
			  GT.CONF.THR,
			  verbose, is.lds, CHR.X, ...){
	if(is.lds) {
		physical <- get("physical")
		open(object)
	}
	if(verbose) message("      Probe stratum ", strata, " of ", length(index.list))
	index <- index.list[[strata]]
##	if(CHR.X) {
##		sample.index <- which(object$gender==2)
##		batches <- split(sample.index, as.character(batch(object))[sample.index])
##	} else {
##		batches <- split(seq_along(batch(object)), as.character(batch(object)))
##	}
	## this message can be confusing if no women are in the dataset
##	if(CHR.X){
##		if(verbose) message("        biallelic cluster medians are estimated using only the women for SNPs on chr. X")
##	}
	batches <- split(seq_along(batch(object)), as.character(batch(object)))
	batchnames <- batchNames(object)
	nr <- length(index)
	nc <- length(batchnames)
	statsA.AA <- statsA.AB <- statsA.BB <- statsB.AA <- statsB.AB <- statsB.BB <- vector("list", nc)
	corrAA <- corrAB <- corrBB <- tau2A.AA <- tau2A.BB <- tau2B.AA <- tau2B.BB <- matrix(NA, nr, nc)
	Ns.AA <- Ns.AB <- Ns.BB <- matrix(NA, nr, nc)
	if(verbose) message("        Begin reading...")
	GG <- as.matrix(calls(object)[index, ])
	CP <- as.matrix(snpCallProbability(object)[index, ])
	AA <- as.matrix(A(object)[index, ])
	BB <- as.matrix(B(object)[index, ])
	FL <- as.matrix(flags(object)[index, ])
	if(verbose) message("        Computing summaries...")
	for(k in seq_along(batches)){
		##note that the genotype frequency for AA would include 'A' on chromosome X for men
		sample.index <- batches[[k]]
		this.batch <- unique(as.character(batch(object)[sample.index]))
		j <- match(this.batch, batchnames)
		G <- GG[, sample.index]
		##NORM <- normal.index[, k]
		xx <- CP[, sample.index]
		highConf <- (1-exp(-xx/1000)) > GT.CONF.THR
		G <- G*highConf
		A <- AA[, sample.index]
		B <- BB[, sample.index]
		## this can be time consuming...do only once
		G.AA <- G==1
		G.AA[G.AA==FALSE] <- NA
		G.AB <- G==2
		G.AB[G.AB==FALSE] <- NA
		G.BB <- G==3
		G.BB[G.BB==FALSE] <- NA
		Ns.AA[, k] <- rowSums(G.AA, na.rm=TRUE)
		Ns.AB[, k] <- rowSums(G.AB, na.rm=TRUE)
		Ns.BB[, k] <- rowSums(G.BB, na.rm=TRUE)
		if(CHR.X){
			gender <- object$gender[sample.index]
			sample.index <- sample.index[gender == 2]
			if(length(sample.index) == 0) next()
		}
		summaryStats <- function(X, INT, FUNS){
			tmp <- matrix(NA, nrow(X), length(FUNS))
			for(j in seq_along(FUNS)){
				FUN <- match.fun(FUNS[j])
				tmp[, j] <- FUN(X*INT, na.rm=TRUE)
			}
			tmp
		}
		stats <- summaryStats(G.AA, A, FUNS=c("rowMedians", "rowMAD"))
		medianA.AA(object)[index, k] <- stats[, 1]
		madA.AA(object)[index, k] <- stats[, 2]

		stats <- summaryStats(G.AB, A, FUNS=c("rowMedians", "rowMAD"))
		medianA.AB(object)[index, k] <- stats[, 1]
		madA.AB(object)[index, k] <- stats[, 2]

		stats <- summaryStats(G.BB, A, FUNS=c("rowMedians", "rowMAD"))
		medianA.BB(object)[index, k] <- stats[, 1]
		madA.BB(object)[index, k] <- stats[, 2]

		stats <- summaryStats(G.AA, B, FUNS=c("rowMedians", "rowMAD"))
		medianB.AA(object)[index, k] <- stats[, 1]
		madB.AA(object)[index, k] <- stats[, 2]

		stats <- summaryStats(G.AB, B, FUNS=c("rowMedians", "rowMAD"))
		medianB.AB(object)[index, k] <- stats[, 1]
		madB.AB(object)[index, k] <- stats[, 2]

		stats <- summaryStats(G.BB, B, FUNS=c("rowMedians", "rowMAD"))
		medianB.BB(object)[index, k] <- stats[, 1]
		madB.BB(object)[index, k] <- stats[, 2]

		## log2 Transform Intensities
		A <- log2(A); B <- log2(B)
		tau2A.AA[, k] <- summaryStats(G.AA, A, FUNS="rowMAD")^2
		tau2A.BB[, k] <- summaryStats(G.BB, A, FUNS="rowMAD")^2
		tau2B.AA[, k] <- summaryStats(G.AA, B, FUNS="rowMAD")^2
		tau2B.BB[, k] <- summaryStats(G.BB, B, FUNS="rowMAD")^2

		corrAA[, k] <- rowCors(A*G.AA, B*G.AA, na.rm=TRUE)
		corrAB[, k] <- rowCors(A*G.AB, B*G.AB, na.rm=TRUE)
		corrBB[, k] <- rowCors(A*G.BB, B*G.BB, na.rm=TRUE)
	}
	if(verbose) message("        Begin writing...")
	N.AA(object)[index,] <- Ns.AA
	N.AB(object)[index,] <- Ns.AB
	N.BB(object)[index,] <- Ns.BB
	corrAA(object)[index,] <- corrAA
	corrAB(object)[index,] <- corrAB
	corrBB(object)[index,] <- corrBB
##	medianA.AA(object)[index,] <- do.call(cbind, lapply(statsA.AA, function(x) x[, 1]))
##	medianA.AB(object)[index,] <- do.call(cbind, lapply(statsA.AB, function(x) x[, 1]))
##	medianA.BB(object)[index,] <- do.call(cbind, lapply(statsA.BB, function(x) x[, 1]))
##	medianB.AA(object)[index,] <- do.call(cbind, lapply(statsB.AA, function(x) x[, 1]))
##	medianB.AB(object)[index,] <- do.call(cbind, lapply(statsB.AB, function(x) x[, 1]))
##	medianB.BB(object)[index,] <- do.call(cbind, lapply(statsB.BB, function(x) x[, 1]))
##
##	madA.AA(object)[index,] <- do.call(cbind, lapply(statsA.AA, function(x) x[, 2]))
##	madA.AB(object)[index,] <- do.call(cbind, lapply(statsA.AB, function(x) x[, 2]))
##	madA.BB(object)[index,] <- do.call(cbind, lapply(statsA.BB, function(x) x[, 2]))
##	madB.AA(object)[index,] <- do.call(cbind, lapply(statsB.AA, function(x) x[, 2]))
##	madB.AB(object)[index,] <- do.call(cbind, lapply(statsB.AB, function(x) x[, 2]))
##	madB.BB(object)[index,] <- do.call(cbind, lapply(statsB.BB, function(x) x[, 2]))
	tau2A.AA(object)[index, ] <- tau2A.AA
##	tau2A.AB(object)[index, ] <- tau2A.AB
	tau2A.BB(object)[index, ] <- tau2A.BB
	tau2B.AA(object)[index, ] <- tau2B.AA
##	tau2B.AB(object)[index, ] <- tau2B.AB
	tau2B.BB(object)[index, ] <- tau2B.BB
	if(is.lds) return(TRUE) else return(object)
}

crlmmCopynumber <- function(object,
			    MIN.SAMPLES=10,
			    SNRMin=5,
			    MIN.OBS=1,
			    DF.PRIOR=50,
			    bias.adj=FALSE,
			    prior.prob=rep(1/4,4),
			    seed=1,
			    verbose=TRUE,
			    GT.CONF.THR=0.80,
			    MIN.NU=2^3,
			    MIN.PHI=2^3,
			    THR.NU.PHI=TRUE,
			    type=c("SNP", "NP", "X.SNP", "X.NP")){
	stopifnot(type %in% c("SNP", "NP", "X.SNP", "X.NP"))
	batch <- batch(object)
	is.snp <- isSnp(object)
	is.autosome <- chromosome(object) < 23
	is.annotated <- !is.na(chromosome(object))
	is.X <- chromosome(object) == 23
	is.lds <- is(calls(object), "ffdf") | is(calls(object), "ff_matrix")
	if(is.lds) stopifnot(isPackageLoaded("ff"))
	samplesPerBatch <- table(as.character(batch(object)))
	if(any(samplesPerBatch < MIN.SAMPLES)){
		warning("The following batches have fewer than ", MIN.SAMPLES, ":")
		message(paste(names(samplesPerBatch)[samplesPerBatch < MIN.SAMPLES], collapse=", "))
		message("Not estimating copy number for the above batches")
	}
	mylabel <- function(marker.type){
		switch(marker.type,
		       SNP="autosomal SNPs",
		       NP="autosomal nonpolymorphic markers",
		       X.SNP="chromosome X SNPs",
		       X.NP="chromosome X nonpolymorphic markers")
	}
	if(verbose) message("Computing summary statistics of the genotype clusters for each batch")
##	for(i in seq_along(type)){
	for(i in c(1, 2, 4)){ ## do not do X.SNP.  Do this during fit.lm3
		marker.type <- type[i]
		if(verbose) message(paste("...", mylabel(marker.type)))
		##if(verbose) message(paste("Computing summary statistics for ", mylabel(marker.type), " genotype clusters for each batch")
		marker.index <- whichMarkers(marker.type,
					     is.snp,
					     is.autosome,
					     is.annotated,
					     is.X)
		object <- genotypeSummary(object=object,
					  GT.CONF.THR=GT.CONF.THR,
					  type=marker.type,
					  verbose=verbose,
					  marker.index=marker.index,
					  is.lds=is.lds)
	}
	if(verbose) message("Imputing unobserved genotype medians and shrinking the variances (within-batch, across loci) ")##SNPs only
	for(i in seq_along(type)){
		marker.type <- type[i]
		if(!marker.type %in% c("SNP", "X.SNP")) next()
		message(paste("...", mylabel(marker.type)))
		marker.index <- whichMarkers(marker.type, is.snp,
					     is.autosome, is.annotated, is.X)
		object <- shrinkSummary(object=object,
					MIN.OBS=MIN.OBS,
					MIN.SAMPLES=MIN.SAMPLES,
					DF.PRIOR=DF.PRIOR,
					type=marker.type,
					verbose=verbose,
					marker.index=marker.index,
					is.lds=is.lds)
	}
	if(verbose) message("Estimating parameters for copy number")##SNPs only
	for(i in seq_along(type)){
		marker.type <- type[i]
		message(paste("...", mylabel(marker.type)))
		CHR.X <- ifelse(marker.type %in% c("X.SNP", "X.NP"), TRUE, FALSE)
		marker.index <- whichMarkers(marker.type, is.snp,
					     is.autosome, is.annotated, is.X)
		object <- estimateCnParameters(object=object,
					       type=marker.type,
					       SNRMin=SNRMin,
					       DF.PRIOR=DF.PRIOR,
					       GT.CONF.THR=GT.CONF.THR,
					       MIN.SAMPLES=MIN.SAMPLES,
					       MIN.OBS=MIN.OBS,
					       MIN.NU=MIN.NU,
					       MIN.PHI=MIN.PHI,
					       THR.NU.PHI=THR.NU.PHI,
					       verbose=verbose,
					       marker.index=marker.index,
					       is.lds=is.lds,
					       CHR.X=CHR.X)
	}
	return(object)
}
crlmmCopynumber2 <- function(){
	.Defunct(msg="The crlmmCopynumber2 function has been deprecated. The function crlmmCopynumber should be used instead.  crlmmCopynumber will support large data using ff provided that the ff package is loaded.")
}
crlmmCopynumberLD <- crlmmCopynumber2


estimateCnParameters <- function(object,
				 type=c("SNP", "NP", "X.SNP", "X.NP"),
				 verbose=TRUE,
				 SNRMin=5,
				 DF.PRIOR=50,
				 GT.CONF.THR=0.95,
				 MIN.SAMPLES=10,
				 MIN.OBS=1,
				 MIN.NU=8,
				 MIN.PHI=8,
				 THR.NU.PHI=TRUE,
				 marker.index,
				 CHR.X,
				 is.lds){
	if(missing(marker.index)){
		batch <- batch(object)
		is.snp <- isSnp(object)
		is.autosome <- chromosome(object) < 23
		is.annotated <- !is.na(chromosome(object))
		is.X <- chromosome(object) == 23
		is.lds <- is(calls(object), "ffdf") | is(calls(object), "ff_matrix")
		if(is.lds) stopifnot(isPackageLoaded("ff"))
		CHR.X <- ifelse(type[[1]] %in% c("X.SNP", "X.NP"), TRUE, FALSE)
		marker.index <- whichMarkers(type[[1]], is.snp, is.autosome, is.annotated, is.X)
	}
	lmFxn <- function(type){
		switch(type,
		       SNP="fit.lm1",
		       NP="fit.lm2",
		       X.SNP="fit.lm3",
		       X.NP="fit.lm4")
	}
	myfun <- lmFxn(type[[1]])
	FUN <- get(myfun)
	if(is.lds){
		index.list <- splitIndicesByLength(marker.index, ocProbesets())
		ocLapply(seq(along=index.list),
			 FUN,
			 index.list=index.list,
			 marker.index=marker.index,
			 object=object,
			 batchSize=ocProbesets(),
			 SNRMin=SNRMin,
			 MIN.SAMPLES=MIN.SAMPLES,
			 MIN.OBS=MIN.OBS,
			 DF=DF.PRIOR,
			 GT.CONF.THR=GT.CONF.THR,
			 THR.NU.PHI=THR.NU.PHI,
			 MIN.NU=MIN.NU,
			 MIN.PHI=MIN.PHI,
			 verbose=verbose,
			 is.lds=is.lds,
			 CHR.X=CHR.X,
			 neededPkgs="crlmm")
	} else {
		##FUN <- match.fun(FUN)
		object <- FUN(strata=1,
			      index.list=list(marker.index),
			      marker.index=marker.index,
			      object=object,
			      batchSize=ocProbesets(),
			      SNRMin=SNRMin,
			      MIN.SAMPLES=MIN.SAMPLES,
			      MIN.OBS=MIN.OBS,
			      DF.PRIOR=DF.PRIOR,
			      GT.CONF.THR=GT.CONF.THR,
			      THR.NU.PHI=THR.NU.PHI,
			      MIN.NU=MIN.NU,
			      MIN.PHI=MIN.PHI,
			      is.lds=is.lds,
			      CHR.X=CHR.X,
			      verbose=verbose)
	}
	message("finished")
	return(object)
}
