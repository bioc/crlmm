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
getFeatureData.Affy <- function(cdfName, copynumber=FALSE){
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
	M[index, "chromosome"] <- snpProbes[, grep("chr", colnames(snpProbes))]
	M[index, "isSnp"] <- 1L
	index <- which(is.na(M[, "isSnp"]))
	M[index, "isSnp"] <- 1L
	if(copynumber){
		index <- match(rownames(cnProbes), rownames(M)) #only snp probes in M get assigned position
		M[index, "position"] <- cnProbes[, grep("pos", colnames(cnProbes))]
		M[index, "chromosome"] <- cnProbes[, grep("chr", colnames(cnProbes))]
		M[index, "isSnp"] <- 0L
	}
	##A few of the snpProbes do not match -- I think it is chromosome Y.
	M[is.na(M[, "isSnp"]), "isSnp"] <- 1L
	return(new("AnnotatedDataFrame", data=data.frame(M)))
}

construct <- function(filenames,
		      cdfName,
		      copynumber=TRUE,
		      sns, verbose=TRUE, batch, fns){
	if(!missing(batch)){
		stopifnot(length(batch) == length(sns))
	}
	if(missing(sns) & missing(filenames)) stop("one of filenames or samplenames (sns) must be provided")
	if(verbose) message("Initializing container for copy number estimation")
	featureData <- getFeatureData.Affy(cdfName, copynumber=copynumber)
	if(!missing(fns)){
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
	open(callSet)
	mixtureParams <- matrix(NA, 4, length(filenames))
	is.snp <- isSnp(callSet)
	snp.index <- which(is.snp)
	FUN <- ifelse(is.lds, "snprma2", "snprma")
	snprmaFxn <- function(FUN,...){
		switch(FUN,
		       snprma=snprma(...),
		       snprma2=snprma2(...))
	}
	snprmaRes <- snprmaFxn(FUN, filenames=filenames,
			       mixtureSampleSize=mixtureSampleSize,
			       fitMixture=TRUE,
			       eps=eps,
			       verbose=verbose,
			       seed=seed,
			       cdfName=cdfName,
			       sns=sns)
	if(verbose) message("Finished preprocessing.")
	if(is.lds){
		open(snprmaRes[["A"]])
		open(snprmaRes[["B"]])
		bb = ocProbesets()*ncol(A)*8
		ffrowapply(A(callSet)[i1:i2, ] <- snprmaRes[["A"]][i1:i2, ], X=snprmaRes[["A"]], BATCHBYTES=bb)
		ffrowapply(B(callSet)[i1:i2, ] <- snprmaRes[["B"]][i1:i2, ], X=snprmaRes[["B"]], BATCHBYTES=bb)
	} else{
		A(callSet)[snp.index, ] <- snprmaRes[["A"]]
		B(callSet)[snp.index, ] <- snprmaRes[["B"]]
	}
	pData(callSet)$SKW <- snprmaRes[["SKW"]]
	pData(callSet)$SNR <- snprmaRes[["SNR"]]
	mixtureParams <- snprmaRes$mixtureParams
	np.index <- which(!is.snp)
	if(verbose) message("Normalizing nonpolymorphic markers")
	FUN <- ifelse(is.lds, "cnrma2", "cnrma")
	## main purpose is to update 'alleleA'
	cnrmaFxn <- function(FUN,...){
		switch(FUN,
		       cnrma=cnrma(...),
		       cnrma2=cnrma2(...))
	}
	## consider passing only A for NPs.
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
	tmp <- crlmmGTfxn(FUN,
			  A=snprmaRes[["A"]],
			  B=snprmaRes[["B"]],
			  SNR=snprmaRes[["SNR"]],
			  mixtureParams=snprmaRes[["mixtureParams"]],
			  cdfName=cdfName,
			  row.names=NULL,
			  col.names=sampleNames(callSet),
			  probs=probs,
			  DF=DF,
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
		ffrowapply(snpCall(callSet)[i1:i2, ] <- tmp[["calls"]][i1:i2, ], X=tmp[["calls"]], BATCHBYTES=bb)
		ffrowapply(snpCallProbability(callSet)[i1:i2, ] <- tmp[["confs"]][i1:i2, ], X=tmp[["confs"]], BATCHBYTES=bb)
		close(tmp[["calls"]])
		close(tmp[["confs"]])
	} else {
		calls(callSet)[snp.index, ] <- tmp[["calls"]]
		snpCallProbability(callSet)[snp.index, ] <- tmp[["confs"]]
	}
	callSet$gender <- tmp$gender
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
	if(verbose) message("Probe stratum ", strata, " of ", length(index.list))
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
		B <- batches[[k]]
		this.batch <- unique(as.character(batch(object)[B]))

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
		##
		##---------------------------------------------------------------------------
		## Impute sufficient statistics for unobserved genotypes (plate-specific)
		##---------------------------------------------------------------------------
		unobservedAA <- NN[, 1] < MIN.OBS
		unobservedAB <- NN[, 2] < MIN.OBS
		unobservedBB <- NN[, 3] < MIN.OBS
		unobserved.index <- vector("list", 3)
		unobserved.index[[1]] <- which(unobservedAA & (NN[, 2] >= MIN.OBS & NN[, 3] >= MIN.OBS))
		unobserved.index[[2]] <- which(unobservedAB & (NN[, 1] >= MIN.OBS & NN[, 3] >= MIN.OBS))
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
	if(verbose) message("Probe stratum ", strata, " of ", length(index.list))
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

fit.lm2 <- function(strata,
		    index.list,
		    object,
		    MIN.SAMPLES,
		    THR.NU.PHI,
		    MIN.NU,
		    MIN.PHI,
		    verbose, is.lds, CHR.X, ...){
	if(is.lds) {physical <- get("physical"); open(object)}
	if(verbose) message("Probe stratum ", strata, " of ", length(index.list))
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
		crosshyb <- max(median(medianA.AA[, k]) - median(medianA.AA.np[, k]), 0)
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


summarizeMaleXGenotypes <- function(marker.index,
				    batches,
				    object,
				    GT.CONF.THR,
				    MIN.OBS,
				    MIN.SAMPLES,
				    verbose,
				    is.lds,
				    DF.PRIOR,...){
	nr <- length(marker.index)
	nc <- length(batchNames(object))
	NN.Mlist <- imputed.medianA <- imputed.medianB <- shrink.madA <- shrink.madB <- vector("list", nc)
	gender <- object$gender
	GG <- as.matrix(calls(object)[marker.index, gender==1])
	CP <- as.matrix(snpCallProbability(object)[marker.index, gender==1])
	AA <- as.matrix(A(object)[marker.index, gender==1])
	BB <- as.matrix(B(object)[marker.index, gender==1])
	for(k in seq_along(batches)){
		B <- batches[[k]]
		this.batch <- unique(as.character(batch(object)[B]))
		gender <- object$gender[B]
		if(sum(gender==1) < MIN.SAMPLES) next()
		sns.batch <- sampleNames(object)[B]
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
		N.AA.M <- rowSums(G.AA, na.rm=TRUE)
		N.AB.M <- rowSums(G.AB, na.rm=TRUE)
		N.BB.M <- rowSums(G.BB, na.rm=TRUE)
		summaryStats <- function(X, INT, FUNS){
			tmp <- matrix(NA, nrow(X), length(FUNS))
			for(j in seq_along(FUNS)){
				FUN <- match.fun(FUNS[j])
				tmp[, j] <- FUN(X*INT, na.rm=TRUE)
			}
			tmp
		}
		statsA.AA <- summaryStats(G.AA, A, FUNS=c("rowMedians", "rowMAD"))
		statsA.AB <- summaryStats(G.AB, A, FUNS=c("rowMedians", "rowMAD"))
		statsA.BB <- summaryStats(G.BB, A, FUNS=c("rowMedians", "rowMAD"))
		statsB.AA <- summaryStats(G.AA, B, FUNS=c("rowMedians", "rowMAD"))
		statsB.AB <- summaryStats(G.AB, B, FUNS=c("rowMedians", "rowMAD"))
		statsB.BB <- summaryStats(G.BB, B, FUNS=c("rowMedians", "rowMAD"))
		medianA <- cbind(statsA.AA[, 1], statsA.AB[, 1], statsA.BB[, 1])
		medianB <- cbind(statsB.AA[, 1], statsB.AB[, 1], statsB.BB[, 1])
		madA <- cbind(statsA.AA[, 1], statsA.AB[, 1], statsA.BB[, 1])
		madB <- cbind(statsB.AA[, 1], statsB.AB[, 1], statsB.BB[, 1])
		rm(statsA.AA, statsA.AB, statsA.BB, statsB.AA, statsB.AB, statsB.BB)

		NN.M <- cbind(N.AA.M, N.AB.M, N.BB.M)
		NN.Mlist[[k]] <- NN.M

		shrink.madA[[k]] <- shrink(madA, NN.M, DF.PRIOR)
		shrink.madB[[k]] <- shrink(madB, NN.M, DF.PRIOR)

		##---------------------------------------------------------------------------
		## SNPs that we'll use for imputing location/scale of unobserved genotypes
		##---------------------------------------------------------------------------
		index.complete <- indexComplete(NN.M[, -2], medianA, medianB, MIN.OBS)

		##---------------------------------------------------------------------------
		## Impute sufficient statistics for unobserved genotypes (plate-specific)
		##---------------------------------------------------------------------------
		res <- imputeCenterX(medianA, medianB, NN.M, index.complete, MIN.OBS)
		imputed.medianA[[k]] <- res[[1]]
		imputed.medianB[[k]] <- res[[2]]
	}
	return(list(madA=shrink.madA,
		    madB=shrink.madB,
		    NN.M=NN.Mlist,
		    medianA=imputed.medianA,
		    medianB=imputed.medianB))
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
	if(verbose) message("Probe stratum ", strata, " of ", length(index.list))
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
		res <- summarizeMaleXGenotypes(marker.index=marker.index, batches=batches,
					       object=object, GT.CONF.THR=GT.CONF.THR,
					       MIN.SAMPLES=MIN.SAMPLES,
					       MIN.OBS=MIN.OBS,
					       verbose=verbose, is.lds=is.lds,
					       DF.PRIOR=DF.PRIOR/2)
		madA.Mlist <- res[["madA"]]
		madB.Mlist <- res[["madB"]]
		medianA.Mlist <- res[["medianA"]]
		medianB.Mlist <- res[["medianB"]]
		NN.Mlist <- res[["NN.M"]]
		rm(res)
		## Need N, median, mad
	}
	if(enough.females){
		N.AA.F <- as.matrix(N.AA(object)[marker.index, ])
		N.AB.F <- as.matrix(N.AB(object)[marker.index, ])
		N.BB.F <- as.matrix(N.BB(object)[marker.index, ])
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
	}
	for(k in seq_along(batches)){
		B <- batches[[k]]
		this.batch <- unique(as.character(batch(object)[B]))
		gender <- object$gender[B]
		enough.men <- sum(gender==1) >= MIN.SAMPLES
		enough.women <- sum(gender==2) >= MIN.SAMPLES
		if(!enough.men & !enough.women) {
			if(verbose) message(paste("fewer than", MIN.SAMPLES, "men and women in batch", this.batch, ". CHR X copy number not available. "))
			next()
		}
		if(enough.women){
			medianA.F <- cbind(medianA.AA[, k], medianA.AB[, k], medianA.BB[, k])
			medianB.F <- cbind(medianB.AA[, k], medianB.AB[, k], medianB.BB[, k])
			madA.F <- cbind(madA.AA[, k], madA.AB[, k], madA.BB[, k])
			madB.F <- cbind(madB.AA[, k], madB.AB[, k], madB.BB[, k])
			NN.F <- cbind(N.AA.F[, k], N.AB.F[, k], N.BB.F[, k])
		}
		if(enough.men){
			madA.M <- madA.Mlist[[k]]
			madB.M <- madB.Mlist[[k]]
			medianA.M <- medianA.Mlist[[k]]
			medianB.M <- medianB.Mlist[[k]]
			NN.M <- NN.Mlist[[k]]
		}
		if(enough.men & enough.women){
			betas <- fit.wls(cbind(NN.M[, c(1,3)], NN.F),
					 sigma=cbind(madA.M[, c(1,3)], madA.F),
					 allele="A",
					 Y=cbind(medianA.M[, c(1,3)], medianA.F),
					 autosome=FALSE)
			nuA[, k] <- betas[1, ]
			phiA[, k] <- betas[2, ]
			phiA2[, k] <- betas[3, ]
			betas <- fit.wls(cbind(NN.M[, c(1,3)], NN.F),
					 sigma=cbind(madB.M[, c(1,3)], madB.F),
					 allele="B",
					 Y=cbind(medianB.M[, c(1,3)], medianB.F),
					 autosome=FALSE)
			nuB[, k] <- betas[1, ]
			phiB[, k] <- betas[2, ]
			phiB2[, k] <- betas[3, ]
		}
		if(enough.men & !enough.women){
			betas <- fit.wls(NN.M[, c(1,3)],
					 sigma=madA.M[, c(1,3)],
					 allele="A",
					 Y=medianA.M[, c(1,3)],
					 autosome=FALSE,
					 X=cbind(1, c(0, 1)))
			nuA[, k] <- betas[1, ]
			phiA[, k] <- betas[2, ]
			betas <- fit.wls(NN.M[, c(1,3)],
					 sigma=madB.M[, c(1,3)],
					 allele="B",
					 Y=medianB.M[, c(1,3)],
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
		phiA2[phiA2 < MIN.PHI] <- MIN.PHI
		phiB[phiB < MIN.PHI] <- MIN.PHI
		phiB2[phiB2 < MIN.PHI] <- MIN.PHI
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
	if(verbose) message("Probe stratum ", strata, " of ", length(index.list))
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
		stopifnot(nrow(tmp) > 100)
		X <- cbind(1, log2(c(tmp[, 1], tmp[, 2])))
		Y <- log2(c(tmp[, 3], tmp[, 4]))
		betahat <- solve(crossprod(X), crossprod(X, Y))
		X.men <- cbind(1, medianA.AA.M[, k])
		Yhat1 <- as.numeric(X.men %*% betahat)
		## put intercept and slope for men in nuB and phiB
		phiB[, k] <- 2^(Yhat1)
		nuB[, k] <- 2^(medianA.AA.M[, k]) - phiB[, k]

		X.fem <- cbind(1, medianA.AA.F[, k])
		Yhat2 <- as.numeric(X.fem %*% betahat)
		phiA[, k] <- 2^(Yhat2)
		nuA[, k] <- 2^(medianA.AA.F[, k]) - 2*phiA[, k]
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
	idx2 <- sample(length(fid), 10^5)
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
	fid <- fid[match(row.names, names(fid))]
	np.index <- match(row.names, rownames(A))
	gns <- names(fid)
	set.seed(seed)
	idx2 <- sample(length(fid), 10^5)
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
	list(muA, muB)
}

indexComplete <- function(NN, medianA, medianB, MIN.OBS){
	Nindex <- which(rowSums(NN > MIN.OBS) == ncol(NN))
	correct.order <- which(medianA[, 1] > medianA[, ncol(medianA)] & medianB[, ncol(medianB)] > medianB[, 1])
	index.complete <- intersect(Nindex, correct.order)
	size <- min(5000, length(index.complete))
	if(size == 5000) index.complete <- sample(index.complete, 5000, replace=TRUE)
	if(length(index.complete) < 100){
		stop("fewer than 100 snps pass criteria for imputing unobserved genotype location/scale")
	}
	return(index.complete)
}

imputeCentersForMonomorphicSnps <- function(medianA, medianB, index.complete, unobserved.index){
	cols <- c(3, 1, 2)
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
	list(medianA=medianA, medianB=medianB)
}


imputeCenterX <- function(muA, muB, Ns, index.complete, MIN.OBS){
	index1 <- which(Ns[, 1] == 0 & Ns[, 3] > MIN.OBS)
	if(length(index1) > 0){
		X <- cbind(1, muA[index.complete, 3], muB[index.complete, 3])
		Y <- cbind(1, muA[index.complete, 1], muB[index.complete, 1])
##		X <- cbind(1, muA[index.complete[[1]], 3], muB[index.complete[[1]], 3])
##		Y <- cbind(1, muA[index.complete[[1]], 1], muB[index.complete[[1]], 1])
		betahat <- solve(crossprod(X), crossprod(X,Y))
		##now with the incomplete SNPs
		X <- cbind(1, muA[index1, 3], muB[index1, 3])
		mus <- X %*% betahat
		muA[index1, 1] <- mus[, 2]
		muB[index1, 1] <- mus[, 3]
	}
	index1 <- which(Ns[, 3] == 0)
	if(length(index1) > 0){
		X <- cbind(1, muA[index.complete, 1], muB[index.complete, 1])
		Y <- cbind(1, muA[index.complete, 3], muB[index.complete, 3])
##		X <- cbind(1, muA[index.complete[[2]], 1], muB[index.complete[[2]], 1])
##		Y <- cbind(1, muA[index.complete[[2]], 3], muB[index.complete[[2]], 3])
		betahat <- solve(crossprod(X), crossprod(X,Y))
		##now with the incomplete SNPs
		X <- cbind(1, muA[index1, 1], muB[index1, 1])
		mus <- X %*% betahat
		muA[index1, 3] <- mus[, 2]
		muB[index1, 3] <- mus[, 3]
	}
	list(muA, muB)
}

##posteriorProbability.snps <- function(object, cnOptions, tmp.objects=list()){
##	I <- isSnp(object)
##	gender <- object$gender
##	CHR <- unique(chromosome(object))
##	##batch <- unique(object$batch)
##	batch <- unique(batch(object))
##	if(CHR == 23){
##		phiAX <- getParam(object, "phiAX", batch)[I]
##		phiBX <- getParam(object, "phiBX", batch)[I]
##	}
##	A <- A(object)[I, ]
##	B <- B(object)[I, ]
##	sig2A <- getParam(object, "sig2A", batch)[I]
##	sig2B <- getParam(object, "sig2B", batch)[I]
##	tau2A <- getParam(object, "tau2A", batch)[I]
##	tau2B <- getParam(object, "tau2B", batch)[I]
##	corrA.BB <- getParam(object, "corrA.BB", batch)[I]
##	corrB.AA <- getParam(object, "corrB.AA", batch)[I]
##	corr <- getParam(object, "corr", batch)[I]
##	nuA <- getParam(object, "nuA", batch)[I]
##	nuB <- getParam(object, "nuB", batch)[I]
##	phiA <- getParam(object, "phiA", batch)[I]
##	phiB <- getParam(object, "phiB", batch)[I]
##	normal <- tmp.objects[["normal"]][I, ]
##	prior.prob <- cnOptions$prior.prob
##	emit <- array(NA, dim=c(nrow(A), ncol(A), 10))##SNPs x sample x 'truth'
##	lA <- log2(A)
##	lB <- log2(B)
##	X <- cbind(lA, lB)
##	counter <- 1##state counter
##	for(CT in 0:3){
##		for(CA in 0:CT){
##			cat(".")
##			CB <- CT-CA
##			A.scale <- sqrt(tau2A*(CA==0) + sig2A*(CA > 0))
##			B.scale <- sqrt(tau2B*(CA==0) + sig2B*(CA > 0))
##			if(CA == 0 & CB == 0) rho <- 0
##			if(CA == 0 & CB > 0) rho <- corrA.BB
##			if(CA > 0 & CB == 0) rho <- corrB.AA
##			if(CA > 0 & CB > 0) rho <- corr
##			if(CHR == 23){
##				##means <- cbind(suppressWarnings(log2(nuA[, p]+CA*phiA[, p] + CB*phiAx[, p])), suppressWarnings(log2(nuB[, p]+CB*phiB[, p] + CA*phiBx[, p])))
##				meanA <- suppressWarnings(log2(nuA+CA*phiA + CB*phiAX))
##				meanB <- suppressWarnings(log2(nuB+CB*phiB + CA*phiBX))
##			} else{
##				##means <- cbind(suppressWarnings(log2(nuA+CA*phiA)), suppressWarnings(log2(nuB+CB*phiB)))
##				meanA <- suppressWarnings(log2(nuA+CA*phiA))
##				meanB <- suppressWarnings(log2(nuB+CB*phiB))
##				covs <- rho*A.scale*B.scale
##				A.scale2 <- A.scale^2
##				B.scale2 <- B.scale^2
##			}
##			Q.x.y <- 1/(1-rho^2)*(((lA - meanA)/A.scale)^2 + ((lB - meanB)/B.scale)^2 - 2*rho*((lA - meanA)*(lB - meanB))/(A.scale*B.scale))
##			f.x.y <- 1/(2*pi*A.scale*B.scale*sqrt(1-rho^2))*exp(-0.5*Q.x.y)
##			emit[, , counter] <- f.x.y
##			counter <- counter+1
##		}
##	}
##	priorProb <- cnOptions$prior.prob
##	homDel <- priorProb[1]*emit[, , 1]
##	hemDel <- priorProb[2]*emit[, , c(2, 3)] # + priorProb[3]*emit[, c(4, 5, 6)] + priorProb[4]*emit[, c(7:10)]
##	norm <- priorProb[3]*emit[, , 4:6]
##	amp <- priorProb[4]*emit[, , 7:10]
##	##sum over the different combinations within each copy number state
##	hemDel <- hemDel[, , 1] + hemDel[, , 2]
##	norm <- norm[, , 1] + norm[, , 2] + norm[, , 3]
##	amp <- amp[, , 1] + amp[, , 2] + amp[ , , 3] + amp[, , 4]
##	total <- homDel + hemDel + norm + amp
##	homDel <- homDel/total
##	hemDel <- hemDel/total
##	norm <- norm/total
##	amp <- amp/total
##	tmp.objects$posteriorProb <- list(hemDel=hemDel, norm=norm, amp=amp)
##	##envir[["posteriorProb"]] <- list(hemDel=hemDel, norm=norm, amp=amp)
##	posteriorProb <- array(NA, dim=c(nrow(A), ncol(A), 4))
##	posteriorProb[, , 1] <- homDel
##	posteriorProb[, , 2] <- hemDel
##	posteriorProb[, , 3] <- norm
##	posteriorProb[, , 4] <- amp
##	return(list(tmp.objects, posteriorProb))
##}
##
##biasAdj <- function(object, cnOptions, tmp.objects){
##	gender <- object$gender
##	CHR <- unique(chromosome(object))
##	I <- isSnp(object)
##	A <- A(object)[I, ]
##	normal <- tmp.objects[["normal"]][I, ]
##	results <- posteriorProbability.snps(object=object, cnOptions=cnOptions, tmp.objects=tmp.objects)
##	posteriorProb <- results[[2]]
##	tmp.objects <- results[[1]]
##	mostLikelyState <- apply(posteriorProb, c(1, 2), function(x) order(x, decreasing=TRUE)[1])
##	if(CHR == 23){
##		##so state index 3 is the most likely state for men and women
##		mostLikelyState[, gender==1] <- mostLikelyState[, gender==1] + 1
##	}
##	proportionSamplesAltered <- rowMeans(mostLikelyState != 3)
##	ii <- proportionSamplesAltered < 0.8 & proportionSamplesAltered > 0.01
##	##  only exclude observations from one tail, depending on
##	##  whether more are up or down
##	moreup <- rowSums(mostLikelyState > 3) > rowSums(mostLikelyState < 3) ##3 is normal
##	## if equal, which points have a high posterior probability of altered copy number.  drop those.
##	NORM <- matrix(FALSE, nrow(A), ncol(A))
##	##NORM[proportionSamplesAltered > 0.8, ] <- FALSE
##	ratioUp <- posteriorProb[, , 4]/posteriorProb[, , 3]
##	NORM[ii & moreup, ] <- ratioUp[moreup & ii] < 1  ##normal more likely
##	ratioDown <- posteriorProb[, , 2]/posteriorProb[, , 3]
##	NORM[ii & !moreup, ] <- ratioDown[!moreup & ii] < 1  ##normal more likely
##	normal <- NORM*normal
##	tmp <- tmp.objects[["normal"]]
##	tmp[I, ] <- normal
##	tmp.objects[["normal"]] <- tmp
##	return(tmp.objects)
##}
##
##bias2 <- function(strata.index,
##		  snpBatches,
##		  index,
##		  object,
##		  normal,
##		  prior.prob,
##		  MIN.SAMPLES,
##		  verbose){
##	open(object)
##	open(normal)
##
##	nps <- snpBatches[[strata.index]]
##	nuA <- lM(object)$nuA[nps, , drop=FALSE]
##	phiA <- lM(object)$phiA[nps, , drop=FALSE]
##	sig2A <- lM(object)$sig2A[nps, , drop=FALSE]
##	AA <- as.matrix(A(object)[nps, ])
##	batches <- split(seq(along=batch(object)), batch(object))
##	batches <- batches[sapply(batches, length) >= MIN.SAMPLES]
##
##	cn.lik <- matrix(NA, length(nps)*ncol(object), 4)
##	argmax.cn <- emit[nps, ]
##	norm <- matrix(1L, length(nps), ncol(object))
##
##	for(k in batches){
##		J <- match(unique(batch(object)[k]), unique(batch(object)))
##		lT <- log2(AA[, k])
##		counter <- 1 ##state counter
##		for(CT in c(0, 1.5, 2, 2.5)){
##			##sds <- sqrt(sig2A[, J]*(CT==0) + sig2A[ , J]*(CT > 0))
##			sds <- sqrt(sig2A[, J])
##			means <- suppressWarnings(log2(nuA[, J]+CT*phiA[, J]))
##			lik <- log(dnorm(lT, mean=means, sd=sds))
##			##emit[[counter]][nps, ] <- tmp
##			cn.lik[, counter] <- as.numeric(lik)
##			counter <- counter+1
##		}
##		outlier <- matrix(rowSums(cn.lik < -10) == 4, length(nps), ncol(object))
##		argmax.cn.lik <- apply(cn.lik, 1, function(x) order(x, decreasing=TRUE)[1])
##		argmax.cn <- matrix(argmax.cn.lik, length(nps), length(k))
##
##		isUp <- argmax.cn > 3
##		prUp <- rowMeans(isUp)
##
##		isDn <- argmax.cn < 3
##		prDn <- rowMeans(isDn)
##
##		index <- which(prUp > 0.05 & prUp > prDn)
##		##if proportion up greater than 5%, trim the high cn est.
##		norm[index, k] <- argmax.cn[index, ] > 3
##
##		index <- which(prDn > 0.05 & prDn > prUp)
##		norm[index, k] <- argmax.cn[index, ] < 3
##		norm[index, k] <- norm[index, k]*!outlier
##	}
##	normal[nps, ] <- norm
##	TRUE
##}


##biasAdjust <- function(object, prior.prob=rep(1/4, 4), MIN.SAMPLES=10, verbose=TRUE){
##	load(file.path(ldPath(), "normal.rda"))
##	autosomeIndex.nps <- (1:nrow(object))[chromosome(object) < 23 & !isSnp(object) & !is.na(chromosome(object))]
##
####	emit <- initializeBigMatrix("emit",
####				     nrow(object),
####				     ncol(object),
####				     vmode="double")
##	if(verbose) message("Bias adjustment for nonpolymorphic loci on chromosomes 1-22.")
##	snpBatches <- splitIndicesByLength(autosomeIndex.nps, ocProbesets())
##	ocLapply(seq(along=snpBatches),
##		 bias2,
##		 index=autosomeIndex.nps,
##		 snpBatches=snpBatches,
##		 object=object,
##		 normal=normal,
##		 prior.prob=prior.prob,
##		 MIN.SAMPLES=MIN.SAMPLES,
##		 verbose=verbose)
##
##	if(verbose) message("Bias adjustment for polymorphic loci on chromosomes 1-22.")
##	autosomeIndex.snps <- (1:nrow(object))[chromosome(object) < 23 & isSnp(object) & !is.na(chromosome(object))]
##	snpBatches <- splitIndicesByLength(autosomeIndex.snps, ocProbesets())
##	ocLapply(seq(along=snpBatches),
##		 bias1,
##		 index=autosomeIndex.snps,
##		 snpBatches=snpBatches,
##		 object=object,
##		 normal=normal,
##		 prior.prob=prior.prob,
##		 emit=emit,
##		 MIN.SAMPLES=MIN.SAMPLES,
##		 verbose=verbose)
##}
##biasAdjNP <- function(object, cnOptions, tmp.objects){
##	##batch <- unique(object$batch)
##	batch <- unique(batch(object))
##	normalNP <- tmp.objects[["normal"]][!isSnp(object), ]
##	CHR <- unique(chromosome(object))
##	A <- A(object)[!isSnp(object), ]
##	sig2A <- getParam(object, "sig2A", batch)
##	gender <- object$gender
##	##Assume that on the log-scale, that the background variance is the same...
##	tau2A <- sig2A
##	nuA <- getParam(object, "nuA", batch)
##	phiA <- getParam(object, "phiA", batch)
##	prior.prob <- cnOptions$prior.prob
##	emit <- array(NA, dim=c(nrow(A), ncol(A), 4))##SNPs x sample x 'truth'
##	lT <- log2(A)
##	I <- isSnp(object)
##	counter <- 1 ##state counter
####	for(CT in 0:3){
####		sds <- sqrt(tau2A[I]*(CT==0) + sig2A[I]*(CT > 0))
####		means <- suppressWarnings(log2(nuA[I]+CT*phiA[I]))
####		tmp <- dnorm(lT, mean=means, sd=sds)
####		emit[, , counter] <- tmp
####		counter <- counter+1
####	}
####	mostLikelyState <- apply(emit, c(1, 2), function(x) order(x, decreasing=TRUE)[1])
##	counter <- 1
##	for(CT in c(0,1,2,2.5)){
##		sds <- sqrt(tau2A[I]*(CT==0) + sig2A[I]*(CT > 0))
##		means <- suppressWarnings(log2(nuA[I]+CT*phiA[I]))
##		tmp <- dnorm(lT, mean=means, sd=sds)
##		emit[, , counter] <- tmp
##		counter <- counter+1
##	}
##	mostLikelyState <- apply(emit, c(1, 2), function(x) order(x, decreasing=TRUE)[1])
##
##	if(CHR == 23){
##		## the state index for male on chromosome 23  is 2
##		## add 1 so that the state index is 3 for 'normal' state
##		mostLikelyState[, gender=="male"] <- mostLikelyState[, gender==1] + 1
##	}
##	tmp3 <- mostLikelyState != 3
##	##Those near 1 have NaNs for nu and phi.  this occurs by NaNs in the muA[,, "A"] or muA[, , "B"] for X chromosome
##	proportionSamplesAltered <- rowMeans(tmp3)##prop normal
##	ii <- proportionSamplesAltered < 0.75
##	moreup <- rowSums(mostLikelyState > 3) > rowSums(mostLikelyState < 3)
##	notUp <-  mostLikelyState[ii & moreup, ] <= 3
##	notDown <- mostLikelyState[ii & !moreup, ] >= 3
##	NORM <- matrix(TRUE, nrow(A), ncol(A))
##	NORM[ii & moreup, ] <- notUp
##	NORM[ii & !moreup, ] <- notDown
##	normalNP <- normalNP*NORM
##
##	##flagAltered <- which(proportionSamplesAltered > 0.5)
##	##envir[["flagAlteredNP"]] <- flagAltered
##	normal <- tmp.objects[["normal"]]
##	normal[!isSnp(object), ] <- normalNP
##	tmp.objects[["normal"]] <- normal
##	return(tmp.objects)
##}




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
	if(is(snp[[1]], "ff")){
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
	if(is(A.snp, "ff")){
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
			return("too few females to estimate within genotype summary statistics on CHR X")
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
		gender <- object$gender
		if(sum(gender == 2) < 3) {
			return("too few females to estimate within genotype summary statistics on CHR X")
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
	summaryFxn <- function(type){
		switch(type,
		       SNP="summarizeSnps",
		       NP="summarizeNps",
		       X.SNP="summarizeSnps",
		       X.NP="summarizeNps")
	}
	FUN <- summaryFxn(type[[1]])
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
	} else {
		FUN <- match.fun(FUN)
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
	if(verbose) message("Probe stratum ", strata, " of ", length(index.list))
	index <- index.list[[strata]]
	if(CHR.X) {
		sample.index <- which(object$gender==2)
		batches <- split(sample.index, as.character(batch(object))[sample.index])
	} else {
		batches <- split(seq_along(batch(object)), as.character(batch(object)))
	}
	batchnames <- batchNames(object)
	nr <- length(index)
	nc <- length(batchnames)
	N.AA <- medianA.AA <- madA.AA <- tau2A.AA <- matrix(NA, nr, nc)
	AA <- as.matrix(A(object)[index, ])
	for(k in seq_along(batches)){
		B <- batches[[k]]
		N.AA[, k] <- length(B)
		this.batch <- unique(as.character(batch(object)[B]))
		j <- match(this.batch, batchnames)
		##NORM <- normal.index[, k]
		A <- AA[, B]
		medianA.AA[, k] <- rowMedians(A, na.rm=TRUE)
		madA.AA[, k] <- rowMAD(A, na.rm=TRUE)
		## log2 Transform Intensities
		A <- log2(A)
		tau2A.AA[, k] <- rowMAD(A, na.rm=TRUE)^2
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
	if(verbose) message("Probe stratum ", strata, " of ", length(index.list))
	index <- index.list[[strata]]
	if(CHR.X) {
		sample.index <- which(object$gender==2)
		batches <- split(sample.index, as.character(batch(object))[sample.index])
	} else {
		batches <- split(seq_along(batch(object)), as.character(batch(object)))
	}
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
	if(verbose) message("        Computing summaries...")
	for(k in seq_along(batches)){
		B <- batches[[k]]
		this.batch <- unique(as.character(batch(object)[B]))
		j <- match(this.batch, batchnames)
		G <- GG[, B]
		##NORM <- normal.index[, k]
		xx <- CP[, B]
		highConf <- (1-exp(-xx/1000)) > GT.CONF.THR
		##G <- G*highConf*NORM
		G <- G*highConf
		A <- AA[, B]
		B <- BB[, B]
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
		summaryStats <- function(X, INT, FUNS){
			tmp <- matrix(NA, nrow(X), length(FUNS))
			for(j in seq_along(FUNS)){
				FUN <- match.fun(FUNS[j])
				tmp[, j] <- FUN(X*INT, na.rm=TRUE)
			}
			tmp
		}
		statsA.AA[[k]] <- summaryStats(G.AA, A, FUNS=c("rowMedians", "rowMAD"))
		statsA.AB[[k]] <- summaryStats(G.AB, A, FUNS=c("rowMedians", "rowMAD"))
		statsA.BB[[k]] <- summaryStats(G.BB, A, FUNS=c("rowMedians", "rowMAD"))
		statsB.AA[[k]] <- summaryStats(G.AA, B, FUNS=c("rowMedians", "rowMAD"))
		statsB.AB[[k]] <- summaryStats(G.AB, B, FUNS=c("rowMedians", "rowMAD"))
		statsB.BB[[k]] <- summaryStats(G.BB, B, FUNS=c("rowMedians", "rowMAD"))
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
	medianA.AA(object)[index,] <- do.call(cbind, lapply(statsA.AA, function(x) x[, 1]))
	medianA.AB(object)[index,] <- do.call(cbind, lapply(statsA.AB, function(x) x[, 1]))
	medianA.BB(object)[index,] <- do.call(cbind, lapply(statsA.BB, function(x) x[, 1]))
	medianB.AA(object)[index,] <- do.call(cbind, lapply(statsB.AA, function(x) x[, 1]))
	medianB.AB(object)[index,] <- do.call(cbind, lapply(statsB.AB, function(x) x[, 1]))
	medianB.BB(object)[index,] <- do.call(cbind, lapply(statsB.BB, function(x) x[, 1]))

	madA.AA(object)[index,] <- do.call(cbind, lapply(statsA.AA, function(x) x[, 2]))
	madA.AB(object)[index,] <- do.call(cbind, lapply(statsA.AB, function(x) x[, 2]))
	madA.BB(object)[index,] <- do.call(cbind, lapply(statsA.BB, function(x) x[, 2]))
	madB.AA(object)[index,] <- do.call(cbind, lapply(statsB.AA, function(x) x[, 2]))
	madB.AB(object)[index,] <- do.call(cbind, lapply(statsB.AB, function(x) x[, 2]))
	madB.BB(object)[index,] <- do.call(cbind, lapply(statsB.BB, function(x) x[, 2]))
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
			    GT.CONF.THR=0.95,
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
		message(paste(samplesPerBatch[samplesPerBatch < MIN.SAMPLES], collapse=", "))
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
	for(i in seq_along(type)){
		## do all types
		marker.type <- type[i]
		if(verbose) message(paste("...", mylabel(marker.type)))
		##if(verbose) message(paste("Computing summary statistics for ", mylabel(marker.type), " genotype clusters for each batch")
		marker.index <- whichMarkers(marker.type, is.snp,
					     is.autosome, is.annotated, is.X)
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
	FUN <- lmFxn(type[[1]])
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
		FUN <- match.fun(FUN)
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
