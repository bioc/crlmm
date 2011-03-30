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
	M[is.na(M[, "chromosome"]), "isSnp"] <- 1L
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
	if(any(is.na(batch))){
		stop("NA's in batch variable")
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
	if(!missing(sns)){
		sampleNames(cnSet) <- sns
		protocolData <- annotatedDataFrameFrom(A(cnSet), byrow=FALSE)
	} else {
		sampleNames(cnSet) <- basename(filenames)
		protocolData <- getProtocolData.Affy(filenames)
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
	if(!is.lds) stop("this function now requires that the ff package be loaded")
	if(missing(cdfName)) stop("must specify cdfName")
	if(!isValidCdfName(cdfName)) stop("cdfName not valid.  see validCdfNames")
	if(missing(sns)) sns <- basename(filenames)
	stopifnot(!missing(batch))
	##---------------------------------------------------------------------------
	##
	## from snprma2.  Goal is to initialize A and B with
	## appropriate dimension for snps+nps
	##
	##---------------------------------------------------------------------------
	if (missing(sns)) sns <- basename(filenames)
	if (missing(cdfName))
		cdfName <- read.celfile.header(filenames[1])[["cdfName"]]
	pkgname <- getCrlmmAnnotationName(cdfName)
	stopifnot(require(pkgname, character.only=TRUE, quietly=!verbose))
	if(verbose) message("Loading annotations and mixture model parameters.")
	obj <- loader("preprocStuff.rda", .crlmmPkgEnv, pkgname)
	pnsa <- getVarInEnv("pnsa")
	pnsb <- getVarInEnv("pnsb")
	gns <- getVarInEnv("gns")
	rm(list=obj, envir=.crlmmPkgEnv)
	rm(obj)
	if(verbose) message("Initializing objects.")
	mixtureParams <- initializeBigMatrix("crlmmMixt-", 4, length(filenames), "double")
	SNR <- initializeBigVector("crlmmSNR-", length(filenames), "double")
	SKW <- initializeBigVector("crlmmSKW-", length(filenames), "double")
	featureData <- getFeatureData(cdfName, copynumber=TRUE)
	nr <- nrow(featureData); nc <- length(sns)
	A <- initializeBigMatrix("crlmmA-", nr, length(filenames), "integer")
	B <- initializeBigMatrix("crlmmB-", nr, length(filenames), "integer")
	rownames(A) <- rownames(B) <- featureNames(featureData)
	sampleBatches <- splitIndicesByNode(seq(along=filenames))
	if(verbose) message("Processing ", length(filenames), " files.")
	ocLapply(sampleBatches, rsprocessCEL, filenames=filenames,
		 fitMixture=TRUE, A=A, B=B, SKW=SKW, SNR=SNR,
		 mixtureParams=mixtureParams, eps=eps, seed=seed,
		 mixtureSampleSize=mixtureSampleSize, pkgname=pkgname,
		 neededPkgs=c("crlmm", pkgname))
	## Now we initialize a CNSet object, cloning A and B
	if(verbose) message("Cloning A and B matrices to store genotype calls and confidence scores.")
	## The clones will be used to store calls and confidence scores
	open(A)
	open(B)
	##user do
	## options("ffbatchbytes"=XXX) to make this efficient
	call <- clone(A)
	callProbability=clone(B)
	close(A)
	close(B)
	cnSet <- new("CNSet",
		     alleleA=A,
		     alleleB=B,
		     call=call,
		     callProbability=callProbability,
		     featureData=featureData,
		     annotation=cdfName,
		     batch=batch)
	if(!missing(sns)){
		sampleNames(cnSet) <- sns
	} else {
		sampleNames(cnSet) <- basename(filenames)
	}
	protocolData <- getProtocolData.Affy(filenames)
	rownames(pData(protocolData)) <- sns
	protocolData(cnSet) <- protocolData
	##protocolData(cnSet) <- protocolData
	pd <- data.frame(matrix(NA, nc, 3), row.names=sns)
	colnames(pd)=c("SKW", "SNR", "gender")
	phenoData(cnSet) <- new("AnnotatedDataFrame", data=pd)
	stopifnot(validObject(cnSet))
	snp.I <- isSnp(cnSet)
	snp.index <- which(snp.I)
	pData(cnSet)$SKW <- SKW
	pData(cnSet)$SNR <- SNR
	np.index <- which(!snp.I)
	if(verbose) message("Normalizing nonpolymorphic markers")
	FUN <- ifelse(is.lds, "cnrma2", "cnrma")
	if(verbose) message("Quantile normalizing nonpolymorphic markers")
	cnrma2(A=A(cnSet),
	       filenames=filenames,
	       row.names=featureNames(cnSet)[np.index],
	       cdfName=cdfName,
	       sns=sampleNames(cnSet),
	       seed=seed,
	       verbose=verbose)
	tmp <- rscrlmmGT2(A=calls(cnSet),
			  B=snpCallProbability(cnSet),
			  SNR=SNR,
			  mixtureParams=mixtureParams,
			  cdfName=cdfName,
			  row.names=NULL,
			  col.names=sampleNames(cnSet),
			  SNRMin=SNRMin,
			  recallMin=recallMin,
			  recallRegMin=recallRegMin,
			  gender=gender,
			  verbose=verbose,
			  returnParams=returnParams,
			  badSNP=badSNP,
			  snp.names=featureNames(cnSet)[snp.index])
	if(verbose) message("Genotyping finished.  Updating container with genotype calls and confidence scores.")
	cnSet$gender <- tmp[["gender"]]
	return(cnSet)
}

rsprocessCEL <- function(i, filenames, fitMixture, A, B, SKW, SNR,
			 mixtureParams, eps, seed, mixtureSampleSize,
			 pkgname){
	obj1 <- loader("preprocStuff.rda", .crlmmPkgEnv, pkgname)
	obj2 <- loader("genotypeStuff.rda", .crlmmPkgEnv, pkgname)
	obj3 <- loader("mixtureStuff.rda", .crlmmPkgEnv, pkgname)
	autosomeIndex <- getVarInEnv("autosomeIndex")
	pnsa <- getVarInEnv("pnsa")
	pnsb <- getVarInEnv("pnsb")
	fid <- getVarInEnv("fid")
	reference <- getVarInEnv("reference")
	aIndex <- getVarInEnv("aIndex")
	bIndex <- getVarInEnv("bIndex")
	SMEDIAN <- getVarInEnv("SMEDIAN")
	theKnots <- getVarInEnv("theKnots")
	gns <- getVarInEnv("gns")
	rm(list=c(obj1, obj2, obj3), envir=.crlmmPkgEnv)
	rm(obj1, obj2, obj3)

	## for mixture
	set.seed(seed)
	idx <- sort(sample(autosomeIndex, mixtureSampleSize))
	##for skewness. no need to do everything
	idx2 <- sample(length(fid), 10^5)

	open(A)
	open(B)
	open(SKW)
	open(mixtureParams)
	open(SNR)

	## RS ADDED
	index <- match(gns, rownames(A))

	for (k in i){
		y <- as.matrix(read.celfile(filenames[k], intensity.means.only=TRUE)[["INTENSITY"]][["MEAN"]][fid])
		x <- log2(y[idx2])
		SKW[k] <- mean((x-mean(x))^3)/(sd(x)^3)
		rm(x)
		y <- normalize.quantiles.use.target(y, target=reference)
		## RS: add index for row assignment
		A[index, k] <- intMedianSummaries(y[aIndex, 1, drop=FALSE], pnsa)
		B[index, k] <- intMedianSummaries(y[bIndex, 1, drop=FALSE], pnsb)
		rm(y)
		if(fitMixture){
			S <- (log2(A[idx,k])+log2(B[idx, k]))/2 - SMEDIAN
			M <- log2(A[idx, k])-log2(B[idx, k])
			tmp <- fitAffySnpMixture56(S, M, theKnots, eps=eps)
			rm(S, M)
			mixtureParams[, k] <- tmp[["coef"]]
			SNR[k] <- tmp[["medF1"]]^2/(tmp[["sigma1"]]^2+tmp[["sigma2"]]^2)
			rm(tmp)
		} else {
			mixtureParams[, k] <- NA
			SNR[k] <- NA
		}
	}
	close(A)
	close(B)
	close(SKW)
	close(mixtureParams)
	close(SNR)
	rm(list=ls())
	gc(verbose=FALSE)
	TRUE
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
				    verbose, is.lds=TRUE){
	##if(is.lds) {physical <- get("physical"); open(object)}
	open(object)
	if(verbose) message("      Probe stratum ", strata, " of ", length(index.list))
	marker.index <- index.list[[strata]]
	batches <- split(seq_along(batch(object)), as.character(batch(object)))
	batches <- batches[sapply(batches, length) >= MIN.SAMPLES]
	batch.names <- names(batches)
	batch.index <- which(batchNames(object) %in% batch.names)
	N.AA <- as.matrix(N.AA(object)[marker.index, batch.index])
	N.AB <- as.matrix(N.AB(object)[marker.index, batch.index])
	N.BB <- as.matrix(N.BB(object)[marker.index, batch.index])
	medianA.AA <- as.matrix(medianA.AA(object)[marker.index, batch.index])
	medianA.AB <- as.matrix(medianA.AB(object)[marker.index, batch.index])
	medianA.BB <- as.matrix(medianA.BB(object)[marker.index, batch.index])
	medianB.AA <- as.matrix(medianB.AA(object)[marker.index, batch.index])
	medianB.AB <- as.matrix(medianB.AB(object)[marker.index, batch.index])
	medianB.BB <- as.matrix(medianB.BB(object)[marker.index, batch.index])
	madA.AA <- as.matrix(madA.AA(object)[marker.index, batch.index])
	madA.AB <- as.matrix(madA.AB(object)[marker.index, batch.index])
	madA.BB <- as.matrix(madA.BB(object)[marker.index, batch.index])
	madB.AA <- as.matrix(madB.AA(object)[marker.index, batch.index])
	madB.AB <- as.matrix(madB.AB(object)[marker.index, batch.index])
	madB.BB <- as.matrix(madB.BB(object)[marker.index, batch.index])
	medianA <- medianB <- shrink.madB <- shrink.madA <- vector("list", length(batch.names))
	shrink.tau2A.AA <- tau2A.AA <- as.matrix(tau2A.AA(object)[marker.index, batch.index])
	shrink.tau2B.BB <- tau2B.BB <- as.matrix(tau2B.BB(object)[marker.index, batch.index])
	shrink.tau2A.BB <- tau2A.BB <- as.matrix(tau2A.BB(object)[marker.index, batch.index])
	shrink.tau2B.AA <- tau2B.AA <- as.matrix(tau2B.AA(object)[marker.index, batch.index])
	shrink.corrAA <- corrAA <- as.matrix(corrAA(object)[marker.index, batch.index])
	shrink.corrAB <- corrAB <- as.matrix(corrAB(object)[marker.index, batch.index])
	shrink.corrBB <- corrBB <- as.matrix(corrBB(object)[marker.index, batch.index])
	flags <- as.matrix(flags(object)[marker.index, batch.index])
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
	flags(object)[marker.index, batch.index] <- flags
	medianA.AA(object)[marker.index, batch.index] <- do.call("cbind", lapply(medianA, function(x) x[, 1]))
	medianA.AB(object)[marker.index, batch.index] <- do.call("cbind", lapply(medianA, function(x) x[, 2]))
	medianA.BB(object)[marker.index, batch.index] <- do.call("cbind", lapply(medianA, function(x) x[, 3]))
	medianB.AA(object)[marker.index, batch.index] <- do.call("cbind", lapply(medianB, function(x) x[, 1]))
	medianB.AB(object)[marker.index, batch.index] <- do.call("cbind", lapply(medianB, function(x) x[, 2]))
	medianB.BB(object)[marker.index, batch.index] <- do.call("cbind", lapply(medianB, function(x) x[, 3]))
	##
	madA.AA(object)[marker.index, batch.index] <- do.call("cbind", lapply(shrink.madA, function(x) x[, 1]))
	madA.AB(object)[marker.index, batch.index] <- do.call("cbind", lapply(shrink.madA, function(x) x[, 2]))
	madA.BB(object)[marker.index, batch.index] <- do.call("cbind", lapply(shrink.madA, function(x) x[, 3]))
	madB.AA(object)[marker.index, batch.index] <- do.call("cbind", lapply(shrink.madB, function(x) x[, 1]))
	madB.AB(object)[marker.index, batch.index] <- do.call("cbind", lapply(shrink.madB, function(x) x[, 2]))
	madB.BB(object)[marker.index, batch.index] <- do.call("cbind", lapply(shrink.madB, function(x) x[, 3]))
	##
	corrAA(object)[marker.index, batch.index] <- shrink.corrAA
	corrAB(object)[marker.index, batch.index] <- shrink.corrAB
	corrBB(object)[marker.index, batch.index] <- shrink.corrBB
	tau2A.AA(object)[marker.index, batch.index] <- shrink.tau2A.AA
	tau2A.BB(object)[marker.index, batch.index] <- shrink.tau2A.BB
	tau2B.AA(object)[marker.index, batch.index] <- shrink.tau2B.AA
	tau2B.BB(object)[marker.index, batch.index] <- shrink.tau2B.BB
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
	open(object$gender)
	if(is.lds) {physical <- get("physical"); open(object)}
	if(verbose) message("      Probe stratum ", strata, " of ", length(index.list))
	snps <- index.list[[strata]]
	batches <- split(seq_along(batch(object)), as.character(batch(object)))
	batches <- batches[sapply(batches, length) >= MIN.SAMPLES]
	##batchnames <- batchNames(object)
	batch.names <- names(batches)
	batch.index <- which(batchNames(object) %in% batch.names)
	if(length(batches) > 1 && "grandMean" %in% batchNames(object))
		batch.index <- c(batch.index, match("grandMean", batchNames(object)))
	N.AA <- as.matrix(N.AA(object)[snps, batch.index])
	N.AB <- as.matrix(N.AB(object)[snps, batch.index])
	N.BB <- as.matrix(N.BB(object)[snps, batch.index])
	medianA.AA <- as.matrix(medianA.AA(object)[snps, batch.index])
	medianA.AB <- as.matrix(medianA.AB(object)[snps, batch.index])
	medianA.BB <- as.matrix(medianA.BB(object)[snps, batch.index])
	medianB.AA <- as.matrix(medianB.AA(object)[snps, batch.index])
	medianB.AB <- as.matrix(medianB.AB(object)[snps, batch.index])
	medianB.BB <- as.matrix(medianB.BB(object)[snps, batch.index])
	madA.AA <- as.matrix(madA.AA(object)[snps, batch.index])
	madA.AB <- as.matrix(madA.AB(object)[snps, batch.index])
	madA.BB <- as.matrix(madA.BB(object)[snps, batch.index])
	madB.AA <- as.matrix(madB.AA(object)[snps, batch.index])
	madB.AB <- as.matrix(madB.AB(object)[snps, batch.index])
	madB.BB <- as.matrix(madB.BB(object)[snps, batch.index])
	tau2A.AA <- as.matrix(tau2A.AA(object)[snps, batch.index])
	tau2B.BB <- as.matrix(tau2B.BB(object)[snps, batch.index])
	tau2A.BB <- as.matrix(tau2A.BB(object)[snps, batch.index])
	tau2B.AA <- as.matrix(tau2B.AA(object)[snps, batch.index])
	corrAA <- as.matrix(corrAA(object)[snps, batch.index])
	corrAB <- as.matrix(corrAB(object)[snps, batch.index])
	corrBB <- as.matrix(corrBB(object)[snps, batch.index])
	nuA <- as.matrix(nuA(object)[snps,  batch.index])
	phiA <- as.matrix(phiA(object)[snps, batch.index])
	nuB <- as.matrix(nuB(object)[snps, batch.index])
	phiB <- as.matrix(phiB(object)[snps, batch.index])
	flags <- as.matrix(flags(object)[snps, batch.index])
	for(k in seq(along=batches)){
		B <- batches[[k]]
		if(length(B) < MIN.SAMPLES) next()
		this.batch <- unique(as.character(batch(object)[B]))
		medianA <- cbind(medianA.AA[, k], medianA.AB[, k], medianA.BB[, k])
		medianB <- cbind(medianB.AA[, k], medianB.AB[, k], medianB.BB[, k])
		madA <- cbind(madA.AA[, k], madA.AB[, k], madA.BB[, k])
		madB <- cbind(madB.AA[, k], madB.AB[, k], madB.BB[, k])
		NN <- cbind(N.AA[, k], N.AB[, k], N.BB[, k])
		## regress on the medians using the standard errors (hence the division by N) as weights
		res <- fit.wls(NN=NN, sigma=madA, allele="A", Y=medianA, autosome=!CHR.X)
		nuA[, k] <- res[1, ]
		phiA[, k] <- res[2, ]
		rm(res)
		res <- fit.wls(NN=NN, sigma=madB, allele="B", Y=medianB, autosome=!CHR.X)##allele="B", Ystar=YB, W=wB, Ns=Ns)
		nuB[, k] <- res[1, ]
		phiB[, k] <- res[2, ]
	}
	##---------------------------------------------------------------------------
	##
	## Grand mean
	##
	##---------------------------------------------------------------------------
##	if(length(batches) > 1 && "grandMean" %in% batchNames(object)){
	##  TODO:  There are NA's in the medianA.AA for the grandMean and 0's in the madA
	##         - both need to be handled prior to estimating a grand intercept and slope
	if(FALSE){
		## then the last column is for the grandMean
		k <- ncol(medianA.AA)
		medianA <- cbind(medianA.AA[, k], medianA.AB[, k], medianA.BB[, k])
		medianB <- cbind(medianB.AA[, k], medianB.AB[, k], medianB.BB[, k])
		madA <- cbind(madA.AA[, k], madA.AB[, k], madA.BB[, k])
		madB <- cbind(madB.AA[, k], madB.AB[, k], madB.BB[, k])
		NN <- cbind(N.AA[, k], N.AB[, k], N.BB[, k])
		index <- which(rowSums(is.na(medianA)) == 0)
		res <- fit.wls(NN=NN[index, ], sigma=madA[index, ], allele="A", Y=medianA[index, ], autosome=!CHR.X)
		nuA[, k] <- res[1, ]
		phiA[, k] <- res[2, ]
		rm(res)
		res <- fit.wls(NN=NN, sigma=madB, allele="B", Y=medianB, autosome=!CHR.X)##allele="B", Ystar=YB, W=wB, Ns=Ns)
		nuB[, k] <- res[1, ]
		phiB[, k] <- res[2, ]
	}
	if(THR.NU.PHI){
		nuA[nuA < MIN.NU] <- MIN.NU
		nuB[nuB < MIN.NU] <- MIN.NU
		phiA[phiA < MIN.PHI] <- MIN.PHI
		phiB[phiB < MIN.PHI] <- MIN.PHI
	}
	nuA(object)[snps, batch.index] <- nuA
	nuB(object)[snps, batch.index] <- nuB
	phiA(object)[snps, batch.index] <- phiA
	phiB(object)[snps, batch.index] <- phiB
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
	open(object$gender)
	if(verbose) message("      Probe stratum ", strata, " of ", length(index.list))
	marker.index <- index.list[[strata]]
	batches <- split(seq_along(batch(object)), as.character(batch(object)))
	batches <- batches[sapply(batches, length) >= MIN.SAMPLES]
	batch.names <- names(batches)
	batch.index <- which(batchNames(object) %in% batch.names)
	##
	ii <- isSnp(object) & chromosome(object) < 23 & !is.na(chromosome(object))
	flags <- as.matrix(flags(object)[ii, batch.index])
	fns <- featureNames(object)[ii]
	fns.noflags <- fns[rowSums(flags, na.rm=T) == 0]
	snp.index <- sample(match(fns.noflags, featureNames(object)), 5000)
	##
	nuA.np <- as.matrix(nuA(object)[marker.index, batch.index])
	phiA.np <- as.matrix(phiA(object)[marker.index, batch.index])
	tau2A.AA <- as.matrix(tau2A.AA(object)[marker.index, batch.index])
	##
	nuA.snp <- as.matrix(nuA(object)[snp.index, batch.index])
	nuB.snp <- as.matrix(nuB(object)[snp.index, batch.index])
	phiA.snp <- as.matrix(phiA(object)[snp.index, batch.index])
	phiB.snp <- as.matrix(phiB(object)[snp.index, batch.index])
	medianA.AA <- as.matrix(medianA.AA(object)[snp.index, batch.index])
	medianB.BB <- as.matrix(medianB.BB(object)[snp.index, batch.index])
	medianA.AA.np <- as.matrix(medianA.AA(object)[marker.index, batch.index])
	for(k in seq_along(batches)){
		B <- batches[[k]]
		this.batch <- unique(as.character(batch(object)[B]))
		X <- cbind(1, log2(c(medianA.AA[, k], medianB.BB[, k])))
		Y <- log2(c(phiA.snp[, k], phiB.snp[, k]))
		not.missing <- rowSums(is.na(X)) == 0 & !is.na(Y)
		X <- X[not.missing, ]
		Y <- Y[not.missing]
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
	nuA(object)[marker.index, batch.index] <- nuA.np
	phiA(object)[marker.index, batch.index] <- phiA.np
	if(is.lds) { close(object); return(TRUE)}
	return(object)
}

summarizeMaleXNps <- function(marker.index,
			      batches,
			      object, MIN.SAMPLES){
	nr <- length(marker.index)
	nc <- length(batchNames(object))
	NN.Mlist <- imputed.medianA <- imputed.medianB <- shrink.madA <- shrink.madB <- vector("list", nc)
	open(object$gender)
	gender <- object$gender[]
	open(A(object))
	AA <- as.matrix(A(object)[marker.index, gender==1])
	madA.AA <- medianA.AA <- matrix(NA, nr, nc)
	colnames(madA.AA) <- colnames(medianA.AA) <- batchNames(object)
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
		medianA.AA[, k] <- rowMedians(AA[, J, drop=FALSE], na.rm=TRUE)
		madA.AA[, k] <- rowMAD(AA[, J, drop=FALSE], na.rm=TRUE)
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
	I <- unlist(batches)
	open(object$gender)
	if(gender == "male"){
		gender.index <- which(object$gender[] == 1)
	} else {
		gender.index <- which(object$gender[] == 2)
	}
	batches <- lapply(batches, function(x, gender.index) intersect(x, gender.index), gender.index)
	batch.names <- names(batches)
	batch.index <- which(batchNames(object) %in% batch.names)
	nr <- length(marker.index)
	nc <- length(batch.index)
	NN.Mlist <- medianA <- medianB <- shrink.madA <- shrink.madB <- vector("list", nc)
	names(NN.Mlist) <- names(medianA) <- names(medianB) <- names(shrink.madA) <- names(shrink.madB) <- batch.names
	open(calls(object))
	open(snpCallProbability(object))
	open(A(object))
	open(B(object))
	GG <- as.matrix(calls(object)[marker.index, ])
	CP <- as.matrix(snpCallProbability(object)[marker.index, ])
	AA <- as.matrix(A(object)[marker.index, ])
	BB <- as.matrix(B(object)[marker.index, ])
	for(k in seq_along(batches)){
		sample.index <- batches[[k]]
		if(length(sample.index) == 0) next()
		this.batch <- unique(as.character(batch(object)[sample.index]))
		##gender <- object$gender[sample.index]
		if(length(sample.index) < MIN.SAMPLES) next()
		sns.batch <- sampleNames(object)[sample.index]
		##subset GG apppriately
		sns <- colnames(GG)
		J <- sns%in%sns.batch
		G <- GG[, J, drop=FALSE]
		xx <- CP[, J, drop=FALSE]
		p <- i2p(xx)
		highConf <- (1-exp(-xx/1000)) > GT.CONF.THR
		G <- G*highConf
		A <- AA[, J, drop=FALSE]
		B <- BB[, J, drop=FALSE]
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
			medianA[[k]] <- res[[1]]
			medianB[[k]] <- res[[2]]

			## RS: For Monomorphic SNPs a mixture model may be better
			## RS: Further, we can improve estimation by borrowing strength across batch
			unobserved.index[[1]] <- which(unobservedAA & unobservedAB)
			unobserved.index[[2]] <- which(unobservedBB & unobservedAB)
			unobserved.index[[3]] <- which(unobservedAA & unobservedBB) ## strange
			res <- imputeCentersForMonomorphicSnps(medianA[[k]], medianB[[k]],
							       index.complete,
							       unobserved.index)
			medianA[[k]] <- res[[1]]; medianB[[k]] <- res[[2]]
		}
		medianA[[k]] <- res[[1]]
		medianB[[k]] <- res[[2]]
	}
	close(calls(object))
	close(snpCallProbability(object))
	close(A(object))
	close(B(object))
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
	##if(is.lds) {physical <- get("physical"); open(object)}
	if(verbose) message("      Probe stratum ", strata, " of ", length(index.list))
	open(object$gender)
	gender <- object$gender[]
	enough.males <- sum(gender==1) >= MIN.SAMPLES
	enough.females <- sum(gender==2) >= MIN.SAMPLES
	if(!enough.males & !enough.females){
		message(paste("fewer than", MIN.SAMPLES, "men and women.  Copy number not estimated for CHR X"))
		return(object)
	}
	marker.index <- index.list[[strata]]
	batches <- split(seq_along(batch(object)), as.character(batch(object)))
	batches <- batches[sapply(batches, length) >= MIN.SAMPLES]
	batch.names <- names(batches)
	batch.index <- which(batchNames(object) %in% batch.names)

	nuA <- as.matrix(nuA(object)[marker.index, batch.index])
	nuB <- as.matrix(nuB(object)[marker.index, batch.index])
	phiA <- as.matrix(phiA(object)[marker.index, batch.index])
	phiB <- as.matrix(phiB(object)[marker.index, batch.index])
	phiA2 <- as.matrix(phiPrimeA(object)[marker.index, batch.index])
	phiB2 <- as.matrix(phiPrimeB(object)[marker.index, batch.index])
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
		if(is.null(sample.index)) next()
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
			if(any(madA.F == 0)){
				message("Zeros in median absolute deviation.  Not estimating CN for chrX for batch ", names(batches)[k])
				next()
			}
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
			if(any(madA.F == 0)){
				message("Zeros in median absolute deviation.  Not estimating CN for chrX for batch ", names(batches)[k])
				next()
			}
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
	nuA(object)[marker.index, batch.index] <- nuA
	nuB(object)[marker.index, batch.index] <- nuB
	phiA(object)[marker.index, batch.index] <- phiA
	phiB(object)[marker.index, batch.index] <- phiB
	phiPrimeA(object)[marker.index, batch.index] <- phiA2
	phiPrimeB(object)[marker.index, batch.index] <- phiB2
	if(is.lds) {close(object); return(TRUE)} else return(object)
}

##nonpolymorphic markers on X
fit.lm4 <- function(strata,
		    index.list,
		    object,
		    MIN.SAMPLES,
		    THR.NU.PHI,
		    MIN.NU,
		    MIN.PHI,
		    verbose, is.lds=TRUE, ...){
	##if(is.lds) {physical <- get("physical"); open(object)}
	## exclude batches that have fewer than MIN.SAMPLES
	open(object$gender)
	gender <- object$gender[]
	enough.males <- sum(gender==1) >= MIN.SAMPLES
	enough.females <- sum(gender==2) >= MIN.SAMPLES
	if(!enough.males & !enough.females){
		message(paste("fewer than", MIN.SAMPLES, "men and women.  Copy number not estimated for CHR X"))
		return(object)
	}
	excludeBatch <- names(table(batch(object)))[table(batch(object)) < MIN.SAMPLES]
	if(length(excludeBatch) > 0){
		bns <- unique(batch(object))
		batch.index <- which(!bns %in% excludeBatch)
		sample.index <- which(!batch(object)%in% excludeBatch)
	} else {
		sample.index <- 1:ncol(object)
		batch.index <- as.integer(as.factor(unique(batch(object))))
	}
	if(verbose) message("      Probe stratum ", strata, " of ", length(index.list))
	marker.index <- index.list[[strata]]
	batches <- split(seq_along(batch(object)), as.character(batch(object)))
	batches <- batches[batch.index]
	nc <- length(batch.index)
	if(enough.males){
		res <- summarizeMaleXNps(marker.index=marker.index,
					 batches=batches,
					 object=object, MIN.SAMPLES=MIN.SAMPLES)
		medianA.AA.M <- res[["medianA.AA"]][, batch.index, drop=FALSE]
		madA.AA.M <- res[["madA.AA"]][, batch.index, drop=FALSE]
	}
	medianA.AA.F <- as.matrix(medianA.AA(object)[marker.index, batch.index, drop=FALSE]) ## median for women
	madA.AA.F <- as.matrix(madA.AA(object)[marker.index, batch.index, drop=FALSE]) ## median for women
	split.gender <- split(gender[sample.index], as.character(batch(object)[sample.index]))
	N.M <- sapply(split.gender, function(x) sum(x==1))
	N.F <- sapply(split.gender, function(x) sum(x==2))
	nuA <- as.matrix(nuA(object)[marker.index, batch.index, drop=FALSE])
	nuB <- as.matrix(nuB(object)[marker.index, batch.index, drop=FALSE])
	phiA <- as.matrix(phiA(object)[marker.index, batch.index, drop=FALSE])
	phiB <- as.matrix(phiB(object)[marker.index, batch.index, drop=FALSE])
	ii <- isSnp(object) & chromosome(object) < 23 & !is.na(chromosome(object))
	fns <- featureNames(object)[ii]
	flags <- as.matrix(flags(object)[ii, batch.index])
	fns.noflags <- fns[rowSums(flags, na.rm=T) == 0]
	snp.index <- sample(match(fns.noflags, featureNames(object)), 10000)
	##
	N.AA <- as.matrix(N.AA(object)[snp.index, batch.index, drop=FALSE])
	N.AB <- as.matrix(N.AB(object)[snp.index, batch.index, drop=FALSE])
	N.BB <- as.matrix(N.BB(object)[snp.index, batch.index, drop=FALSE])
	enoughAA <- rowSums(N.AA < 5) == 0
	enoughAB <- rowSums(N.AB < 5) == 0
	enoughBB <- rowSums(N.BB < 5) == 0
	snp.index <- snp.index[enoughAA & enoughAB & enoughBB]
	if(length(snp.index) < 100){
		message("too few snps pass criteria for estimating parameters for NP markers on chr X")
		return(object)
	}
	nuA.snp <- as.matrix(nuA(object)[snp.index, batch.index, drop=FALSE])
	nuA.snp.notmissing <- rowSums(is.na(nuA.snp)) == 0
	nuA.snp.notnegative <- rowSums(as.matrix(nuA(object)[snp.index, batch.index, drop=FALSE]) < 20) == 0
	snp.index <- snp.index[nuA.snp.notmissing & nuA.snp.notnegative]
	medianA.AA.snp <- as.matrix(medianA.AA(object)[snp.index, batch.index])
	medianB.BB.snp <- as.matrix(medianB.BB(object)[snp.index, batch.index])
	nuA.snp <- as.matrix(nuA(object)[snp.index, batch.index])
	nuB.snp <- as.matrix(nuB(object)[snp.index, batch.index])
	phiA.snp <- as.matrix(phiA(object)[snp.index, batch.index])
	phiB.snp <- as.matrix(phiB(object)[snp.index, batch.index])
##	pseudoAR <- position(object)[snp.index] < 2709520 | (position(object)[snp.index] > 154584237 & position(object)[snp.index] < 154913754)
##	pseudoAR[is.na(pseudoAR)] <- FALSE
	for(k in seq_along(batches)){
		B <- batches[[k]]
		this.batch <- unique(as.character(batch(object)[B]))
		gender <- object$gender[B]
		enough.men <- N.M[[this.batch]] >= MIN.SAMPLES
		enough.women <- N.F[[this.batch]] >= MIN.SAMPLES
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
		notmissing.index <- which(rowSums(is.na(X)) == 0 & !is.na(Y))
		X <- X[notmissing.index, ]
		Y <- Y[notmissing.index]
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
	nuA(object)[marker.index, batch.index] <- nuA
	phiA(object)[marker.index, batch.index] <- phiA
	nuB(object)[marker.index, batch.index] <- nuB
	phiB(object)[marker.index, batch.index] <- phiB
	##if(is.lds) {close(object); return(TRUE)} else return(object)
	TRUE
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
	open(A)
	##idx2 <- sample(length(fid), 10^5)
	for (k in i){
		y <- as.matrix(read.celfile(filenames[k], intensity.means.only=TRUE)[["INTENSITY"]][["MEAN"]][fid])
		A[np.index, k] <- as.integer(normalize.quantiles.use.target(y, target=reference))
	}
	close(A)
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
			  type="SNP",
			  MIN.OBS=1,
			  MIN.SAMPLES=10,
			  DF.PRIOR=50,
			  verbose=TRUE,
			  marker.index,
			  is.lds=TRUE){
	stopifnot(type[[1]] == "SNP")
	CHR.X <- FALSE ## this is no longer needed
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
	TRUE
}

genotypeSummary <- function(object,
			    GT.CONF.THR=0.80,
			    type=c("SNP", "NP", "X.SNP", "X.NP"), ##"X.snps", "X.nps"),
			    verbose=TRUE,
			    marker.index,
			    is.lds){
	if(type == "X.SNP" | type=="X.NP"){
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
	TRUE
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
	open(A(object))
	open(object$gender)
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
		I.A <- AA[, sample.index, drop=FALSE]
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
			  verbose, is.lds=TRUE, CHR.X, ...){
##	if(is.lds) {
##		physical <- get("physical")
##		open(object)
##	}
	if(verbose) message("      Probe stratum ", strata, " of ", length(index.list))
	index <- index.list[[strata]]
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
	summaryStats <- function(X, INT, FUNS){
		tmp <- matrix(NA, nrow(X), length(FUNS))
		for(j in seq_along(FUNS)){
			FUN <- match.fun(FUNS[j])
			tmp[, j] <- FUN(X*INT, na.rm=TRUE)
		}
		tmp
	}
	if(verbose) message("        Computing summaries...")
	for(k in seq_along(batches)){
		##note that the genotype frequency for AA would include 'A' on chromosome X for men
		sample.index <- batches[[k]]
		this.batch <- unique(as.character(batch(object)[sample.index]))
		j <- match(this.batch, batchnames)
		G <- GG[, sample.index, drop=FALSE]
		xx <- CP[, sample.index, drop=FALSE]
		highConf <- (1-exp(-xx/1000)) > GT.CONF.THR
		## Some markers may have genotype confidences scores that are ALL below the threshold
		## For these markers, provide statistical summaries based on all the samples and flag
		## Provide summaries for these markers and flag to indicate the problem
		## RS: check whether we need the next to lines and, if so, effect downstream
		ii <- which(rowSums(highConf) == 0)
		if(length(ii) > 0) highConf[ii, ] <- TRUE
		G <- G*highConf
		## table(rowSums(G==0))
		##G <- G*highConf*NORM
		A <- AA[, sample.index, drop=FALSE]
		B <- BB[, sample.index, drop=FALSE]
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
		stats <- summaryStats(G.AA, A, FUNS=c("rowMedians", "rowMAD"))
		medianA.AA(object)[index, k] <- stats[, 1]
		##missing.index <- which(rowSums(is.na(medianA.AA(object)[index, ,drop=FALSE])) > 0)
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
		rm(A, B, G.AA, G.AB, G.BB, xx, highConf, G)
		gc()
	}
	##---------------------------------------------------------------------------
	## grand mean
	##---------------------------------------------------------------------------
	if(length(batches) > 1 && "grandMean" %in% batchNames(object)){
		##k <- k+1
		k <- match("grandMean", batchNames(object))
		if(verbose) message("        computing grand means...")
		##G <- GG[, B]
		##NORM <- normal.index[, k]
		##xx <- CP[, B]
		##highConf <- (1-exp(-xx/1000)) > GT.CONF.THR
		highConf <- (1-exp(-CP/1000)) > GT.CONF.THR
		##G <- G*highConf
		## Some markers may have genotype confidences scores that are ALL below the threshold
		## For these markers, provide statistical summaries based on all the samples and flag
		## Provide summaries for these markers and flag to indicate the problem
		ii <- which(rowSums(highConf) == 0)
		if(length(ii) > 0) highConf[ii, ] <- TRUE
		GG <- GG*highConf
		## table(rowSums(G==0))
		##G <- G*highConf*NORM
##		A <- AA[, B]
##		B <- BB[, B]
		## this can be time consuming...do only once
##		G.AA <- G==1
		G.AA <- GG==1
		G.AA[G.AA==FALSE] <- NA
		G.AB <- GG==2
		G.AB[G.AB==FALSE] <- NA
		G.BB <- GG==3
		G.BB[G.BB==FALSE] <- NA
		Ns.AA[, k] <- rowSums(G.AA, na.rm=TRUE)
		Ns.AB[, k] <- rowSums(G.AB, na.rm=TRUE)
		Ns.BB[, k] <- rowSums(G.BB, na.rm=TRUE)
		## Calculate row medians and MADs
		stats <- summaryStats(G.AA, AA, FUNS=c("rowMedians", "rowMAD"))
		medianA.AA(object)[index, k] <- stats[, 1]
		madA.AA(object)[index, k] <- stats[, 2]
		stats <- summaryStats(G.AB, AA, FUNS=c("rowMedians", "rowMAD"))
		medianA.AB(object)[index, k] <- stats[, 1]
		madA.AB(object)[index, k] <- stats[, 2]
		stats <- summaryStats(G.BB, AA, FUNS=c("rowMedians", "rowMAD"))
		medianA.BB(object)[index, k] <- stats[, 1]
		madA.BB(object)[index, k] <- stats[, 2]
		stats <- summaryStats(G.AA, BB, FUNS=c("rowMedians", "rowMAD"))
		medianB.AA(object)[index, k] <- stats[, 1]
		madB.AA(object)[index, k] <- stats[, 2]
		stats <- summaryStats(G.AB, BB, FUNS=c("rowMedians", "rowMAD"))
		medianB.AB(object)[index, k] <- stats[, 1]
		madB.AB(object)[index, k] <- stats[, 2]
		stats <- summaryStats(G.BB, BB, FUNS=c("rowMedians", "rowMAD"))
		medianB.BB(object)[index, k] <- stats[, 1]
		madB.BB(object)[index, k] <- stats[, 2]
		## log2 Transform Intensities
		AA <- log2(AA); BB <- log2(BB)
		tau2A.AA[, k] <- summaryStats(G.AA, AA, FUNS="rowMAD")^2
		tau2A.BB[, k] <- summaryStats(G.BB, AA, FUNS="rowMAD")^2
		tau2B.AA[, k] <- summaryStats(G.AA, BB, FUNS="rowMAD")^2
		tau2B.BB[, k] <- summaryStats(G.BB, BB, FUNS="rowMAD")^2
		corrAA[, k] <- rowCors(AA*G.AA, BB*G.AA, na.rm=TRUE)
		corrAB[, k] <- rowCors(AA*G.AB, BB*G.AB, na.rm=TRUE)
		corrBB[, k] <- rowCors(AA*G.BB, BB*G.BB, na.rm=TRUE)
		##
		## TODO:   fill in NAs -- use code from shrinkGenotypeSummaries
		##
		rm(GG, CP, AA, BB, FL, stats)
		gc()
	}
	if(verbose) message("        Begin writing...")
	N.AA(object)[index,] <- Ns.AA
	N.AB(object)[index,] <- Ns.AB
	N.BB(object)[index,] <- Ns.BB
	corrAA(object)[index,] <- corrAA
	corrAB(object)[index,] <- corrAB
	corrBB(object)[index,] <- corrBB
	tau2A.AA(object)[index, ] <- tau2A.AA
	tau2A.BB(object)[index, ] <- tau2A.BB
	tau2B.AA(object)[index, ] <- tau2B.AA
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
	if(GT.CONF.THR >= 1 | GT.CONF.THR < 0) stop("GT.CONF.THR must be in [0,1)")
	batch <- batch(object)
	is.snp <- isSnp(object)
	is.autosome <- chromosome(object) < 23
	is.annotated <- !is.na(chromosome(object))
	is.X <- chromosome(object) == 23
	is.lds <- is(calls(object), "ffdf") | is(calls(object), "ff_matrix")
	if(is.lds) stopifnot(isPackageLoaded("ff"))
	samplesPerBatch <- table(as.character(batch(object)))
	open(object)
	if(any(samplesPerBatch < MIN.SAMPLES)){
		tmp <- paste(names(samplesPerBatch)[samplesPerBatch < MIN.SAMPLES], collapse=", ")
		message("The following batches have fewer than ", MIN.SAMPLES, " samples: ",  tmp)
		message("Not estimating copy number parameters for batch ", tmp)
	}
	mylabel <- function(marker.type){
		switch(marker.type,
		       SNP="autosomal SNPs",
		       NP="autosomal nonpolymorphic markers",
		       X.SNP="chromosome X SNPs",
		       X.NP="chromosome X nonpolymorphic markers")
	}
	if(verbose) message("Computing summary statistics of the genotype clusters for each batch")
	I <- which(type %in% c("SNP", "NP", "X.NP"))
	for(j in seq_along(I)){ ## do not do X.SNP.  Do this during fit.lm3
		i <- I[j]
		marker.type <- type[i]
		if(verbose) message(paste("...", mylabel(marker.type)))
		##if(verbose) message(paste("Computing summary statistics for ", mylabel(marker.type), " genotype clusters for each batch")
		marker.index <- whichMarkers(marker.type, is.snp,
					     is.autosome, is.annotated, is.X)
		genotypeSummary(object=object,
				GT.CONF.THR=GT.CONF.THR,
				type=marker.type,
				verbose=verbose,
				marker.index=marker.index,
				is.lds=is.lds)
	}
	if(verbose) message("Imputing unobserved genotype medians and shrinking the variances (within-batch, across loci) ")##SNPs only
	if("SNP" %in% type){
		marker.index <- whichMarkers("SNP", is.snp,
					     is.autosome, is.annotated, is.X)
		shrinkSummary(object=object,
			      MIN.OBS=MIN.OBS,
			      MIN.SAMPLES=MIN.SAMPLES,
			      DF.PRIOR=DF.PRIOR,
			      type="SNP",
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
	close(object)
	TRUE
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


imputeAA.AB <- function(index, N.AA, N.AB, N.BB,
			medianA.AA, medianA.AB, medianA.BB,
			medianB.AA, medianB.AB, medianB.BB){
	gt.to.impute <- 1:2
	imputed <- rep(FALSE, length(index))
	for(i in index){
		Ns <- cbind(N.AA[i, ], N.AB[i, ], N.BB[i, ])
		medianA <- cbind(medianA.AA[i, ], medianA.AB[i, ], medianA.BB[i, ])
		medianB <- cbind(medianB.AA[i, ], medianB.AB[i, ], medianB.BB[i, ])
		colnames(medianA) <- colnames(medianB) <- c("AA", "AB", "BB")

		##Find batches with sufficient data
		K <- which(rowSums(Ns < 3) == 0)
		if(length(K) < 1) next() ## nothing we can do.  use information from other snps
		X <- cbind(1, medianA[K, -gt.to.impute, drop=FALSE], medianB[K, -gt.to.impute, drop=FALSE])
		Y <- cbind(medianA[K, gt.to.impute, drop=FALSE], medianB[K, gt.to.impute, drop=FALSE])
		tmp <- tryCatch(betahat <- solve(crossprod(X), crossprod(X,Y)), error=function(e) NULL)
		if(is.null(tmp)) {
			##cat(".");
			next()
		}
		## Get data from observed genotypes in batch with insufficient data for genotypes 'gt.to.impute'
		L <- which(rowSums(Ns < 3) == 2)
		if(length(L) == 0) next()
		Z <- cbind(1, medianA[L, -gt.to.impute, drop=FALSE], medianB[L, -gt.to.impute, drop=FALSE])
		imputedVals <- Z %*% betahat
		medianA.AA[i, L] <- imputedVals[, 1]
		medianA.AB[i, L] <- imputedVals[, 2]
		medianB.AA[i, L] <- imputedVals[, 3]
		medianB.AB[i, L] <- imputedVals[, 4]
		imputed[i] <- TRUE
	}
	return(list(medianA.AA=medianA.AA, medianA.AB=medianA.AB, medianA.BB=medianA.BB,
		    medianB.AA=medianB.AA, medianB.AB=medianB.AB, medianB.BB=medianB.BB,
		    imputed=imputed))
}

imputeAB.BB <- function(index, N.AA, N.AB, N.BB,
			medianA.AA, medianA.AB, medianA.BB,
			medianB.AA, medianB.AB, medianB.BB){
	gt.to.impute <- 2:3
	imputed <- rep(FALSE, length(index))
	for(j in seq_along(index)){
		i <- index[j]
		Ns <- cbind(N.AA[i, ], N.AB[i, ], N.BB[i, ])
		medianA <- cbind(medianA.AA[i, ], medianA.AB[i, ], medianA.BB[i, ])
		medianB <- cbind(medianB.AA[i, ], medianB.AB[i, ], medianB.BB[i, ])
		colnames(medianA) <- colnames(medianB) <- c("AA", "AB", "BB")

		##Find batches with sufficient data
		K <- which(rowSums(Ns < 3) == 0)
		if(length(K) < 1) next() ## nothing we can do.  use information from other snps
		X <- cbind(1, medianA[K, -gt.to.impute, drop=FALSE], medianB[K, -gt.to.impute, drop=FALSE])
		Y <- cbind(medianA[K, gt.to.impute, drop=FALSE], medianB[K, gt.to.impute, drop=FALSE])
		tmp <- tryCatch(betahat <- solve(crossprod(X), crossprod(X,Y)), error=function(e) NULL)
		if(is.null(tmp)) {
##			cat(".");
			next()
		}

		## Get data from observed genotypes in batch with insufficient data for genotypes 'gt.to.impute'
		L <- which(rowSums(Ns < 3) == 2)
		if(length(L) == 0) next()
		Z <- cbind(1, medianA[L, -gt.to.impute, drop=FALSE], medianB[L, -gt.to.impute, drop=FALSE])
		imputedVals <- Z %*% betahat
		medianA.AB[i, L] <- imputedVals[, 1]
		medianA.BB[i, L] <- imputedVals[, 2]
		medianB.AB[i, L] <- imputedVals[, 3]
		medianB.BB[i, L] <- imputedVals[, 4]
		imputed[j] <- TRUE
	}
	return(list(medianA.AA=medianA.AA, medianA.AB=medianA.AB, medianA.BB=medianA.BB,
		    medianB.AA=medianB.AA, medianB.AB=medianB.AB, medianB.BB=medianB.BB,
		    imputed=imputed))
}

imputeAA <- function(index, N.AA, N.AB, N.BB,
		     medianA.AA, medianA.AB, medianA.BB,
		     medianB.AA, medianB.AB, medianB.BB){
	gt.to.impute <- 1
	imputed <- rep(FALSE, length(index))
	for(j in seq_along(index)){
		i <- index[j]
		Ns <- cbind(N.AA[i, ], N.AB[i, ], N.BB[i, ])
		medianA <- cbind(medianA.AA[i, ], medianA.AB[i, ], medianA.BB[i, ])
		medianB <- cbind(medianB.AA[i, ], medianB.AB[i, ], medianB.BB[i, ])
		colnames(medianA) <- colnames(medianB) <- c("AA", "AB", "BB")

		##Find batches with sufficient data
		K <- which(rowSums(Ns < 3) == 0)
		if(length(K) < 1) next() ## nothing we can do.  use information from other snps
		X <- cbind(1, medianA[K, -gt.to.impute, drop=FALSE], medianB[K, -gt.to.impute, drop=FALSE])
		Y <- cbind(medianA[K, gt.to.impute, drop=FALSE], medianB[K, gt.to.impute, drop=FALSE])
		tmp <- tryCatch(betahat <- solve(crossprod(X), crossprod(X,Y)), error=function(e) NULL)
		if(is.null(tmp)) {
			##cat(".");
			next()
		}
		## Get data from observed genotypes in batch with insufficient data for genotypes 'gt.to.impute'
		L <- which(rowSums(Ns < 3) == 1)
		if(length(L) == 0) next()
		Z <- cbind(1, medianA[L, -gt.to.impute, drop=FALSE], medianB[L, -gt.to.impute, drop=FALSE])
		imputedVals <- Z %*% betahat
		medianA.AA[i, L] <- imputedVals[, 1]
		medianB.AA[i, L] <- imputedVals[, 2]
		imputed[j] <- TRUE
	}
	return(list(medianA.AA=medianA.AA, medianA.AB=medianA.AB, medianA.BB=medianA.BB,
		    medianB.AA=medianB.AA, medianB.AB=medianB.AB, medianB.BB=medianB.BB,
		    imputed=imputed))
}

imputeBB <- function(index, N.AA, N.AB, N.BB,
		     medianA.AA, medianA.AB, medianA.BB,
		     medianB.AA, medianB.AB, medianB.BB){
	gt.to.impute <- 3
	imputed <- rep(FALSE, length(index))
	for(j in seq_along(index)){
		i <- index[j]
		Ns <- cbind(N.AA[i, ], N.AB[i, ], N.BB[i, ])
		medianA <- cbind(medianA.AA[i, ], medianA.AB[i, ], medianA.BB[i, ])
		medianB <- cbind(medianB.AA[i, ], medianB.AB[i, ], medianB.BB[i, ])
		colnames(medianA) <- colnames(medianB) <- c("AA", "AB", "BB")

		##Find batches with sufficient data
		K <- which(rowSums(Ns < 3) == 0)
		if(length(K) < 1) next() ## nothing we can do.  use information from other snps
		X <- cbind(1, medianA[K, -gt.to.impute, drop=FALSE], medianB[K, -gt.to.impute, drop=FALSE])
		Y <- cbind(medianA[K, gt.to.impute, drop=FALSE], medianB[K, gt.to.impute, drop=FALSE])
		colnames(Y) <- c("A.BB", "B.BB")
		tmp <- tryCatch(betahat <- solve(crossprod(X), crossprod(X,Y)), error=function(e) NULL)
		if(is.null(tmp)) {
			##cat(".");
			next()
		}
##		else{
##			R <- Y-crossprod(X, betahat)
##			RSS <- t(R)%*%R
##		}
		## Get data from observed genotypes in batch with insufficient data for genotypes 'gt.to.impute'
		L <- which(rowSums(Ns < 3) == 1)
		if(length(L) == 0) next()
		Z <- cbind(1, medianA[L, -gt.to.impute, drop=FALSE], medianB[L, -gt.to.impute, drop=FALSE])
		imputedVals <- Z %*% betahat
		medianA.BB[i, L] <- imputedVals[, 1]
		medianB.BB[i, L] <- imputedVals[, 2]
		imputed[j] <- TRUE
	}
	return(list(medianA.AA=medianA.AA, medianA.AB=medianA.AB, medianA.BB=medianA.BB,
		    medianB.AA=medianB.AA, medianB.AB=medianB.AB, medianB.BB=medianB.BB,
		    imputed=imputed))
}

imputeAcrossBatch <- function(N.AA, N.AB, N.BB,
			      medianA.AA, medianA.AB, medianA.BB,
			      medianB.AA, medianB.AB, medianB.BB){
	N.missing <- (N.AA < 3) + (N.AB < 3) + (N.BB < 3)
	## find all indices in which one or more batches need to have 2 genotypes imputed
	missingAA.AB <- (N.AA < 3) & (N.AB < 3)
	missingAB.BB <- (N.AB < 3) & (N.BB < 3)
	missingAA <- (N.AA < 3) & (N.AB >= 3)
	missingBB <- (N.BB < 3) & (N.AB >= 3)
	index <- list(AA.AB=which(rowSums(missingAA.AB) > 0),
		      AB.BB=which(rowSums(missingAB.BB) > 0),
		      AA=which(rowSums(missingAA) > 0),
		      BB=which(rowSums(missingBB) > 0))
	imputeNone <- which(rowSums(N.missing == 0) > 0)
	## only works if there are batches with complete data
	index <- lapply(index, intersect, y=imputeNone)
	##indices.to.update <- rep(1:4, each=sapply(index, length))
	updated <- vector("list", 4)
	names(updated) <- c("AA.AB", "AB.BB", "AA", "BB")
	if(length(index[["AA.AB"]] > 0)){
		res <- imputeAA.AB(index[["AA.AB"]],
				   N.AA,
				   N.AB,
				   N.BB,
				   medianA.AA,
				   medianA.AB,
				   medianA.BB,
				   medianB.AA, medianB.AB, medianB.BB)
		updated$AA.AB <- res$imputed
	}
	if(length(index[["AB.BB"]] > 0)){
		res <- imputeAB.BB(index[["AB.BB"]],
				   N.AA,
				   N.AB,
				   N.BB,
				   res[["medianA.AA"]],
				   res[["medianA.AB"]],
				   res[["medianA.BB"]],
				   res[["medianB.AA"]],
				   res[["medianB.AB"]],
				   res[["medianB.BB"]])
		updated$AB.BB <- res$imputed
	}
	if(length(index[["AA"]] > 0)){
		res <- imputeAA(index[["AA"]],
				N.AA,
				N.AB,
				N.BB,
				res[["medianA.AA"]],
				res[["medianA.AB"]],
				res[["medianA.BB"]],
				res[["medianB.AA"]], res[["medianB.AB"]], res[["medianB.BB"]])
		updated$AA <- res$imputed
	}
	if(length(index[["BB"]] > 0)){
		res <- imputeBB(index[["BB"]],
				N.AA,
				N.AB,
				N.BB,
				res[["medianA.AA"]],
				res[["medianA.AB"]],
				res[["medianA.BB"]],
				res[["medianB.AA"]],
				res[["medianB.AB"]],
				res[["medianB.BB"]])
		updated$BB <- res$imputed
	}
	updated.indices <- unlist(updated)
	return(list(res, updated))
}


calculatePosteriorMean <- function(object, type=c("SNP", "NP", "X.SNP", "X.NP"), verbose=TRUE,
				   prior.prob=c(1/7, 1/7, 3/7, 1/7, 1/7),
				   CN=0:4, scale.sd=1){
	stopifnot(type %in% c("SNP", "NP", "X.SNP", "X.NP"))
	stopifnot(sum(prior.prob)==1)
	stopifnot(length(CN) == length(prior.prob))
	batch <- batch(object)
	is.snp <- isSnp(object)
	is.autosome <- chromosome(object) < 23
	is.annotated <- !is.na(chromosome(object))
	is.X <- chromosome(object) == 23
	is.lds <- is(calls(object), "ffdf") | is(calls(object), "ff_matrix")
	if(is.lds) stopifnot(isPackageLoaded("ff"))
	samplesPerBatch <- table(as.character(batch(object)))
	if(!"posteriorMean" %in% assayDataElementNames(object)){
		message("adding <posteriorMean> slot to assayData")
		pM <- matrix(NA, nrow(object), ncol(object), dimnames=list(featureNames(object), sampleNames(object)))
		tmp <- assayDataNew(alleleA=A(object),
				    alleleB=B(object),
				    call=calls(object),
				    callProbability=snpCallProbability(object),
				    posteriorMean=pM)
		assayData(object) <- tmp
	}
	## add new assay data element for posterior probabilities
	mylabel <- function(marker.type){
		switch(marker.type,
		       SNP="autosomal SNPs",
		       NP="autosomal nonpolymorphic markers",
		       X.SNP="chromosome X SNPs",
		       X.NP="chromosome X nonpolymorphic markers")
	}
	if(type=="SNP"){
		if(verbose) message(paste("...", mylabel(type)))
		marker.index <- whichMarkers("SNP",
					     is.snp,
					     is.autosome,
					     is.annotated,
					     is.X)
		marker.list <- splitIndicesByLength(marker.index, ocProbesets())
		emit <- ocLapply(seq_along(marker.list),
				 posteriorMean.snp,
				 object=object,
				 index.list=marker.list,
				 verbose=verbose,
				 prior.prob=prior.prob,
				 CN=CN,
				 scale.sd=scale.sd)
		#for(i in seq_along(marker.list)){
		#index <- marker.list[[i]]

	#}
	} else stop("type not available")
	if(length(emit)==1) emit <- emit[[1]] else stop("need to rbind elements of emit list?")
	##tmp <- do.call("rbind", emit)
	match.index <- match(rownames(emit), featureNames(object))
	S <- length(prior.prob)
	pM <- matrix(0, length(match.index), ncol(object), dimnames=list(featureNames(object)[match.index], sampleNames(object)))
	for(i in 1:S) pM <- pM + CN[i]*emit[, , i]
	posteriorMean(object)[match.index, ] <- pM
	return(object)
}

posteriorMean.snp <- function(stratum, object, index.list, CN,
			      prior.prob=c(1/7, 1/7, 3/7, 1/7, 1/7), is.lds=TRUE, verbose=TRUE,
			      scale.sd=1){
	if(length(scale.sd) == 1) rep(scale.sd,2)
	if(verbose) message("Probe stratum ", stratum, " of ", length(index.list))
	index <- index.list[[stratum]]
	test <- tryCatch(open(A(object)), error=function(e) NULL)
	if(!is.null(test)){
		open(B(object))
		open(tau2A.AA(object))
		open(tau2B.BB(object))
		open(tau2A.BB(object))
		open(tau2B.AA(object))
		open(corrAA(object))
		open(corrAB(object))
		open(corrBB(object))
		open(nuA(object))
		open(nuB(object))
		open(phiA(object))
		open(phiB(object))
	}
	a <- log2(as.matrix(A(object)[index, ]))
	b <- log2(as.matrix(B(object)[index, ]))
	NN <- Ns(object, i=index)[, , 1]
	sig2A <- as.matrix(tau2A.AA(object)[index,])
	sig2B <- as.matrix(tau2B.BB(object)[index,])
	tau2A <- as.matrix(tau2A.BB(object)[index,])
	tau2B <- as.matrix(tau2B.AA(object)[index,])
	corrAA <- as.matrix(corrAA(object)[index, ])
	corrAB <- as.matrix(corrAB(object)[index, ])
	corrBB <- as.matrix(corrBB(object)[index, ])
	nuA <- as.matrix(nuA(object)[index, ])
	phiA <- as.matrix(phiA(object)[index, ])
	nuB <- as.matrix(nuB(object)[index, ])
	phiB <- as.matrix(phiB(object)[index, ])
	if(!is.null(test)){
		close(A(object))
		close(B(object))
		close(tau2A.AA(object))
		close(tau2B.BB(object))
		close(tau2A.BB(object))
		close(tau2B.AA(object))
		close(corrAA(object))
		close(corrAB(object))
		close(corrBB(object))
		close(nuA(object))
		close(nuB(object))
		close(phiA(object))
		close(phiB(object))
	}
	S <- length(prior.prob)
	emit <- array(NA, dim=c(nrow(a), ncol(a), S))##SNPs x sample x 'truth'
	sample.index <- split(1:ncol(object), batch(object))
	sum.mymatrix <- function(...){
		x <- list(...)
		return(x[[1]] + do.call(sum, x[-1]))
	}
	numberGenotypes <- function(CT){
		stopifnot(length(CT)==1)
		copynumber <- paste("cn", CT, sep="")
		switch(copynumber,
		       cn0=1,
		       cn1=2,
		       cn2=3,
		       cn3=4,
		       cn4=4,
		       cn5=6, NULL)
	}
	##emit <- vector("list", length(sample.index))
	for(j in seq_along(sample.index)){
		cat("batch ", j, "\n")
		J <- sample.index[[j]]
		probs <- array(NA, dim=c(nrow(a), length(J), S))
		for(k in seq_along(CN)){
			CT <- CN[k]
			## 5: AAAAA, AAAAB, AAABB, AABBB, ABBBB, BBBBB
			##CN=4
			## AAAA, AAAB, AABB, ABBB, BBBB:  L = 4
			##CN=3
			##  AAA, AAB, ABB, BBB ; L = 4
			## CN=2
			##  AA, AB, BB; L=3
			## CN = 1: A, B; L=2
			## CN = 0: null; L=1
			L <- numberGenotypes(CT)
			stopifnot(!is.null(L))
			f.x.y <- vector("list", L)
			for(i in seq_along(f.x.y)){
				CA <- (0:CT)[i]
				CB <- CT-CA
				A.scale <- sqrt(tau2A*(CA==0) + sig2A*(CA > 0))
				B.scale <- sqrt(tau2B*(CB==0) + sig2B*(CB > 0))
				if(CA == 0 | CB == 0){
					A.scale <- A.scale*scale.sd[1]
					B.scale <- B.scale*scale.sd[1]
				} else { ## one or both greater than zero
					A.scale <- A.scale*scale.sd[2]
					B.scale <- B.scale*scale.sd[2]
				}
				A.scale <- as.numeric(A.scale)
				B.scale <- as.numeric(B.scale)
				A.scale2 <- A.scale^2
				B.scale2 <- B.scale^2
				if(CA == 0 & CB == 0) rho <- 0
				if(CA == 0 & CB > 0) rho <- corrBB
				if(CA > 0 & CB == 0) rho <- corrAA
				if(CA > 0 & CB > 0) rho <- corrAB
				rho <- as.numeric(rho)
				meanA <- as.numeric(log2(nuA+CA*phiA))
				meanB <- as.numeric(log2(nuB+CB*phiB))
				covs <- rho*A.scale*B.scale
				x <- a[, J]
				y <- b[, J]
				Q.x.y <- 1/(1-rho^2)*((x - meanA)^2/A.scale2 + (y - meanB)^2/B.scale2 - (2*rho*((x - meanA)*(y - meanB)))/(A.scale*B.scale))
				f.x.y[[i]] <- 1/(2*pi*A.scale*B.scale*sqrt(1-rho^2))*exp(-0.5*Q.x.y)
				class(f.x.y[[i]]) <- "mymatrix"
			}
			probs[, , k] <- do.call("sum", f.x.y)
			##if none of the states are likely (outlier), assign NA
			##		emit[, , counter] <- f.x.y
			##		counter <- counter+1
		}
		emit[, J, ] <- probs
	}
	for(i in 1:S) emit[, , i] <- emit[, , i]*prior.prob[i]
	total <- matrix(0, nrow(emit), ncol(emit))
	for(i in 1:S) total <- total+emit[, , i]
	## how to handle outliers...
	##  - use priors (posterior mean will likely be near 2).  smoothing needs to take into account the uncertainty
	##  - need uncertainty estimates for posterior means
	for(i in 1:S) emit[, , i] <- emit[, , i]/total
	dimnames(emit) <- list(featureNames(object)[index], sampleNames(object), paste("states", 1:S, sep="_"))
	## for one marker/one sample, the emission probs must sum to 1
	return(emit)
}




rscrlmmGT2 <- function(A, B, SNR, mixtureParams, cdfName, row.names=NULL,
                     col.names=NULL, probs=c(1/3, 1/3, 1/3), DF=6,
                     SNRMin=5, recallMin=10, recallRegMin=1000,
                     gender=NULL, desctrucitve=FALSE, verbose=TRUE,
                     returnParams=FALSE, badSNP=.7, snp.names){
	pkgname <- getCrlmmAnnotationName(cdfName)
	stopifnot(require(pkgname, character.only=TRUE, quietly=!verbose))
	open(SNR)
	open(A)
	open(B)
	open(mixtureParams)
	## expect objects to be ff
	keepIndex <- which( SNR[] > SNRMin)
	if(length(keepIndex)==0) stop("No arrays above quality threshold!")
	if(is.null(rownames(A))){
		loader("preprocStuff.rda", .crlmmPkgEnv, pkgname)
		gns <- getVarInEnv("gns", .crlmmPkgEnv)
		stopifnot(nrow(A) == length(gns))
		index <- seq(length=nrow(A))
	}
	if(!missing(snp.names)){
		stopifnot(!is.null(rownames(A)))
		##verify that A has only snps.  otherwise, calling function must pass rownames
		index <- match(snp.names, rownames(A))
	}
	snpBatches <- splitIndicesByLength(index, ocProbesets())
	NR <- length(unlist(snpBatches))
	if(verbose) cat("Calling", NR, "SNPs for recalibration... ")
	NC <- ncol(A)
	##
	if(verbose) message("Loading annotations.")
	obj1 <- loader("genotypeStuff.rda", .crlmmPkgEnv, pkgname)
	obj2 <- loader("mixtureStuff.rda", .crlmmPkgEnv, pkgname)
	## this is toget rid of the 'no visible binding' notes
	## variable definitions
	XIndex <- getVarInEnv("XIndex")
	autosomeIndex <- getVarInEnv("autosomeIndex")
	YIndex <- getVarInEnv("YIndex")
	SMEDIAN <- getVarInEnv("SMEDIAN")
	theKnots <- getVarInEnv("theKnots")
	regionInfo <- getVarInEnv("regionInfo")
	params <- getVarInEnv("params")
	rm(list=c(obj1, obj2), envir=.crlmmPkgEnv)
	rm(obj1, obj2)
	##
	## IF gender not provide, we predict
	## FIXME: XIndex may be greater than ocProbesets()
	if(is.null(gender)){
		if(verbose) message("Determining gender.")
		##    XMedian <- apply(log2(A[XIndex,, drop=FALSE])+log2(B[XIndex,, drop=FALSE]), 2, median)/2
		XMedian <- ocLapply(splitIndicesByNode(1:NC), predictGender, theA=A, theB=B, XIndex=XIndex, neededPkgs="crlmm")
		XMedian <- unlist(XMedian)
		if(sum(SNR[] > SNRMin)==1){
			gender <- which.min(c(abs(XMedian-8.9), abs(XMedian-9.5)))
		}else{
			gender <- kmeans(XMedian, c(min(XMedian[SNR[]>SNRMin]), max(XMedian[SNR[]>SNRMin])))[["cluster"]]
		}
	}
	##
	Indexes <- list(autosomeIndex, XIndex, YIndex)
	cIndexes <- list(keepIndex,
			 keepIndex[which(gender[keepIndex]==2)],
			 keepIndex[which(gender[keepIndex]==1)])
	if(verbose) cat("Calling", NR, "SNPs for recalibration... ")
	## call C
	fIndex <- which(gender==2)
	mIndex <- which(gender==1)
	## different here
	## use gtypeCallerR in batches
	##snpBatches <- splitIndicesByLength(1:nrow(A), ocProbesets())
	newparamsBatch <- vector("list", length(snpBatches))
	process1 <- function(idxBatch, snpBatches, autosomeIndex, XIndex,
			     YIndex, A, B, mixtureParams, fIndex, mIndex,
			     params, cIndexes, SMEDIAN, theKnots, DF, probs, batchSize){
		open(A)
		open(B)
		open(mixtureParams)
		snps <- snpBatches[[idxBatch]]
		rSnps <- range(snps)
		last <- (idxBatch-1)*batchSize
		IndexesBatch <- list(autosomeIndex[autosomeIndex %in% snps]-last,
				     XIndex[XIndex %in% snps]-last,
				     YIndex[YIndex %in% snps]-last)
		IndexesBatch <- lapply(IndexesBatch, as.integer)
		tmpA <- as.matrix(A[snps,])
		tmpB <- as.matrix(B[snps,])
		## newparamsBatch[[idxBatch]]
		tmp <- gtypeCallerR(tmpA, tmpB, fIndex, mIndex,
				    params[["centers"]][snps,],
				    params[["scales"]][snps,],
				    params[["N"]][snps,],
				    IndexesBatch, cIndexes,
				    sapply(IndexesBatch, length),
				    sapply(cIndexes, length), SMEDIAN,
				    theKnots, mixtureParams[], DF, probs, 0.025)
		rm(snps, rSnps, IndexesBatch, tmpA, tmpB, last)
		gc(verbose=FALSE)
		close(A)
		close(B)
		close(mixtureParams)
		tmp
	}
	##
	newparamsBatch <- ocLapply(seq(along=snpBatches), process1,
				   snpBatches=snpBatches,
				   autosomeIndex=autosomeIndex, XIndex=XIndex,
				   YIndex=YIndex, A=A, B=B,
				   mixtureParams=mixtureParams, fIndex=fIndex,
				   mIndex=mIndex, params=params,
				   cIndexes=cIndexes, SMEDIAN=SMEDIAN,
				   theKnots=theKnots, DF=DF, probs=probs,
				   batchSize=ocProbesets())
	newparams <- vector("list", 3)
	names(newparams) <- c("centers", "scales", "N")
	newparams[["centers"]] <- do.call("rbind", lapply(newparamsBatch, "[[", 1))
	newparams[["scales"]] <- do.call("rbind", lapply(newparamsBatch, "[[", 2))
	newparams[["N"]] <- do.call("rbind", lapply(newparamsBatch, "[[", 3))
	rm(newparamsBatch); gc(verbose=FALSE)
	if(verbose) message("Done.")
	if(verbose) message("Estimating recalibration parameters.")
	d <- newparams[["centers"]] - params$centers
	##
	##regression
	Index <- intersect(which(pmin(newparams[["N"]][, 1],
				      newparams[["N"]][, 2],
				      newparams[["N"]][, 3]) > recallMin &
				 !apply(regionInfo, 1, any)),
			   autosomeIndex)
	if(length(Index) < recallRegMin){
		warning("Recalibration not possible. Possible cause: small sample size.")
		newparams <- params
		dev <- vector("numeric", nrow(newparams[["centers"]]))
		SS <- matrix(Inf, 3, 3)
		DD <- 0
	}else{
		data4reg <- as.data.frame(newparams[["centers"]][Index,])
		names(data4reg) <- c("AA", "AB", "BB")
		regParams <- cbind(  coef(lm(AA~AB*BB, data=data4reg)),
				   c(coef(lm(AB~AA+BB, data=data4reg)), 0),
				   coef(lm(BB~AA*AB, data=data4reg)))
		rownames(regParams) <- c("intercept", "X", "Y", "XY")
		rm(data4reg)
		##
		minN <- 3
		newparams[["centers"]][newparams[["N"]] < minN] <- NA
		Index <- setdiff(which(rowSums(is.na(newparams[["centers"]]))==1), YIndex)
		if(verbose) message("Filling out empty centers", appendLF=FALSE)
		for(i in Index){
			if(verbose) if(i%%10000==0) message(".", appendLF=FALSE)
			mu <- newparams[["centers"]][i, ]
			j <- which(is.na(mu))
			newparams[["centers"]][i, j] <- c(1, mu[-j], prod(mu[-j]))%*%regParams[, j]
			rm(mu, j)
		}
		##
		##remaing NAs are made like originals
		if(length(YIndex)>0){
			noMoveIndex <- union(setdiff(which(rowSums(is.na(newparams[["centers"]]))>0), YIndex),
					     YIndex[rowSums(is.na(newparams[["centers"]][YIndex, ])>1)])
		}
		snps2ignore <- which(rowSums(is.na(newparams[["centers"]])) > 0)
		snps2keep <- setdiff(autosomeIndex, snps2ignore)
		rm(snps2ignore)
		newparams[["centers"]][is.na(newparams[["centers"]])] <- params[["centers"]][is.na(newparams[["centers"]])]
		if(verbose) cat("\n")
		##
		if(verbose) message("Calculating and standardizing size of shift... ", appendLF=FALSE)
		GG <- DD <- newparams[["centers"]] - params[["centers"]]
		DD <- sweep(DD, 2, colMeans(DD[autosomeIndex, ]))
		SS <- cov(DD[autosomeIndex, ])
		SSI <- solve(SS)
		dev <- vector("numeric", nrow(DD))
		if(length(YIndex)){
			dev[-YIndex] <- apply(DD[-YIndex, ], 1, function(x) x%*%SSI%*%x)
			dev[-YIndex] <- 1/sqrt( (2*pi)^3*det(SS))*exp(-0.5*dev[-YIndex])
			##Now Y (only two params)
			SSY <- SS[c(1, 3), c(1, 3)]
			SSI <- solve(SSY)
			dev[YIndex] <- apply(DD[YIndex, c(1, 3)], 1, function(x) x%*%SSI%*%x)
			dev[YIndex] <- 1/sqrt( (2*pi)^2*det(SSY))*exp(-0.5*dev[YIndex])
		} else {
			dev=apply(DD,1,function(x) x%*%SSI%*%x)
			dev=1/sqrt( (2*pi)^3*det(SS))*exp(-0.5*dev)
		}
	}
	if (verbose) message("OK")
	##
	## BC: must keep SD
	params[-2] <- newparams[-2]
	rm(newparams)
	gc(verbose=FALSE)
	##
	if(verbose) message("Calling ", NR, " SNPs... ", appendLF=FALSE)
	##
	## ###################
	## ## MOVE TO C#######
	##
	## running in batches
	process2 <- function(idxBatch, snpBatches, autosomeIndex, XIndex,
			     YIndex, A, B, mixtureParams, fIndex, mIndex,
			     params, cIndexes, SMEDIAN, theKnots, DF, probs,
			     regionInfo, batchSize){
		open(A)
		open(B)
		open(mixtureParams)
		snps <- snpBatches[[idxBatch]]
		tmpA <- as.matrix(A[snps,])
		tmpB <- as.matrix(B[snps,])
		rSnps <- range(snps)
		last <- (idxBatch-1)*batchSize
		IndexesBatch <- list(autosomeIndex[autosomeIndex %in% snps]-last,
				     XIndex[XIndex %in% snps]-last,
				     YIndex[YIndex %in% snps]-last)
		IndexesBatch <- lapply(IndexesBatch, as.integer)
		ImNull <- gtypeCallerR2(tmpA, tmpB, fIndex, mIndex,
					params[["centers"]][snps,],
					params[["scales"]][snps,],
					params[["N"]][snps,],
					IndexesBatch, cIndexes,
					sapply(IndexesBatch, length),
					sapply(cIndexes, length),
					SMEDIAN, theKnots, mixtureParams[],
					DF, probs, 0.025,
					which(regionInfo[snps, 2]),
					which(regionInfo[snps, 1]))
		A[snps,] <- tmpA
		B[snps,] <- tmpB
		rm(tmpA, tmpB, snps, rSnps, IndexesBatch, ImNull, last)
		gc(verbose=FALSE)
		close(A)
		close(B)
		close(mixtureParams)
	}
	##
	ocLapply(seq(along=snpBatches), process2, snpBatches=snpBatches,
		 autosomeIndex=autosomeIndex, XIndex=XIndex, YIndex=YIndex,
		 A=A, B=B, mixtureParams=mixtureParams, fIndex=fIndex,
		 mIndex=mIndex, params=params, cIndexes=cIndexes,
		 SMEDIAN=SMEDIAN, theKnots=theKnots, DF=DF, probs=probs,
		 regionInfo=regionInfo, batchSize=ocProbesets())
	##  END MOVE TO C#######
	## ##################
	##
	dev <- dev/(dev+1/383)
	if(!is.null(row.names)){ rownames(A) <- rownames(B) <- names(dev) <- row.names}
	if(!is.null(col.names)){ colnames(A) <- colnames(B) <- col.names}
	##
	if(length(Index) >= recallRegMin){
		tmp4batchQC <- DD[autosomeIndex,]/(params[["N"]][autosomeIndex,]+1)
		tmpSnpQc <- dev[autosomeIndex]
		SS <- cov(tmp4batchQC[tmpSnpQc < badSNP,])
		batchQC <- mean(diag(SS))
	}else{
		SS <- matrix(0, 3, 3)
		batchQC <- Inf
	}
	##
	if(verbose) message("Done.")
	if (returnParams){
		return(list(calls=A, confs=B, SNPQC=dev, batchQC=batchQC, params=params, DD=DD, covDD=SS, gender=gender, pkgname=pkgname))
	}else{
		return(list(calls=A, confs=B, SNPQC=dev, batchQC=batchQC, DD=DD, covDD=SS, gender=gender, pkgname=pkgname))
	}
}

crlmm2.2 <- function(filenames, row.names=TRUE, col.names=TRUE,
                   probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5, gender=NULL,
                   save.it=FALSE, load.it=FALSE, intensityFile,
                   mixtureSampleSize=10^5, eps=0.1, verbose=TRUE,
                   cdfName, sns, recallMin=10, recallRegMin=1000,
                   returnParams=FALSE, badSNP=.7){
  if ((load.it || save.it) && missing(intensityFile))
    stop("'intensityFile' is missing, and you chose either load.it or save.it")
  if (missing(sns)) sns <- basename(filenames)
  if (!missing(intensityFile))
    if (load.it & !file.exists(intensityFile)){
      load.it <- FALSE
      message("File ", intensityFile, " does not exist.")
      message("Not loading it, but running SNPRMA from scratch.")
    }
  if (!load.it){
    res <- snprma2(filenames, fitMixture=TRUE,
                   mixtureSampleSize=mixtureSampleSize, verbose=verbose,
                   eps=eps, cdfName=cdfName, sns=sns)
    open(res[["A"]])
    open(res[["B"]])
    open(res[["SNR"]])
    open(res[["mixtureParams"]])
    if(save.it){
      t0 <- proc.time()
      save(res, file=intensityFile)
      t0 <- proc.time()-t0
      if (verbose) message("Used ", t0[3], " seconds to save ", intensityFile, ".")
    }
  }else{
    if (verbose) message("Loading ", intensityFile, ".")
    obj <- load(intensityFile)
    if (verbose) message("Done.")
    if (obj != "res")
      stop("Object in ", intensityFile, " seems to be invalid.")
  }
  if(row.names) row.names=res$gns else row.names=NULL
  if(col.names) col.names=res$sns else col.names=NULL
  res2 <- rscrlmmGT2(res[["A"]], res[["B"]], res[["SNR"]],
                   res[["mixtureParams"]], res[["cdfName"]],
                   gender=gender, row.names=row.names,
                   col.names=col.names, recallMin=recallMin,
                   recallRegMin=1000, SNRMin=SNRMin,
                   returnParams=returnParams, badSNP=badSNP,
                   verbose=verbose)

  res2[["SNR"]] <- res[["SNR"]]
  res2[["SKW"]] <- res[["SKW"]]
  return(list2SnpSet(res2, returnParams=returnParams))
}
