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

construct <- function(filenames, cdfName, copynumber=FALSE,
		      sns, verbose=TRUE, batch, fns){
	if(!missing(batch)){
		stopifnot(length(batch) == length(sns))
	}
	if(missing(sns) & missing(filenames)) stop("one of filenames or samplenames (sns) must be provided")
	if(verbose) message("Initializing container for assay data elements alleleA, alleleB, call, callProbability")
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

##genotype <- function(filenames,
##		     cdfName,
##		     batch,
##		     mixtureSampleSize=10^5,
##		     eps=0.1,
##		     verbose=TRUE,
##		     seed=1,
##		     sns,
##		     copynumber=TRUE,
##		     probs=rep(1/3, 3),
##		     DF=6,
##		     SNRMin=5,
##		     recallMin=10,
##		     recallRegMin=1000,
##		     gender=NULL,
##		     returnParams=TRUE,
##		     badSNP=0.7){
##	if(missing(cdfName)) stop("must specify cdfName")
##	if(!isValidCdfName(cdfName)) stop("cdfName not valid.  see validCdfNames")
##	if(missing(sns)) sns <- basename(filenames)
##	## callSet contains potentially very big matrices
##	## More big matrices are created within snprma, that will then be removed.
##	callSet <- construct(filenames=filenames,
##			     cdfName=cdfName,
##			     copynumber=copynumber,
##			     sns=sns,
##			     verbose=verbose,
##			     batch=batch)
##	##lM(callSet) <- initializeParamObject(list(featureNames(callSet), unique(protocolData(callSet)$batch)))
##	mixtureParams <- matrix(NA, 4, length(filenames))
##	snp.index <- which(isSnp(callSet)==1)
##	sample.strata <- splitIndicesByLength(1:ncol(callSet), ocSamples())
##	iter <- 1
####	for(j in batches){
##	j <- unlist(sample.strata)
##		if(verbose) message("Batch ", iter, " of ", length(sample.strata))
##		snprmaRes <- snprma(filenames=filenames[j],
##				    mixtureSampleSize=mixtureSampleSize,
##				    fitMixture=TRUE,
##				    eps=eps,
##				    verbose=verbose,
##				    seed=seed,
##				    cdfName=cdfName,
##				    sns=sns[j])
##		stopifnot(identical(featureNames(callSet)[snp.index], snprmaRes$gns))
##		pData(callSet)$SKW[j] <- snprmaRes$SKW
##		pData(callSet)$SNR[j] <- snprmaRes$SNR
##		suppressWarnings(A(callSet)[snp.index, j] <- snprmaRes[["A"]])
##		suppressWarnings(B(callSet)[snp.index, j] <- snprmaRes[["B"]])
##		mixtureParams[, j] <- snprmaRes$mixtureParams
##		rm(snprmaRes); ##gc()
##		if(copynumber){
##			np.index <- which(isSnp(callSet) == 0)
##			cnrmaRes <- cnrma(filenames=filenames[j],
##					  cdfName=cdfName,
##					  row.names=featureNames(callSet)[np.index],
##					  sns=sns[j],
##					  seed=seed,
##					  verbose=verbose)
##			stopifnot(identical(featureNames(callSet)[np.index], rownames(cnrmaRes)))
##			A(callSet)[np.index, j] <- cnrmaRes
##			rm(cnrmaRes); ##gc()
##		}
##		## as.matrix needed when ffdf is used
##		tmp <- crlmmGT(A=as.matrix(A(callSet)[snp.index, j]),
##			       B=as.matrix(B(callSet)[snp.index, j]),
##			       SNR=callSet$SNR[j],
##			       mixtureParams=mixtureParams[, j],
##			       cdfName=annotation(callSet),
##			       row.names=featureNames(callSet)[snp.index],
##			       col.names=sampleNames(callSet)[j],
##			       probs=probs,
##			       DF=DF,
##			       SNRMin=SNRMin,
##			       recallMin=recallMin,
##			       recallRegMin=recallRegMin,
##			       gender=gender,
##			       verbose=verbose,
##			       returnParams=returnParams,
##			       badSNP=badSNP)
##		snpCall(callSet)[snp.index, j] <- tmp[["calls"]]
##		snpCallProbability(callSet)[snp.index, j] <- tmp[["confs"]]
##		callSet$gender[j] <- tmp$gender
##		iter <- iter+1
####	}
##	return(callSet)
##}



##genotypeLargeData
genotype <- function(filenames,
		       cdfName,
		       batch,
		       mixtureSampleSize=10^5,
		       eps=0.1,
		       verbose=TRUE,
		       seed=1,
		       sns,
		       copynumber=TRUE,
		       probs=rep(1/3, 3),
		       DF=6,
		       SNRMin=5,
		       recallMin=10,
		       recallRegMin=1000,
		       gender=NULL,
		       returnParams=TRUE,
		       badSNP=0.7){
	is.lds <- ifelse(isPackageLoaded("ff"), TRUE, FALSE)
	if(!copynumber){
		FUN <- ifelse(is.lds, "crlmm2", "crlmm")
		callSet <- FUN(filenames=filenames,
			       cdfName=cdfName,
			       mixtureSampleSize=mixtureSampleSize,
			       eps=eps,
			       verbose=verbose,
			       sns=sns,
			       probs=probs,
			       DF=DF,
			       SNRMin=SNRMin,
			       recallMin=recallMin,
			       recallRegMin=recallRegMin,
			       gender=gender,
			       returnParams=returnParams,
			       badSNP=badSNP)
		return(callSet)
	}
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
##	stopifnot(identical(featureNames(callSet)[snp.index], snprmaRes$gns))
	pData(callSet)$SKW <- snprmaRes$SKW
	pData(callSet)$SNR <- snprmaRes$SNR
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
	} else {
		calls(callSet)[snp.index, ] <- tmp[["calls"]]
		snpCallProbability(callSet)[snp.index, ] <- tmp[["confs"]]
	}
	callSet$gender <- tmp$gender
	close(callSet)
	return(callSet)
}
genotypeLD <- genotype2 <- genotype

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

applyByGenotype <- function(x, FUN, G){
	FUN <- match.fun(FUN)
	tmp <- matrix(NA, nrow(x), 3)
	for(j in 1:3){
		GT <- G==j
		GT[GT == FALSE] <- NA
		gt.x <- GT*x
		tmp[, j] <- FUN(gt.x, na.rm=TRUE)
	}
	tmp
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

corByGenotype <- function(A, B, G, Ns, which.cluster=c(1,2,3)[1]){##, DF.PRIOR){
	x <- A * (G == which.cluster)
	x[x==0] <- NA
	y <- B * (G == which.cluster)
	res <- as.matrix(rowCors(x, y, na.rm=TRUE))
##	cors <- shrink(res, Ns[, which.cluster], DF.PRIOR)
	res
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


##fit.wls <- function(allele, Ystar, W, Ns, autosome=TRUE){
fit.wls <- function(NN, sigma, allele, Y, autosome, X){
	##		Np <- NN
##		Np[Np < 1] <- 1
##		vA2 <- vA^2/Np
##		vB2 <- vB^2/Np
##		wA <- sqrt(1/vA2)
##		wB <- sqrt(1/vB2)
##		YA <- muA*wA
##		YB <- muB*wB
	Np <- NN
	Np[Np < 1] <- 1
	W <- (sigma/sqrt(Np))^-1
	Ystar <- Y*W
	complete <- which(rowSums(is.na(W)) == 0 & rowSums(is.na(Ystar)) == 0)
##	if(any(!is.finite(W))){## | any(!is.finite(V))){
##		i <- which(rowSums(!is.finite(W)) > 0)
##		stop("Possible zeros in the within-genotype estimates of the spread (vA, vB). ")
##	}
##	NOHET <- mean(NN[, 2], na.rm=TRUE) < 0.05
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
##	if(NOHET) X <- X[-2, ] ##more than 1 X chromosome, but all homozygous
	##How to quickly generate Xstar, Xstar = diag(W) %*% X
	##Xstar <- apply(W, 1, generateX, X)
	ww <- rep(1, ncol(Ystar))
	for(i in complete){
		betahat[, i] <- dqrlsWrapper(W[i, ] * X, Ystar[i, ], ww)
		##ssr <- sum((Ystar[i, ] - matrix(Xstar[, i], nrow(X), ncol(X)) %*% matrix(betahat[, i], ncol(X), 1))^2)
	}
	return(betahat)
}

##nuphiAlleleX <- function(allele, Ystar, W, Ns, chrX=FALSE){
##	complete <- rowSums(is.na(W)) == 0
##	if(any(!is.finite(W))){## | any(!is.finite(V))){
##		i <- which(rowSums(!is.finite(W)) > 0)
##		stop("Possible zeros in the within-genotype estimates of the spread (vA, vB). ")
##	}
##	##NOHET <- mean(Ns[, 2], na.rm=TRUE) < 0.05
##	if(missing(allele)) stop("must specify allele")
##	if(allele == "A") X <- cbind(1, c(1, 0, 2, 1, 0), c(0, 1, 0, 1, 2))
##	if(allele == "B") X <- cbind(1, c(0, 1, 0, 1, 2), c(1, 0, 2, 1, 0))
##	##if(NOHET) X <- X[-2, ] ##more than 1 X chromosome, but all homozygous
##	##How to quickly generate Xstar, Xstar = diag(W) %*% X
##	Xstar <- apply(W, 1, generateX, X)
##	IXTX <- apply(Xstar, 2, generateIXTX, nrow=nrow(X))
##
##	betahat <- matrix(NA, 3, nrow(Ystar))
##	ses <- matrix(NA, 3, nrow(Ystar))
##	##betahat <- matrix(NA, 2, nrow(Ystar))
##	##ses <- matrix(NA, 2, nrow(Ystar))
##	for(i in 1:nrow(Ystar)){
##		betahat[, i] <- crossprod(matrix(IXTX[, i], ncol(X), ncol(X)), crossprod(matrix(Xstar[, i], nrow=nrow(X)), Ystar[i, ]))
##		ssr <- sum((Ystar[i, ] - matrix(Xstar[, i], nrow(X), ncol(X)) %*% matrix(betahat[, i], ncol(X), 1))^2)
##		ses[, i] <- sqrt(diag(matrix(IXTX[, i], ncol(X), ncol(X)) * ssr))
##	}
##	nu <- betahat[1, ]
##	phi <- betahat[2, ]
##	phi2 <- betahat[3, ]
##	return(list(nu, phi, phi2))
##}


##nuphiAllele <- function(object, allele, Ystar, W, tmp.objects, cnOptions){
##	I <- isSnp(object)
##	Ystar <- Ystar[I, ]
##	rownames(Ystar) <- featureNames(object)[isSnp(object)]
##	complete <- rowSums(is.na(W)) == 0 & I
##	W <- W[I, ]
##	if(any(!is.finite(W))){## | any(!is.finite(V))){
##		i <- which(rowSums(!is.finite(W)) > 0)
##		stop("Possible zeros in the within-genotype estimates of the spread (vA, vB). ")
##	}
##	Ns <- tmp.objects[["Ns"]]
##	Ns <- Ns[I, ]
##	CHR <- unique(chromosome(object))
##	batch <- unique(batch(object))
##	nuA <- getParam(object, "nuA", batch)
##	nuB <- getParam(object, "nuB", batch)
##	phiA <- getParam(object, "phiA", batch)
##	phiB <- getParam(object, "phiB", batch)
##	if(CHR==23){
##		phiAX <- getParam(object, "phiAX", batch)
##		phiBX <- getParam(object, "phiBX", batch)
##	}
##	NOHET <- mean(Ns[, "AB"], na.rm=TRUE) < 0.05
##	if(missing(allele)) stop("must specify allele")
##	if(CHR == 23){
##		if(length(grep("AB", colnames(W))) > 0){
##			if(allele == "A") X <- cbind(1, c(1, 0, 2, 1, 0), c(0, 1, 0, 1, 2))
##			if(allele == "B") X <- cbind(1, c(0, 1, 0, 1, 2), c(1, 0, 2, 1, 0))
##		} else{
##		if(allele == "A") X <- cbind(1, c(1, 0, 2, 0), c(0, 1, 0, 2))
##			if(allele == "B") X <- cbind(1, c(0, 1, 0, 2), c(1, 0, 2, 0))
##		}
##		betahat <- matrix(NA, 3, nrow(Ystar))
##	} else {##autosome
##		if(allele == "A") X <- cbind(1, 2:0) else X <- cbind(1, 0:2)
##		if(NOHET) X <- X[-2, ] ##more than 1 X chromosome, but all homozygous
##		betahat <- matrix(NA, 2, nrow(Ystar))
##	}
##	ww <- rep(1, ncol(Ystar))
##	II <- which(rowSums(is.nan(Ystar)) == 0)
##	for(i in II){
##		betahat[, i] <- dqrlsWrapper(W[i, ] * X, Ystar[i, ], ww)
##	}
##	if(allele == "A"){
##		nuA[complete] <- betahat[1, ]
##		phiA[complete] <- betahat[2, ]
##		object <- pr(object, "nuA", batch, nuA)
##		object <- pr(object, "phiA", batch, phiA)
##		if(CHR == 23){
##			phiAX[complete] <- betahat[3, ]
##			object <- pr(object, "phiAX", batch, phiAX)
##		}
##	}
##	if(allele=="B"){
##		nuB[complete] <- betahat[1, ]
##		phiB[complete] <- betahat[2, ]
##		object <- pr(object, "nuB", batch, nuB)
##		object <- pr(object, "phiB", batch, phiB)
##		if(CHR == 23){
##			phiBX[complete] <- betahat[3, ]
##			object <- pr(object, "phiBX", batch, phiBX)
##		}
##	}
##	return(object)
##}



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


##replace with release/crlmm/R/cnrma-functions
##crlmmCopynumber <- function(object,
##			    chromosome=1:23,
##			    which.batches,
##			    MIN.SAMPLES=10,
##			    SNRMin=5,
##			    MIN.OBS=3,
##			    DF.PRIOR=50,
##			    bias.adj=FALSE,
##			    prior.prob=rep(1/4,4),
##			    seed=1,
##			    verbose=TRUE,
##			    GT.CONF.THR=0.99,
##			    PHI.THR=2^6,
##			    nHOM.THR=5,
##			    MIN.NU=2^3,
##			    MIN.PHI=2^3,
##			    THR.NU.PHI=TRUE,
##			    thresholdCopynumber=TRUE,
##			    weighted.lm=TRUE){
##	ffIsLoaded <- class(calls(object))[[1]] == "ff"
##	if(ffIsLoaded){
##		open(object)
##		open(object$SKW)
##		open(object$SNR)
##	}
##	stopifnot("batch" %in% varLabels(protocolData(object)))
##	stopifnot("chromosome" %in% fvarLabels(object))
##	stopifnot("position" %in% fvarLabels(object))
##	stopifnot("isSnp" %in% fvarLabels(object))
##	##batch <- object$batch
##	batch <- batch(object)
##	if(ffIsLoaded) {
##		open(object$SNR)
##		SNR <- object$SNR[, ]
##	} else SNR <-  object$SNR
##	batches <- split((1:ncol(object))[SNR > SNRMin], batch[SNR > SNRMin])
##	if(any(sapply(batches, length) < MIN.SAMPLES)) message("Excluding batches with fewer than ", MIN.SAMPLES, " samples")
##	batches <- batches[sapply(batches, length) >= MIN.SAMPLES]
##	if(missing(which.batches)) which.batches <- seq(along=batches)
##	for(i in chromosome){
##		if(verbose) cat("Chromosome ", i, "\n")
##		if(i >= 24) next()
##		ii <- which.batches[1]
##		for(j in batches[which.batches]){
##			if(verbose) message("Batch ", ii, " of ", length(which.batches))
##			row.index <- which(chromosome(object) == i)
##			##Note that ffdf assayDataElements are data.frames after subsetting(not matrices)
####			ca <- as.matrix(CA(object)[row.index, j])
####			cb <- as.matrix(CB(object)[row.index, j])
####			dimnames(ca) <- dimnames(cb) <- list(featureNames(object)[row.index], sampleNames(object)[j])
##			tmp <- new("CNSet",
##				   call=as.matrix(calls(object)[row.index, j]),
##				   callProbability=as.matrix(snpCallProbability(object)[row.index, j]),
##				   alleleA=as.matrix(A(object)[row.index, j]),
##				   alleleB=as.matrix(B(object)[row.index, j]),
####				   CA=ca, CB=cb,
##				   phenoData=phenoData(object)[j, ],
##				   annotation=annotation(object))
##			featureData(tmp) <- addFeatureAnnotation(tmp)
##			featureData(tmp) <- lm.parameters(tmp, batch=unique(batch[j]))
##			tmp$batch <- batch(object)[j]
##			tmp <- computeCopynumber(tmp,
##						 MIN.OBS=MIN.OBS,
##						 DF.PRIOR=DF.PRIOR,
##						 bias.adj=bias.adj,
##						 prior.prob=prior.prob,
##						 seed=seed,
##						 verbose=verbose,
##						 GT.CONF.THR=GT.CONF.THR,
##						 PHI.THR=PHI.THR,
##						 nHOM.THR=nHOM.THR,
##						 MIN.NU=MIN.NU,
##						 MIN.PHI=MIN.PHI,
##						 THR.NU.PHI=THR.NU.PHI,
##						 thresholdCopynumber=thresholdCopynumber)
##			fData(tmp) <- fData(tmp)[, -(1:3)]
####			CA(tmp) <- matrix(as.integer(CA(tmp)*100), nrow=nrow(tmp), ncol=ncol(tmp),
####					  dimnames=list(featureNames(tmp), sampleNames(tmp)))
####			CB(tmp) <- matrix(as.integer(CB(tmp)*100), nrow=nrow(tmp), ncol=ncol(tmp),
####					  dimnames=list(featureNames(tmp), sampleNames(tmp)))
####			CA(object)[row.index, j] <- CA(tmp)
####			CB(object)[row.index, j] <- CB(tmp)
##			##ad-hocery
##			batchName <- unique(batch(object)[j])
##			fvarLabels(tmp)[15:17] <- paste(c("corrAB", "corrBB", "corrAA"), batchName, sep=".")
##			fvarLabels(tmp)[13:14] <- paste(c("phiPrimeA", "phiPrimeB"), batchName, sep=".")
##			fvarLabels(tmp) <- gsub("_", ".", fvarLabels(tmp))
##			##fvarLabels(tmp) <- gsub("\\.[1-9]", "", fvarLabels(tmp))
##			if(ffIsLoaded){
##				physical <- get("physical")
##				fData(tmp) <- fData(tmp)[, fvarLabels(tmp) %in% names(physical(lM(object)))]
##				jj <- match(fvarLabels(tmp), names(lM(object)))
##				lM(object)[row.index, jj] <- fData(tmp)
##			} else {
##				nms <- paste(names(lM(object)), batchName, sep=".")
##				fData(tmp) <- fData(tmp)[, fvarLabels(tmp) %in% nms]
##				for(k in seq_along(fvarLabels(tmp))){
##					kk <- match(fvarLabels(tmp)[k], paste(names(lM(object)), batchName, sep="."))
##					column <- match(batchName, colnames(lM(object)[[k]]))
##					lM(object)[[k]][row.index, column] <- fData(tmp)[, k]
##				}
##			}
##			rm(tmp); ##gc()
##			ii <- ii+1
##		}
##	}
##	return(object)
##}

crlmmCopynumberLD <- function(object,
			     which.batches,
			    MIN.SAMPLES=10,
			    SNRMin=5,
			    MIN.OBS=1,
			    DF.PRIOR=50,
			    bias.adj=FALSE,
			    prior.prob=rep(1/4,4),
			    seed=1,
			    verbose=TRUE,
			    GT.CONF.THR=0.99,
			    PHI.THR=2^6,
			    nHOM.THR=5,
			    MIN.NU=2^3,
			    MIN.PHI=2^3,
			    THR.NU.PHI=TRUE,
			    thresholdCopynumber=TRUE){
	stopifnot("batch" %in% varLabels(protocolData(object)))
	stopifnot("chromosome" %in% fvarLabels(object))
	stopifnot("position" %in% fvarLabels(object))
	stopifnot("isSnp" %in% fvarLabels(object))
	##lM(object) <- initializeParamObject(list(featureNames(object), unique(protocolData(object)$batch)))
	batch <- batch(object)
	XIndex.snps <- which(chromosome(object) == 23 & isSnp(object) & !is.na(chromosome(object)))
	##YIndex.snps <- (1:nrow(object))[chromosome(object) == 24 & isSnp(object)]
	XIndex.nps <- which(chromosome(object) == 23 & !isSnp(object) & !is.na(chromosome(object)))
	##autosomeIndex.snps <- (1:nrow(object))[chromosome(object) < 23 & isSnp(object) & !is.na(chromosome(object))]
	autosomeIndex.snps <- which(chromosome(object) < 23 & isSnp(object) & !is.na(chromosome(object)))
	autosomeIndex.nps <- which(chromosome(object) < 23 & !isSnp(object) & !is.na(chromosome(object)))
	##autosomeIndex.nps <- (1:nrow(object))[chromosome(object) < 23 & !isSnp(object) & !is.na(chromosome(object))]

	## Do chromosome X in batches
	Ns <- initializeBigMatrix("Ns", nrow(object), 5)
	colnames(Ns) <- c("A", "B", "AA", "AB", "BB")
	if(!file.exists(file.path(ldPath(), "normal.rda"))){
		normal <- initializeBigMatrix("normal", nrow(object), ncol(object), vmode="integer")
		normal[,] <- 1L
		save(normal, file=file.path(ldPath(), "normal.rda"))
	} else load(file.path(ldPath(), "normal.rda"))
	if(!file.exists(file.path(ldPath(), "snpflags.rda"))){
		snpflags <- initializeBigMatrix("snpflags", nrow(object), length(unique(batch(object))), vmode="integer")
		snpflags[,] <- 0L
		save(snpflags, file=file.path(ldPath(), "snpflags.rda"))
	} else{
		load(file.path(ldPath(), "snpflags.rda"))
	}
	if(verbose) message("Estimating allele-specific copy number at autosomal SNPs")
	snpBatches <- splitIndicesByLength(autosomeIndex.snps, ocProbesets())
	ocLapply(seq(along=snpBatches),
		 fit.lm1,
		 marker.index=autosomeIndex.snps,
		 object=object,
		 Ns=Ns,
		 normal=normal,
		 snpflags=snpflags,
		 snpBatches=snpBatches,
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
		 neededPkgs="crlmm")
	## autosomal NPs
	snpBatches <- splitIndicesByLength(autosomeIndex.nps, ocProbesets())
	if(verbose) message("Estimating total copy number at nonpolymorphic loci")
	ocLapply(seq(along=snpBatches),
			fit.lm2,
			marker.index=autosomeIndex.nps,
			object=object,
			Ns=Ns,
			normal=normal,
			snpflags=snpflags,
			snpBatches=snpBatches,
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
		 neededPkgs="crlmm")
	snpBatches <- splitIndicesByLength(XIndex.snps, ocProbesets())
	if(verbose) message("Estimating allele-specific copy number at polymorphic loci on chromosome X")
	ocLapply(seq(along=snpBatches),
		 fit.lm3,
		 marker.index=XIndex.snps,
		 object=object,
		 Ns=Ns,
		 normal=normal,
		 snpflags=snpflags,
		 snpBatches=snpBatches,
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
		 neededPkgs="crlmm")
	if(verbose) message("Estimating total copy number for nonpolymorphic loci on chromosome X")
	snpBatches <- splitIndicesByLength(XIndex.nps, ocProbesets())
	ocLapply(seq(along=snpBatches),
		 fit.lm4,
		 marker.index=XIndex.nps,
		 object=object,
		 Ns=Ns,
		 normal=normal,
		 snpflags=snpflags,
		 snpBatches=snpBatches,
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
		 neededPkgs="crlmm")
	return(object)
}
crlmmCopynumber2 <- crlmmCopynumberLD




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

		##---------------------------------------------------------------------------
		## SNPs that we'll use for imputing location/scale of unobserved genotypes
		##---------------------------------------------------------------------------
		index.complete <- indexComplete(NN, medianA[[k]], medianB[[k]], MIN.OBS)

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

	madA.AA(object)[marker.index, ] <- do.call("cbind", lapply(shrink.madA, function(x) x[, 1]))
	madA.AB(object)[marker.index, ] <- do.call("cbind", lapply(shrink.madA, function(x) x[, 2]))
	madA.BB(object)[marker.index, ] <- do.call("cbind", lapply(shrink.madA, function(x) x[, 3]))
	madB.AA(object)[marker.index, ] <- do.call("cbind", lapply(shrink.madB, function(x) x[, 1]))
	madB.AB(object)[marker.index, ] <- do.call("cbind", lapply(shrink.madB, function(x) x[, 2]))
	madB.BB(object)[marker.index, ] <- do.call("cbind", lapply(shrink.madB, function(x) x[, 3]))

	corrAA(object)[marker.index, ] <- shrink.corrAA
	corrAB(object)[marker.index, ] <- shrink.corrAB
	corrBB(object)[marker.index, ] <- shrink.corrBB
	tau2A.AA(object)[marker.index,] <- shrink.tau2A.AA
	tau2A.BB(object)[marker.index,] <- shrink.tau2A.BB
	tau2B.AA(object)[marker.index,] <- shrink.tau2B.AA
	tau2B.BB(object)[marker.index,] <- shrink.tau2B.BB
	if(is.lds) return(TRUE) else retrun(object)
}



fit.lm1 <- function(strata,
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
		##res <- fit.wls(allele="A", Ystar=YA, W=wA, Ns=Ns)
		##nuA[, J] <- res[[1]]
		##phiA[, J] <- res[[2]]
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
##	cA[cA < 0.05] <- 0.05
##	cB[cB < 0.05] <- 0.05
##	cA[cA > 5] <-  5
##	cB[cB > 5] <- 5
##	cA <- matrix(as.integer(cA*100), nrow(cA), ncol(cA))
##	cB <- matrix(as.integer(cB*100), nrow(cB), ncol(cB))
##	CA(object)[snps, ] <- cA
##	CB(object)[snps, ] <- cB
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
	marker.index <- index.list[[strata]]
	batches <- split(seq_along(batch(object)), as.character(batch(object)))
	batches <- batches[sapply(batches, length) >= MIN.SAMPLES]


	ii <- isSnp(object) & chromosome(object) < 23 & !is.na(chromosome(object))
	flags <- as.matrix(flags(object)[ii, ])
	fns <- featureNames(object)[ii]
	fns.noflags <- fns[rowSums(flags, na.rm=T) == 0]
	snp.index <- sample(match(fns.noflags, featureNames(object)), 5000)

	##flags <- as.matrix(snpflags[,])
	##noflags <- rowSums(flags, na.rm=TRUE) == 0  ##NA's for unevaluated batches

	nuA.np <- as.matrix(nuA(object)[marker.index, ])
	phiA.np <- as.matrix(phiA(object)[marker.index, ])
	tau2A.AA <- as.matrix(tau2A.AA(object)[marker.index, ])

	##nuA.np <- phiA.np <- sig2A.np <- matrix(NA, length(marker.index), length(unique(batch(object))))
	## for imputation, we need the corresponding parameters of the snps
	##NN <- min(10e3, length(which(ii & noflags)))
	##snp.ind <- sample(which(ii & noflags), NN)
	nuA.snp <- as.matrix(nuA(object)[snp.index, ])
	nuB.snp <- as.matrix(nuB(object)[snp.index, ])
	phiA.snp <- as.matrix(phiA(object)[snp.index, ])
	phiB.snp <- as.matrix(phiB(object)[snp.index, ])
	medianA.AA <- as.matrix(medianA.AA(object)[snp.index,])
	medianB.BB <- as.matrix(medianB.BB(object)[snp.index,])



	medianA.AA.np <- as.matrix(medianA.AA(object)[marker.index,])

##	nnuA.snp <- as.matrix(physical(lM(object))$nuA[snp.ind,])
##	pphiA.snp <- as.matrix(physical(lM(object))$phiA[snp.ind,])
##	nnuB.snp <- as.matrix(physical(lM(object))$nuB[snp.ind,])
##	pphiB.snp <- as.matrix(physical(lM(object))$phiB[snp.ind,])

##	AA.snp <- as.matrix(A(object)[snp.ind, ])
##	BB.snp <- as.matrix(B(object)[snp.ind, ])
##	NNORM.snp <- as.matrix(normal[snp.ind, ])
##	NORM.np <- as.matrix(normal[snps, ])
##	AA.np <- as.matrix(A(object)[marker.index, ])
##	GG <- as.matrix(calls(object)[snp.ind, ])
##	CP <- as.matrix(snpCallProbability(object)[snp.ind, ])
	for(k in seq_along(batches)){
		B <- batches[[k]]
		this.batch <- unique(as.character(batch(object)[B]))
##		phiA.snp <- pphiA.snp[, J]
##		phiB.snp <- pphiB.snp[, J]
##		A.snp <- AA.snp[, k]
##		B.snp <- BB.snp[, k]
##		NORM.snp <- NNORM.snp[, k]
##		G <- GG[, k]
##		xx <- CP[, k]
##		highConf <- (1-exp(-xx/1000)) > GT.CONF.THR
##		G <- G*highConf*NORM.snp
##		G[G==0] <- NA
		##nonpolymorphic
##		A.np <- AA.np[, k]
##		Ns <- applyByGenotype(matrix(1, nrow(G), ncol(G)), rowSums, G)
##		muA <- applyByGenotype(A.snp, rowMedians, G)
##		muB <- applyByGenotype(B.snp, rowMedians, G)
##		muA <- muA[, 1]
##		muB <- muB[, 3]
##		X <- cbind(1, log2(c(muA, muB)))
		X <- cbind(1, log2(c(medianA.AA[, k], medianB.BB[, k])))
		Y <- log2(c(phiA.snp, phiB.snp))
		betahat <- solve(crossprod(X), crossprod(X, Y))
		##
##		mus <- rowMedians(A.np * NORM.np[, k], na.rm=TRUE)
##		averaging across markers, is there a difference in the
##		typical AA intensity for SNPs and the AA intensity for
##		nonpolymorphic loci
##		crosshyb <- max(median(muA) - median(mus), 0)
		crosshyb <- max(median(medianA.AA[, k]) - median(medianA.AA.np[, k]), 0)
##		X <- cbind(1, log2(mus+crosshyb))
		X <- cbind(1, log2(medianA.AA.np[, k] + crosshyb))
		logPhiT <- X %*% betahat
		phiA.np[, k] <- 2^(logPhiT)
		nuA.np[, k] <- medianA.AA.np[,k]-2*phiA.np[, k]
##		cA[, k] <- 1/phiA.np[, J] * (A.np - nuA.np[, J])
##		sig2A.np[, J] <- rowMAD(log2(A.np*NORM.np[, k]), na.rm=TRUE)
##		rm(NORM.snp, highConf, xx, G, Ns, A.np, X, Y, betahat, mus, logPhiT)
##		gc()
	}
	if(THR.NU.PHI){
		nuA.np[nuA.np < MIN.NU] <- MIN.NU
		phiA.np[phiA.np < MIN.PHI] <- MIN.PHI
	}
##	cA[cA < 0.05] <- 0.05
##	cA[cA > 5] <-  5
##	cA <- matrix(as.integer(cA*100), nrow(cA), ncol(cA))
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

##		A <- log2(A); B <- log2(B)
##		tau2A.AA <- summaryStats(G.AA, A, FUNS="rowMAD")^2
##		tau2A.BB <- summaryStats(G.BB, A, FUNS="rowMAD")^2
##		tau2B.AA <- summaryStats(G.AA, B, FUNS="rowMAD")^2
##		tau2B.BB <- summaryStats(G.BB, B, FUNS="rowMAD")^2
		##tau2A <- cbind(tau2A.AA, tau2A.BB)
		##tau2B <- cbind(tau2B.AA, tau2B.BB)
		NN.M <- cbind(N.AA.M, N.AB.M, N.BB.M)
		NN.Mlist[[k]] <- NN.M

		shrink.madA[[k]] <- shrink(madA, NN.M, DF.PRIOR)
		shrink.madB[[k]] <- shrink(madB, NN.M, DF.PRIOR)

##		shrink.tau2A.BB[, k] <- shrink(tau2A.BB[, k, drop=FALSE], NN.M[, 3], DF.PRIOR)[, drop=FALSE]
##		shrink.tau2B.AA[, k] <- shrink(tau2B.AA[, k, drop=FALSE], NN.M[, 1], DF.PRIOR)[, drop=FALSE]
##		shrink.tau2A.AA[, k] <- shrink(tau2A.AA[, k, drop=FALSE], NN.M[, 1], DF.PRIOR)[, drop=FALSE]
##		shrink.tau2B.BB[, k] <- shrink(tau2B.BB[, k, drop=FALSE], NN.M[, 3], DF.PRIOR)[, drop=FALSE]

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
	TRUE
}

fit.lm4 <- function(strata,
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


# steps: quantile normalize hapmap: create 1m_reference_cn.rda object
##cnrma <- function(filenames, cdfName, row.names, sns, seed=1, verbose=FALSE){
##	if(missing(cdfName)) stop("must specify cdfName")
##	pkgname <- getCrlmmAnnotationName(cdfName)
##	require(pkgname, character.only=TRUE) || stop("Package ", pkgname, " not available")
##	if (missing(sns)) sns <- basename(filenames)
##        loader("npProbesFid.rda", .crlmmPkgEnv, pkgname)
##	fid <- getVarInEnv("npProbesFid")
##	fid <- fid[match(row.names, names(fid))]
##	set.seed(seed)
##	idx2 <- sample(length(fid), 10^5) ##for skewness. no need to do everything
##	SKW <- vector("numeric", length(filenames))
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
##	##res3 <- list(NP=NP, SKW=SKW)
##	cat("\n")
##	return(NP)
##}

##getFlags <- function(object, PHI.THR){
##	batch <- unique(object$batch)
##	nuA <- getParam(object, "nuA", batch)
##	nuB <- getParam(object, "nuB", batch)
##	phiA <- getParam(object, "phiA", batch)
##	phiB <- getParam(object, "phiB", batch)
##	negativeNus <- nuA < 1 | nuB < 1
##	negativePhis <- phiA < PHI.THR | phiB < PHI.THR
##	negativeCoef <- negativeNus | negativePhis
##	notfinitePhi <- !is.finite(phiA) | !is.finite(phiB)
##	flags <- negativeCoef | notfinitePhi
##	return(flags)
##}


##instantiateObjects <- function(object, cnOptions){
##	Ns <- matrix(NA, nrow(object), 5)
##	colnames(Ns) <- c("A", "B", "AA", "AB", "BB")
##	vA <- vB <- muB <- muA <- Ns
##	normal <- matrix(TRUE, nrow(object), ncol(object))
##	dimnames(normal) <- list(featureNames(object), sampleNames(object))
##	tmp.objects <- list(vA=vA,
##			    vB=vB,
##			    muB=muB,
##			    muA=muA,
##			    Ns=Ns,
##			    normal=normal)
##        return(tmp.objects)
##}
##
##thresholdCopynumber <- function(object){
##	ca <- CA(object)
##	cb <- CB(object)
##	ca[ca < 0.05] <- 0.05
##	ca[ca > 5] <- 5
##	cb[cb < 0.05] <- 0.05
##	cb[cb > 5] <- 5
##	CA(object) <- ca
##	CB(object) <- cb
##	return(object)
##}
##
####linear model parameters
##lm.parameters <- function(object, batch){##cnOptions){
##	fD <- fData(object)
##	##batch <- object$batch
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
##
##nonpolymorphic <- function(object, cnOptions, tmp.objects){
##	batch <- unique(object$batch)
##	CHR <- unique(chromosome(object))
##	goodSnps <- function(object, PHI.THR, tmp.objects, nAA.THR, nBB.THR){
##		Ns <- tmp.objects[["Ns"]]
##		##Ns <- get("Ns", envir)
##		flags <- getFlags(object, PHI.THR)
##		fewAA <- Ns[, "AA"] < nAA.THR
##		fewBB <- Ns[, "BB"] < nBB.THR
##		flagsA <- flags | fewAA
##		flagsB <- flags | fewBB
##		flags <- list(A=flagsA, B=flagsB)
##		return(flags)
##	}
##	nAA.THR <- cnOptions$nHOM.THR
##	nBB.THR <- cnOptions$nHOM.THR
##	PHI.THR <- cnOptions$PHI.THR
##	snpflags <- goodSnps(object, PHI.THR, tmp.objects, nAA.THR, nBB.THR)
##	flagsA <- snpflags$A
##	flagsB <- snpflags$B
####	if(all(flagsA) | all(flagsB)) stop("all snps are flagged")
##	nuA <- getParam(object, "nuA", batch)
##	nuB <- getParam(object, "nuB", batch)
##	phiA <- getParam(object, "phiA", batch)
##	phiB <- getParam(object, "phiB", batch)
##	sns <- sampleNames(object)
##	muA <- tmp.objects[["muA"]]
##	muB <- tmp.objects[["muB"]]
##	A <- A(object)
##	B <- B(object)
####	CA <- CA(object)
####	CB <- CB(object)
##	if(CHR == 23){
##		phiAX <- getParam(object, "phiAX", batch)
##		phiBX <- getParam(object, "phiBX", batch)
##	}
##	##---------------------------------------------------------------------------
##	## Train on unflagged SNPs
##	##---------------------------------------------------------------------------
##	##Might be best to train using the X chromosome, since for the
##	##X phi and nu have been adjusted for cross-hybridization
##	##plateInd <- plate == uplate[p]
##	##muA <- muA[!flagsA, p, c("A", "AA")]
##	##muB <- muB[!flagsB, p, c("B", "BB")]
##	muA <- muA[!flagsA, "AA"]
##	muB <- muB[!flagsB, "BB"]
##	X <- cbind(1, log2(c(muA, muB)))
##	Y <- log2(c(phiA[!flagsA], phiB[!flagsB]))
##	if(nrow(X) > 5000){
##		ix <- sample(1:nrow(X), 5000)
##	} else {
##		ix <- 1:nrow(X)
##	}
##	betahat <- solve(crossprod(X[ix, ]), crossprod(X[ix, ], Y[ix]))
##	normal <- tmp.objects[["normal"]][!isSnp(object), ]
##	if(CHR == 23){
##		##normalNP <- envir[["normalNP"]]
##		##normalNP <- normalNP[, plate==uplate[p]]
##		##nuT <- envir[["nuT"]]
##		##phiT <- envir[["phiT"]]
##
##		##cnvs <- envir[["cnvs"]]
##                ##loader("cnProbes.rda", pkgname=pkgname, envir=.crlmmPkgEnv)
##                ##cnProbes <- get("cnProbes", envir=.crlmmPkgEnv)
##		##cnProbes <- cnProbes[match(cnvs, rownames(cnProbes)), ]
##
##		##For build Hg18
##		##http://genome.ucsc.edu/cgi-bin/hgGateway
##		##pseudo-autosomal regions on X
##		##chrX:1-2,709,520 and chrX:154584237-154913754, respectively
##		##par:pseudo-autosomal regions
##		pseudoAR <- position(object) < 2709520 | (position(object) > 154584237 & position(object) < 154913754)
##		##pseudoAR <- cnProbes[, "position"] < 2709520 | (cnProbes[, "position"] > 154584237 & cnProbes[, "position"] < 154913754)
##		##in case some of the cnProbes are not annotated
##		pseudoAR[is.na(pseudoAR)] <- FALSE
##		pseudoAR <- pseudoAR[!isSnp(object)]
##		##gender <- envir[["gender"]]
##		gender <- object$gender
##		obj1 <- object[!isSnp(object), ]
##		A.male <- A(obj1[, gender==1])
##		mu1 <- rowMedians(A.male, na.rm=TRUE)
##		##mu1 <- rowMedians(NP[, gender=="male"], na.rm=TRUE)
##		##mu2 <- rowMedians(NP[, gender=="female"], na.rm=TRUE)
##		A.female <- A(obj1[, gender==2])
##		mu2 <- rowMedians(A.female, na.rm=TRUE)
##		mus <- log2(cbind(mu1, mu2))
##		X.men <- cbind(1, mus[, 1])
##		X.fem <- cbind(1, mus[, 2])
##
##		Yhat1 <- as.numeric(X.men %*% betahat)
##		Yhat2 <- as.numeric(X.fem %*% betahat)
##		phi1 <- 2^(Yhat1)
##		phi2 <- 2^(Yhat2)
##		nu1 <- 2^(mus[, 1]) - phi1
##		nu2 <- 2^(mus[, 2]) - 2*phi2
##
##		if(any(pseudoAR)){
##			nu1[pseudoAR] <- 2^(mus[pseudoAR, 1]) - 2*phi1[pseudoAR]
##		}
####		CT1 <- 1/phi1*(A.male-nu1)
####		CT2 <- 1/phi2*(A.female-nu2)
####		CA <- CA(obj1)
####		CA[, gender==1] <- CT1
####		CA[, gender==2] <- CT2
####		CA(object)[!isSnp(object), ] <- CA
##		##only using females to compute the variance
##		##normalNP[, gender=="male"] <- NA
##		normal[, gender==1] <- NA
##		sig2A <- getParam(object, "sig2A", batch)
##		normal.f <- normal[, object$gender==2]
##		sig2A[!isSnp(object)] <- rowMAD(log2(A.female*normal.f), na.rm=TRUE)^2
##		sig2A[!isSnp(object) & is.na(sig2A)] <- median(sig2A[!isSnp(object)], na.rm=TRUE)
##		##sig2T[, p] <- rowMAD(log2(NP*normalNP), na.rm=TRUE)^2
##		object <- pr(object, "sig2A", batch, sig2A)
##
##		nuA[!isSnp(object)] <- nu2
##		phiA[!isSnp(object)] <- phi2
##
##		THR.NU.PHI <- cnOptions$THR.NU.PHI
##		if(THR.NU.PHI){
##			verbose <- cnOptions$verbose
##			##Assign values to object
##			object <- pr(object, "nuA", batch, nuA)
##			object <- pr(object, "phiA", batch, phiA)
##			##if(verbose) message("Thresholding nu and phi")
##			object <- thresholdModelParams(object, cnOptions)
##		} else {
##			object <- pr(object, "nuA", batch, nuA)
##			object <- pr(object, "phiA", batch, phiA)
##		}
##	} else {
##		A <- A(object)[!isSnp(object), ]
##		mus <- rowMedians(A * normal, na.rm=TRUE)
##		crosshyb <- max(median(muA) - median(mus), 0)
##		X <- cbind(1, log2(mus+crosshyb))
##		logPhiT <- X %*% betahat
##		phiA[!isSnp(object)] <- 2^(logPhiT)
##		nuA[!isSnp(object)] <- mus-2*phiA[!isSnp(object)]
##
##		THR.NU.PHI <- cnOptions$THR.NU.PHI
##		if(THR.NU.PHI){
##			verbose <- cnOptions$verbose
##			##Assign values to object
##			object <- pr(object, "nuA", batch, nuA)
##			object <- pr(object, "phiA", batch, phiA)
##			##if(verbose) message("Thresholding nu and phi")
##			object <- thresholdModelParams(object, cnOptions)
##			##reassign values (now thresholded at MIN.NU and MIN.PHI
##			nuA <- getParam(object, "nuA", batch)
##			phiA <- getParam(object, "phiA", batch)
##		}
##		##CA(object)[!isSnp(object), ] <- 1/phiA[!isSnp(object)]*(A - nuA[!isSnp(object)])
##		sig2A <- getParam(object, "sig2A", batch)
##		sig2A[!isSnp(object)] <- rowMAD(log2(A*normal), na.rm=TRUE)^2
##		object <- pr(object, "sig2A", batch, sig2A)
##		##added
##		object <- pr(object, "nuA", batch, nuA)
##		object <- pr(object, "phiA", batch, phiA)
##	}
##	return(object)
##}
##
##sufficient statistics on the intensity scale
##withinGenotypeMoments <- function(object, cnOptions, tmp.objects){
##	normal <- tmp.objects[["normal"]]
##	## muA, muB: robust estimates of the within-genotype center (intensity scale)
##	muA <- tmp.objects[["muA"]]
##	muB <- tmp.objects[["muB"]]
##	## vA, vB: robust estimates of the within-genotype variance (intensity scale)
##	vA <- tmp.objects[["vA"]]
##	vB <- tmp.objects[["vB"]]
##	Ns <- tmp.objects[["Ns"]]
##	G <- snpCall(object)
##	GT.CONF.THR <- cnOptions$GT.CONF.THR
##	CHR <- unique(chromosome(object))
##	A <- A(object)
##	B <- B(object)
####	highConf <- (1-exp(-confs(object)/1000)) > GT.CONF.THR
##	xx <- snpCallProbability(object)
##	highConf <- (1-exp(-xx/1000)) > GT.CONF.THR
##	##highConf <- confs(object) > GT.CONF.THR
##	##highConf <- highConf > GT.CONF.THR
##	if(CHR == 23){
##		gender <- object$gender
####		gender <- envir[["gender"]]
##		IX <- matrix(gender, nrow(G), ncol(G), byrow=TRUE)
####		IX <- IX == "female"
##		IX <- IX == 2  ##2=female, 1=male
##	} else IX <- matrix(TRUE, nrow(G), ncol(G))
##	index <- GT.B <- GT.A <- vector("list", 3)
##	names(index) <- names(GT.B) <- names(GT.A) <- c("AA", "AB", "BB")
##	##--------------------------------------------------
##	##within-genotype sufficient statistics
##	##--------------------------------------------------
##	##GT.B <- GT.A <- list()
##	snpIndicator <- matrix(isSnp(object), nrow(object), ncol(object)) ##RS: added
##	for(j in 1:3){
##		GT <- G==j & highConf & IX & snpIndicator
##		GT <- GT * normal
##		Ns[, j+2] <- rowSums(GT, na.rm=TRUE)
##		GT[GT == FALSE] <- NA
##		GT.A[[j]] <- GT*A
##		GT.B[[j]] <- GT*B
##		index[[j]] <- which(Ns[, j+2] > 0 & isSnp(object)) ##RS: added
##		muA[, j+2] <- rowMedians(GT.A[[j]], na.rm=TRUE)
##		muB[, j+2] <- rowMedians(GT.B[[j]], na.rm=TRUE)
##		vA[, j+2] <- rowMAD(GT.A[[j]], na.rm=TRUE)
##		vB[, j+2] <- rowMAD(GT.B[[j]], na.rm=TRUE)
##
##		##Shrink towards the typical variance
##		DF <- Ns[, j+2]-1
##		DF[DF < 1] <- 1
##		v0A <- median(vA[, j+2], na.rm=TRUE)
##		v0B <- median(vB[, j+2], na.rm=TRUE)
##		if(v0A == 0) v0A <- NA
##		if(v0B == 0) v0B <- NA
##		DF.PRIOR <- cnOptions$DF.PRIOR
##		vA[, j+2] <- (vA[, j+2]*DF + v0A*DF.PRIOR)/(DF.PRIOR+DF)
##		vA[is.na(vA[, j+2]), j+2] <- v0A
##		vB[, j+2] <- (vB[, j+2]*DF + v0B*DF.PRIOR)/(DF.PRIOR+DF)
##		vB[is.na(vB[, j+2]), j+2] <- v0B
##	}
##	if(CHR == 23){
##		k <- 1
##		for(j in c(1,3)){
##			GT <- G==j & highConf & !IX
##			Ns[, k] <- rowSums(GT)
##			GT[GT == FALSE] <- NA
##			muA[, k] <- rowMedians(GT*A, na.rm=TRUE)
##			muB[, k] <- rowMedians(GT*B, na.rm=TRUE)
##			vA[, k] <- rowMAD(GT*A, na.rm=TRUE)
##			vB[, k] <- rowMAD(GT*B, na.rm=TRUE)
##
##			DF <- Ns[, k]-1
##			DF[DF < 1] <- 1
##			v0A <- median(vA[, k], na.rm=TRUE)
##			v0B <- median(vB[, k], na.rm=TRUE)
##			vA[, k] <- (vA[, k]*DF + v0A*DF.PRIOR)/(DF.PRIOR+DF)
##			vA[is.na(vA[, k]), k] <- v0A
##			vB[, k] <- (vB[, k]*DF + v0B*DF.PRIOR)/(DF.PRIOR+DF)
##			vB[is.na(vB[, k]), k] <- v0B
##			k <- k+1
##		}
##	}
##	tmp.objects[["Ns"]] <- Ns
##	tmp.objects[["vA"]] <- vA
##	tmp.objects[["vB"]] <- vB
##	tmp.objects[["muA"]] <- muA
##	tmp.objects[["muB"]] <- muB
##	tmp.objects$index <- index
##	tmp.objects$GT.A <- GT.A
##	tmp.objects$GT.B <- GT.B
##	return(tmp.objects)
##}

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
		if(sum(unobserved.index[[j]]) == 0) next()
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

##oneBatch <- function(object, cnOptions, tmp.objects){
##	muA <- tmp.objects[["muA"]]
##	muB <- tmp.objects[["muB"]]
##	Ns <- tmp.objects[["Ns"]]
##	CHR <- unique(chromosome(object))
##	##---------------------------------------------------------------------------
##	## Impute sufficient statistics for unobserved genotypes (plate-specific)
##	##---------------------------------------------------------------------------
##	MIN.OBS <- cnOptions$MIN.OBS
##	index.AA <- which(Ns[, "AA"] >= MIN.OBS)
##	index.AB <- which(Ns[, "AB"] >= MIN.OBS)
##	index.BB <- which(Ns[, "BB"] >= MIN.OBS)
##	correct.orderA <- muA[, "AA"] > muA[, "BB"]
##	correct.orderB <- muB[, "BB"] > muB[, "AA"]
##	##For chr X, this will ignore the males
##	nobs <- rowSums(Ns[, 3:5] >= MIN.OBS, na.rm=TRUE) == 3
##	index.complete <- which(correct.orderA & correct.orderB & nobs) ##be selective here
##	size <- min(5000, length(index.complete))
##	if(size == 5000) index.complete <- sample(index.complete, 5000)
##	if(length(index.complete) < 200){
##		stop("fewer than 200 snps pass criteria for predicting the sufficient statistics")
##	}
##	index <- tmp.objects[["index"]]
##	index[[1]] <- which(Ns[, "AA"] == 0 & (Ns[, "AB"] >= MIN.OBS & Ns[, "BB"] >= MIN.OBS))
##	index[[2]] <- which(Ns[, "AB"] == 0 & (Ns[, "AA"] >= MIN.OBS & Ns[, "BB"] >= MIN.OBS))
##	index[[3]] <- which(Ns[, "BB"] == 0 & (Ns[, "AB"] >= MIN.OBS & Ns[, "AA"] >= MIN.OBS))
##	mnA <- muA[, 3:5]
##	mnB <- muB[, 3:5]
##	for(j in 1:3){
##		if(length(index[[j]]) == 0) next()
##		X <- cbind(1, mnA[index.complete,  -j, drop=FALSE], mnB[index.complete,  -j, drop=FALSE])
##		Y <- cbind(mnA[index.complete, j], mnB[index.complete,  j])
##		betahat <- solve(crossprod(X), crossprod(X,Y))
##		X <- cbind(1, mnA[index[[j]],  -j, drop=FALSE],  mnB[index[[j]],  -j, drop=FALSE])
##		mus <- X %*% betahat
##		muA[index[[j]], j+2] <- mus[, 1]
##		muB[index[[j]], j+2] <- mus[, 2]
##	}
##	nobsA <- Ns[, "A"] > MIN.OBS
##	nobsB <- Ns[, "B"] > MIN.OBS
##	notMissing <- !(is.na(muA[, "A"]) | is.na(muA[, "B"]) | is.na(muB[, "A"]) | is.na(muB[, "B"]))
##	complete <- list()
##	complete[[1]] <- which(correct.orderA & correct.orderB & nobsA & notMissing) ##be selective here
##	complete[[2]] <- which(correct.orderA & correct.orderB & nobsB & notMissing) ##be selective here
##	size <- min(5000, length(complete[[1]]))
##	if(size > 5000) complete <- lapply(complete, function(x) sample(x, size))
##	if(CHR == 23){
##		index <- list()
##		index[[1]] <- which(Ns[, "A"] == 0)
##		index[[2]] <- which(Ns[, "B"] == 0)
##		cols <- 2:1
##		for(j in 1:2){
##			if(length(index[[j]]) == 0) next()
##			X <- cbind(1, muA[complete[[j]], cols[j]], muB[complete[[j]], cols[j]])
##			Y <- cbind(muA[complete[[j]], j], muB[complete[[j]], j])
##			betahat <- solve(crossprod(X), crossprod(X,Y))
##			X <- cbind(1, muA[index[[j]], cols[j]],  muB[index[[j]], cols[j]])
##			mus <- X %*% betahat
##			muA[index[[j]], j] <- mus[, 1]
##			muB[index[[j]], j] <- mus[, 2]
##		}
##	}
##	##missing two genotypes
##	noAA <- Ns[, "AA"] < MIN.OBS
##	noAB <- Ns[, "AB"] < MIN.OBS
##	noBB <- Ns[, "BB"] < MIN.OBS
##	index[[1]] <- noAA & noAB
##	index[[2]] <- noBB & noAB
##	index[[3]] <- noAA & noBB
####	snpflags <- envir[["snpflags"]]
##	snpflags <- index[[1]] | index[[2]] | index[[3]]
####	snpflags[, p] <- index[[1]] | index[[2]] | index[[3]]
##
##	##---------------------------------------------------------------------------
##	## Two genotype clusters not observed -- would sequence help? (didn't seem to)
##	## 1. extract index of complete data
##	## 2. Regress  mu_missing ~ sequence + mu_observed
##	## 3. solve for nu assuming the median is 2
##	##---------------------------------------------------------------------------
##	cols <- c(3, 1, 2)
##	for(j in 1:3){
##		if(sum(index[[j]]) == 0) next()
##		k <- cols[j]
##		X <- cbind(1, mnA[index.complete, k], mnB[index.complete, k])
##		Y <- cbind(mnA[index.complete,  -k],
##			   mnB[index.complete,  -k])
##		betahat <- solve(crossprod(X), crossprod(X,Y))
##		X <- cbind(1, mnA[index[[j]],  k], mnB[index[[j]],  k])
##		mus <- X %*% betahat
##		muA[index[[j]], -c(1, 2, k+2)] <- mus[, 1:2]
##		muB[index[[j]], -c(1, 2, k+2)] <- mus[, 3:4]
##	}
##	negA <- rowSums(muA < 0) > 0
##	negB <- rowSums(muB < 0) > 0
##	snpflags <- snpflags| negA | negB | rowSums(is.na(muA[, 3:5]), na.rm=TRUE) > 0
##	tmp.objects$snpflags <- snpflags
##	tmp.objects[["muA"]] <- muA
##	tmp.objects[["muB"]] <- muB
##	tmp.objects[["snpflags"]] <- snpflags
##	return(tmp.objects)
##}
##
##Estimate tau2, sigma2, and correlation (updates the object)
##locationAndScale <- function(object, cnOptions, tmp.objects){
##	DF.PRIOR <- cnOptions$DF.PRIOR
##	Ns <- tmp.objects[["Ns"]]
##	index <- tmp.objects[["index"]]
##	index.AA <- index[[1]]
##	index.AB <- index[[2]]
##	index.BB <- index[[3]]
##	rm(index); ##gc()
##
##	GT.A <- tmp.objects[["GT.A"]]
##	GT.B <- tmp.objects[["GT.B"]]
##	AA.A <- GT.A[[1]]
##	AB.A <- GT.A[[2]]
##	BB.A <- GT.A[[3]]
##
##	AA.B <- GT.B[[1]]
##	AB.B <- GT.B[[2]]
##	BB.B <- GT.B[[3]]
##	x <- BB.A[index.BB, ]
##	##batch <- unique(object$batch)
##	batch <- unique(batch(object))
##	tau2A <- getParam(object, "tau2A", batch)
##	tau2A[index.BB] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2
##	DF <- Ns[, "BB"]-1
##	DF[DF < 1] <- 1
##	med <- median(tau2A, na.rm=TRUE)
##	tau2A <- (tau2A * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
##	tau2A[is.na(tau2A) & isSnp(object)] <- med
##	object <- pr(object, "tau2A", batch, tau2A)
##
##	sig2B <- getParam(object, "sig2B", batch)
##	x <- BB.B[index.BB, ]
##	sig2B[index.BB] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2
##	med <- median(sig2B, na.rm=TRUE)
##	sig2B <- (sig2B * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
##	sig2B[is.na(sig2B) & isSnp(object)] <- med
##	object <- pr(object, "sig2B", batch, sig2B)
##
##	tau2B <- getParam(object, "tau2B", batch)
##	x <- AA.B[index.AA, ]
##	tau2B[index.AA] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2
##	DF <- Ns[, "AA"]-1
##	DF[DF < 1] <- 1
##	med <- median(tau2B, na.rm=TRUE)
##	tau2B <- (tau2B * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
##	tau2B[is.na(tau2B) & isSnp(object)] <- med
##	object <- pr(object, "tau2B", batch, tau2B)
##
##	sig2A <- getParam(object, "sig2A", batch)
##	x <- AA.A[index.AA, ]
##	sig2A[index.AA] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2##var(log(IA)|AA)
##	med <- median(sig2A, na.rm=TRUE)
##	sig2A <- (sig2A * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
##	sig2A[is.na(sig2A) & isSnp(object)] <- med
##	object <- pr(object, "sig2A", batch, sig2A)
##
##	if(length(index.AB) > 0){ ##all homozygous is possible
##		corr <- getParam(object, "corr", batch)
##		x <- AB.A[index.AB, ]
##		y <- AB.B[index.AB, ]
##		corr[index.AB] <- rowCors(x, y, na.rm=TRUE)
##		corr[corr < 0] <- 0
##		DF <- Ns[, "AB"]-1
##		DF[DF<1] <- 1
##		med <- median(corr, na.rm=TRUE)
##		corr <- (corr*DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
##		corr[is.na(corr) & isSnp(object)] <- med
##		object <- pr(object, "corr", batch, corr)
##	}
##	corrB.AA <- getParam(object, "corrB.AA", batch)
##	backgroundB <- AA.B[index.AA, ]
##	signalA <- AA.A[index.AA, ]
##	corrB.AA[index.AA] <- rowCors(backgroundB, signalA, na.rm=TRUE)
##	DF <- Ns[, "AA"]-1
##	DF[DF < 1] <- 1
##	med <- median(corrB.AA, na.rm=TRUE)
##	corrB.AA <- (corrB.AA*DF + med*DF.PRIOR)/(DF.PRIOR + DF)
##	corrB.AA[is.na(corrB.AA) & isSnp(object)] <- med
##	object <- pr(object, "corrB.AA", batch, corrB.AA)
##
##	corrA.BB <- getParam(object, "corrA.BB", batch)
##	backgroundA <- BB.A[index.BB, ]
##	signalB <- BB.B[index.BB, ]
##	corrA.BB[index.BB] <- rowCors(backgroundA, signalB, na.rm=TRUE)
##	DF <- Ns[, "BB"]-1
##	DF[DF < 1] <- 1
##	med <- median(corrA.BB, na.rm=TRUE)
##	corrA.BB <- (corrA.BB*DF + med*DF.PRIOR)/(DF.PRIOR + DF)
##	corrA.BB[is.na(corrA.BB) & isSnp(object)] <- med
##	object <- pr(object, "corrA.BB", batch, corrA.BB)
##	return(object)
##}
##
##coefs <- function(object, cnOptions, tmp.objects){
##	##batch <- unique(object$batch)
##	batch <- unique(batch(object))
##	CHR <- unique(chromosome(object))
##	muA <- tmp.objects[["muA"]]
##	muB <- tmp.objects[["muB"]]
##	vA <- tmp.objects[["vA"]]
##	vB <- tmp.objects[["vB"]]
##	Ns <- tmp.objects[["Ns"]]
##	if(CHR != 23){
##		IA <- muA[, 3:5]
##		IB <- muB[, 3:5]
##		vA <- vA[, 3:5]
##		vB <- vB[, 3:5]
##		Np <- Ns[, 3:5]
##	} else {
##		NOHET <- is.na(median(vA[, "AB"], na.rm=TRUE))
##		if(NOHET){
##			IA <- muA[, -4]
##			IB <- muB[, -4]
##			vA <- vA[, -4]
##			vB <- vB[, -4]
##			Np <- Ns[, -4]
##		} else{
##			IA <- muA
##			IB <- muB
##			vA <- vA
##			vB <- vB
##			Np <- Ns
##		}
##
##	}
##	Np[Np < 1] <- 1
##	vA2 <- vA^2/Np
##	vB2 <- vB^2/Np
##	wA <- sqrt(1/vA2)
##	wB <- sqrt(1/vB2)
##	YA <- IA*wA
##	YB <- IB*wB
##	##update lm.coefficients stored in object
##	object <- nuphiAllele(object, allele="A", Ystar=YA, W=wA, tmp.objects=tmp.objects, cnOptions=cnOptions)
##	object <- nuphiAllele(object, allele="B", Ystar=YB, W=wB, tmp.objects=tmp.objects, cnOptions=cnOptions)
##	##---------------------------------------------------------------------------
##	##Estimate crosshyb using X chromosome and sequence information
##	##---------------------------------------------------------------------------
##	##browser()
##	####data(sequences, package="genomewidesnp6Crlmm")
##	##snpflags <- envir[["snpflags"]]
##	##muA <- envir[["muA"]][, p, 3:5]
##	##muB <- envir[["muB"]][, p, 3:5]
##	##Y <- envir[["phiAx"]]
##	##load("sequences.rda")
##	##seqA <- sequences[, "A", ][, 1]
##	##seqA <- seqA[match(snps, names(seqA))]
##	##X <- cbind(1, sequenceDesignMatrix(seqA))
##	##X <- cbind(X, nuA[, p], phiA[, p], nuB[, p], phiB[, p])
##	##missing <- rowSums(is.na(X)) > 0
##	##betahat <- solve(crossprod(X[!missing, ]), crossprod(X[!missing, ], Y[!missing]))
##	return(object)
##}
##
##polymorphic <- function(object, cnOptions, tmp.objects){
##	##batch <- unique(object$batch)
##	batch <- unique(batch(object))
##	CHR <- unique(chromosome(object))
##	vA <- tmp.objects[["vA"]]
##	vB <- tmp.objects[["vB"]]
##	Ns <- tmp.objects[["Ns"]]
##
##	nuA <- getParam(object, "nuA", batch)
##	nuB <- getParam(object, "nuB", batch)
####	nuA.se <- getParam(object, "nuA.se", batch)
####	nuB.se <- getParam(object, "nuB.se", batch)
##
##	phiA <- getParam(object, "phiA", batch)
##	phiB <- getParam(object, "phiB", batch)
####	phiA.se <- getParam(object, "phiA.se", batch)
####	phiB.se <- getParam(object, "phiB.se", batch)
##	A <- A(object)
##	B <- B(object)
##
##	NOHET <- mean(Ns[, "AB"], na.rm=TRUE) < 0.05
##	##---------------------------------------------------------------------------
##	## Estimate CA, CB
##	##---------------------------------------------------------------------------
##	if(CHR == 23){
##		phiAX <- getParam(object, "phiAX", batch)  ##nonspecific hybridization coef
##		phiBX <- getParam(object, "phiBX", batch)  ##nonspecific hybridization coef
##		phistar <- phiBX/phiA
####		tmp <- (B-nuB - phistar*A + phistar*nuA)/phiB
####		copyB <- tmp/(1-phistar*phiAX/phiB)
####		copyA <- (A-nuA-phiAX*copyB)/phiA
####		CB(object) <- copyB  ## multiplies by 100 and converts to integer
####		CA(object) <- copyA
##	} else{
####		CA(object) <- matrix((1/phiA*(A-nuA)), nrow(A), ncol(A))
####		CB(object) <- matrix((1/phiB*(B-nuB)), nrow(B), ncol(B))
##
##	}
##	return(object)
##}
##
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

##bias1 <- function(strata.index,
##		  snpBatches,
##		  index,
##		  object,
##		  normal,
##		  emit,
##		  prior.prob,
##		  MIN.SAMPLES,
##		  verbose){
##
##}

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


##biasAdjNP <- function(plateIndex, envir, priorProb){
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


##getParams <- function(object, batch){
##	##batch <- unique(object$batch)
##	batch <- unique(batch(object))
##	if(length(batch) > 1) stop("batch variable not unique")
##	nuA <- as.numeric(fData(object)[, match(paste("nuA", batch, sep="_"), fvarLabels(object))])
##	nuB <- as.numeric(fData(object)[, match(paste("nuB", batch, sep="_"), fvarLabels(object))])
##	phiA <- as.numeric(fData(object)[, match(paste("phiA", batch, sep="_"), fvarLabels(object))])
##	phiB <- as.numeric(fData(object)[, match(paste("phiB", batch, sep="_"), fvarLabels(object))])
##	tau2A <- as.numeric(fData(object)[, match(paste("tau2A", batch, sep="_"), fvarLabels(object))])
##	tau2B <- as.numeric(fData(object)[, match(paste("tau2B", batch, sep="_"), fvarLabels(object))])
##	sig2A <- as.numeric(fData(object)[, match(paste("sig2A", batch, sep="_"), fvarLabels(object))])
##	sig2B <- as.numeric(fData(object)[, match(paste("sig2B", batch, sep="_"), fvarLabels(object))])
##	corrA.BB <- as.numeric(fData(object)[, match(paste("corrA.BB", batch, sep="_"), fvarLabels(object))])
##	corrB.AA <- as.numeric(fData(object)[, match(paste("corrB.AA", batch, sep="_"), fvarLabels(object))])
##	corr <- as.numeric(fData(object)[, match(paste("corr", batch, sep="_"), fvarLabels(object))])
##	params <- list(nuA=nuA,
##		       nuB=nuB,
##		       phiA=phiA,
##		       phiB=phiB,
##		       tau2A=tau2A,
##		       tau2B=tau2B,
##		       sig2A=sig2A,
##		       sig2B=sig2B,
##		       corrA.BB=corrA.BB,
##		       corrB.AA=corrB.AA,
##		       corr=corr)
##	return(params)
##}
##
##
#### Constrain nu and phi to positive values
##thresholdModelParams <- function(object, cnOptions){
##	MIN.NU <- cnOptions$MIN.NU
##	MIN.PHI <- cnOptions$MIN.PHI
##	batch <- unique(object$batch)
##	nuA <- getParam(object, "nuA", batch)
##	nuA[nuA < MIN.NU] <- MIN.NU
##	object <- pr(object, "nuA", batch, nuA)
##	nuB <- getParam(object, "nuB", batch)
##	if(!all(is.na(nuB))){
##		nuB[nuB < MIN.NU] <- MIN.NU
##		object <- pr(object, "nuB", batch, nuB)
##	}
##	phiA <- getParam(object, "phiA", batch)
##	phiA[phiA < MIN.PHI] <- MIN.PHI
##	object <- pr(object, "phiA", batch, phiA)
##	phiB <- getParam(object, "phiB", batch)
##	if(!all(is.na(phiB))){
##		phiB[phiB < MIN.PHI] <- MIN.PHI
##		object <- pr(object, "phiB", batch, phiB)
##	}
##	phiAX <- as.numeric(getParam(object, "phiAX", batch))
##	if(!all(is.na(phiAX))){
##		phiAX[phiAX < MIN.PHI] <- MIN.PHI
##		object <- pr(object, "phiAX", batch, phiAX)
##	}
##	phiBX <- as.numeric(getParam(object, "phiBX", batch))
##	if(!all(is.na(phiBX))){
##		phiBX[phiBX < MIN.PHI] <- MIN.PHI
##		object <- pr(object, "phiBX", batch, phiBX)
##	}
##	return(object)
##}
##
##cnCNSet <- function(object, cnOptions){
##	PLATE <- unique(batch(object))
##	verbose <- cnOptions$verbose
##	tmp.objects <- instantiateObjects(object, cnOptions)
##	bias.adj <- cnOptions$bias.adj
##	if(bias.adj & ncol(object) <= 15){
##		warning(paste("bias.adj is TRUE, but too few samples to perform this step"))
##		cnOptions$bias.adj <- bias.adj <- FALSE
##	}
##	if(bias.adj){
##		if(verbose) message("Dropping samples with low posterior prob. of normal copy number (samples dropped is locus-specific)")
##		tmp.objects <- biasAdjNP(object, cnOptions, tmp.objects)
##		tmp.objects <- biasAdj(object, cnOptions, tmp.objects)
##		if(verbose) message("Recomputing location and scale parameters")
##	}
##	##update tmp.objects
##	tmp.objects <- withinGenotypeMoments(object,
##					     cnOptions=cnOptions,
##					     tmp.objects=tmp.objects)
##	object <- locationAndScale(object, cnOptions, tmp.objects)
##	tmp.objects <- oneBatch(object,
##				cnOptions=cnOptions,
##				tmp.objects=tmp.objects)
##	##coefs calls nuphiAllele.
##	object <- coefs(object, cnOptions, tmp.objects)
##	##nuA=getParam(object, "nuA", PLATE)
##	THR.NU.PHI <- cnOptions$THR.NU.PHI
##	if(THR.NU.PHI){
##		verbose <- cnOptions$verbose
##		##if(verbose) message("Thresholding nu and phi")
##		object <- thresholdModelParams(object, cnOptions)
##	}
##	##if(verbose) message("\nAllele specific copy number")
##	object <- polymorphic(object, cnOptions, tmp.objects)
##	if(any(!isSnp(object))){  ##there are nonpolymorphic probes
##		##if(verbose) message("\nCopy number for nonpolymorphic probes...")
##		object <- nonpolymorphic(object, cnOptions, tmp.objects)
##	}
##	##---------------------------------------------------------------------------
##	##Note: the replacement method multiples by 100
####	CA(object)[, batch==PLATE] <- CA(object)
####	CB(object)[, batch==PLATE] <- CB(object)
##	##---------------------------------------------------------------------------
##	##update-the plate-specific parameters for copy number
##	object <- pr(object, "nuA", PLATE, getParam(object, "nuA", PLATE))
##	object <- pr(object, "nuA.se", PLATE, getParam(object, "nuA.se", PLATE))
##	object <- pr(object, "nuB", PLATE, getParam(object, "nuB", PLATE))
##	object <- pr(object, "nuB.se", PLATE, getParam(object, "nuB.se", PLATE))
##	object <- pr(object, "phiA", PLATE, getParam(object, "phiA", PLATE))
##	object <- pr(object, "phiA.se", PLATE, getParam(object, "phiA.se", PLATE))
##	object <- pr(object, "phiB", PLATE, getParam(object, "phiB", PLATE))
##	object <- pr(object, "phiB.se", PLATE, getParam(object, "phiB.se", PLATE))
##	object <- pr(object, "tau2A", PLATE, getParam(object, "tau2A", PLATE))
##	object <- pr(object, "tau2B", PLATE, getParam(object, "tau2B", PLATE))
##	object <- pr(object, "sig2A", PLATE, getParam(object, "sig2A", PLATE))
##	object <- pr(object, "sig2B", PLATE, getParam(object, "sig2B", PLATE))
##	object <- pr(object, "phiAX", PLATE, as.numeric(getParam(object, "phiAX", PLATE)))
##	object <- pr(object, "phiBX", PLATE, as.numeric(getParam(object, "phiBX", PLATE)))
##	object <- pr(object, "corr", PLATE, getParam(object, "corr", PLATE))
##	object <- pr(object, "corrA.BB", PLATE, getParam(object, "corrA.BB", PLATE))
##	object <- pr(object, "corrB.AA", PLATE, getParam(object, "corrB.AA", PLATE))
##	##object <- object[order(chromosome(object), position(object)), ]
####	if(cnOptions[["thresholdCopynumber"]]){
####		object <- thresholdCopynumber(object)
####	}
##	return(object)
##}
##


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
constructIlluminaAssayData <- function(np, snp, object, storage.mode="environment", order.index){
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
	A <- rbind(snp[[1]], np[[1]], deparse.level=0)[order.index, ]
	B <- rbind(snp[[2]], np[[2]], deparse.level=0)[order.index, ]
	gt <- stripnames(calls(object))
	emptyMatrix <- matrix(integer(), nrow(np[[1]]), ncol(A))
	gt <- rbind(gt, emptyMatrix, deparse.level=0)[order.index,]
	pr <- stripnames(snpCallProbability(object))
	pr <- rbind(pr, emptyMatrix, deparse.level=0)[order.index, ]
	aD <- assayDataNew(storage.mode,
			   alleleA=A,
			   alleleB=B,
			   call=gt,
			   callProbability=pr)
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
	new.order <- order(fD$chromosome, fD$position)
	fD <- fD[new.order, ]
	aD <- constructIlluminaAssayData(cnAB, res, crlmmResult, order.index=new.order)
	##protocolData(crlmmResult)$batch <- vector("integer", ncol(crlmmResult))
	batch <- vector("integer", ncol(crlmmResult))
	container <- new("CNSet",
			 call=aD[["call"]],
			 callProbability=aD[["callProbability"]],
			 alleleA=aD[["alleleA"]],
			 alleleB=aD[["alleleB"]],
			 phenoData=phenoData(crlmmResult),
			 protocolData=protocolData(crlmmResult),
			 featureData=fD,
			 batch=batch,
			 annotation="human370v1c")
##	lM(container) <- initializeParamObject(list(featureNames(container), unique(protocolData(container)$batch)))
##	lM(container) <- initializeLmFrom(container)
	container
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
			  type=c("SNP", "NP", "X.SNP", "X.NP"), ##"X.snps", "X.nps"),
			  MIN.OBS=1,
			  MIN.SAMPLES=10,
			  DF.PRIOR=50,
			  verbose=TRUE){
	if(type == "X.SNP" | type=="X.NP"){
		gender <- object$gender
		if(sum(gender == 2) < 3) {
			return("too few females to estimate within genotype summary statistics on CHR X")
		}
		CHR.X <- TRUE
	} else CHR.X <- FALSE
	batch <- batch(object)
	is.snp <- isSnp(object)
	is.autosome <- chromosome(object) < 23
	is.annotated <- !is.na(chromosome(object))
	is.X <- chromosome(object) == 23
	is.lds <- is(calls(object), "ffdf") | is(calls(object), "ff_matrix")
	if(is.lds) require(ff)
	whichMarkers <- function(type, is.snp, is.autosome, is.annotated, is.X){
		switch(type,
		       SNP=which(is.snp & is.autosome & is.annotated),
		       NP=which(!is.snp & is.autosome),
		       X.SNP=which(is.snp & is.X),
		       X.NP=which(!is.snp & is.X),
		       stop("'type' must be one of 'SNP', 'NP', 'X.SNP', or 'X.NP'")
	       )
	}
	marker.index <- whichMarkers(type[[1]], is.snp, is.autosome, is.annotated, is.X)
	summaryFxn <- function(type){
		switch(type,
		       SNP="shrinkGenotypeSummaries",
		       X.SNP="shrinkGenotypeSummaries", ## this shrinks for the females only
##		       NP="summarizeNps",
##		       X.SNP="summarizeSnps",
##		       X.NP="summarizeNps")
		       stop())
	}
	FUN <- summaryFxn(type[[1]])
	if(is.lds){
		index.list <- splitIndicesByLength(marker.index, ocProbesets())
		ocLapply(seq(along=index.list),
			 FUN,
			 index.list=index.list,
			 object=object,
			 verbose=verbose,
			 MIN.OBS=MIN.OBS,
			 MIN.SAMPLES=MIN.SAMPLES,
			 DF.PRIOR=DF.PRIOR,
			 is.lds=is.lds,
			 neededPkgs="crlmm")
	} else {
		FUN <- match.fun(FUN)
		object <- FUN(strata.index=1,
			      index.list=list(marker.index),
			      object=object,
			      MIN.OBS=MIN.OBS,
			      MIN.SAMPLES=MIN.SAMPLES,
			      DF.PRIOR=DF.PRIOR,
			      verbose=verbose,
			      is.lds=is.lds)
	}
	return(object)
}

genotypeSummary <- function(object,
			    GT.CONF.THR=0.95,
			    type=c("SNP", "NP", "X.SNP", "X.NP"), ##"X.snps", "X.nps"),
			    verbose=TRUE){
	if(type == "X.SNP" | type=="X.NP"){
		gender <- object$gender
		if(sum(gender == 2) < 3) {
			return("too few females to estimate within genotype summary statistics on CHR X")
		}
		CHR.X <- TRUE
	} else CHR.X <- FALSE
	batch <- batch(object)
	is.snp <- isSnp(object)
	is.autosome <- chromosome(object) < 23
	is.annotated <- !is.na(chromosome(object))
	is.X <- chromosome(object) == 23
	is.lds <- is(calls(object), "ffdf") | is(calls(object), "ff_matrix")
	if(is.lds) require(ff)
	whichMarkers <- function(type, is.snp, is.autosome, is.annotated, is.X){
		switch(type,
		       SNP=which(is.snp & is.autosome & is.annotated),
		       NP=which(!is.snp & is.autosome),
		       X.SNP=which(is.snp & is.X),
		       X.NP=which(!is.snp & is.X),
		       stop("'type' must be one of 'SNP', 'NP', 'X.SNP', or 'X.NP'")
	       )
	}
	marker.index <- whichMarkers(type[[1]], is.snp, is.autosome, is.annotated, is.X)
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
##			 marker.index=marker.index,
			 object=object,
			 batchSize=ocProbesets(),
			 GT.CONF.THR=GT.CONF.THR,
			 verbose=verbose,
			 is.lds=is.lds,
			 CHR.X=CHR.X,
			 neededPkgs="crlmm")
	} else {
		FUN <- match.fun(FUN)
		object <- FUN(strata.index=1,
			      index.list=list(marker.index),
##			      marker.index=marker.index,
			      object=object,
			      batchSize=ocProbesets(),
			      GT.CONF.THR=GT.CONF.THR,
			      verbose=verbose,
			      CHR.X=CHR.X,
			      is.lds=is.lds)
	}
	return(object)
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
			  batchSize,
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
	GG <- as.matrix(calls(object)[index, ])
	CP <- as.matrix(snpCallProbability(object)[index, ])
	AA <- as.matrix(A(object)[index, ])
	BB <- as.matrix(B(object)[index, ])
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
			    GT.CONF.THR=0.99,
			    PHI.THR=2^6,
			    nHOM.THR=5,
			    MIN.NU=2^3,
			    MIN.PHI=2^3,
			    THR.NU.PHI=TRUE,
			    type=c("SNP", "NP", "X.SNP", "X.NP")){
	if(type == "X.SNP" | type=="X.NP"){
		gender <- object$gender
		if(sum(gender == 2) < 3) {
			warning("too few females to estimate within genotype summary statistics on CHR X")
			return(object)
		}
		CHR.X <- TRUE
	} else CHR.X <- FALSE
	batch <- batch(object)
	is.snp <- isSnp(object)
	is.autosome <- chromosome(object) < 23
	is.annotated <- !is.na(chromosome(object))
	is.X <- chromosome(object) == 23
	is.lds <- is(calls(object), "ffdf") | is(calls(object), "ff_matrix")
	if(is.lds) require(ff)
	samplesPerBatch <- table(as.character(batch(object)))
	if(any(samplesPerBatch < MIN.SAMPLES)){
		warning("The following batches have fewer than ", MIN.SAMPLES, ":")
		message(paste(samplesPerBatch[samplesPerBatch < MIN.SAMPLES], collapse=", "))
		message("Not estimating copy number for the above batches")
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
	marker.index <- whichMarkers(type[[1]], is.snp, is.autosome, is.annotated, is.X)
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
			      batchSize=batchSize,
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
