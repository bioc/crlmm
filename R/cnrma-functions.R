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

nuphiAllele <- function(p, allele, Ystar, W, envir){
	Ns <- envir[["Ns"]]
	CHR <- envir[["chrom"]]
	nuA <- envir[["nuA"]]
	nuB <- envir[["nuB"]]
	nuA.se <- envir[["nuA.se"]]
	nuB.se <- envir[["nuB.se"]]
	phiA <- envir[["phiA"]]
	phiB <- envir[["phiB"]]
	phiA.se <- envir[["phiA.se"]]
	phiB.se <- envir[["phiB.se"]]
	if(CHR==23){
		phiAx <- envir[["phiAx"]]
		phiBx <- envir[["phiBx"]]
	}
	complete <- rowSums(is.na(W)) == 0
	NOHET <- mean(Ns[, p, "AB"], na.rm=TRUE) < 0.05
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
	if(any(!is.finite(W))){## | any(!is.finite(V))){
		i <- which(rowSums(!is.finite(W)) > 0)
		browser()
		stop("Inf values in W or V")
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
		nuA[complete, p] <- betahat[1, ]
		phiA[complete, p] <- betahat[2, ]
		nuA.se[complete, p] <- ses[1, ]
		phiA.se[complete, p] <- ses[2, ]
		envir[["nuA"]] <- nuA
		envir[["phiA"]] <- phiA
		envir[["nuA.se"]] <- nuA.se
		envir[["phiA.se"]] <- phiA.se
		if(CHR == 23){
			phiAx[complete, p] <- betahat[3, ]
			envir[["phiAx"]] <- phiAx
		}
	}
	if(allele=="B"){
		nuB[complete, p] <- betahat[1, ]
		phiB[complete, p] <- betahat[2, ]
		nuB.se[complete, p] <- ses[1, ]
		phiB.se[complete, p] <- ses[2, ]
		envir[["nuB"]] <- nuB
		envir[["phiB"]] <- phiB
		envir[["nuB.se"]] <- nuB.se
		envir[["phiB.se"]] <- phiB.se
		if(CHR == 23){
			phiBx[complete, p] <- betahat[3, ]
			envir[["phiBx"]] <- phiBx
		}
	}
}

celDates <- function(celfiles){
	if(!all(file.exists(celfiles))) stop("1 or more cel file does not exist")
	celdates <- vector("character", length(celfiles))
	celtimes <- vector("character", length(celfiles))
	for(i in seq(along=celfiles)){
		if(i %% 100 == 0) cat(".")
		tmp <- read.celfile.header(celfiles[i], info="full")$DatHeader
		tmp <- strsplit(tmp, "\ +")
		celdates[i] <- tmp[[1]][6]
		celtimes[i] <- tmp[[1]][7]
	}
	tmp <- paste(celdates, celtimes)
	celdts <- strptime(tmp, "%m/%d/%y %H:%M:%S")
	return(celdts)
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
	XMedian <- apply(log2(A[XIndex,, drop=FALSE])+log2(B[XIndex,, drop=FALSE]), 2, median)/2
	SNR <- res$SNR
	if(sum(SNR>SNRMin)==1){
		gender <- which.min(c(abs(XMedian-8.9), abs(XMedian-9.5)))
	}else{
		gender <- kmeans(XMedian, c(min(XMedian[SNR>SNRMin]), max(XMedian[SNR>SNRMin])))[["cluster"]]
	}
	return(gender)
}

combineIntensities <- function(res, cnrmaResult, cdfName){
	rownames(res$B) <- rownames(res$A) <- res$gns
	colnames(res$B) <- colnames(res$A) <- res$sns
	if(!is.null(cnrmaResult)){
		NP <- cnrmaResult$NP
		blank <- matrix(NA, nrow(NP), ncol(NP))
		dimnames(blank) <- dimnames(NP)
		A <- rbind(res$A, NP)
		B <- rbind(res$B, blank)
	} else {
		A <- res$A
		B <- res$B
	}
	aD <- assayDataNew("lockedEnvironment",
			   A=A,
			   B=B)
	ABset <- new("ABset",
		     assayData=aD,
		     featureData=annotatedDataFrameFrom(A, byrow=TRUE),
		     phenoData=annotatedDataFrameFrom(A, byrow=FALSE),
		     annotation=cdfName)
	ABset$SNR <- res$SNR
	ABset$gender <- predictGender(res=res, cdfName=cdfName)
	return(ABset)
}

harmonizeSnpSet <- function(crlmmResult, ABset){
	blank <- matrix(NA, length(cnNames(ABset)), ncol(ABset))
	rownames(blank) <- cnNames(ABset)
	colnames(blank) <- sampleNames(ABset)
	crlmmCalls <- rbind(calls(crlmmResult), blank)
	crlmmConf <- rbind(confs(crlmmResult), blank)
	fD <- as.matrix(fData(crlmmResult))
	fD2 <- matrix(NA, nrow(blank), ncol(fD))
	rownames(fD2) <- rownames(blank)
	fD <- rbind(fD, fD2)
	aD <- assayDataNew("lockedEnvironment",
			   calls=crlmmCalls,
			   callProbability=crlmmConf)
	##Make crlmmResult the same dimension as ABset
	fD <- new("AnnotatedDataFrame",
		  data=data.frame(fD),
		  varMetadata=fvarMetadata(crlmmResult))
	crlmmResult <- new("SnpSet",
			   call=crlmmCalls,
			   callProbability=crlmmConf,
			   featureData=fD,
			   phenoData=phenoData(crlmmResult),
			   protocolData=protocolData(ABset),
			   annotation=annotation(ABset))
	stopifnot(all.equal(dimnames(crlmmResult), dimnames(ABset)))
	crlmmResult
}

harmonizeDimnamesTo <- function(object1, object2){
	#object2 should be a subset of object 1
	object2 <- object2[featureNames(object2) %in% featureNames(object1), ]
	object1 <- object1[match(featureNames(object2), featureNames(object1)), ]
	object1 <- object1[, match(sampleNames(object2), sampleNames(object1))]
	stopifnot(all.equal(featureNames(object1), featureNames(object2)))
	stopifnot(all.equal(sampleNames(object1), sampleNames(object2)))
	return(object1)
}

crlmmIlluminaWrapper <- function(sampleSheet, outdir="./", cdfName,
				 save.intermediate=FALSE,
				 splitByChr=TRUE,...){
	if(file.exists(file.path(outdir, "RG.rda"))) load(file.path(outdir, "RG.rda"))
	else {
		RG <- readIdatFiles(sampleSheet=get("samplesheet5"),
                                    arrayInfoColNames=list(barcode=NULL, position="SentrixPosition"),
                                    saveDate=TRUE, path=get("path"))
		J <- match(c("1_A", "3_A", "5_A", "7_A"), sampleNames(RG))
		RG <- RG[, -J]
		if(save.intermediate) save(RG, file=file.path(outdir, "RG.rda"))  ##935M for 91 samples...better not to save this
	}	
	if(!file.exists(file.path(outdir, "res.rda"))){
		crlmmOut <- crlmmIllumina(RG=RG, cdfName=cdfName,
                                          sns=pData(RG)$ID,
                                          returnParams=TRUE,
                                          save.it=TRUE,
                                          intensityFile=file.path(outdir, "res.rda"))
		if(save.intermediate) save(crlmmOut, file=file.path(outdir, "crlmmOut.rda"))				
	} else{
		message("Loading...")		
		load(file.path(outdir, "res.rda"))
		load(file.path(outdir, "crlmmOut.rda"))		
	}
	ABset <- combineIntensities(get("res"), NULL, cdfName=cdfName)
	protocolData(ABset)[["ScanDate"]] <- as.character(pData(RG)$ScanDate)
	crlmmResult <- harmonizeSnpSet(crlmmOut, ABset)
	stopifnot(all.equal(dimnames(crlmmOut), dimnames(ABset)))
	crlmmList <- list(ABset,
			  crlmmResult)
	crlmmList <- as(crlmmList, "CrlmmSetList")
	if(splitByChr){
		message("Saving by chromosome")
		splitByChromosome(crlmmList, cdfName=cdfName, outdir=outdir)
	} else{
		message("Saving crlmmList object to ", outdir)
		save(crlmmList, file=file.path(outdir, "crlmmList.rda"))
	}
	message("CrlmmSetList objects saved to ", outdir)
}
	
	
	
crlmmWrapper <- function(filenames, outdir="./", cdfName="genomewidesnp6",
			 save.it=FALSE,
			 splitByChr=TRUE, ...){
	##no visible binding for res
	if(!file.exists(file.path(outdir, "crlmmResult.rda"))){
		crlmmResult <- crlmm(filenames=filenames, cdfName=cdfName, save.it=TRUE, ...)
		if(save.it) save(crlmmResult, file=file.path(outdir, "crlmmResult.rda"))
	} else {
		message("Loading crlmmResult...")
		load(file.path(outdir, "crlmmResult.rda"))
	}
	if(!file.exists(file.path(outdir, "cnrmaResult.rda"))){
		message("Quantile normalizing the copy number probes")		
		cnrmaResult <- cnrma(filenames=filenames, cdfName=cdfName)
		if(save.it) save(cnrmaResult, file=file.path(outdir, "cnrmaResult.rda"))
	} else{
		message("Loading cnrmaResult...")		
		load(file.path(outdir, "cnrmaResult.rda"))
	}
	load(file.path(outdir, "intensities.rda"))
	ABset <- combineIntensities(get("res"), cnrmaResult, cdfName)
	protocolData(ABset)[["ScanDate"]] <- as.character(celDates(filenames))	
	crlmmResult <- harmonizeSnpSet(crlmmResult, ABset)
	stopifnot(all.equal(dimnames(crlmmResult), dimnames(ABset)))
	crlmmResults <- list(ABset,
			     crlmmResult)
	crlmmResults <- as(crlmmResults, "CrlmmSetList")
	
	if(splitByChr){
		message("Saving by chromosome")
		splitByChromosome(crlmmResults, cdfName=cdfName, outdir=outdir)
	} else{
		message("Saving crlmmList object to ", outdir)
		save(crlmmResults, file=file.path(outdir, "crlmmResults.rda"))
	}
	if(!save.it){
		message("Cleaning up")
		unlink(file.path(outdir, "intensities.rda"))
	}
	return()
}



cnrma <- function(filenames, cdfName="genomewidesnp6", sns, seed=1, verbose=FALSE){
	pkgname <- getCrlmmAnnotationName(cdfName)
	require(pkgname, character.only=TRUE) || stop("Package ", pkgname, " not available")
	if (missing(sns)) sns <- basename(filenames)
        loader("npProbesFid.rda", .crlmmPkgEnv, pkgname)
	fid <- getVarInEnv("npProbesFid")
	set.seed(seed)
	idx2 <- sample(length(fid), 10^5) ##for skewness. no need to do everything
	SKW <- vector("numeric", length(filenames))
	NP <- matrix(NA, length(fid), length(filenames))
	verbose <- TRUE
	if(verbose){
		message("Processing ", length(filenames), " files.")
		if (getRversion() > '2.7.0') pb <- txtProgressBar(min=0, max=length(filenames), style=3)
	}
        loader("1m_reference_cn.rda", .crlmmPkgEnv, pkgname)
	reference <- getVarInEnv("reference")
	for(i in seq(along=filenames)){
		y <- as.matrix(read.celfile(filenames[i], intensity.means.only=TRUE)[["INTENSITY"]][["MEAN"]][fid])
		x <- log2(y[idx2])
		SKW[i] <- mean((x-mean(x))^3)/(sd(x)^3)
		rm(x)
		NP[, i] <- as.integer(normalize.quantiles.use.target(y, target=reference))
		if (verbose)
			if (getRversion() > '2.7.0') setTxtProgressBar(pb, i)
			else cat(".")
	}
	dimnames(NP) <- list(names(fid), sns)
	##dimnames(NP) <- list(map[, "man_fsetid"], sns)
	res3 <- list(NP=NP, SKW=SKW)
	return(res3)
}

getFlags <- function(phi.thr, envir){
	nuA <- get("nuA", envir)
	nuB <- get("nuB", envir)
	phiA <- get("phiA", envir)
	phiB <- get("phiB", envir)
	
	negativeNus <- nuA < 1 | nuB < 1
	negativePhis <- phiA < phi.thr | phiB < phi.thr
	negativeCoef <- negativeNus | negativePhis

	notfinitePhi <- !is.finite(phiA) | !is.finite(phiB)
	flags <- negativeCoef | notfinitePhi
	return(flags)
}

goodSnps <- function(phi.thr, envir, fewAA=20, fewBB=20){
	Ns <- get("Ns", envir)
	flags <- getFlags(phi.thr=phi.thr, envir)
	fewAA <- Ns[, , "AA"] < fewAA
	fewBB <- Ns[, , "BB"] < fewBB
	flagsA <- flags | fewAA
	flagsB <- flags | fewBB
	flags <- list(A=flagsA, B=flagsB)
	return(flags)
}

##Needs to allow for NULL NP
instantiateObjects <- function(calls, conf, NP, plate, envir,
			       chrom,
			       A, B,
			       gender, SNRmin=5, SNR,
                               pkgname,
			       locusSet){
	envir[["chrom"]] <- chrom
	A <- A[, SNR > SNRmin]
	B <- B[, SNR > SNRmin]
	calls <- calls[, SNR > SNRmin]
	conf <- conf[, SNR > SNRmin]
	if(!is.null(NP))
		NP <- NP[, SNR > SNRmin]
	plate <- plate[SNR > SNRmin]
	uplate <- unique(plate)
	SNR <- SNR[SNR > SNRmin]
	
	envir[["uplate"]] <- uplate
	envir[["plate"]] <- plate	
	envir[["NP"]] <- NP
	envir[["A"]] <- A
	envir[["B"]] <- B
	envir[["calls"]] <- calls
	envir[["conf"]] <- conf
	snps <- rownames(calls)
	if(!is.null(NP)){	
		cnvs <- rownames(NP)
	} else cnvs <- NULL
	sns <- basename(colnames(calls))
	if(!is.null(NP))
		stopifnot(identical(colnames(calls), colnames(NP)))
	envir[["sns"]] <- sns
	envir[["snps"]] <- snps
	envir[["cnvs"]] <- cnvs

	if(chrom == 23){
		if(is.null(gender)){
			message("Estimating gender")
			XMedian <- apply(log2(A[, , drop=FALSE]) + log2(B[, , drop=FALSE]), 2, median)/2
			gender <- kmeans(XMedian, c(min(XMedian[SNR>SNRmin]), max(XMedian[SNR>SNRmin])))[["cluster"]]			
			gender[gender==2] <- "female"
			gender[gender=="1"] <- "male"
			envir[["gender"]] <- gender
		} else envir[["gender"]] <- gender
		phiAx <- matrix(NA, nrow(calls), length(uplate))
		envir[["phiAx"]] <- phiAx
		envir[["phiBx"]] <- phiAx
	}
	CA <- CB <- matrix(NA, nrow(calls), ncol(calls))
	envir[["CA"]] <- CA
	envir[["CB"]] <- CB
	
	Ns <- array(NA, dim=c(nrow(calls), length(uplate), 5))
	dimnames(Ns)[[3]] <- c("A", "B", "AA", "AB", "BB")

	envir[["Ns"]] <- envir[["muB"]] <- envir[["muA"]] <- Ns
	envir[["vB"]] <- envir[["vA"]] <- Ns

	if(!is.null(NP)){
		CT.sds <- CT <- matrix(NA, nrow(NP), length(sns))
		nuT <- matrix(NA, nrow(NP), length(uplate))
		phiT <- nuT
		normalNP <- matrix(TRUE, nrow(NP), ncol(NP))
		sig2T <- matrix(NA, nrow(NP), length(uplate))		
	} else{
		sig2T <- normalNP <- nuT <- phiT <- CT.sds <- CT <- NULL
	}
	envir[["CT"]] <- CT
	envir[["CT.sds"]] <- CT.sds
	envir[["nuT"]] <- nuT
	envir[["phiT"]] <- phiT
	plates.completed <- rep(FALSE, length(uplate))
	envir[["plates.completed"]] <- plates.completed
	steps <- rep(FALSE, 4)
	names(steps) <- c("suffStats", "coef", "snp-cn", "np-cn")
	envir[["steps"]] <- steps
	snpflags <- matrix(FALSE, length(snps), length(uplate))
	npflags <- matrix(FALSE, length(cnvs), length(uplate))
	envir[["snpflags"]] <- snpflags
	envir[["npflags"]] <- npflags
	tau2A <- matrix(NA, nrow(calls), length(uplate))
	envir[["tau2A"]] <- tau2A
	envir[["tau2B"]] <- tau2A
	envir[["sig2A"]] <- tau2A
	envir[["sig2B"]] <- tau2A

	envir[["sig2T"]] <- sig2T
	envir[["corr"]] <- tau2A
	envir[["corrA.BB"]] <- tau2A
	envir[["corrB.AA"]] <- tau2A
	envir[["nuB"]] <- envir[["nuA"]] <- tau2A
	envir[["phiB"]] <- envir[["phiA"]] <- tau2A
	envir[["nuB.se"]] <- envir[["nuA.se"]] <- tau2A
	envir[["phiB.se"]] <- envir[["phiA.se"]] <- tau2A
	normal <- matrix(TRUE, nrow(A), ncol(A))
	envir[["normal"]] <- normal
	envir[["normalNP"]] <- normalNP
}

computeCopynumber <- function(object,
			      CHR,
			      bias.adj=FALSE,
			      batch,
			      SNRmin=5,
			      cdfName="genomewidesnp6", ...){
	if(ncol(object) < 10)
		stop("Must have at least 10 samples in each batch to estimate model parameters....preferably closer to 90 samples per batch")

	##require(oligoClasses)
	if(missing(CHR)) stop("Must specify CHR")
	if(CHR == 24) stop("Nothing available yet for chromosome Y")
	if(missing(batch)) stop("Must specify batch")
	if(length(batch) != ncol(object[[1]])) stop("Batch must be the same length as the number of samples")
	##the AB intensities
	Nset <- object[[1]]
	##indices of polymorphic loci
	index <- snpIndex(Nset)
	ABset <- Nset[index, ]
	##indices of nonpolymorphic loci
	NPset <- Nset[-index, ]
	##genotypes/confidences
	snpset <- object[[2]][index,]
	##previous version of compute copy number
	envir <- new.env()
	message("Fitting model for copy number estimation...")

	.computeCopynumber(chrom=CHR,
			   A=A(ABset),
			   B=B(ABset),
			   calls=calls(snpset),
			   conf=confs(snpset),
			   NP=A(NPset),
			   plate=batch,
			   envir=envir,
			   SNR=ABset$SNR,
			   bias.adj=FALSE,
			   SNRmin=SNRmin,
			   cdfName=cdfName,
			   ...)

	if(bias.adj){
		message("Running bias adjustment...")
		.computeCopynumber(chrom=CHR,
				   A=A(ABset),
				   B=B(ABset),
				   calls=calls(snpset),
				   conf=confs(snpset),
				   NP=A(NPset),
				   plate=batch,
				   envir=envir,
				   SNR=ABset$SNR,
				   bias.adj=TRUE,
				   SNRmin=SNRmin,
				   cdfName=cdfName,
				   ...)
	}
	message("Organizing results...")			
	locusSet <- list2locusSet(envir, ABset=ABset, NPset=NPset, CHR=CHR, cdfName=cdfName)
	if(anyDuplicated(position(locusSet))){
		message("duplicated physical positions removed from CopyNumberSet object")
		locusSet <- locusSet[!duplicated(position(locusSet)), ]
	}
	message("Thresholding model parameters")
	locusSet <- thresholdModelParams(locusSet)
	object[[3]] <- locusSet
	message("harmonizing dimnames of the elements in CrlmmSetList...")
	object <- .harmonizeDimnames(object)
	message("Reordering features by chromosome and physical position...")
	object <- object[order(chromosome(object), position(object)), ]	
	return(object)
}

cnIllumina <- function(object,
		       CHR,
		       bias.adj=FALSE,
		       batch,
		       SNRmin=5,
		       cdfName="genomewidesnp6", ...){
	if(missing(CHR)) stop("Must specify CHR")
	if(missing(batch)) stop("Must specify batch")
	if(length(batch) != ncol(object[[1]])) stop("Batch must be the same length as the number of samples")
	ABset <- object[[1]]
	snpset <- object[[2]]
	envir <- new.env()
	NP <- NULL
	message("Fitting model for copy number estimation...")
	.computeCopynumber(chrom=CHR,
			   A=A(ABset),
			   B=B(ABset),
			   calls=calls(snpset),
			   conf=confs(snpset),
			   NP=NULL,
			   plate=batch,
			   envir=envir,
			   SNR=ABset$SNR,
			   bias.adj=FALSE,
			   SNRmin=SNRmin,
			   cdfName=cdfName,
			   ...)

	if(bias.adj){
		message("Running bias adjustment...")
		.computeCopynumber(chrom=CHR,
				   A=A(ABset),
				   B=B(ABset),
				   calls=calls(snpset),
				   conf=confs(snpset),
				   NP=NULL,
				   plate=batch,
				   envir=envir,
				   SNR=ABset$SNR,
				   bias.adj=TRUE,
				   SNRmin=SNRmin,
				   cdfName=cdfName,
				   ...)
	}
	message("Organizing results...")			
	locusSet <- list2locusSet(envir, ABset=ABset, NPset=NULL, CHR=CHR, cdfName=cdfName)
	return(locusSet)
}

##getCopyNumberEnvironment <- function(crlmmSetList, cnSet){
##	envir <- new.env()
##	envir[["A"]] <- A(crlmmSetList)[snpIndex(crlmmSetList), ]
##	envir[["B"]] <- A(crlmmSetList)[snpIndex(crlmmSetList), ]
##	envir[["CA"]] <- CA(cnSet)[snpIndex(cnSet), ]/100
##	envir[["CB"]] <- CB(cnSet)[snpIndex(cnSet), ]/100		
##	envir[["NP"]] <- A(crlmmSetList)[cnIndex(crlmmSetList), ]
##	envir[["calls"]] <- calls(crlmmSetList)[snpIndex(crlmmSetList), ]
##	envir[["chrom"]] <- unique(chromosome(cnSet))
##	envir[["cnvs"]] <- cnNames(crlmmSetList)
##	envir[["conf"]] <- conf(crlmmSetList)
##	envir[["corr"]] <- fData(cnSet)$corr
##	envir[["corrA.BB"]] <- fData(cnSet)$corrA.BB
##	envir[["corrB.AA"]] <- fData(cnSet)$corrB.AA
##	envir[["CT"]] <- CA(cnSet)[cnIndex(cnSet), ]/100
##	##envir[["CT.sds"]] <- 
##	##envir[["GT.A"]]
##	##envir[["GT.B"]]
##	##envir[["index"]]
##	##envir[["muA"]]
##	##envir[["muB"]]
##	##envir[["normal"]]
##	##envir[["normalNP"]]
##	##envir[["npflags"]]
##	##envir[["Ns"]]
##	nuA <- fData(cnSet)[, grep("nuA", fvarLabels(cnSet))]
##	nuB <- fData(cnSet)[, grep("nuB", fvarLabels(cnSet))]
##	phiA <- fData(cnSet)[, grep("phiA", fvarLabels(cnSet))]
##	phiB <- fData(cnSet)[, grep("phiB", fvarLabels(cnSet))]
##	envir[["nuA"]] <- nuA[snpIndex(cnSet), ]
##	envir[["nuB"]] <- nuB[snpIndex(cnSet), ]
##	envir[["phiA"]] <- phiA[snpIndex(cnSet), ]
##	envir[["phiB"]] <- phiB[snpIndex(cnSet), ]
##	envir[["nuT"]] <- nuA[cnIndex(cnSet), ]
##	envir[["phiT"]] <- phiA[cnIndex(cnSet), ]
##	envir[["plate"]] <- cnSet$batch
##	
##}

updateNuPhi <- function(crlmmSetList, cnSet){
	##TODO: remove the use of environments.
	##repopulate the environment
	crlmmSetList <- crlmmSetList[, match(sampleNames(cnSet), sampleNames(crlmmSetList))]
	crlmmSetList <- crlmmSetList[match(featureNames(cnSet), featureNames(crlmmSetList)), ]
	##envir <- getCopyNumberEnvironment(crlmmSetList, cnSet)
	Nset <- crlmmSetList[[1]]
	##indices of polymorphic loci
	index <- snpIndex(Nset)
	ABset <- Nset[index, ]
	##indices of nonpolymorphic loci
	NPset <- Nset[-index, ]
	##genotypes/confidences
	snpset <- crlmmSetList[[2]][index,]
	##previous version of compute copy number
	envir <- new.env()	
	message("Running bias adjustment...")
##	.computeCopynumber(chrom=CHR,
##			   A=A(ABset),
##			   B=B(ABset),
##			   calls=calls(snpset),
##			   conf=confs(snpset),
##			   NP=A(NPset),
##			   plate=batch,
##			   envir=envir,
##			   SNR=ABset$SNR,
##			   bias.adj=TRUE,
##			   SNRmin=SNRmin,
##			   ...)	
}

list2locusSet <- function(envir, ABset, NPset, CHR, cdfName="genomewidesnp6"){
	if(missing(CHR)) stop("Must specify chromosome")
	pkgname <- paste(cdfName, "Crlmm", sep="")	
	path <- system.file("extdata", package=pkgname)
	loader("cnProbes.rda", pkgname=pkgname, envir=.crlmmPkgEnv)
	cnProbes <- get("cnProbes", envir=.crlmmPkgEnv)
	loader("snpProbes.rda", pkgname=pkgname, envir=.crlmmPkgEnv)
	snpProbes <- get("snpProbes", envir=.crlmmPkgEnv)	
	##require(oligoClasses) || stop("oligoClasses package not available")
	if(length(ls(envir)) == 0) stop("environment empty")
	batch <- envir[["plate"]]
	##SNP copy number	
	CA <- envir[["CA"]]
	dimnames(CA) <- list(envir[["snps"]], envir[["sns"]])	
	CB <- envir[["CB"]]
	dimnames(CB) <- dimnames(CA)

	##NP copy number
	if(!is.null(NPset)){
		CT <- envir[["CT"]]
		rownames(CT) <- envir[["cnvs"]]
		colnames(CT) <- envir[["sns"]]
		sig2T <- envir[["sig2T"]]
		rownames(sig2T) <- rownames(CT)
		nuT <- envir[["nuT"]]
		colnames(nuT) <- paste("nuT", unique(batch), sep="_")
		phiT <- envir[["phiT"]]
		colnames(phiT) <- paste("phiT", unique(batch), sep="_")
		naMatrix <- matrix(NA, nrow(CT), ncol(CT))
		dimnames(naMatrix) <- dimnames(CT)
	} else{
		sig2T <- nuT <- phiT <- naMatrix <- CT <- NULL
	}
	CA <- rbind(CA, CT)
	CB <- rbind(CB, naMatrix)	

	##SNP parameters
	tau2A <- envir[["tau2A"]]
	colnames(tau2A) <- paste("tau2A", unique(batch), sep="_")
	tau2B <- envir[["tau2B"]]
	colnames(tau2B) <- paste("tau2B", unique(batch), sep="_")
	sig2A <- envir[["sig2A"]]
	colnames(sig2A) <- paste("sig2A", unique(batch), sep="_")
	sig2B <- envir[["sig2B"]]
	colnames(sig2B) <- paste("sig2B", unique(batch), sep="_")
	nuA <- envir[["nuA"]]
	colnames(nuA) <- paste("nuA", unique(batch), sep="_")
	nuB <- envir[["nuB"]]
	colnames(nuB) <- paste("nuB", unique(batch), sep="_")
	phiA <- envir[["phiA"]]
	colnames(phiA) <- paste("phiA", unique(batch), sep="_")
	phiB <- envir[["phiB"]]
	colnames(phiB) <- paste("phiB", unique(batch), sep="_")
	corr <- envir[["corr"]]
	colnames(corr) <- paste("corr", unique(batch), sep="_")
	corrA.BB <- envir[["corrA.BB"]]
	colnames(corrA.BB) <- paste("corrA.BB", unique(batch), sep="_")
	corrB.AA <- envir[["corrB.AA"]]
	colnames(corrB.AA) <- paste("corrB.AA", unique(batch), sep="_")


	##Combine SNP and NP parameters
	if(!is.null(NPset)){
		naMatrixParams <- matrix(NA, nrow(CT), length(unique(batch)))
		dimnames(naMatrixParams) <- list(rownames(CT), unique(batch))
	} else{
		naMatrixParams <- NULL
	}
	tau2A <- rbind(tau2A, naMatrixParams)
	tau2B <- rbind(tau2B, naMatrixParams)
	sig2A <- rbind(sig2A, sig2T)
	sig2B <- rbind(sig2B, naMatrixParams)
	corr <- rbind(corr, naMatrixParams)
	corrA.BB <- rbind(corrA.BB, naMatrixParams)
	corrB.AA <- rbind(corrB.AA, naMatrixParams)
	nuA <- rbind(nuA, nuT)
	phiA <- rbind(phiA, phiT)
	nuB <- rbind(nuB, naMatrixParams)
	phiB <- rbind(phiB, naMatrixParams)
	rownames(tau2A) <- rownames(tau2B) <- rownames(sig2A) <- rownames(sig2B) <- rownames(CA)
	rownames(corr) <- rownames(corrA.BB) <- rownames(corrB.AA) <- rownames(CA)
	rownames(nuA) <- rownames(phiA) <- rownames(nuB) <- rownames(phiB) <- rownames(CA)	
	##phenodata
	phenodata <- phenoData(ABset)
	phenodata <- phenodata[match(envir[["sns"]], sampleNames(phenodata)), ]
	if(!("batch" %in% varLabels(phenodata))) phenodata$batch <- envir[["plate"]]

	##Feature Data
	position.snp <- snpProbes[match(envir[["snps"]], rownames(snpProbes)), "position"]
	names(position.snp) <- envir[["snps"]]
	if(!is.null(NPset)){
		position.np <- cnProbes[match(envir[["cnvs"]], rownames(cnProbes)), "position"]
		names(position.np) <- envir[["cnvs"]]
	} else position.np <- NULL
	position <- c(position.snp, position.np)
	if(!(identical(names(position), rownames(CA)))){
		position <- position[match(rownames(CA), names(position))]
	}
	if(sum(duplicated(names(position))) > 0){
		warning("Removing rows with NA identifiers...")
		##RS: fix this
		I <- which(!is.na(names(position)))
	}  else I <- seq(along=names(position))
	fd <- data.frame(cbind(CHR,
			       position[I],
			       tau2A[I,, drop=FALSE],
			       tau2B[I,, drop=FALSE],
			       sig2A[I,, drop=FALSE],
			       sig2B[I,, drop=FALSE],
			       nuA[I,, drop=FALSE],
			       nuB[I,, drop=FALSE],
			       phiA[I,, drop=FALSE],
			       phiB[I,, drop=FALSE],
			       corr[I,, drop=FALSE],
			       corrA.BB[I,, drop=FALSE],
			       corrB.AA[I,, drop=FALSE]))
	colnames(fd)[1:2] <- c("chromosome", "position")
	rownames(fd) <- rownames(CA)[I]
	fD <- new("AnnotatedDataFrame",
		  data=fd,
		  varMetadata=data.frame(labelDescription=colnames(fd)))	
	assayData <- assayDataNew("lockedEnvironment",
				  CA=CA[I, ],
				  CB=CB[I, ])
	cnset <- new("CopyNumberSet",
		      assayData=assayData,
		      featureData=fD,
		      phenoData=phenodata,
		      annotation="genomewidesnp6")
	cnset <- thresholdCopyNumberSet(cnset)
	return(cnset)
}

thresholdCopyNumberSet <- function(object){
	ca <- CA(object)
	cb <- CB(object)
	ca[ca < 0.05] <- 0.05
	ca[ca > 5] <- 5
	cb[cb < 0.05] <- 0.05
	cb[cb > 5] <- 5
	ca <- matrix(as.integer(ca*100), nrow(ca), ncol(ca))
	cb <- matrix(as.integer(cb*100), nrow(cb), ncol(cb))
	rownames(ca) <- rownames(cb) <- featureNames(object)
	colnames(ca) <- colnames(cb) <- sampleNames(object)
	CA(object) <- ca
	CB(object) <- cb
	return(object)
}


.computeCopynumber <- function(chrom,
			       A,
			       B,
			       calls,
			       conf,
			       NP,
			       plate,
			       MIN.OBS=5,
			       envir,
			       P,
			       DF.PRIOR=50,
			       CONF.THR=0.99,
			       bias.adj=FALSE,
			       priorProb,
			       gender=NULL,
			       SNR,
			       SNRmin, seed=123,
			       cdfName,
			       verbose=TRUE, ...){
	if(missing(cdfName)) stop("cdfName must be provided")
	require(paste(cdfName, "Crlmm", sep=""), character.only=TRUE) || stop(paste("cdf ", cdfName, "Crlmm", " not available.", sep=""))
	if(!missing(plate)){
		if(length(plate) != ncol(A))
			stop("plate must the same length as the number of columns of A")
	}
	message("Using ", DF.PRIOR, " df for inverse chi squares.")		
	set.seed(seed)
	if(length(ls(envir)) == 0) {
		instantiateObjects(calls=calls,
				   conf=conf,
				   NP=NP,
				   plate=plate,
				   envir=envir,
				   chrom=chrom,
				   A=A, B=B,
				   gender=gender,
				   SNR=SNR,
				   SNRmin=SNRmin,
				   pkgname=cdfName)
	}
	plate <- envir[["plate"]]
	uplate <- envir[["uplate"]]	
	calls <- envir[["calls"]]
	A <- envir[["A"]]
	B <- envir[["B"]]
	conf <- envir[["conf"]]
	NP <- envir[["NP"]]
	if(bias.adj){
		##assign uniform priors for total copy number states
		if(missing(priorProb)) priorProb <- rep(1/4, 4)
		envir[["steps"]] <- rep(FALSE, 4)
	}
	##will be updating these objects
	if(verbose) message("Sufficient statistics")
	if(missing(P)) P <- seq(along=uplate)
	steps <- envir[["steps"]]
	if(!steps[1]){
		for(p in P){
			cat(".")
			if(sum(plate == uplate[p]) < 10) next()
			J <- plate==uplate[p]
			if(!is.null(NP)){
				npMatrix <- NP[, J]
			} else npMatrix <- NULL
			oneBatch(plateIndex=p,
				 A=A[, J],
				 B=B[, J],
				 calls=calls[, J],
				 conf=conf[, J],
				 gender=NULL,
				 npMatrix,
				 plate[J],
				 MIN.OBS=1,
				 envir=envir,
				 DF.PRIOR=DF.PRIOR,
				 CONF.THR=CONF.THR,
				 bias.adj=bias.adj,
				 priorProb=priorProb,...)
		}
		steps[1] <- TRUE
		envir[["steps"]] <- steps
	}
	if(!steps[2]){
		message("\nEstimating coefficients")	
		for(p in P){
			cat(".")
			coefs(plateIndex=p, conf=conf[, plate==uplate[p]],
			      envir=envir, CONF.THR=CONF.THR, MIN.OBS=MIN.OBS)
		}
		steps[2] <- TRUE
		envir[["steps"]] <- steps		
	}
	if(!steps[3]){
		message("\nAllele specific copy number")	
		for(p in P){
			cat(".")
			polymorphic(plateIndex=p,
				    A=A[, plate==uplate[p]],
				    B=B[, plate==uplate[p]],
				    envir=envir)
		}
		steps[3] <- TRUE
		envir[["steps"]] <- steps						
	}
	if(!steps[4]){
		if(!is.null(NP)){
			message("\nCopy number for nonpolymorphic probes...")	
			for(p in P){
				cat(".")
				nonpolymorphic(plateIndex=p,
					       NP=NP[, plate==uplate[p]],
					       envir=envir)
			}
		}
		steps[4] <- TRUE
		envir[["step"]] <- steps
	}
}

nonpolymorphic <- function(plateIndex, NP, envir, CONF.THR=0.99, DF.PRIOR=50, pkgname="genomewidesnp6Crlmm"){
	p <- plateIndex
	CHR <- envir[["chrom"]]
	plate <- envir[["plate"]]
	uplate <- envir[["uplate"]]
	plates.completed <- envir[["plates.completed"]]
	if(!plates.completed[p]) return()
	snpflags <- goodSnps(phi.thr=2^6, envir=envir, fewAA=10, fewBB=10)
	flagsA <- snpflags$A[, p]
	flagsB <- snpflags$B[, p]
	if(all(flagsA) | all(flagsB)) stop("all snps are flagged")
	nuA <- envir[["nuA"]][, p]
	nuB <- envir[["nuB"]][, p]
	phiA <- envir[["phiA"]][, p]
	phiB <- envir[["phiB"]][, p]
	uplate <- envir[["uplate"]]
	sns <- envir[["sns"]]
	muA <- envir[["muA"]][, p, ]
	muB <- envir[["muB"]][, p, ]
	nuT <- envir[["nuT"]]
	phiT <- envir[["phiT"]]
	sig2T <- envir[["sig2T"]]
	A <- envir[["A"]][, plate==uplate[p]]
	B <- envir[["B"]][, plate==uplate[p]]
	CA <- envir[["A"]][, plate==uplate[p]]
	CB <- envir[["B"]][, plate==uplate[p]]
	if(CHR == 23){
		phiAx <- envir[["phiAx"]][, p]
		phiBx <- envir[["phiBx"]][, p]
	}
	##---------------------------------------------------------------------------
	## Train on unflagged SNPs
	##---------------------------------------------------------------------------
	##Might be best to train using the X chromosome, since for the
	##X phi and nu have been adjusted for cross-hybridization
	plateInd <- plate == uplate[p]
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
	if(CHR == 23){
		normalNP <- envir[["normalNP"]]
		normalNP <- normalNP[, plate==uplate[p]]
		nuT <- envir[["nuT"]]
		phiT <- envir[["phiT"]]
		
		cnvs <- envir[["cnvs"]]
                loader("cnProbes.rda", pkgname=pkgname, envir=.crlmmPkgEnv)
                cnProbes <- get("cnProbes", envir=.crlmmPkgEnv)
		cnProbes <- cnProbes[match(cnvs, rownames(cnProbes)), ]

		##For build Hg18
		##http://genome.ucsc.edu/cgi-bin/hgGateway
		##pseudo-autosomal regions on X
		##chrX:1-2,709,520 and chrX:154584237-154913754, respectively
		##par:pseudo-autosomal regions
		pseudoAR <- cnProbes[, "position"] < 2709520 | (cnProbes[, "position"] > 154584237 & cnProbes[, "position"] < 154913754)
		gender <- envir[["gender"]]
		mu1 <- rowMedians(NP[, gender=="male"], na.rm=TRUE)
		mu2 <- rowMedians(NP[, gender=="female"], na.rm=TRUE)
		mus <- log2(cbind(mu1, mu2))
		X.men <- cbind(1, mus[, 1])
		X.fem <- cbind(1, mus[, 2])
		
		Yhat1 <- as.numeric(X.men %*% betahat)
		Yhat2 <- as.numeric(X.fem %*% betahat)
		phi1 <- 2^(Yhat1)
		phi2 <- 2^(Yhat2)
		nu1 <- 2^(mus[, 1]) - phi1
		nu2 <- 2^(mus[, 2]) - 2*phi2		
		nu1[pseudoAR] <- 2^(mus[pseudoAR, 1]) - 2*phi1[pseudoAR]
		CT1 <- 1/phi1*(NP[, gender=="male"]-nu1)
		CT2 <- 1/phi2*(NP[, gender=="female"]-nu2)
		CT1 <- matrix(as.integer(100*CT1), nrow(CT1), ncol(CT1))
		CT2 <- matrix(as.integer(100*CT2), nrow(CT2), ncol(CT2))
		CT <- envir[["CT"]]
		CT[, plate==uplate[p] & gender=="male"] <- CT1
		CT[, plate==uplate[p] & gender=="female"] <- CT2
		envir[["CT"]] <- CT

		##only using females to compute the variance
		normalNP[, gender=="male"] <- NA		
		sig2T[, p] <- rowMAD(log2(NP*normalNP), na.rm=TRUE)^2
		nuT[, p] <- nu2
		phiT[, p] <- phi2
		envir[["sig2T"]] <- sig2T
		envir[["CT"]] <- CT
		envir[["phiT"]] <- nuT
		envir[["nuT"]] <- phiT
	} else {
		normalNP <- envir[["normalNP"]]
		normalNP <- normalNP[, plate==uplate[p]]
		mus <- rowMedians(NP * normalNP, na.rm=TRUE)
		crosshyb <- median(muA) - median(mus)
		X <- cbind(1, log2(mus+crosshyb))
		##X <- cbind(1, log2(mus))
		logPhiT <- X %*% betahat
		phiT[, p] <- 2^(logPhiT)
		nuT[, p] <- mus - 2*phiT[, p]
		T <- 1/phiT[, p]*(NP - nuT[, p])
		CT <- envir[["CT"]]
		CT[, plate==uplate[p]] <- matrix(as.integer(100*T), nrow(T), ncol(T))

		##Variance for prediction region
		##sig2T[, plate==uplate[[p]]] <- rowMAD(log2(NP*normalNP), na.rm=TRUE)^2
		sig2T[, p] <- rowMAD(log2(NP*normalNP), na.rm=TRUE)^2	
		envir[["sig2T"]] <- sig2T
		envir[["CT"]] <- CT
		envir[["phiT"]] <- phiT
		envir[["nuT"]] <- nuT
	}
	##---------------------------------------------------------------------------
	## For NA SNPs, treat as nonpolymorphic
	##---------------------------------------------------------------------------
}

##sufficient statistics on the intensity scale
withinGenotypeMoments <- function(p, A, B, calls, conf, CONF.THR, DF.PRIOR, envir){
	CHR <- envir[["chrom"]]
	Ns <- envir[["Ns"]]
	muA <- envir[["muA"]]
	muB <- envir[["muB"]]
	vA <- envir[["vA"]]
	vB <- envir[["vB"]]
	plate <- envir[["plate"]]
	uplate <- envir[["uplate"]]
	normal <- envir[["normal"]][, plate==uplate[p]]
	G <- calls; rm(calls); gc()


	highConf <- 1-exp(-conf/1000)
	highConf <- highConf > CONF.THR
	if(CHR == 23){
		gender <- envir[["gender"]]
		IX <- matrix(gender, nrow(G), ncol(G), byrow=TRUE)
		IX <- IX == "female"
	} else IX <- matrix(TRUE, nrow(G), ncol(G))
	index <- GT.B <- GT.A <- vector("list", 3)
	names(index) <- names(GT.B) <- names(GT.A) <- c("AA", "AB", "BB")
	##--------------------------------------------------
	##within-genotype sufficient statistics
	##--------------------------------------------------
	GT.B <- GT.A <- list()
	for(j in 1:3){
		##GT <- G==j & highConf & IX & normal
		GT <- G==j & highConf & IX
		GT <- GT * normal
		Ns[, p, j+2] <- rowSums(GT, na.rm=TRUE)		
		GT[GT == FALSE] <- NA
		GT.A[[j]] <- GT*A
		GT.B[[j]] <- GT*B
		index[[j]] <- which(Ns[, p, j+2] > 0)
		muA[, p, j+2] <- rowMedians(GT.A[[j]], na.rm=TRUE)
		muB[, p, j+2] <- rowMedians(GT.B[[j]], na.rm=TRUE)
		vA[, p, j+2] <- rowMAD(GT.A[[j]], na.rm=TRUE)
		vB[, p, j+2] <- rowMAD(GT.B[[j]], na.rm=TRUE)

		##Shrink towards the typical variance
		DF <- Ns[, p, j+2]-1
		DF[DF < 1] <- 1
		v0A <- median(vA[, p, j+2], na.rm=TRUE)
		v0B <- median(vB[, p, j+2], na.rm=TRUE)
		if(v0A == 0) v0A <- NA
		if(v0B == 0) v0B <- NA
		vA[, p, j+2] <- (vA[, p, j+2]*DF + v0A*DF.PRIOR)/(DF.PRIOR+DF)
		vA[is.na(vA[, p, j+2]), p, j+2] <- v0A
		vB[, p, j+2] <- (vB[, p, j+2]*DF + v0B*DF.PRIOR)/(DF.PRIOR+DF)
		vB[is.na(vB[, p, j+2]), p, j+2] <- v0B
	}
	if(CHR == 23){
		k <- 1
		for(j in c(1,3)){
			GT <- G==j & highConf & !IX 
			Ns[, p, k] <- rowSums(GT)
			GT[GT == FALSE] <- NA
			muA[, p, k] <- rowMedians(GT*A, na.rm=TRUE)
			muB[, p, k] <- rowMedians(GT*B, na.rm=TRUE)
			vA[, p, k] <- rowMAD(GT*A, na.rm=TRUE)
			vB[, p, k] <- rowMAD(GT*B, na.rm=TRUE)
			
			DF <- Ns[, p, k]-1
			DF[DF < 1] <- 1
			v0A <- median(vA[, p, k], na.rm=TRUE)
			v0B <- median(vB[, p, k], na.rm=TRUE)
			vA[, p, k] <- (vA[, p, k]*DF + v0A*DF.PRIOR)/(DF.PRIOR+DF)
			vA[is.na(vA[, p, k]), p, k] <- v0A
			vB[, p, k] <- (vB[, p, k]*DF + v0B*DF.PRIOR)/(DF.PRIOR+DF)
			vB[is.na(vB[, p, k]), p, k] <- v0B			
			k <- k+1
		}
	}
	envir[["GT.A"]] <- GT.A
	envir[["GT.B"]] <- GT.B
	envir[["Ns"]] <- Ns
	envir[["index"]] <- index
	envir[["muA"]] <- muA
	envir[["muB"]] <- muB
	envir[["vA"]] <- vA
	envir[["vB"]] <- vB
}
	

oneBatch <- function(plateIndex,
		     A,
		     B,
		     calls,
		     conf,
		     gender,
		     NP,
		     plate,
		     MIN.OBS=1,
		     envir,
		     DF.PRIOR=50,
		     CONF.THR=0.99,
		     trim=0,
		     bias.adj=FALSE, priorProb, ...){
	p <- plateIndex
	CHR <- envir[["chrom"]]
	if(bias.adj){
		nuA <- envir[["nuA"]]
		if(all(is.na(nuA))){
			message("Background and signal coefficients have not yet been estimated -- can not do bias correction yet")
			stop("Must run computeCopynumber a second time with bias.adj=TRUE to do the adjustment")
		}
		message("running bias adjustment")		
		##adjustment for nonpolymorphic probes
		if(!is.null(NP))
			biasAdjNP(plateIndex=p, envir=envir, priorProb=priorProb)
		##adjustment for SNPs
		biasAdj(plateIndex=p, envir=envir, priorProb=priorProb)
		message("Recomputing location and scale parameters")		
	}
	withinGenotypeMoments(p=p,
			      A=A,
			      B=B,
			      calls=calls,
			      conf=conf,
			      CONF.THR=CONF.THR,
			      DF.PRIOR=DF.PRIOR,
			      envir=envir)
	GT.A <- envir[["GT.A"]]
	GT.B <- envir[["GT.B"]]
	index <- envir[["index"]]
	locationAndScale(p=p, GT.A=GT.A, GT.B=GT.B, index=index, envir=envir, DF.PRIOR=DF.PRIOR)
	muA <- envir[["muA"]]
	muB <- envir[["muB"]]
	Ns <- envir[["Ns"]]

	##---------------------------------------------------------------------------
	## Impute sufficient statistics for unobserved genotypes (plate-specific)
	##---------------------------------------------------------------------------
	index.AA <- which(Ns[, p, "AA"] >= 3)
	index.AB <- which(Ns[, p, "AB"] >= 3)
	index.BB <- which(Ns[, p, "BB"] >= 3)
	correct.orderA <- muA[, p, "AA"] > muA[, p, "BB"]
	correct.orderB <- muB[, p, "BB"] > muB[, p, "AA"]
	##For chr X, this will ignore the males 
	nobs <- rowSums(Ns[, p, 3:5] >= MIN.OBS, na.rm=TRUE) == 3
	index.complete <- which(correct.orderA & correct.orderB & nobs) ##be selective here
	size <- min(5000, length(index.complete))
	if(size == 5000) index.complete <- sample(index.complete, 5000)
	if(length(index.complete) < 200){
		warning("fewer than 200 snps pass criteria for predicting the sufficient statistics")
		stop()
	}
	index[[1]] <- which(Ns[, p, "AA"] == 0 & (Ns[, p, "AB"] >= MIN.OBS & Ns[, p, "BB"] >= MIN.OBS))
	index[[2]] <- which(Ns[, p, "AB"] == 0 & (Ns[, p, "AA"] >= MIN.OBS & Ns[, p, "BB"] >= MIN.OBS))
	index[[3]] <- which(Ns[, p, "BB"] == 0 & (Ns[, p, "AB"] >= MIN.OBS & Ns[, p, "AA"] >= MIN.OBS))
	mnA <- muA[, p, 3:5]
	mnB <- muB[, p, 3:5]
	for(j in 1:3){
		if(length(index[[j]]) == 0) next()
		X <- cbind(1, mnA[index.complete,  -j], mnB[index.complete,  -j])
		Y <- cbind(mnA[index.complete, j], mnB[index.complete,  j])
		betahat <- solve(crossprod(X), crossprod(X,Y))
		X <- cbind(1, mnA[index[[j]],  -j],  mnB[index[[j]],  -j])
		mus <- X %*% betahat
		muA[index[[j]], p, j+2] <- mus[, 1]
		muB[index[[j]], p, j+2] <- mus[, 2]
	}
	nobsA <- Ns[, p, "A"] > 10
	nobsB <- Ns[, p, "B"] > 10
	notMissing <- !(is.na(muA[, p, "A"]) | is.na(muA[, p, "B"]) | is.na(muB[, p, "A"]) | is.na(muB[, p, "B"]))
	complete <- list()
	complete[[1]] <- which(correct.orderA & correct.orderB & nobsA & notMissing) ##be selective here
	complete[[2]] <- which(correct.orderA & correct.orderB & nobsB & notMissing) ##be selective here	
	size <- min(5000, length(complete[[1]]))
	if(size == 5000) complete <- lapply(complete, function(x) sample(x, size))
	if(CHR == 23){
		index <- list()
		index[[1]] <- which(Ns[, p, "A"] == 0)
		index[[2]] <- which(Ns[, p, "B"] == 0)
		cols <- 2:1
		for(j in 1:2){
			if(length(index[[j]]) == 0) next()
			X <- cbind(1, muA[complete[[j]], p, cols[j]], muB[complete[[j]], p, cols[j]])
			Y <- cbind(muA[complete[[j]], p, j], muB[complete[[j]], p, j])
			betahat <- solve(crossprod(X), crossprod(X,Y))
			X <- cbind(1, muA[index[[j]], p, cols[j]],  muB[index[[j]], p, cols[j]])
			mus <- X %*% betahat
			muA[index[[j]], p, j] <- mus[, 1]
			muB[index[[j]], p, j] <- mus[, 2]
		}
	}
	##missing two genotypes
	noAA <- Ns[, p, "AA"] < MIN.OBS
	noAB <- Ns[, p, "AB"] < MIN.OBS
	noBB <- Ns[, p, "BB"] < MIN.OBS
	index[[1]] <- noAA & noAB
	index[[2]] <- noBB & noAB
	index[[3]] <- noAA & noBB
	snpflags <- envir[["snpflags"]]
	snpflags[, p] <- index[[1]] | index[[2]] | index[[3]]

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
		muA[index[[j]], p, -c(1, 2, k+2)] <- mus[, 1:2]
		muB[index[[j]], p, -c(1, 2, k+2)] <- mus[, 3:4]
	}
	negA <- rowSums(muA[, p, ] < 0) > 0
	negB <- rowSums(muB[, p, ] < 0) > 0	
	snpflags[, p] <- snpflags[, p] | negA | negB | rowSums(is.na(muA[, p, 3:5]), na.rm=TRUE) > 0
	envir[["snpflags"]] <- snpflags
	dn.Ns <- dimnames(Ns)
	Ns <- array(as.integer(Ns), dim=dim(Ns))
	dimnames(Ns)[[3]] <- dn.Ns[[3]]
	envir[["Ns"]] <- Ns
	envir[["muA"]] <- muA
	envir[["muB"]] <- muB
	plates.completed <- envir[["plates.completed"]]
	plates.completed[p] <- TRUE
	envir[["plates.completed"]] <- plates.completed
}

##Estimate tau2, sigma2, and correlation
locationAndScale <- function(p, GT.A, GT.B, index, envir, DF.PRIOR){
	tau2A <- envir[["tau2A"]]
	tau2B <- envir[["tau2B"]]
	sig2A <- envir[["sig2A"]]
	sig2B <- envir[["sig2B"]]

	corr <- envir[["corr"]]
	corrA.BB <- envir[["corrA.BB"]]
	corrB.AA <- envir[["corrB.AA"]]
	Ns <- get("Ns", envir)
	
	index.AA <- index[[1]]
	index.AB <- index[[2]]
	index.BB <- index[[3]]
	rm(index); gc()

	AA.A <- GT.A[[1]]
	AB.A <- GT.A[[2]]
	BB.A <- GT.A[[3]]
	
	AA.B <- GT.B[[1]]
	AB.B <- GT.B[[2]]
	BB.B <- GT.B[[3]]	
	x <- BB.A[index.BB, ]
	tau2A[index.BB, p] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2
	DF <- Ns[, p, "BB"]-1
	DF[DF < 1] <- 1
	med <- median(tau2A[, p], na.rm=TRUE)
	tau2A[, p] <- (tau2A[, p] * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
	tau2A[is.na(tau2A[, p]), p] <- med
	
	x <- BB.B[index.BB, ]
	sig2B[index.BB, p] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2	
	med <- median(sig2B[, p], na.rm=TRUE)
	sig2B[, p] <- (sig2B[, p] * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
	sig2B[is.na(sig2B[, p]), p] <- med
	
	x <- AA.B[index.AA, ]
	tau2B[index.AA, p] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2		
	DF <- Ns[, p, "AA"]-1
	DF[DF < 1] <- 1
	med <- median(tau2B[, p], na.rm=TRUE)
	tau2B[, p] <- (tau2B[, p] * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
	tau2B[is.na(tau2B[, p]), p] <- med
	
	x <- AA.A[index.AA, ]
	sig2A[index.AA, p] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2##var(log(IA)|AA)	
	med <- median(sig2A[, p], na.rm=TRUE)
	sig2A[, p] <- (sig2A[, p]*DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
	sig2A[is.na(sig2A[, p]), p] <- med	

	if(length(index.AB) > 0){ ##all homozygous is possible
		x <- AB.A[index.AB, ]
		y <- AB.B[index.AB, ]
		corr[index.AB, p] <- rowCors(x, y, na.rm=TRUE)
		corr[corr < 0] <- 0
		DF <- Ns[, p, "AB"]-1
		DF[DF<1] <- 1
		med <- median(corr[, p], na.rm=TRUE)
		corr[, p] <- (corr[, p]*DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
		corr[is.na(corr[, p]), p] <- med
	}
	backgroundB <- AA.B[index.AA, ]
	signalA <- AA.A[index.AA, ]
	corrB.AA[index.AA, p] <- rowCors(backgroundB, signalA, na.rm=TRUE)
	DF <- Ns[, p, "AA"]-1
	DF[DF < 1] <- 1
	med <- median(corrB.AA[, p], na.rm=TRUE)
	corrB.AA[, p] <- (corrB.AA[, p]*DF + med*DF.PRIOR)/(DF.PRIOR + DF)
	corrB.AA[is.na(corrB.AA[, p]), p] <- med

	backgroundA <- BB.A[index.BB, ]
	signalB <- BB.B[index.BB, ]
	corrA.BB[index.BB, p] <- rowCors(backgroundA, signalB, na.rm=TRUE)
	DF <- Ns[, p, "BB"]-1
	DF[DF < 1] <- 1
	med <- median(corrA.BB[, p], na.rm=TRUE)
	corrA.BB[, p] <- (corrA.BB[, p]*DF + med*DF.PRIOR)/(DF.PRIOR + DF)
	corrA.BB[is.na(corrA.BB[, p]), p] <- med

	envir[["tau2A"]] <- tau2A
	envir[["tau2B"]] <- tau2B
	envir[["sig2A"]] <- sig2A
	envir[["sig2B"]] <- sig2B
	envir[["corr"]] <- corr
	envir[["corrB.AA"]] <- corrB.AA
	envir[["corrA.BB"]] <- corrA.BB	
}

coefs <- function(plateIndex, conf, MIN.OBS=3, envir, CONF.THR=0.99){
	p <- plateIndex
	plates.completed <- envir[["plates.completed"]]
	if(!plates.completed[p]) return()
	CHR <- envir[["chrom"]]
	plate <- envir[["plate"]]
	muA <- envir[["muA"]]
	muB <- envir[["muB"]]
	vA <- envir[["vA"]]
	vB <- envir[["vB"]]
	Ns <- envir[["Ns"]]
	uplate <- envir[["uplate"]]
	if(CHR != 23){
		IA <- muA[, p, 3:5]
		IB <- muB[, p, 3:5]
		vA <- vA[, p, 3:5]
		vB <- vB[, p, 3:5]
		Np <- Ns[, p, 3:5]
	} else {
		NOHET <- is.na(median(vA[, p, "AB"], na.rm=TRUE))
		if(NOHET){
			IA <- muA[, p, -4]
			IB <- muB[, p, -4]
			vA <- vA[, p, -4]
			vB <- vB[, p, -4]
			Np <- Ns[, p, -4]
		} else{
			IA <- muA[, p, ]
			IB <- muB[, p, ]
			vA <- vA[, p, ]
			vB <- vB[, p, ]
			Np <- Ns[, p, ]
		}
	}
	##---------------------------------------------------------------------------
	## Estimate nu and phi
	##---------------------------------------------------------------------------
##	if(NOHET){
##		##only homozygous
##		Np <- Np[, -2]
##		Np[Np < 1] <- 1
##		IA <- IA[, c(1, 3)]
##		IB <- IB[, c(1, 3)]		
##		vA <- vA[, c(3,5)]
##		vB <- vB[, c(3,5)]
##	}else 	Np[Np < 1] <- 1
	Np[Np < 1] <- 1
	vA2 <- vA^2/Np
	vB2 <- vB^2/Np
	wA <- sqrt(1/vA2)
	wB <- sqrt(1/vB2)
	YA <- IA*wA
	YB <- IB*wB
	nuphiAllele(p=p, allele="A", Ystar=YA, W=wA, envir=envir)
	nuphiAllele(p=p, allele="B", Ystar=YB, W=wB, envir=envir)

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
}

polymorphic <- function(plateIndex, A, B, envir){
	p <- plateIndex
	plates.completed <- envir[["plates.completed"]]
	if(!plates.completed[p]) return()
	CHR <- envir[["chrom"]]
	plate <- envir[["plate"]]
	uplate <- envir[["uplate"]]
	vA <- envir[["vA"]]
	vB <- envir[["vB"]]
	nuA <- envir[["nuA"]][, p]
	nuB <- envir[["nuB"]][, p]
	nuA.se <- envir[["nuA.se"]]
	nuB.se <- envir[["nuB.se"]]
	phiA <- envir[["phiA"]][, p]
	phiB <- envir[["phiB"]][, p]
	phiA.se <- envir[["phiA.se"]]
	phiB.se <- envir[["phiB.se"]]
	Ns <- get("Ns", envir)
	CA <- get("CA", envir)
	CB <- get("CB", envir)
	NOHET <- mean(Ns[, p, "AB"], na.rm=TRUE) < 0.05
	##---------------------------------------------------------------------------
	## Estimate CA, CB
	##---------------------------------------------------------------------------
	if(CHR == 23){
		phiAx <- as.matrix(envir[["phiAx"]])
		phiBx <- as.matrix(envir[["phiBx"]])
		phiAx <- phiAx[, p]  ##nonspecific hybridization slope
		phiBx <- phiBx[, p]  ##nonspecific hybridization slope
		phistar <- phiBx/phiA  
		tmp <- (B-nuB - phistar*A + phistar*nuA)/phiB
		copyB <- tmp/(1-phistar*phiAx/phiB)
		copyA <- (A-nuA-phiAx*copyB)/phiA
		CB[, plate==uplate[p]] <- matrix(as.integer(100*copyB), nrow(copyB), ncol(copyB))
		CA[, plate==uplate[p]] <- matrix(as.integer(100*copyA), nrow(copyA), ncol(copyA))
	} else{
		CA[, plate==uplate[p]] <- matrix(as.integer(100*1/phiA*(A-nuA)), nrow(A), ncol(A))
		CB[, plate==uplate[p]] <- matrix(as.integer(100*1/phiB*(B-nuB)), nrow(A), ncol(A))
	}
	assign("CA", CA, envir)
	assign("CB", CB, envir)
}




biasAdj <- function(plateIndex, envir, priorProb, PROP=0.75){
	gender <- envir[["gender"]]
	CHR <- envir[["chrom"]]
	if(CHR == 23){
		phiAx <- envir[["phiAx"]]
		phiBx <- envir[["phiBx"]]
	}
	A <- envir[["A"]]
	B <- envir[["B"]]
	sig2A <- envir[["sig2A"]]
	sig2B <- envir[["sig2B"]]
	tau2A <- envir[["tau2A"]]
	tau2B <- envir[["tau2B"]]
	corrA.BB <- envir[["corrA.BB"]]
	corrB.AA <- envir[["corrB.AA"]]
	corr <- envir[["corr"]]
	nuA <- envir[["nuA"]]
	nuB <- envir[["nuB"]]
	phiA <- envir[["phiA"]]
	phiB <- envir[["phiB"]]
	normal <- envir[["normal"]]
	p <- plateIndex
	plate <- envir[["plate"]]
	if(missing(priorProb)) priorProb <- rep(1/4, 4) ##uniform
	emit <- array(NA, dim=c(nrow(A), ncol(A), 10))##SNPs x sample x 'truth'	
	lA <- log2(A)
	lB <- log2(B)	
	X <- cbind(lA, lB)
	counter <- 1##state counter								
	for(CT in 0:3){
		for(CA in 0:CT){
			cat(".")
			CB <- CT-CA
			A.scale <- sqrt(tau2A[, p]*(CA==0) + sig2A[, p]*(CA > 0))
			B.scale <- sqrt(tau2B[, p]*(CA==0) + sig2B[, p]*(CA > 0))
			if(CA == 0 & CB == 0) rho <- 0
			if(CA == 0 & CB > 0) rho <- corrA.BB[, p]
			if(CA > 0 & CB == 0) rho <- corrB.AA[, p]
			if(CA > 0 & CB > 0) rho <- corr[, p]
			if(CHR == 23){
				##means <- cbind(suppressWarnings(log2(nuA[, p]+CA*phiA[, p] + CB*phiAx[, p])), suppressWarnings(log2(nuB[, p]+CB*phiB[, p] + CA*phiBx[, p])))
				meanA <- suppressWarnings(log2(nuA[, p]+CA*phiA[, p] + CB*phiAx[, p]))
				meanB <- suppressWarnings(log2(nuB[, p]+CB*phiB[, p] + CA*phiBx[, p]))				
			} else{
				##means <- cbind(suppressWarnings(log2(nuA[, p]+CA*phiA[, p])), suppressWarnings(log2(nuB[, p]+CB*phiB[, p])))
				meanA <- suppressWarnings(log2(nuA[, p]+CA*phiA[, p]))
				meanB <- suppressWarnings(log2(nuB[, p]+CB*phiB[, p]))
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
	envir[["posteriorProb"]] <- list(hemDel=hemDel, norm=norm, amp=amp)
	posteriorProb <- array(NA, dim=c(nrow(A), ncol(A), 4))
	posteriorProb[, , 1] <- homDel
	posteriorProb[, , 2] <- hemDel
	posteriorProb[, , 3] <- norm
	posteriorProb[, , 4] <- amp
	mostLikelyState <- apply(posteriorProb, c(1, 2), function(x) order(x, decreasing=TRUE)[1])
	if(CHR == 23){
		##so state index 3 is the most likely state for men and women
		mostLikelyState[, gender=="male"] <- mostLikelyState[, gender=="male"] + 1
	}
	proportionSamplesAltered <- rowMeans(mostLikelyState != 3)
	ii <- proportionSamplesAltered < 0.8 & proportionSamplesAltered > 0.01
	
	##  only exclude observations from one tail, depending on
	##  whether more are up or down
	
	moreup <- rowSums(mostLikelyState > 3) > rowSums(mostLikelyState < 3) ##3 is normal
	NORM <- matrix(FALSE, nrow(A), ncol(A))
	NORM[proportionSamplesAltered > 0.8, ] <- FALSE
	ratioUp <- posteriorProb[, , 4]/posteriorProb[, , 3]
	NORM[ii & moreup, ] <- ratioUp[moreup & ii] < 1  ##normal more likely
	ratioDown <- posteriorProb[, , 2]/posteriorProb[, , 3]
	NORM[ii & !moreup, ] <- ratioDown[!moreup & ii] < 1  ##normal more likely
	normal <- NORM*normal
	
	flagAltered <- which(proportionSamplesAltered > 0.5)
	envir[["flagAltered"]] <- flagAltered
	envir[["normal"]] <- normal
}


posteriorNonpolymorphic <- function(plateIndex, envir, priorProb, cnStates=0:6){
	p <- plateIndex
	CHR <- envir[["chrom"]]
	if(missing(priorProb)) priorProb <- rep(1/length(cnStates), length(cnStates)) ##uniform	
	plate <- envir[["plate"]]
	uplate <- envir[["plate"]]
	NP <- envir[["NP"]][, plate==uplate[p]]
	nuT <- envir[["nuT"]][, p]
	phiT <- envir[["phiT"]][, p]
	sig2T <- envir[["sig2T"]][, p]
	##Assuming background variance for np probes is the same on the log-scale
	emit <- array(NA, dim=c(nrow(NP), ncol(NP), length(cnStates)))##SNPs x sample x 'truth'
	lT <- log2(NP)
	sds <- sqrt(sig2T)
	counter <- 1##state counter	
	for(CT in cnStates){
		cat(".")
		if(CHR == 23) browser()
		means <- suppressWarnings(log2(nuT + CT*phiT))
		emit[, , counter] <- dnorm(lT, mean=means, sd=sds)
		counter <- counter+1
	}
	for(j in seq(along=cnStates)){
		emit[, , j] <- priorProb[j]*emit[, , j]
	}
	homDel <- emit[, , 1]
	hemDel <- emit[, , 2]
	norm <- emit[, , 3]
	amp <- emit[, , 4]
	amp4 <- emit[, , 5]
	amp5 <- emit[, , 6]
	amp6 <- emit[, , 7]
	total <- homDel+hemDel+norm+amp+amp4+amp5+amp6
	weights <- array(NA, dim=c(nrow(NP), ncol(NP), length(cnStates)))
	weights[, , 1] <- homDel/total
	weights[, , 2] <- hemDel/total
	weights[, , 3] <- norm/total
	weights[, , 4] <- amp/total
	weights[, , 5] <- amp4/total
	weights[, , 6] <- amp5/total
	weights[, , 7] <- amp6/total
	##posterior mode
	posteriorMode <- apply(weights, c(1, 2), function(x) order(x, decreasing=TRUE)[1])
	posteriorMode <- posteriorMode-1
	##sns <- envir[["sns"]]
	##colnames(posteriorMode) <- sns
	##envir[["np.posteriorMode"]] <- posteriorMode
	##envir[["np.weights"]] <- weights
	posteriorMeans <- 0*homDel/total + 1*hemDel/total + 2*norm/total + 3*amp/total + 4*amp4/total + 5*amp5/total + 6*amp6/total
	##colnames(posteriorMeans) <- sns
	##envir[["np.posteriorMeans"]] <- posteriorMeans
	return(posteriorMode)
}

posteriorWrapper <- function(envir){
	snp.PM <- matrix(NA, length(envir[["snps"]]), length(envir[["sns"]]))
	np.PM <- matrix(NA, length(envir[["cnvs"]]), length(envir[["sns"]]))
	plate <- envir[["plate"]]
	uplate <- envir[["uplate"]]
	for(p in seq(along=uplate)){
		tmp <- expectedC(plateIndex=p, envir=envir)
		snp.PM[, plate==uplate[p]] <- tmp
		##snp.pm <- env[["posteriorMode"]]
		##trace(posteriorNonpolymorphic, browser)
		tmp <- posteriorNonpolymorphic(plateIndex=p, envir=envir)
		np.PM[, plate==uplate[p]] <- tmp##env[["np.posteriorMode"]]
		##pMode <- rbind(snp.pm, np.pm)
		##rownames(pMode) <- c(env[["snps"]], env[["cnvs"]])
		##dn <- dimnames(pMode)
		##pMode <- matrix(as.integer(pMode), nrow(pMode), ncol(pMode))
	}
	PM <- rbind(snp.PM, np.PM)
	PM <- matrix(as.integer(PM), nrow(PM), ncol(PM))
	dns <- list(c(envir[["snps"]], envir[["cnvs"]]), envir[["sns"]])
	dimnames(PM) <- dns
	return(PM)
}





##for polymorphic probes
expectedC <- function(plateIndex, envir, priorProb, cnStates=0:6){
	p <- plateIndex
	CHR <- envir[["chrom"]]
	if(missing(priorProb)) priorProb <- rep(1/length(cnStates), length(cnStates)) ##uniform	
	plate <- envir[["plate"]]
	uplate <- envir[["uplate"]]
	A <- envir[["A"]]
	B <- envir[["B"]]
	A <- A[, plate==uplate[p]]
	B <- B[, plate==uplate[p]]
	calls <- envir[["calls"]]	
	calls <- calls[, plate==unique(plate)[p]]
	probA <- sqrt(rowMeans(calls == 1, na.rm=TRUE))
	probB <- sqrt(rowMeans(calls == 3, na.rm=TRUE))
	sig2A <- envir[["sig2A"]]
	sig2B <- envir[["sig2B"]]
	tau2A <- envir[["tau2A"]]
	tau2B <- envir[["tau2B"]]
	corrA.BB <- envir[["corrA.BB"]]
	corrB.AA <- envir[["corrB.AA"]]
	corr <- envir[["corr"]]
	nuA <- envir[["nuA"]]
	nuB <- envir[["nuB"]]
	phiA <- envir[["phiA"]]
	phiB <- envir[["phiB"]]
	emit <- array(NA, dim=c(nrow(A), ncol(A), 28))##SNPs x sample x 'truth'
	##AAAA, AAAB, AABB, ABBB, BBBB
	##AAAAA, AAAAB, AAABB, AABBB, ABBBB, BBBBB
	##AAAAAA, AAAAAB, AAAABB, AAABBB, AABBBB, ABBBBB, BBBBBB
	lA <- log2(A)
	lB <- log2(B)	
	X <- cbind(lA, lB)	
	counter <- 1##state counter
	for(CT in cnStates){
		cat(".")
		for(CA in 0:CT){
			CB <- CT-CA
			A.scale <- sqrt(tau2A[, p]*(CA==0) + sig2A[, p]*(CA > 0))
			B.scale <- sqrt(tau2B[, p]*(CB==0) + sig2B[, p]*(CB > 0))
			scale <- c(A.scale, B.scale)
			if(CA == 0 & CB == 0) rho <- 0
			if(CA == 0 & CB > 0) rho <- corrA.BB[, p]
			if(CA > 0 & CB == 0) rho <- corrB.AA[, p]
			if(CA > 0 & CB > 0) rho <- corr[, p]
			if(CHR == 23) browser()
			means <- cbind(suppressWarnings(log2(nuA[, p]+CA*phiA[, p])), suppressWarnings(log2(nuB[, p]+CB*phiB[, p])))
			covs <- rho*A.scale*B.scale
			A.scale2 <- A.scale^2
			B.scale2 <- B.scale^2			
			##ensure positive definite			
			##Sigma <- as.matrix(nearPD(matrix(c(A.scale^2, covs,
			##covs, B.scale^2), 2, 2))[[1]])
			m <- 1##snp counter				
			for(i in 1:nrow(A)){
				Sigma <- matrix(c(A.scale2[i], covs[i], covs[i], B.scale2[i]), 2,2)
				xx <- matrix(X[i, ], ncol=2)
				tmp <- dmvnorm(xx, mean=means[i, ], sigma=Sigma) 				
				##Using HWE: P(CA=ca, CB=cb|CT=c)				
				ptmp <- (probA[i]^CA)*(probB[i]^CB)*tmp
				emit[m, , counter] <- ptmp
				m <- m+1				
			}
			counter <- counter+1			
		}
	}
	##priorProb=P(CT=c)
	homDel <- priorProb[1]*emit[, , 1]
	hemDel <- priorProb[2]*emit[, , c(2, 3)] # + priorProb[3]*emit[, c(4, 5, 6)] + priorProb[4]*emit[, c(7:10)]
	norm <- priorProb[3]*emit[, , 4:6]
	amp <- priorProb[4]*emit[, , 7:10]
	amp4 <- priorProb[5]*emit[, , 11:15]
	amp5 <- priorProb[6]*emit[, , 16:21]
	amp6 <- priorProb[7]*emit[, , 22:28]	
	##sum over the different combinations within each copy number state
	hemDel <- apply(hemDel, c(1,2), sum)
	norm <- apply(norm, c(1, 2), sum)
	amp <- apply(amp, c(1,2), sum)
	amp4 <- apply(amp4, c(1,2), sum)
	amp5 <- apply(amp5, c(1,2), sum)
	amp6 <- apply(amp6, c(1,2), sum)
	total <- homDel+hemDel+norm+amp+amp4+amp5+amp6
	weights <- array(NA, dim=c(nrow(homDel), ncol(A), 7))
	weights[, , 1] <- homDel/total
	weights[, , 2] <- hemDel/total
	weights[, , 3] <- norm/total
	weights[, , 4] <- amp/total
	weights[, , 5] <- amp4/total
	weights[, , 6] <- amp5/total
	weights[, , 7] <- amp6/total
	##posterior mode
	posteriorMode <- apply(weights, c(1, 2), function(x) order(x, decreasing=TRUE)[1])
	posteriorMode <- posteriorMode-1
	##This is for one plate.  Need to instantiate a much bigger
	##object in the environment
	
	##envir[["posteriorMode"]] <- posteriorMode
	##weights <- list(homDel/total, hemDel/total, norm/total, amp/total, amp4/total, amp5/total, amp6/total)
	##names(weights) <- c(cnStates)
	##envir[["weights"]] <- weights
	posteriorMeans <- 0*homDel/total + 1*hemDel/total + 2*norm/total + 3*amp/total + 4*amp4/total + 5*amp5/total + 6*amp6/total
	##sns <- envir[["sns"]]
	##colnames(posteriorMeans) <- sns
	##envir[["posteriorMeans"]] <- posteriorMeans
	return(posteriorMode)
}

biasAdjNP <- function(plateIndex, envir, priorProb){
	p <- plateIndex
	normalNP <- envir[["normalNP"]]
	CHR <- envir[["chrom"]]
	NP <- envir[["NP"]]
	plate <- envir[["plate"]]
	uplate <- envir[["uplate"]]
	sig2T <- envir[["sig2T"]]
	gender <- envir[["gender"]]
	normalNP <- normalNP[, plate==uplate[p]]	
	NP <- NP[, plate==uplate[p]]
	sig2T <- sig2T[, p]


	##Assume that on the log-scale, that the background variance is the same...
	tau2T <- sig2T	
	nuT <- envir[["nuT"]]
	nuT <- nuT[, p]
	phiT <- envir[["phiT"]]
	phiT <- phiT[, p]
	
	if(missing(priorProb)) priorProb <- rep(1/4, 4) ##uniform
	emit <- array(NA, dim=c(nrow(NP), ncol(NP), 4))##SNPs x sample x 'truth'	
	lT <- log2(NP)
	counter <- 1 ##state counter
	for(CT in 0:3){
		sds <- sqrt(tau2T*(CT==0) + sig2T*(CT > 0))
		means <- suppressWarnings(log2(nuT+CT*phiT))
		tmp <- dnorm(lT, mean=means, sd=sds)
		emit[, , counter] <- tmp
		counter <- counter+1
	}
	mostLikelyState <- apply(emit, c(1, 2), function(x) order(x, decreasing=TRUE)[1])
	if(CHR == 23){
		## the state index for male on chromosome 23  is 2
		## add 1 so that the state index is 3 for 'normal' state
		mostLikelyState[, gender=="male"] <- mostLikelyState[, gender=="male"] + 1
	}
	tmp3 <- mostLikelyState != 3
	##Those near 1 have NaNs for nu and phi.  this occurs by NaNs in the muA[,, "A"] or muA[, , "B"] for X chromosome
	proportionSamplesAltered <- rowMeans(tmp3)##prop normal
	ii <- proportionSamplesAltered < 0.75
	moreup <- rowSums(mostLikelyState > 3) > rowSums(mostLikelyState < 3)
	notUp <-  mostLikelyState[ii & moreup, ] <= 3
	notDown <- mostLikelyState[ii & !moreup, ] >= 3
	NORM <- matrix(TRUE, nrow(NP), ncol(NP))
	NORM[ii & moreup, ] <- notUp
	NORM[ii & !moreup, ] <- notDown
	normalNP <- normalNP*NORM

	flagAltered <- which(proportionSamplesAltered > 0.5)
	envir[["flagAlteredNP"]] <- flagAltered
	tmp <- envir[["normalNP"]]
	tmp[, plate==uplate[p]] <- normalNP
	envir[["normalNP"]] <- tmp
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

computeEmission <- function(crlmmResults, copyNumberStates=0:5, MIN=2^3,
			    EMIT.THR,
			    scaleSds=TRUE){
	##threshold small nu's and phis
	cnset <- thresholdModelParams(crlmmResults[[3]], MIN=MIN)
	index <- order(chromosome(cnset), position(cnset))
	
	if(any(diff(index) > 1)) stop("must be ordered by chromosome and physical position")
	emissionProbs <- array(NA, dim=c(nrow(cnset), ncol(cnset), length(copyNumberStates)))
	dimnames(emissionProbs) <- list(featureNames(crlmmResults),
					sampleNames(crlmmResults),
					paste("copy.number_", copyNumberStates, sep=""))	
	b <- batch(cnset)
	for(i in seq(along=unique(b))){
		if(i == 1) cat("Computing emission probabilities \n")
		message("batch ", unique(b)[i], "...")
		emissionProbs[, b == unique(b)[i], ] <- .getEmission(crlmmResults, cnset, batch=unique(b)[i],
				copyNumberStates=copyNumberStates,
				scaleSds=scaleSds)
	}
	if(missing(EMIT.THR)){
		EMIT.THR <- quantile(emissionProbs, probs=0.25, na.rm=TRUE)
		message("Thresholding emission probabilities at a small negative value (", round(EMIT.THR, 1), ") to reduce influence of outliers.")
		emissionProbs[emissionProbs < EMIT.THR] <- EMIT.THR
	}
	emissionProbs
}

thresholdModelParams <- function(object, MIN=2^3){
	nuA <- fData(object)[, grep("nuA", fvarLabels(object))]
	nuB <- fData(object)[, grep("nuB", fvarLabels(object))]
	phiA <- fData(object)[, grep("phiA", fvarLabels(object))]
	phiB <- fData(object)[, grep("phiB", fvarLabels(object))]

	nuA[nuA < MIN] <- MIN
	nuB[nuB < MIN] <- MIN
	phiA[phiA < MIN] <- MIN
	phiB[phiB < MIN] <- MIN
	fData(object)[, grep("nuA", fvarLabels(object))] <- nuA
	fData(object)[, grep("nuB", fvarLabels(object))] <- nuB
	fData(object)[, grep("phiA", fvarLabels(object))] <- phiA
	fData(object)[, grep("phiB", fvarLabels(object))] <- phiB
	object
}

.getEmission <- function(crlmmResults, cnset, batch, copyNumberStates, scaleSds=TRUE){
	if(length(batch) > 1) stop("batch variable not unique")
	crlmmResults <- crlmmResults[, cnset$batch==batch]
	cnset <- cnset[, cnset$batch == batch]

##	a <- A(crlmmResults)
##	b <- B(crlmmResults)	
##	sds.a <- apply(log2(a), 2, sd, na.rm=TRUE)
##	sds.b <- apply(log2(b), 2, sd, na.rm=TRUE)
	if(scaleSds){
		a <- CA(crlmmResults)
		b <- CB(crlmmResults)
		sds.a <- apply(a, 2, sd, na.rm=TRUE)
		sds.b <- apply(b, 2, sd, na.rm=TRUE)	
	
		sds.a <- log2(sds.a/median(sds.a))
		sds.b <- log2(sds.b/median(sds.b))
		
		sds.a <- matrix(sds.a, nrow(cnset), ncol(cnset), byrow=TRUE)
		sds.b <- matrix(sds.b, nrow(cnset), ncol(cnset), byrow=TRUE)

	} else sds.a <- sds.b <- matrix(0, nrow(cnset), ncol(cnset))
##	ca <- CA(cnset)
##	sds.ca <- apply(ca, 2, sd, na.rm=T)
##	sds.ca <- sds.ca/median(sds.ca)
##	sds.scale <- sds/median(sds)  #scale snp-specific variance by measure of the relative sample noise
	
	emissionProbs <- array(NA, dim=c(nrow(crlmmResults[[1]]),
				   ncol(crlmmResults[[1]]), length(copyNumberStates)))
	snpset <- cnset[snpIndex(cnset), ]
	params <- getParams(snpset, batch=batch)
	##attach(params)
	corr <- params[["corr"]]
	corrA.BB <- params[["corrA.BB"]]
	corrB.AA <- params[["corrB.AA"]]	
	nuA <- params[["nuA"]]
	nuB <- params[["nuB"]]
	phiA <- params[["phiA"]]
	phiB <- params[["phiB"]]
	sig2A <- params[["sig2A"]]
	sig2B <- params[["sig2B"]]
	tau2A <- params[["tau2A"]]
	tau2B <- params[["tau2B"]]
	a <- as.numeric(log2(A(crlmmResults[snpIndex(crlmmResults), ])))
	b <- as.numeric(log2(B(crlmmResults[snpIndex(crlmmResults), ])))
	for(k in seq(along=copyNumberStates)){
		##cat(k, " ")
		CN <- copyNumberStates[k]
		f.x.y <- matrix(0, nrow(snpset), ncol(snpset))
		for(CA in 0:CN){
			CB <- CN-CA
			sigmaA <- sqrt(tau2A*(CA==0) + sig2A*(CA > 0))
			sigmaB <- sqrt(tau2B*(CB==0) + sig2B*(CB > 0))
			if(CA == 0 & CB > 0) r <- corrA.BB
			if(CA > 0 & CB == 0) r <- corrB.AA
			if(CA > 0 & CB > 0) r <- corr
			if(CA == 0 & CB == 0) r <- 0
			muA <- log2(nuA+CA*phiA)
			muB <- log2(nuB+CB*phiB)

			sigmaA <- matrix(sigmaA, nrow=length(sigmaA), ncol=ncol(cnset), byrow=FALSE)
			sigmaB <- matrix(sigmaB, nrow=length(sigmaB), ncol=ncol(cnset), byrow=FALSE)
			sigmaA <- sigmaA+sds.a[snpIndex(crlmmResults), ]
			sigmaB <- sigmaB+sds.b[snpIndex(crlmmResults), ]			

			##might want to allow the variance to be sample-specific
			##TODO:
			## rho, sd.A, and sd.B are locus-specific
			## Some samples are more noisy than others.
			##
			##  - scale the variances by a sample-specific estimate of the variances
			## var(I_A, ijp) = sigma_A_ip * sigma_A_jp

			meanA <- as.numeric(matrix(muA, nrow(snpset), ncol(snpset)))
			meanB <- as.numeric(matrix(muB, nrow(snpset), ncol(snpset)))			
			rho <- as.numeric(matrix(r, nrow(snpset), ncol(snpset)))
			sd.A <- as.numeric(matrix(sigmaA, nrow(snpset), ncol(snpset)))
			sd.B <- as.numeric(matrix(sigmaB, nrow(snpset), ncol(snpset)))
			Q.x.y <- 1/(1-rho^2)*(((a - meanA)/sd.A)^2 + ((b - meanB)/sd.B)^2 - 2*rho*((a - meanA)*(b - meanB))/(sd.A*sd.B))
			##for CN states > 1, assume that any of the possible genotypes are equally likely a priori...just take the sum
			##for instance, for state copy number 2 there are three combinations: AA, AB, BB
			##   -- two of the three combinations should be near zero.
			## TODO: copy-neutral LOH would put near-zero mass on CA > 0, CB > 0
			f.x.y <- f.x.y + matrix(1/(2*pi*sd.A*sd.B*sqrt(1-rho^2))*exp(-0.5*Q.x.y), nrow(snpset), ncol(snpset))
		}
		emissionProbs[snpIndex(crlmmResults), , k] <- log(f.x.y)
	}


	##**************************************************
	##
	##Emission probabilities for nonpolymorphic probes
	##	
	##**************************************************	
	cnset <- cnset[cnIndex(cnset), ]
	params <- getParams(cnset, batch=batch)
	nuA <- params[["nuA"]]
	phiA <- params[["phiA"]]
	sig2A <- params[["sig2A"]]
	a <- as.numeric(log2(A(crlmmResults[cnIndex(crlmmResults), ])))
	for(k in seq(along=copyNumberStates)){
		CT <- copyNumberStates[k]
		mus.matrix=matrix(log2(nuA + CT*phiA), nrow(cnset), ncol(cnset))
		mus <- as.numeric(matrix(log2(nuA + CT*phiA), nrow(cnset), ncol(cnset)))
		##Again, should make sds sample-specific
		sds.matrix <- matrix(sqrt(sig2A), nrow(cnset), ncol(cnset))

		sds.matrix <- sds.matrix + sds.a[cnIndex(crlmmResults), ]
		sds <- as.numeric(sds.matrix)
		suppressWarnings(tmp <- matrix(dnorm(a, mean=mus, sd=sds), nrow(cnset), ncol(cnset)))
		emissionProbs[cnIndex(crlmmResults), , k] <- log(tmp)
	}
	emissionProbs
}
