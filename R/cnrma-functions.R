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
	return(new("AnnotatedDataFrame", data=data.frame(M)))
	##list(snpIndex, npIndex, fns)
	##crlmmOpts$snpRange <- range(snpIndex)
	##crlmmOpts$npRange <- range(npIndex)
}
construct <- function(filenames, cdfName, copynumber=FALSE, sns){
	if(missing(sns)) sns <- basename(filenames)
	protocolData <- getProtocolData.Affy(filenames)
	featureData <- getFeatureData.Affy(cdfName, copynumber=copynumber)
	nr <- nrow(featureData); nc <- length(filenames)
	callSet <- new("SnpSuperSet", 
			 alleleA=initializeBigMatrix(name="A", nr, nc),
			 alleleB=initializeBigMatrix(name="B", nr, nc),
			 call=initializeBigMatrix(name="call", nr, nc),
			 callProbability=initializeBigMatrix(name="callPr", nr,nc),
			 featureData=featureData,
			 annotation=cdfName)
	protocolData(callSet) <- protocolData
	sampleNames(callSet) <- sns
	return(callSet)
}
setReplaceMethod("calls", "SnpSuperSet", function(object, value) assayDataElementReplace(object, "call", value))
setReplaceMethod("confs", "SnpSuperSet", function(object, value) assayDataElementReplace(object, "callProbability", value))
setMethod("confs", "SnpSuperSet", function(object) assayDataElement(object, "callProbability"))
crlmm.batch <- function(object,
			batchSize,
			mixtureParams,
			probs=rep(1/3,3),
			DF=6,
			SNRMin=5,
			recallMin=10,
			recallRegMin=1000,
			gender=NULL,
			desctrucitve=FALSE,
			verbose=TRUE,
			returnParams=FALSE,
			badSNP=0.7){
	##Call in batches to reduce ram
	BS <- batchSize
	nc <- ncol(object)
	if(nc > BS){
		N <- ceiling(nc/BS)
		S <- ceiling(nc/N)
		colindex <- split(1:nc, rep(1:nc, each=S, length.out=nc))
	} else {
		colindex <- list(1:nc)
	}
	if(length(colindex) > 1)
		message("Calling genotypes in batches of size ", length(colindex[[1]]), " to reduce required RAM")
	row.index <- which(isSnp(object)==1 | is.na(isSnp(object)))
	for(i in seq(along=colindex)){
		col.index <- colindex[[i]]
		tmp <- crlmmGT(A=A(object)[row.index, col.index],
			       B=B(object)[row.index, col.index],
			       SNR=object$SNR[col.index],
			       mixtureParams=mixtureParams[, col.index],
			       cdfName=annotation(object),
			       row.names=featureNames(object)[row.index],
			       col.names=sampleNames(object)[col.index],
			       probs=probs,
			       DF=DF,
			       SNRMin=SNRMin,
			       recallMin=recallMin,
			       recallRegMin=recallRegMin,
			       gender=gender,
			       desctrucitve=desctrucitve,
			       verbose=verbose,
			       returnParams=returnParams,
			       badSNP=badSNP)
		##ensure matrix is passed
		calls(object)[row.index, col.index] <- tmp[["calls"]]
		confs(object)[row.index, col.index] <- tmp[["confs"]]
		object$gender[col.index] <- tmp$gender
		rm(tmp); gc()
	}
	return(object)
}

genotype <- function(filenames, cdfName, mixtureSampleSize=10^5,
		     fitMixture=TRUE,
		     eps=0.1, verbose=TRUE,
		     seed=1, sns, copynumber=FALSE,
		     batchSize=1000,
		     probs=rep(1/3, 3),
		     DF=6,
		     SNRMin=5,
		     recallMin=10,
		     recallRegMin=1000,
		     gender=NULL,
		     returnParams=TRUE,
		     badSNP=0.7){
	if(missing(cdfName)) stop("must specify cdfName")
	if(!isValidCdfName(cdfName)) stop("cdfName not valid.  see validCdfNames")
	if(missing(sns)) sns <- basename(filenames)
	## callSet contains potentially very big matrices
	## More big matrices are created within snprma, that will then be removed.
	snprmaRes <- snprma(filenames=filenames,
			    mixtureSampleSize=mixtureSampleSize,
			    fitMixture=fitMixture,
			    eps=eps,
			    verbose=verbose,
			    seed=seed,
			    cdfName=cdfName,
			    sns=sns)
	message("Initializing container for assay data elements alleleA, alleleB, call, callProbability")
	callSet <- construct(filenames=filenames,
			     cdfName=cdfName,
			     copynumber=copynumber,
			     sns=sns)
	sampleStats <- data.frame(SKW=snprmaRes$SKW,
				  SNR=snprmaRes$SNR)
	pD <- new("AnnotatedDataFrame",
		  data=sampleStats,
		  varMetadata=data.frame(labelDescription=colnames(sampleStats)))
	sampleNames(pD) <- sampleNames(callSet)
	phenoData(callSet) <- pD
	if(!copynumber){
		## A and B are the right size
		A(callSet) <- snprmaRes[["A"]]  ## should work regardless of ff, matrix
		B(callSet) <- snprmaRes[["B"]]
	} else {
		## A and B are not big enough to hold the nonpolymorphic markers
		index <- which(fData(callSet)[, "isSnp"] == 1)
		## Inefficient.  
		## write one column at a time (????).  (we don't want to bring in the whole matrix if its huge)
		if(isPackageLoaded("ff")){
			for(j in 1:ncol(callSet)){
				A(callSet)[index, j] <- snprmaRes[["A"]][, j]  
				B(callSet)[index, j] <- snprmaRes[["B"]][, j]
			}
			delete(snprmaRes[["A"]])##removes the file on disk
			delete(snprmaRes[["B"]])##removes the file on disk
		} else {
			A(callSet)[index, ] <- snprmaRes[["A"]]
			B(callSet)[index, ] <- snprmaRes[["B"]]
		}
	}
	gc()
	callSet <- crlmm.batch(object=callSet,
				batchSize=batchSize,
				mixtureParams=snprmaRes$mixtureParams,
				probs=probs,
				DF=DF,
				SNRMin=SNRMin,
				recallMin=recallMin,
				recallRegMin=recallRegMin,
				gender=gender,
				verbose=verbose,
				returnParams=returnParams,
				badSNP=badSNP)
	return(callSet)
}





##---------------------------------------------------------------------------
##---------------------------------------------------------------------------
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

combineIntensities <- function(res, NP=NULL, callSet){
	rownames(res$B) <- rownames(res$A) <- res$gns
	colnames(res$B) <- colnames(res$A) <- res$sns
	if(!is.null(NP)){
		blank <- matrix(NA, nrow(NP), ncol(NP))
		dimnames(blank) <- dimnames(NP)
		A <- rbind(res$A, NP)
		B <- rbind(res$B, blank)
	} else {
		A <- res$A
		B <- res$B
	}
	dimnames(B) <- dimnames(A)
	index.snps <- match(res$gns, rownames(A))
	callsConfs <- calls <- matrix(NA, nrow(A), ncol(A), dimnames=dimnames(A))
	
	calls[index.snps, ] <- calls(callSet)
	callsConfs[index.snps, ] <- assayData(callSet)[["callProbability"]]
	fd <- data.frame(matrix(NA, nrow(calls), length(fvarLabels(callSet))))
	fd[index.snps, ] <- fData(callSet)
	rownames(fd) <- rownames(A)
	colnames(fd) <- fvarLabels(callSet)
	fD <- new("AnnotatedDataFrame",
		  data=data.frame(fd),
		  varMetadata=data.frame(labelDescription=colnames(fd), row.names=colnames(fd)))
	superSet <- new("CNSet",
			CA=matrix(NA, nrow(A), ncol(A), dimnames=dimnames(A)),
			CB=matrix(NA, nrow(A), ncol(A), dimnames=dimnames(A)),
			alleleA=A,
			alleleB=B,
			call=calls,
			callProbability=callsConfs,
			featureData=fD,
			phenoData=phenoData(callSet),
			experimentData=experimentData(callSet),
			protocolData=protocolData(callSet),
			annotation=annotation(callSet))
	return(superSet)
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

crlmmCopynumber <- function(filenames, cnOptions, object, ...){
	if(!missing(object)){
		stopifnot(class(object) == "CNSet")
		createIntermediateObjects <- FALSE
	} else {
		createIntermediateObjects <- TRUE
		crlmmResults <- crlmmWrapper(filenames, cnOptions, ...)
		snprmaResult <- crlmmResults[["snprmaResult"]]
		cnrmaResult <- crlmmResults[["cnrmaResult"]]
		callSet <- crlmmResults[["callSet"]]
		rm(crlmmResults); gc()
		annotation(callSet) <- cnOptions[["cdfName"]]
		stopifnot(identical(featureNames(callSet), snprmaResult[["gns"]]))
		path <- system.file("extdata", package=paste(annotation(callSet), "Crlmm", sep=""))	
		load(file.path(path, "snpProbes.rda"))
		snpProbes <- get("snpProbes")
		load(file.path(path, "cnProbes.rda"))
		cnProbes <- get("cnProbes")	
		k <- grep("chr", colnames(snpProbes))
		if(length(k) < 1) stop("chr or chromosome not in colnames(snpProbes)")
	}
	set.seed(cnOptions[["seed"]])  ##for reproducibility
	chromosome <- cnOptions[["chromosome"]]
	SNRmin <- cnOptions[["SNRmin"]]
	for(CHR in chromosome){
		cat("Chromosome ", CHR, "\n")
		if(createIntermediateObjects){
			snps <- rownames(snpProbes)[snpProbes[, k] == CHR]
			cnps <- rownames(cnProbes)[cnProbes[, k] == CHR]
			index.snps <- match(snps, featureNames(callSet))
			index.nps <- match(cnps, rownames(cnrmaResult[["NP"]]))
			if(!is.null(cnrmaResult)){
				npI <- cnrmaResult$NP[index.nps, ]
			} else npI <- NULL
			snpI <- list(A=snprmaResult$A[index.snps, ],
				     B=snprmaResult$B[index.snps, ],
				     sns=snprmaResult$sns,
				     gns=snprmaResult$gns[index.snps],
				     SNR=snprmaResult$SNR,
				     SKW=snprmaResult$SKW,
				     mixtureParams=snprmaResult$mixtureParams,
				     cdfName=snprmaResult$cdfName)
			cnOptions[["batch"]] <- cnOptions[["batch"]][snpI[["SNR"]]  >= SNRmin]
			cnSet <- combineIntensities(res=snpI,
						    NP=npI,
						    callSet=callSet[index.snps, ])
			if(any(cnSet$SNR > SNRmin)){
				message(paste("Excluding samples with SNR < ", SNRmin))
				cnSet <- cnSet[, cnSet$SNR >= SNRmin]
			}
			rm(snpI, npI, snps, cnps, index.snps, index.nps); gc()
			pData(cnSet)$batch <- cnOptions[["batch"]]
			featureData(cnSet) <- lm.parameters(cnSet, cnOptions)
		} else {
			cnSet <- object
		}
		if(CHR != 24){		
			cnSet <- computeCopynumber(cnSet, cnOptions)
		} else{
			message("Copy number estimates not available for chromosome Y.  Saving only the 'callSet' object for this chromosome")
			alleleSet <- cnSet
			save(alleleSet, file=file.path(cnOptions[["outdir"]], paste("alleleSet_", CHR, ".rda", sep="")))
			rm(cnSet, alleleSet); gc()
			next()
		}
		if(length(chromosome) == 1){
			if(cnOptions[["save.cnset"]]){
				save(cnSet, file=file.path(cnOptions[["outdir"]], paste("cnSet_", CHR, ".rda", sep="")))
			}
		}
		if(length(chromosome) > 1){
			save(cnSet, file=file.path(cnOptions[["outdir"]], paste("cnSet_", CHR, ".rda", sep="")))
			rm(cnSet); gc()
		} else {
			return(cnSet)
		}
	}
	saved.objects <- list.files(cnOptions[["outdir"]], pattern="cnSet", full.names=TRUE)		
	return(saved.objects)
}



crlmmWrapper <- function(filenames, cnOptions, ...){
	cdfName <- cnOptions[["cdfName"]]
	load.it <- cnOptions[["load.it"]]
	save.it <- cnOptions[["save.it"]]
	splitByChr <- cnOptions[["splitByChr"]]
	crlmmFile <- cnOptions[["crlmmFile"]]
	intensityFile=cnOptions[["intensityFile"]]
	rgFile=cnOptions[["rgFile"]]
	##use.ff=cnOptions[["use.ff"]]
	outdir <- cnOptions[["outdir"]]
	if(missing(cdfName)) stop("cdfName is missing -- a valid cdfName is required.  See crlmm:::validCdfNames()")
	platform <- whichPlatform(cdfName)
	if(!(platform %in% c("affymetrix", "illumina"))){
		stop("Only 'affymetrix' and 'illumina' platforms are supported at this time.")
	} else {
		message("Checking whether annotation package for the ", platform, " platform is available")
		if(!isValidCdfName(cdfName)){
			cat("FALSE\n")
			stop(cdfName, " is not a valid entry.  See crlmm:::validCdfNames(platform)")
		} else cat("TRUE\n")
	}
	if(missing(intensityFile)) stop("must specify 'intensityFile'.")
	if(missing(crlmmFile)) stop("must specify 'crlmmFile'.")
	if(platform == "illumina"){
		if(missing(rgFile)){
			##stop("must specify 'rgFile'.")
			rgFile <- file.path(dirname(crlmmFile), "rgFile.rda")
			message("rgFile not specified.  Using ", rgFile)
		}
		if(!load.it){
			RG <- readIdatFiles(...)
			if(save.it) save(RG, file=rgFile)
		}
		if(load.it & !file.exists(rgFile)){
			message("load.it is TRUE, but rgFile not present.  Attempting to read the idatFiles.")
			RG <- readIdatFiles(...)
			if(save.it) save(RG, file=rgFile)
		}
		if(load.it & file.exists(rgFile)){
			message("Loading RG file")
			load(rgFile)
			RG <- get("RG")
		}
	}
	if(!(file.exists(dirname(crlmmFile)))) stop(dirname(crlmmFile), " does not exist.")
	if(!(file.exists(dirname(intensityFile)))) stop(dirname(intensityFile), " does not exist.")
	##---------------------------------------------------------------------------
	## FIX
	outfiles <- file.path(dirname(crlmmFile), paste("crlmmSetList_", 1:24, ".rda", sep=""))
	if(load.it & all(file.exists(outfiles))){
		load(outfiles[1])
		crlmmSetList <- get("crlmmSetList")
		if(!all(sampleNames(crlmmSetList) == basename(filenames))){
			stop("load.it is TRUE, but sampleNames(crlmmSetList != basename(filenames))")
		} else{
			return("load.it is TRUE and 'crlmmSetList_<CHR>.rda' objects found. Nothing to do...")
		}
	}
	if(load.it){
		if(!file.exists(crlmmFile)){
			message("load.it is TRUE, but ", crlmmFile, " does not exist.  Rerunning the genotype calling algorithm") 
			load.it <- FALSE
		}
	}
	if(platform == "affymetrix"){
		if(!file.exists(crlmmFile) | !load.it){
			callSet <- crlmm(filenames=filenames,
					     cdfName=cdfName,
					     save.it=TRUE,
					     load.it=load.it,
					     intensityFile=intensityFile)
			message("Quantile normalizing the copy number probes...")		
			cnrmaResult <- cnrma(filenames=filenames, cdfName=cdfName, outdir=outdir)
			if(save.it){
				message("Saving callSet and cnrmaResult to", crlmmFile)
				save(callSet, cnrmaResult, file=crlmmFile)
			}
		} else {
			message("Loading ", crlmmFile, "...")
			load(intensityFile)				
			load(crlmmFile)
			callSet <- get("callSet")
			cnrmaResult <- get("cnrmaResult")
		}
		scanDates <- data.frame(ScanDate=sapply(filenames, celfileDate))
		protocolData(callSet) <- new("AnnotatedDataFrame",
					     data=scanDates,
					     varMetadata=data.frame(labelDescription=colnames(scanDates),
					     row.names=colnames(scanDates)))
	}
	if(platform == "illumina"){
		if(!file.exists(crlmmFile) | !load.it){		
			callSet <- crlmmIllumina(RG=RG,
						 cdfName=cdfName,
						 sns=sampleNames(RG),
						 returnParams=TRUE,
						 save.it=TRUE,
						 intensityFile=intensityFile)
			if(save.it) save(callSet, file=crlmmFile)
		} else {
			message("Loading ", crlmmFile, "...")
			load(crlmmFile)
			callSet <- get("callSet")
		}
		protocolData(callSet) <- protocolData(RG)
	}
	if(platform=="affymetrix") {
		protocolData(callSet)[["ScanDate"]] <- as.character(celDates(filenames))
		sampleNames(protocolData(callSet)) <- sampleNames(callSet)
	}	
	load(intensityFile)
	snprmaResult <- get("res")
	if(platform=="illumina"){
		if(exists("cnAB")){
			np.A <- as.integer(cnAB$A)
			np.B <- as.integer(cnAB$B)
			np <- ifelse(np.A > np.B, np.A, np.B)
			np <- matrix(np, nrow(cnAB$A), ncol(cnAB$A))
			rownames(np) <- cnAB$gns
			colnames(np) <- cnAB$sns
			cnAB$NP <- np
			##sampleNames(callSet) <- res$sns
			sampleNames(callSet) <- cnAB$sns
			cnrmaResult <- get("cnAB")
		} else cnrmaResult <- NULL
	}
	if(platform=="affymetrix"){
		if(exists("cnrmaResult")){
			cnrmaResult <- get("cnrmaResult")
		} else cnrmaResult <- NULL
	}
	crlmmResults <- list(snprmaResult=snprmaResult,
			     cnrmaResult=cnrmaResult,
			     callSet=callSet)

	if(!save.it){
		message("Cleaning up")		
		unlink(intensityFile)
	}
	return(crlmmResults)
}

validCdfNames <- function(){
	c("genomewidesnp6",
	  "genomewidesnp5",
	  "human370v1c",
	  "human370quadv3c",
	  "human550v3b",
	  "human650v3a",
	  "human610quadv1b",
	  "human660quadv1a",
	  "human1mduov3b")
}

isValidCdfName <- function(cdfName){
	chipList <- validCdfNames()
	result <- cdfName %in% chipList	
	if(!(result)){
		warning("cdfName must be one of the following: ",
			chipList)
	}
	return(result)
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


# steps: quantile normalize hapmap: create 1m_reference_cn.rda object
cnrma <- function(filenames, cdfName, sns, seed=1, verbose=FALSE, outdir){
	if(missing(cdfName)) stop("must specify cdfName")
	pkgname <- getCrlmmAnnotationName(cdfName)
	require(pkgname, character.only=TRUE) || stop("Package ", pkgname, " not available")
	if (missing(sns)) sns <- basename(filenames)
        loader("npProbesFid.rda", .crlmmPkgEnv, pkgname)
	fid <- getVarInEnv("npProbesFid")
	set.seed(seed)
	idx2 <- sample(length(fid), 10^5) ##for skewness. no need to do everything
	SKW <- vector("numeric", length(filenames))
##	if(bigmemory){
##		NP <- filebacked.big.matrix(length(pnsa), length(filenames),
##					    type="integer",
##					    init=as.integer(0),
##					    backingpath=outdir,
##					    backingfile="NP.bin",
##					    descriptorfile="NP.desc")
##	} else{
		NP <- matrix(NA, length(fid), length(filenames))
##	}
	verbose <- TRUE
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
	##if(!is.matrix(reference)) stop("target distribution for quantile normalization not available.")
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
	cat("\n")
	return(res3)
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


instantiateObjects <- function(object, cnOptions){
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

##preprocessOptions <- function(crlmmFile="snpsetObject.rda",
##			      intensityFile="normalizedIntensities.rda",
##			      rgFile="rgFile.rda"){
##
##}

cnOptions <- function(outdir="./",
		      cdfName,
		      crlmmFile,
		      intensityFile,
		      rgFile="rgFile.rda",
		      save.it=TRUE,
		      save.cnset=TRUE,
		      load.it=TRUE,
		      splitByChr=TRUE,
		      MIN.OBS=3,
		      MIN.SAMPLES=10,
		      batch=NULL,
		      DF.PRIOR=50,
		      bias.adj=FALSE,
		      prior.prob=rep(1/4, 4),
		      SNRmin=4,
		      chromosome=1:24,
		      seed=123,
		      verbose=TRUE,
		      GT.CONF.THR=0.99,
		      PHI.THR=2^6,##used in nonpolymorphic fxn for training
		      nHOM.THR=5, ##used in nonpolymorphic fxn for training
		      MIN.NU=2^3,
		      MIN.PHI=2^3,
		      THR.NU.PHI=TRUE,
		      thresholdCopynumber=TRUE,
		      unlink=TRUE,
		      ##hiddenMarkovModel=FALSE,
		      ##circularBinarySegmentation=FALSE,
##		      cbsOpts=NULL,
		      ##hmmOpts=NULL,
		      ...){
	if(missing(cdfName)) stop("must specify cdfName")
	if(!file.exists(outdir)){
		message(outdir, " does not exist.  Trying to create it.")
		dir.create(outdir, recursive=TRUE)
	}
	stopifnot(isValidCdfName(cdfName))
##	if(hiddenMarkovModel){
##		hmmOpts <- hmmOptions(...)
##	}
	if(missing(crlmmFile)){
		crlmmFile <- file.path(outdir, "snpsetObject.rda")
	}
	if(missing(intensityFile)){
		intensityFile <- file.path(outdir, "normalizedIntensities.rda")
	}
	if(is.null(batch))
		stop("must specify batch -- should be the same length as the number of files")
	list(outdir=outdir,
	     cdfName=cdfName,
	     crlmmFile=crlmmFile,
	     intensityFile=intensityFile,
	     rgFile=file.path(outdir, rgFile),
	     save.it=save.it,
	     save.cnset=save.cnset,
	     load.it=load.it,
	     splitByChr=splitByChr,
	     MIN.OBS=MIN.OBS,
	     MIN.SAMPLES=MIN.SAMPLES,
	     batch=batch,
	     DF.PRIOR=DF.PRIOR,
	     GT.CONF.THR=GT.CONF.THR,
	     prior.prob=prior.prob,
	     bias.adj=bias.adj,
	     SNRmin=SNRmin,
	     chromosome=chromosome,
	     seed=seed,
	     verbose=verbose,
	     PHI.THR=PHI.THR,
	     nHOM.THR=nHOM.THR,
	     MIN.NU=MIN.NU,
	     MIN.PHI=MIN.PHI,
	     THR.NU.PHI=THR.NU.PHI,
	     thresholdCopynumber=thresholdCopynumber,
	     unlink=unlink)
##	     hiddenMarkovModel=hiddenMarkovModel,
##	     circularBinarySegmentation=circularBinarySegmentation,
##	     cbsOpts=cbsOpts,
##	     hmmOpts=hmmOpts) ##remove SnpSuperSet object
}

##linear model parameters
lm.parameters <- function(object, cnOptions){
	fD <- fData(object)
	batch <- object$batch
	uplate <- unique(batch)
	parameterNames <- c(paste("tau2A", uplate, sep="_"),
			    paste("tau2B", uplate, sep="_"),
			    paste("sig2A", uplate, sep="_"),
			    paste("sig2B", uplate, sep="_"),
			    paste("nuA", uplate, sep="_"),
			    paste("nuA.se", uplate, sep="_"),			    
			    paste("nuB", uplate, sep="_"),
			    paste("nuB.se", uplate, sep="_"),			    			    
			    paste("phiA", uplate, sep="_"),
			    paste("phiA.se", uplate, sep="_"),			    
			    paste("phiB", uplate, sep="_"),
			    paste("phiB.se", uplate, sep="_"),			    
			    paste("phiAX", uplate, sep="_"),
			    paste("phiBX", uplate, sep="_"),			    
			    paste("corr", uplate, sep="_"),
			    paste("corrA.BB", uplate, sep="_"),
			    paste("corrB.AA", uplate, sep="_"))
	pMatrix <- data.frame(matrix(numeric(0),
				     nrow(object),
				     length(parameterNames)),
				     row.names=featureNames(object))
	colnames(pMatrix) <- parameterNames
	fD2 <- cbind(fD, pMatrix)
	new("AnnotatedDataFrame", data=fD2,
	    varMetadata=data.frame(labelDescription=colnames(fD2),
	    row.names=colnames(fD2)))
}

nonpolymorphic <- function(object, cnOptions, tmp.objects){
	batch <- unique(object$batch)
	CHR <- unique(chromosome(object))
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
			verbose <- cnOptions$verbose
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
			verbose <- cnOptions$verbose
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
		stop("fewer than 200 snps pass criteria for predicting the sufficient statistics")
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
	nobsA <- Ns[, "A"] > 10
	nobsB <- Ns[, "B"] > 10
	notMissing <- !(is.na(muA[, "A"]) | is.na(muA[, "B"]) | is.na(muB[, "A"]) | is.na(muB[, "B"]))
	complete <- list()
	complete[[1]] <- which(correct.orderA & correct.orderB & nobsA & notMissing) ##be selective here
	complete[[2]] <- which(correct.orderA & correct.orderB & nobsB & notMissing) ##be selective here	
	size <- min(5000, length(complete[[1]]))
	if(size == 5000) complete <- lapply(complete, function(x) sample(x, size))
	if(CHR == 23){
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

computeCopynumber.SnpSuperSet <- function(object, cnOptions){
##	use.ff <- cnOptions[["use.ff"]]
##	if(!use.ff){
##		object <- as(object, "CrlmmSet")
##	} else	object <- as(object, "CrlmmSetFF")
	bias.adj <- cnOptions[["bias.adj"]]
	##must be FALSE to initialize parameters
	cnOptions[["bias.adj"]] <- FALSE
	## Add linear model parameters to the CrlmmSet object
	featureData(object) <- lm.parameters(object, cnOptions)
	if(!isValidCdfName(annotation(object))) stop(annotation(object), " not supported.")
	object <- as(object, "CNSet")
	object <- computeCopynumber.CNSet(object, cnOptions)
	if(bias.adj==TRUE){## run a second time
		object <- computeCopynumber.CNSet(object, cnOptions)
	}
	return(object)
}


computeCopynumber.CNSet <- function(object, cnOptions){
	CHR <- unique(chromosome(object))
	batch <- object$batch
	if(length(batch) != ncol(object)) stop("Batch must be the same length as the number of samples")
	MIN.SAMPLES <- cnOptions$MIN.SAMPLES
	verbose <- cnOptions$verbose
	for(i in seq(along=unique(batch))){
		PLATE <- unique(batch)[i]
		if(sum(batch == PLATE) < MIN.SAMPLES) {
			message("Skipping plate ", PLATE)
			next()
		}		
		object.batch <- object[, batch==PLATE]
		tmp.objects <- instantiateObjects(object.batch,
						  cnOptions)
		bias.adj <- cnOptions$bias.adj
		if(bias.adj & ncol(object) <= 15){
			warning(paste("bias.adj is TRUE, but too few samples to perform this step"))
			cnOptions$bias.adj <- bias.adj <- FALSE
		}
		if(bias.adj){
			if(verbose) message("Dropping samples with low posterior prob. of normal copy number (samples dropped is locus-specific)")
			tmp.objects <- biasAdjNP(object.batch, cnOptions, tmp.objects)
			tmp.objects <- biasAdj(object.batch, cnOptions, tmp.objects)
			if(verbose) message("Recomputing location and scale parameters")
		}
		##update tmp.objects
		tmp.objects <- withinGenotypeMoments(object.batch,
						     cnOptions=cnOptions,
						     tmp.objects=tmp.objects)
		object.batch <- locationAndScale(object.batch, cnOptions, tmp.objects)
		tmp.objects <- oneBatch(object.batch,
					cnOptions=cnOptions,
					tmp.objects=tmp.objects)
		##coefs calls nuphiAllele.
		object.batch <- coefs(object.batch, cnOptions, tmp.objects)
		##nuA=getParam(object.batch, "nuA", PLATE)
		THR.NU.PHI <- cnOptions$THR.NU.PHI
		if(THR.NU.PHI){
			verbose <- cnOptions$verbose
			if(verbose) message("Thresholding nu and phi")
			object.batch <- thresholdModelParams(object.batch, cnOptions)
		}		
		if(verbose) message("\nAllele specific copy number")	
		object.batch <- polymorphic(object.batch, cnOptions, tmp.objects)
		if(any(!isSnp(object))){ ## there are nonpolymorphic probes
			if(verbose) message("\nCopy number for nonpolymorphic probes...")	
			object.batch <- nonpolymorphic(object.batch, cnOptions, tmp.objects)
		}
		##---------------------------------------------------------------------------
		##Note: the replacement method multiples by 100
		CA(object)[, batch==PLATE] <- CA(object.batch)
		CB(object)[, batch==PLATE] <- CB(object.batch)
		##---------------------------------------------------------------------------
		##update-the plate-specific parameters for copy number
		object <- pr(object, "nuA", PLATE, getParam(object.batch, "nuA", PLATE))
		object <- pr(object, "nuA.se", PLATE, getParam(object.batch, "nuA.se", PLATE))
		object <- pr(object, "nuB", PLATE, getParam(object.batch, "nuB", PLATE))
		object <- pr(object, "nuB.se", PLATE, getParam(object.batch, "nuB.se", PLATE))
		object <- pr(object, "phiA", PLATE, getParam(object.batch, "phiA", PLATE))
		object <- pr(object, "phiA.se", PLATE, getParam(object.batch, "phiA.se", PLATE))
		object <- pr(object, "phiB", PLATE, getParam(object.batch, "phiB", PLATE))
		object <- pr(object, "phiB.se", PLATE, getParam(object.batch, "phiB.se", PLATE))
		object <- pr(object, "tau2A", PLATE, getParam(object.batch, "tau2A", PLATE))
		object <- pr(object, "tau2B", PLATE, getParam(object.batch, "tau2B", PLATE))				
		object <- pr(object, "sig2A", PLATE, getParam(object.batch, "sig2A", PLATE))
		object <- pr(object, "sig2B", PLATE, getParam(object.batch, "sig2B", PLATE))		
		object <- pr(object, "phiAX", PLATE, as.numeric(getParam(object.batch, "phiAX", PLATE)))
		object <- pr(object, "phiBX", PLATE, as.numeric(getParam(object.batch, "phiBX", PLATE)))
		object <- pr(object, "corr", PLATE, getParam(object.batch, "corr", PLATE))
		object <- pr(object, "corrA.BB", PLATE, getParam(object.batch, "corrA.BB", PLATE))
		object <- pr(object, "corrB.AA", PLATE, getParam(object.batch, "corrB.AA", PLATE))		
		rm(object.batch, tmp.objects); gc();
	}
	object <- object[order(chromosome(object), position(object)), ]
	if(cnOptions[["thresholdCopynumber"]]){
		object <- thresholdCopynumber(object)
	}
	return(object)
}







