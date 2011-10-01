crlmmGT2 <- function(A, B, SNR, mixtureParams, cdfName, row.names=NULL,
                     col.names=NULL, probs=c(1/3, 1/3, 1/3), DF=6,
                     SNRMin=5, recallMin=10, recallRegMin=1000,
                     gender=NULL, desctrucitve=FALSE, verbose=TRUE,
                     returnParams=FALSE, badSNP=.7){
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
	snp.names <- snpNames(pkgname)
	stopifnot(!is.null(rownames(A)))
	index <- match(snp.names, rownames(A))
##	if(!missing(snp.names)){
##
##		##verify that A has only snps.  otherwise, calling function must pass rownames
##		index <- match(snp.names, rownames(A))
##	}
	snpBatches <- splitIndicesByLength(index, ocProbesets())
	NR <- length(unlist(snpBatches))
	if(verbose) message("Calling ", NR, " SNPs for recalibration... ")
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
		gender <- imputeGender(A, B, XIndex, YIndex, SNR, SNRMin)
	}
	##
	Indexes <- list(autosomeIndex, XIndex, YIndex)
	cIndexes <- list(keepIndex,
			 keepIndex[which(gender[keepIndex]==2)],
			 keepIndex[which(gender[keepIndex]==1)])
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
	##if(verbose) message("Calling process1")
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
