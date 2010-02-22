crlmm <- function(filenames, row.names=TRUE, col.names=TRUE,
                  probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5, gender=NULL,
                  save.it=FALSE, load.it=FALSE,
		  snprmaFile="snp_rmaResult.rda",
		  callsFile="genotypes.rda",
		  confsFile="confs.rda",
		  AFile="A.rda",
		  BFile="B.rda",
                  mixtureSampleSize=10^5, eps=0.1, verbose=TRUE,
                  cdfName, sns, recallMin=10, recallRegMin=1000,
                  returnParams=FALSE, badSNP=0.7,
		  crlmmBatchSize=1000){
	BS <- crlmmBatchSize
	if(load.it & file.exists(snprmaFile) & file.exists(callsFile)) return() ##nothing to do
	if (load.it & !file.exists(snprmaFile)){
		##stop("'intensityFile' is missing, and you chose either load.it or save.it")
		message("'snprmaFile' does not exist and you chose to load.it.  Rerunning snprma...")
		load.it <- FALSE
	} 
	if (missing(sns)) sns <- basename(filenames)
	if (!load.it){
		res <- snprma(filenames,
			      fitMixture=TRUE,
			      mixtureSampleSize=mixtureSampleSize,
			      verbose=verbose,
			      eps=eps,
			      cdfName=cdfName,
			      sns=sns,
			      AFile=AFile,
			      BFile=BFile)
		save(res, file=snprmaFile)##file.path(dirname(snprmaFile), "res.rda"))
	} else {
		message("Loading snprmaFile...")
		load(snprmaFile)
		res <- get("res")
	}
	SKW <- res[["SKW"]]
	SNR <- res[["SNR"]]	
	load(AFile)
	load(BFile)
	if(isPackageLoaded("ff")) {open(A); open(B)}
	if(row.names) row.names <- res$gns
	if(col.names) col.names <- res$sns	
##	}else{
##		if (verbose) message("Loading ", snprmaFile, ".")
##		obj <- load(snprmaFile)
##		if (verbose) message("Done.")
##		if (obj != "res")
##			stop("Object in ", snprmaFile, " seems to be invalid.")
##	}
	##path <- system.file("extdata", package=paste(cdfName, "Crlmm", sep=""))
	##load(file.path(path, "snpProbes.rda"))
	gns <- res[["gns"]]
	sns <- colnames(A)
	nc <- ncol(A)	
	##
	##Ensures that if ff package is loaded, ordinary matrices are still passed to crlmmGT
	if(nc > BS){
		## genotype the samples in batches

		##number crlmm batches
		N <- ceiling(nc/BS)
		##samples per batch
		S <- ceiling(nc/N)
		colindex <- split(1:nc, rep(1:nc, each=S, length.out=nc))
	} else {
		colindex <- list(1:nc)
	}
	if(length(colindex) > 1)
		message("Calling genotypes in batches of size ", length(colindex[[1]]), " to reduce required RAM")	
	if(isPackageLoaded("ff")){
		confs <- initializeBigMatrix(nrow(A), ncol(A))
		calls <- initializeBigMatrix(nrow(A), ncol(A))
	} else{
		confs <- matrix(NA, nrow(A), ncol(A))
		calls <- matrix(NA, nrow(A), ncol(A))
	}	
	sex <- vector("list", length(colindex))
	for(i in seq(along=colindex)){
		aMatrix <- A[1:length(gns), colindex[[i]]]
		bMatrix <- B[1:length(gns), colindex[[i]]]
		snr <- res[["SNR"]][colindex[[i]]]
		mix <- res[["mixtureParameters"]][, colindex[[i]]]
		##	res2 <- crlmmGT(res[["A"]], res[["B"]], res[["SNR"]],
		res2 <- tryCatch(crlmmGT(aMatrix,
					 bMatrix,
					 snr,
					 mix,
					 res[["cdfName"]],
					 gender=gender,
					 row.names=row.names,
					 col.names=col.names[colindex[[i]]],
					 recallMin=recallMin,
					 recallRegMin=1000,
					 SNRMin=SNRMin,
					 returnParams=returnParams,
					 badSNP=badSNP,
					 verbose=verbose), error=function(e) NULL)
		if(is.null(res2)) {
			unlink(callsFile)
			unlink(confsFile)
			stop("error in call to crlmmGT")
		}
		rm(aMatrix, bMatrix, snr, mix); gc()
		calls[1:length(gns), colindex[[i]]] <- res2[["calls"]]
		confs[1:length(gns), colindex[[i]]] <- res2[["confs"]]
		sex[[i]] <- res2[["gender"]]
		rm(res2); gc()
	}
	if(isPackageLoaded("ff")) close(B)
	rm(B); gc()
	if(isPackageLoaded("ff")) close(calls)
	save(calls, file=callsFile)
	rm(calls); gc()
	if(isPackageLoaded("ff")) close(confs)
	save(confs, file=confsFile)
	rm(confs); gc()
	res$gender <- unlist(sex)
	save(res, file=snprmaFile)
	if(isPackageLoaded("ff")) close(A)
	rm(list=ls()); gc()
	return()
}

crlmmGT <- function(A, B, SNR, mixtureParams, cdfName, row.names=NULL,
                    col.names=NULL, probs=c(1/3, 1/3, 1/3), DF=6,
                    SNRMin=5, recallMin=10, recallRegMin=1000,
                    gender=NULL, desctrucitve=FALSE, verbose=TRUE,
                    returnParams=FALSE, badSNP=.7){
  keepIndex <- which(SNR>SNRMin)
  if(length(keepIndex)==0) stop("No arrays above quality threshold!")
  NC <- ncol(A)
  NR <- nrow(A)
  pkgname <- getCrlmmAnnotationName(cdfName)
  if(!require(pkgname, character.only=TRUE, quietly=!verbose)){
    suggCall <- paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep="")
    msg <- paste("If", pkgname,
                 "is installed on an alternative location, please load it manually by using",
                 suggCall)
    message(strwrap(msg))
    stop("Package ", pkgname, " could not be found.")
    rm(suggCall, msg)
  }
  if(verbose) message("Loading annotations.")
  loader("genotypeStuff.rda", .crlmmPkgEnv, pkgname)
  loader("mixtureStuff.rda", .crlmmPkgEnv, pkgname)

  ## this is toget rid of the 'no visible binding' notes
  ## variable definitions
  XIndex <- getVarInEnv("XIndex")
  autosomeIndex <- getVarInEnv("autosomeIndex")
  YIndex <- getVarInEnv("YIndex")
  SMEDIAN <- getVarInEnv("SMEDIAN")
  theKnots <- getVarInEnv("theKnots")
  regionInfo <- getVarInEnv("regionInfo")
  params <- getVarInEnv("params")
  
  ##IF gender not provide, we predict
  if(is.null(gender)){
    if(verbose) message("Determining gender.")
    XMedian <- apply(log2(A[XIndex,, drop=FALSE])+log2(B[XIndex,, drop=FALSE]), 2, median)/2
    if(sum(SNR>SNRMin)==1){
      gender <- which.min(c(abs(XMedian-8.9), abs(XMedian-9.5)))
    }else{
      gender <- kmeans(XMedian, c(min(XMedian[SNR>SNRMin]), max(XMedian[SNR>SNRMin])))[["cluster"]]
    }
  }
  Indexes <- list(autosomeIndex, XIndex, YIndex)
  cIndexes <- list(keepIndex, 
                   keepIndex[which(gender[keepIndex]==2)], 
                   keepIndex[which(gender[keepIndex]==1)])
  
  if(verbose) cat("Calling", NR, "SNPs for recalibration... ")

  ## call C
  fIndex <- which(gender==2)
  mIndex <- which(gender==1)
  newparams <- gtypeCallerR(A, B, fIndex, mIndex,
                            params[["centers"]], params[["scales"]], params[["N"]],
                            Indexes, cIndexes,
                            sapply(Indexes, length), sapply(cIndexes, length),
                            SMEDIAN, theKnots,
                            mixtureParams, DF, probs, 0.025)
  gc(verbose=FALSE)
  names(newparams) <- c("centers", "scales", "N")
  
  if(verbose) message("Done.")
  if(verbose) message("Estimating recalibration parameters.")
  d <- newparams[["centers"]] - params$centers

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
  
    minN <- 3
    newparams[["centers"]][newparams[["N"]] < minN] <- NA
    Index <- setdiff(which(rowSums(is.na(newparams[["centers"]]))==1), YIndex)
    if(verbose) cat("Filling out empty centers")
    for(i in Index){
      if(verbose) if(i%%10000==0)cat(".")
      mu <- newparams[["centers"]][i, ]
      j <- which(is.na(mu))
      newparams[["centers"]][i, j] <- c(1, mu[-j], prod(mu[-j]))%*%regParams[, j]
    }
    
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
  
    if(verbose) message("Calculating and standardizing size of shift.")
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
    
  ## BC: must keep SD
  params[-2] <- newparams[-2]
  
  rm(newparams);gc(verbose=FALSE)  
  if(verbose) cat("Calling", NR, "SNPs... ")
  ## ###################
  ## ## MOVE TO C#######
  ImNull <- gtypeCallerR2(A, B, fIndex, mIndex, params[["centers"]],
                          params[["scales"]], params[["N"]], Indexes,
                          cIndexes, sapply(Indexes, length),
                          sapply(cIndexes, length), SMEDIAN, theKnots,
                          mixtureParams, DF, probs, 0.025,
                          which(regionInfo[,2]),
                          which(regionInfo[,1]))
  gc(verbose=FALSE)
  ##  END MOVE TO C#######
  ## ##################
  
  dev <- dev/(dev+1/383)
  if(!is.null(row.names)){ rownames(A) <- rownames(B) <- names(dev) <- row.names}
  if(!is.null(col.names)){ colnames(A) <- colnames(B) <- col.names}

  if(length(Index) >= recallRegMin){
   tmp4batchQC <- DD[autosomeIndex,]/(params[["N"]][autosomeIndex,]+1)
   tmpSnpQc <- dev[autosomeIndex]
   SS <- cov(tmp4batchQC[tmpSnpQc < badSNP,])
   batchQC <- mean(diag(SS))
  }else{
    SS <- matrix(0, 3, 3)
    batchQC <- Inf
  }
  
  if(verbose) message("Done.")
  if (returnParams){
    return(list(calls=A, confs=B, SNPQC=dev, batchQC=batchQC, params=params, DD=DD, covDD=SS, gender=gender, pkgname=pkgname))
  }else{
    return(list(calls=A, confs=B, SNPQC=dev, batchQC=batchQC, DD=DD, covDD=SS, gender=gender, pkgname=pkgname))
  }
}


gtypeCallerR <- function(A, B, fIndex, mIndex, theCenters, theScales,
                         theNs, Indexes, cIndexes, nIndexes,
                         ncIndexes, SMEDIAN, knots, params, dft,
                         probs, trim){

  stopifnot(!missing(A), !missing(B), dim(A)==dim(B),
            nrow(A)==nrow(theCenters), nrow(A)==nrow(theScales),
            nrow(A) == nrow(theNs), length(knots)==3, nrow(params)==4,
            ncol(params)==ncol(A), length(trim)==1, length(probs)==3)

  ## make code robust
  ## check types before passing to C
  
  .Call("gtypeCallerPart1", A, B,
        as.integer(fIndex), as.integer(mIndex),
        as.numeric(theCenters), as.numeric(theScales),
        as.integer(theNs), lapply(Indexes, as.integer), lapply(cIndexes, as.integer), as.integer(nIndexes), as.integer(ncIndexes),
        as.numeric(SMEDIAN), as.numeric(knots), as.numeric(params),
        as.integer(dft), as.numeric(probs), as.numeric(trim),
        PACKAGE="crlmm")
  
}

gtypeCallerR2 <- function(A, B, fIndex, mIndex, theCenters, theScales,
                         theNs, Indexes, cIndexes, nIndexes,
                         ncIndexes, SMEDIAN, knots, params, dft,
                         probs, trim, noTraining, noInfo){

  stopifnot(!missing(A), !missing(B), dim(A)==dim(B),
            nrow(A)==nrow(theCenters), nrow(A)==nrow(theScales),
            nrow(A) == nrow(theNs), length(knots)==3, nrow(params)==4,
            ncol(params)==ncol(A), length(trim)==1, length(probs)==3)

  .Call("gtypeCallerPart2", A, B,
        as.integer(fIndex), as.integer(mIndex),
        as.numeric(theCenters), as.numeric(theScales),
        as.integer(theNs), Indexes, cIndexes, nIndexes, ncIndexes,
        as.numeric(SMEDIAN), as.numeric(knots), as.numeric(params),
        as.integer(dft), as.numeric(probs), as.numeric(trim),
        as.integer(noTraining), as.integer(noInfo), PACKAGE="crlmm")
  
}
