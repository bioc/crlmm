setMethod("crlmm", "AffymetrixAlleleSet", function(object, filenames){
	message("Initializing AlleleSet")
	obj <- construct(object, filenames)
	obj <- snprma(obj, filenames)
	ops <- crlmmOptions(obj)
	if(ops$copynumber & ops$nonpolymorphic)
		obj <- cnrma(obj, filenames)
	message("Initializing CallSet")
	object <- as(obj, "CallSet")
	rm(obj); gc()
	callSet <- crlmm.batch(object)

})
setMethod("crlmm", "IlluminaAlleleSet", function(object, filenames){
	RG <- construct(object, filenames)
	RG <- readIdatFiles(RG, filenames)
	XY <- RGtoXY(RG)
	rm(RG); gc()
	storageMode(XY) <- "environment"
	stripNorm <- crlmmOptions(XY)$readOpts[["stripNorm"]]
	useTarget <- crlmmOptions(XY)$readOpts[["useTarget"]]
	verbose <- crlmmOptions(XY)$verbose
	if(stripNorm)
		XY = stripNormalize(XY, useTarget=useTarget, verbose=verbose)
	## See the coercion step ... this puts X -> A and Y -> B
	alleleSet <- as(XY, "IlluminaAlleleSet")
	##Do i need to pass zero.  Is the only thing updated the mixtureParams and SNR?
	alleleSet <- preprocessInfinium2(alleleSet, zero=Z(XY))
	rm(XY); gc()
	##Initialize call set
	callSet <- as(alleleSet, "CallSet")
	callSet <- crlmm.batch(callSet)
	return(callSet)
})

##setMethod("getFeatureData", "SmallDataOptions", function(object){
##        if(object$platform != "Illumina") stop("only for Illumina platforms")
##        message("reading first idat file to extract feature data")
##        fileExt <- object$illuminaOpts[["fileExt"]]
##        filenames <- object$filenames
##        grnfile = paste(filenames[1], fileExt$green, sep=fileExt$sep)
##        if(!file.exists(grnfile)){
##                stop(paste(grnfile, " does not exist. Check fileExt argument"))
##        }
##        G <- readIDAT(grnfile)
##        idsG = rownames(G$Quants)
##        nr <- length(idsG)
##        fD <- new("AnnotatedDataFrame", data=data.frame(row.names=idsG))##, varMetadata=data.frame(labelDescript
##        return(fD)
##})

##crlmm.IlluminaAlleleSet <- function(object, filenames){
##	message("Initializing CallSet")
##	object <- as(object, "IlluminaCallSet")
##	rm(obj); gc()
##	##Call in batches to reduce ram
##	nc <- ncol(object)
##	object$gender <- crlmmOptions(object)$crlmmOpts[["gender"]]
##	if(is.null(object$gender)) object$gender <- rep(NA, nc)
##	cOps <- crlmmOptions(object)$crlmmOpts
##	BS <- cOps$batchSize
##	gc()
##	if(nc > BS){
##		N <- ceiling(nc/BS)
##		S <- ceiling(nc/N)
##		colindex <- split(1:nc, rep(1:nc, each=S, length.out=nc))
##	} else {
##		colindex <- list(1:nc)
##	}
##	if(length(colindex) > 1)
##		message("Calling genotypes in batches of size ", length(colindex[[1]]), " to reduce required RAM")
##	row.index <- which(isSnp(object)==1 | is.na(isSnp(object)))
##	if(length(colindex) > 1){
##		for(i in seq(along=colindex)){
##			col.index <- colindex[[i]]
##			tmp <- crlmm.AffymetrixCallSet(object[row.index, col.index])
##			calls(object)[row.index, col.index] <- calls(tmp)
##			confs(object)[row.index, col.index] <- confs(tmp)
##			object$gender[col.index] <- tmp$gender
##			rm(tmp); gc()
##		}
##	} else {
##		col.index <- colindex[[1]]
##		tmp <- crlmm.AffymetrixCallSet(object[row.index, col.index])  ##ensure matrix is passed
##		## enables writing to file
##		calls(object)[row.index, col.index] <- calls(tmp)
##		confs(object)[row.index, col.index] <- confs(tmp)
##		object$gender <- tmp$gender
##	}
##	##Return a CallSet object
##	class(object) <- "CallSet"
##	return(object)	
##}


crlmm.batch <- function(object){
	##Call in batches to reduce ram
	nc <- ncol(object)
	object$gender <- crlmmOptions(object)$crlmmOpts[["gender"]]
	cOps <- crlmmOptions(object)$crlmmOpts
	BS <- cOps$batchSize
	gc()
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
	if(length(colindex) > 1){
		for(i in seq(along=colindex)){
			col.index <- colindex[[i]]
			tmp <- crlmm.CallSet(object[row.index, col.index])
			calls(object)[row.index, col.index] <- calls(tmp)
			confs(object)[row.index, col.index] <- confs(tmp)
			object$gender[col.index] <- tmp$gender
			rm(tmp); gc()
		}
	} else {
		col.index <- colindex[[1]]
		tmp <- crlmm.CallSet(object[row.index, col.index])  ##ensure matrix is passed
		## enables writing to file
		calls(object)[row.index, col.index] <- calls(tmp)
		confs(object)[row.index, col.index] <- confs(tmp)
		object$gender <- tmp$gender
	}
	##Return a CallSet object
	class(object) <- "CallSet"
	return(object)
}

crlmm.CallSet <- function(object){
##		  filenames, row.names=TRUE, col.names=TRUE,
##                  probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5, gender=NULL,
##                  save.it=FALSE, load.it=FALSE,
##		  snprmaFile="snp_rmaResult.rda",
##		  callsFile="genotypes.rda",
##		  confsFile="confs.rda",
##		  AFile="A.rda",
##		  BFile="B.rda",
##                  mixtureSampleSize=10^5, eps=0.1, verbose=TRUE,
##                  cdfName, sns, recallMin=10, recallRegMin=1000,
##                  returnParams=FALSE, badSNP=0.7,
##		  crlmmBatchSize=1000){
	ops <- crlmmOptions(object)
	cOps <- ops$crlmmOpts
	cdfName <- annotation(object)
	verbose <- ops$verbose
	probs <- cOps$probs
	mixtureParams <- cOps$mixtureParams
	DF <- cOps$DF
	SNRMin <- cOps$SNRMin
	recallMin <- cOps$recallMin
	recallRegMin <- cOps$recallRegMin
	badSNP <- cOps$badSNP
	gender <- object$gender
	res2 <- tryCatch(crlmmGT(A=as.matrix(A(object)),
				 B=as.matrix(B(object)),
				 SNR=object$SNR,
				 mixtureParams=mixtureParams,
				 cdfName=cdfName,
				 row.names=featureNames(object),
				 col.names=sampleNames(object),
				 gender=gender,
				 recallMin=recallMin,
				 recallRegMin=recallRegMin,
				 SNRMin=SNRMin,
				 returnParams=TRUE,
				 badSNP=badSNP,
				 verbose=verbose), error=function(e) NULL)
	if(is.null(res2)) stop("error in call to crlmmGT")
	stopifnot(identical(featureNames(object), rownames(res2[["calls"]])))
	stopifnot(identical(sampleNames(object), colnames(res2[["calls"]])))
	calls(object) <- res2[["calls"]]
	##P <- res2[["confs"]]
	##X <- -1000*log(1-P)
	confs(object) <- res2[["confs"]]
	object$gender <- res2[["gender"]]
	rm(res2); gc()
	return(object)
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
