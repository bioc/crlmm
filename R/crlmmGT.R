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
  if(!require(pkgname, character.only=TRUE)){
    suggCall <- paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep="")
    msg <- paste("If", pkgname, "is installed on an alternative location, please load it manually by using", suggCall)
    message(strwrap(msg))
    stop("Package ", pkgname, " could not be found.")
    rm(suggCall, msg)
  }

  if(verbose) message("Loading annotations.")
  data(genotypeStuff, mixtureStuff, package=pkgname, envir=.crlmmPkgEnv)

  ## this is toget rid of the 'no visible binding' notes
  ## variable definitions
  XIndex <- getVarInEnv("XIndex")
  autosomeIndex <- getVarInEnv("autosomeIndex")
  YIndex <- getVarInEnv("YIndex")
  SMEDIAN <- getVarInEnv("SMEDIAN")
  theKnots <- getVarInEnv("theKnots")
  regionInfo <- getVarInEnv("regionInfo")
  
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
  
  if(verbose) cat("Calling", NR, "SNPs for recalibration")

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
    warning("Recallibration not possible.")
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
    newparams[["centers"]][is.na(newparams[["centers"]])] <- params[["centers"]][is.na(newparams[["centers"]])]
    if(verbose) cat("\n")
  
    if(verbose) message("Calculating and standardizing size of shift.")
    DD <- newparams[["centers"]] - params[["centers"]]
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
  if(verbose) cat("Calling", NR, "SNPs")
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
    batchQC <- Inf
  }
  
  if(verbose) message("Done.")
  if (returnParams){
    return(list(calls=A, confs=B, SNPQC=dev, batchQC=batchQC, params=params, DD=DD))
  }else{
    return(list(calls=A, confs=B, SNPQC=dev, batchQC=batchQC, DD=DD))
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





##################
##################
### THIS IS TEMPORARY NOT OFFICIALLY USED
##################
####################
crlmmGTTNoN <- function(A, B, SNR, mixtureParams, cdfName,
                         row.names=NULL, col.names=NULL, probs=c(1/3,
                         1/3, 1/3), DF=6, SNRMin=6, gender=NULL,
                         desctrucitve=FALSE, verbose=TRUE){
  keepIndex <- which(SNR>SNRMin)
  if(length(keepIndex)==0) stop("No arrays above quality threshold!")
  
  NC <- ncol(A)
  NR <- nrow(A)
  
  if(verbose) message("Loading annotations.")
  data(genotypeStuff, mixtureStuff, package=pkgname, envir=.crlmmPkgEnv)

  ## this is toget rid of the 'no visible binding' notes
  ## variable definitions
  XIndex <- getVarInEnv("XIndex")
  autosomeIndex <- getVarInEnv("autosomeIndex")
  YIndex <- getVarInEnv("YIndex")
  SMEDIAN <- getVarInEnv("SMEDIAN")
  theKnots <- getVarInEnv("theKnots")
  regionInfo <- getVarInEnv("regionInfo")
  
  ##IF gender not provide, we predict
  if(is.null(gender)){
    if(verbose) message("Determining gender.")
    XMedian <- apply(log2(A[XIndex,, drop=FALSE])+log2(B[XIndex,, drop=FALSE]), 2, median)/2
    if(sum(SNR>SNRMin)==1) gender <- which.min(c(abs(XMedian-8.9), abs(XMedian-9.5))) else  gender <- kmeans(XMedian, c(min(XMedian[SNR>SNRMin]), max(XMedian[SNR>SNRMin])))[["cluster"]]
  }
  
  Indexes <- list(autosomeIndex, XIndex, YIndex)
  cIndexes <- list(keepIndex, 
                   keepIndex[which(gender[keepIndex]==2)], 
                   keepIndex[which(gender[keepIndex]==1)])
  
  if(verbose) cat("Calling", NR, "SNPs for recalibration")

  ## call C
  fIndex <- which(gender==2)
  mIndex <- which(gender==1)
  t0 <- proc.time()
  newparams <- gtypeCallerRTNoN(A, B, fIndex, mIndex,
                            params[["centers"]], params[["scales"]], params[["N"]],
                            Indexes, cIndexes,
                            sapply(Indexes, length), sapply(cIndexes, length),
                            SMEDIAN, theKnots,
                            mixtureParams, DF, probs, 0.025)
  t0 <- proc.time()-t0
  message("Part 1 took ", t0[3], " seconds.")
  names(newparams) <- c("centers", "scales", "N")
  
  if(verbose) message("Done.")
  if(verbose) message("Estimating recalibration parameters.")
  d <- newparams[["centers"]] - params$centers

  ##regression 
  MIN <- 10
  Index <- intersect(which(pmin(newparams[["N"]][, 1], newparams[["N"]][, 2], newparams[["N"]][, 3])>MIN & !apply(regionInfo, 1, any)), autosomeIndex)
  data4reg <- as.data.frame(newparams[["centers"]][Index,])
  names(data4reg) <- c("AA", "AB", "BB")
  regParams <- cbind(  coef(lm(AA~AB*BB, data=data4reg)),
                     c(coef(lm(AB~AA+BB, data=data4reg)), 0), 
                       coef(lm(BB~AA*AB, data=data4reg)))
  rownames(regParams) <- c("intercept", "X", "Y", "XY")
  rm(data4reg)
  
  minN <- 3
  newparams[["centers"]][newparams[["N"]]<minN] <- NA
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
  newparams[["centers"]][is.na(newparams[["centers"]])] <- params[["centers"]][is.na(newparams[["centers"]])]
  if(verbose) cat("\n")
  
  if(verbose) message("Calculating and standardizing size of shift.")
  DD <- newparams[["centers"]] - params[["centers"]]
  MM <- colMeans(DD[autosomeIndex, ])
  DD <- sweep(DD, 2, MM)
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

  ## BC: must keep SD
  params[-2] <- newparams[-2]
  rm(newparams);gc(verbose=FALSE)  
  if(verbose) cat("Calling", NR, "SNPs")
  ## ###################
  ## ## MOVE TO C#######
  t0 <- proc.time()
  ImNull <- gtypeCallerR2TNoN(A, B, fIndex, mIndex, params[["centers"]],
                          params[["scales"]], params[["N"]], Indexes,
                          cIndexes, sapply(Indexes, length),
                          sapply(cIndexes, length), SMEDIAN, theKnots,
                          mixtureParams, DF, probs, 0.025,
                          which(regionInfo[,2]),
                          which(regionInfo[,1]))
  t0 <- proc.time()-t0
  message("Part 2 took ", t0[3], " seconds.")
  ##  END MOVE TO C#######
  ## ##################
  
  dev <- dev/(dev+1/383)
  if(!is.null(row.names)){ rownames(A) <- rownames(B) <- names(dev) <- row.names}
  if(!is.null(col.names)){ colnames(A) <- colnames(B) <- col.names}
  
  if(verbose) message("Done.")
  return(list(calls=A, confs=B, SNPQC=dev, batchQC=mean(diag(SS))))
}


gtypeCallerRTNoN <- function(A, B, fIndex, mIndex, theCenters, theScales,
                         theNs, Indexes, cIndexes, nIndexes,
                         ncIndexes, SMEDIAN, knots, params, dft,
                         probs, trim){

  stopifnot(!missing(A), !missing(B), dim(A)==dim(B),
            nrow(A)==nrow(theCenters), nrow(A)==nrow(theScales),
            nrow(A) == nrow(theNs), length(knots)==3, nrow(params)==4,
            ncol(params)==ncol(A), length(trim)==1, length(probs)==3)
  
  .Call("gtypeCallerPart1TNoN", A, B, fIndex, mIndex, theCenters,
        theScales, theNs, Indexes, cIndexes, nIndexes, ncIndexes,
        SMEDIAN, knots, params, dft, probs, trim, PACKAGE="crlmm")
  
}

gtypeCallerR2TNoN <- function(A, B, fIndex, mIndex, theCenters, theScales,
                         theNs, Indexes, cIndexes, nIndexes,
                         ncIndexes, SMEDIAN, knots, params, dft,
                         probs, trim, noTraining, noInfo){

  stopifnot(!missing(A), !missing(B), dim(A)==dim(B),
            nrow(A)==nrow(theCenters), nrow(A)==nrow(theScales),
            nrow(A) == nrow(theNs), length(knots)==3, nrow(params)==4,
            ncol(params)==ncol(A), length(trim)==1, length(probs)==3)
  
  .Call("gtypeCallerPart2TNoN", A, B, fIndex, mIndex, theCenters,
        theScales, theNs, Indexes, cIndexes, nIndexes, ncIndexes,
        SMEDIAN, knots, params, dft, probs, trim, noTraining, noInfo, PACKAGE="crlmm")
  
}

crlmmGTNormalNoN <- function(A, B, SNR, mixtureParams, cdfName,
                         row.names=NULL, col.names=NULL, probs=c(1/3,
                         1/3, 1/3), DF=6, SNRMin=6, gender=NULL,
                         desctrucitve=FALSE, verbose=TRUE){
  keepIndex <- which(SNR>SNRMin)
  if(length(keepIndex)==0) stop("No arrays above quality threshold!")
  
  NC <- ncol(A)
  NR <- nrow(A)
  
  if(verbose) message("Loading annotations.")
  data(genotypeStuff, mixtureStuff, package=pkgname, envir=.crlmmPkgEnv)

  ## this is toget rid of the 'no visible binding' notes
  ## variable definitions
  XIndex <- getVarInEnv("XIndex")
  autosomeIndex <- getVarInEnv("autosomeIndex")
  YIndex <- getVarInEnv("YIndex")
  SMEDIAN <- getVarInEnv("SMEDIAN")
  theKnots <- getVarInEnv("theKnots")
  regionInfo <- getVarInEnv("regionInfo")
  
  ##IF gender not provide, we predict
  if(is.null(gender)){
    if(verbose) message("Determining gender.")
    XMedian <- apply(log2(A[XIndex,, drop=FALSE])+log2(B[XIndex,, drop=FALSE]), 2, median)/2
    if(sum(SNR>SNRMin)==1) gender <- which.min(c(abs(XMedian-8.9), abs(XMedian-9.5))) else  gender <- kmeans(XMedian, c(min(XMedian[SNR>SNRMin]), max(XMedian[SNR>SNRMin])))[["cluster"]]
  }
  
  Indexes <- list(autosomeIndex, XIndex, YIndex)
  cIndexes <- list(keepIndex, 
                   keepIndex[which(gender[keepIndex]==2)], 
                   keepIndex[which(gender[keepIndex]==1)])
  
  if(verbose) cat("Calling", NR, "SNPs for recalibration")

  ## call C
  fIndex <- which(gender==2)
  mIndex <- which(gender==1)
  t0 <- proc.time()
  newparams <- gtypeCallerRNormalNoN(A, B, fIndex, mIndex,
                            params[["centers"]], params[["scales"]], params[["N"]],
                            Indexes, cIndexes,
                            sapply(Indexes, length), sapply(cIndexes, length),
                            SMEDIAN, theKnots,
                            mixtureParams, DF, probs, 0.025)
  t0 <- proc.time()-t0
  message("Part 1 took ", t0[3], " seconds.")
  names(newparams) <- c("centers", "scales", "N")
  
  if(verbose) message("Done.")
  if(verbose) message("Estimating recalibration parameters.")
  d <- newparams[["centers"]] - params$centers

  ##regression 
  MIN <- 10
  Index <- intersect(which(pmin(newparams[["N"]][, 1], newparams[["N"]][, 2], newparams[["N"]][, 3])>MIN & !apply(regionInfo, 1, any)), autosomeIndex)
  data4reg <- as.data.frame(newparams[["centers"]][Index,])
  names(data4reg) <- c("AA", "AB", "BB")
  regParams <- cbind(  coef(lm(AA~AB*BB, data=data4reg)),
                     c(coef(lm(AB~AA+BB, data=data4reg)), 0), 
                       coef(lm(BB~AA*AB, data=data4reg)))
  rownames(regParams) <- c("intercept", "X", "Y", "XY")
  rm(data4reg)
  
  minN <- 3
  newparams[["centers"]][newparams[["N"]]<minN] <- NA
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
  newparams[["centers"]][is.na(newparams[["centers"]])] <- params[["centers"]][is.na(newparams[["centers"]])]
  if(verbose) cat("\n")
  
  if(verbose) message("Calculating and standardizing size of shift.")
  DD <- newparams[["centers"]] - params[["centers"]]
  MM <- colMeans(DD[autosomeIndex, ])
  DD <- sweep(DD, 2, MM)
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

  ## BC: must keep SD
  params[-2] <- newparams[-2]
  rm(newparams);gc(verbose=FALSE)  
  if(verbose) cat("Calling", NR, "SNPs")
  ## ###################
  ## ## MOVE TO C#######
  t0 <- proc.time()
  ImNull <- gtypeCallerR2NormalNoN(A, B, fIndex, mIndex, params[["centers"]],
                          params[["scales"]], params[["N"]], Indexes,
                          cIndexes, sapply(Indexes, length),
                          sapply(cIndexes, length), SMEDIAN, theKnots,
                          mixtureParams, DF, probs, 0.025,
                          which(regionInfo[,2]),
                          which(regionInfo[,1]))
  t0 <- proc.time()-t0
  message("Part 2 took ", t0[3], " seconds.")
  ##  END MOVE TO C#######
  ## ##################
  
  dev <- dev/(dev+1/383)
  if(!is.null(row.names)){ rownames(A) <- rownames(B) <- names(dev) <- row.names}
  if(!is.null(col.names)){ colnames(A) <- colnames(B) <- col.names}
  
  if(verbose) message("Done.")
  return(list(calls=A, confs=B, SNPQC=dev, batchQC=mean(diag(SS))))
}


gtypeCallerRNormalNoN <- function(A, B, fIndex, mIndex, theCenters, theScales,
                         theNs, Indexes, cIndexes, nIndexes,
                         ncIndexes, SMEDIAN, knots, params, dft,
                         probs, trim){

  stopifnot(!missing(A), !missing(B), dim(A)==dim(B),
            nrow(A)==nrow(theCenters), nrow(A)==nrow(theScales),
            nrow(A) == nrow(theNs), length(knots)==3, nrow(params)==4,
            ncol(params)==ncol(A), length(trim)==1, length(probs)==3)
  
  .Call("gtypeCallerPart1NormalNoN", A, B, fIndex, mIndex, theCenters,
        theScales, theNs, Indexes, cIndexes, nIndexes, ncIndexes,
        SMEDIAN, knots, params, dft, probs, trim, PACKAGE="crlmm")
  
}

gtypeCallerR2NormalNoN <- function(A, B, fIndex, mIndex, theCenters, theScales,
                         theNs, Indexes, cIndexes, nIndexes,
                         ncIndexes, SMEDIAN, knots, params, dft,
                         probs, trim, noTraining, noInfo){

  stopifnot(!missing(A), !missing(B), dim(A)==dim(B),
            nrow(A)==nrow(theCenters), nrow(A)==nrow(theScales),
            nrow(A) == nrow(theNs), length(knots)==3, nrow(params)==4,
            ncol(params)==ncol(A), length(trim)==1, length(probs)==3)
  
  .Call("gtypeCallerPart2NormalNoN", A, B, fIndex, mIndex, theCenters,
        theScales, theNs, Indexes, cIndexes, nIndexes, ncIndexes,
        SMEDIAN, knots, params, dft, probs, trim, noTraining, noInfo, PACKAGE="crlmm")
  
}


