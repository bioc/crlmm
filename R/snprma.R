snprma <- function(filenames, mixtureSampleSize=10^5, fitMixture=FALSE, eps=0.1, verbose=TRUE, seed=1, cdfName, sns){
  if (missing(sns)) sns <- basename(filenames)
  ##ADD CHECK TO SEE IF LOADED
  if (missing(cdfName))
    cdfName <- read.celfile.header(filenames[1])$cdfName
##  stuffDir <- changeToCrlmmAnnotationName(cdfName)
  pkgname <- getCrlmmAnnotationName(cdfName)
  if(!require(pkgname, character.only=TRUE, quietly=!verbose)){
    suggCall <- paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep="")
    msg <- paste("If", pkgname, "is installed on an alternative location, please load it manually by using", suggCall)
    message(strwrap(msg))
    stop("Package ", pkgname, " could not be found.")
    rm(suggCall, msg)
  }
  
  if(verbose) message("Loading annotations and mixture model parameters.")
  data(preprocStuff, genotypeStuff, mixtureStuff, package=pkgname, envir=.crlmmPkgEnv)
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

  ##We will read each cel file, summarize, and run EM one by one
  ##We will save parameters of EM to use later
  mixtureParams <- matrix(0, 4, length(filenames))
  SNR <- vector("numeric", length(filenames))
  SKW <- vector("numeric", length(filenames))

  ## This is the sample for the fitting of splines
  ## BC: I like better the idea of the user passing the seed,
  ##     because this might intefere with other analyses
  ##     (like what happened to GCRMA)
  set.seed(seed)
  
  idx <- sort(sample(autosomeIndex, mixtureSampleSize))

  ##S will hold (A+B)/2 and M will hold A-B
  ##NOTE: We actually dont need to save S. Only for pics etc...
  ##f is the correction. we save to avoid recomputing
  A <- matrix(as.integer(0), length(pnsa), length(filenames))
  B <- matrix(as.integer(0), length(pnsb), length(filenames))
  
  if(verbose){
    message("Processing ", length(filenames), " files.")
    if (getRversion() > '2.7.0') pb <- txtProgressBar(min=0, max=length(filenames), style=3)
  }
  ##We start looping throug cel files
  idx2 <- sample(length(fid), 10^5) ##for skewness. no need to do everything
  for(i in seq(along=filenames)){
    y <- as.matrix(read.celfile(filenames[i], intensity.means.only=TRUE)[["INTENSITY"]][["MEAN"]][fid])
    x <- log2(y[idx2])
    SKW[i] <- mean((x-mean(x))^3)/(sd(x)^3)
    rm(x)
    y <- normalize.quantiles.use.target(y, target=reference)
    A[, i] <- intMedianSummaries(y[aIndex, 1, drop=FALSE], pnsa)
    B[, i] <- intMedianSummaries(y[bIndex, 1, drop=FALSE], pnsb)
    
    ##Now to fit the EM
    if(fitMixture){
      S <- (log2(A[idx, i])+log2(B[idx, i]))/2 - SMEDIAN
      M <- log2(A[idx, i])-log2(B[idx, i])
      
      ##we need to test the choice of eps.. it is not the max diff between funcs
      tmp <- fitAffySnpMixture56(S, M, theKnots, eps=eps)
      
      mixtureParams[, i] <- tmp[["coef"]]
      SNR[i] <- tmp[["medF1"]]^2/(tmp[["sigma1"]]^2+tmp[["sigma2"]]^2)
    }
    if (verbose)
      if (getRversion() > '2.7.0') setTxtProgressBar(pb, i)
      else cat(".")
  }
  if (verbose)
    if (getRversion() > '2.7.0') close(pb)
    else cat("\n")
  if (!fitMixture) SNR <- mixtureParams <- NA
  ## gns comes from preprocStuff.rda
  list(A=A, B=B, sns=sns, gns=gns, SNR=SNR, SKW=SKW, mixtureParams=mixtureParams, cdfName=cdfName)
}
