crlmm <- function(filenames, row.names=TRUE, col.names=TRUE,
                  probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5, gender=NULL,
                  save.it=FALSE, load.it=FALSE, intensityFile,
                  mixtureSampleSize=10^5, eps=0.1, verbose=TRUE,
                  cdfName, sns, recallMin=10, recallRegMin=1000,
                  returnParams=FALSE, badSNP=.7){
  if ((load.it | save.it) & missing(intensityFile))
    stop("'intensityFile' is missing, and you chose either load.it or save.it")
  
  if (missing(sns)) sns <- basename(filenames)
  if (!missing(intensityFile))
    if (load.it & !file.exists(intensityFile)){
      load.it <- FALSE
      message("File ", intensityFile, " does not exist.")
      message("Not loading it, but running SNPRMA from scratch.")
    }
  if (!load.it){
    res <- snprma(filenames, fitMixture=TRUE,
                  mixtureSampleSize=mixtureSampleSize, verbose=verbose,
                  eps=eps, cdfName=cdfName, sns=sns)
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

  res2 <- crlmmGT(res[["A"]], res[["B"]], res[["SNR"]],
                  res[["mixtureParams"]], res[["cdfName"]],
                  gender=gender, row.names=row.names,
                  col.names=col.names, recallMin=recallMin,
                  recallRegMin=1000, SNRMin=SNRMin,
                  returnParams=returnParams, badSNP=badSNP,
                  verbose=verbose)

  res2[["SNR"]] <- res[["SNR"]]
  res2[["SKW"]] <- res[["SKW"]]
  return(res2)
}


###############################
####### THIS IS TEMPORARY NOT OFFICIALLY USED
#####################################

crlmmTNoN <- function(filenames, row.names=TRUE, col.names=TRUE,
                  probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5, gender=NULL,
                  save.it=FALSE, load.it=FALSE,
                  intensityFile="tmpcrlmmintensities.rda",
                  desctrucitve=FALSE, mixtureSampleSize=10^5, eps=0.1,
                  verbose=TRUE){
  if (load.it & !file.exists(intensityFile)){
    load.it <- FALSE
    message("File ", intensityFile, " does not exist.")
    message("Not loading it, but running SNPRMA from scratch.")
  }
  if (!load.it){
    res <- snprma(filenames, fitMixture=TRUE,
                  mixtureSampleSize=mixtureSampleSize, verbose=verbose,
                  eps=eps)
    if(save.it) save(res, file=intensityFile)
  }else{
    message("Loading ", intensityFile, ".")
    obj <- load(intensityFile)
    message("Done.")
    if (obj != "res")
      stop("Object in ", intensityFile, " seems to be invalid.")
  }
  if(row.names) row.names=res$gns else row.names=NULL
  if(col.names) col.names=res$sns else col.names=NULL
  res2 <- crlmmGTTNoN(res[["A"]], res[["B"]], res[["SNR"]],
                  res[["mixtureParams"]], res[["cdfName"]],
                  gender=gender, row.names=row.names,
                  col.names=col.names)
  res2[["SNR"]] <- res[["SNR"]]
  return(res2)
}

crlmmNormalNoN <- function(filenames, row.names=TRUE, col.names=TRUE,
                  probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5, gender=NULL,
                  save.it=FALSE, load.it=FALSE,
                  intensityFile="tmpcrlmmintensities.rda",
                  desctrucitve=FALSE, mixtureSampleSize=10^5, eps=0.1,
                  verbose=TRUE){
  if (load.it & !file.exists(intensityFile)){
    load.it <- FALSE
    message("File ", intensityFile, " does not exist.")
    message("Not loading it, but running SNPRMA from scratch.")
  }
  if (!load.it){
    res <- snprma(filenames, fitMixture=TRUE,
                  mixtureSampleSize=mixtureSampleSize, verbose=verbose,
                  eps=eps)
    if(save.it) save(res, file=intensityFile)
  }else{
    message("Loading ", intensityFile, ".")
    obj <- load(intensityFile)
    message("Done.")
    if (obj != "res")
      stop("Object in ", intensityFile, " seems to be invalid.")
  }
  if(row.names) row.names=res$gns else row.names=NULL
  if(col.names) col.names=res$sns else col.names=NULL
  res2 <- crlmmGTNormalNoN(res[["A"]], res[["B"]], res[["SNR"]],
                  res[["mixtureParams"]], res[["cdfName"]],
                  gender=gender, row.names=row.names,
                  col.names=col.names)
  res2[["SNR"]] <- res[["SNR"]]
  return(res2)
}
