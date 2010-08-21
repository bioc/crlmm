setMethod("cnIndex", "eSet", function(object){
	index <- match(cnNames(object), featureNames(object), nomatch=0)
	index[index != 0]
  })

setMethod("cnNames", "eSet", function(object) {
	path <- system.file("extdata", package=paste(annotation(object), "Crlmm", sep=""))
	load(file.path(path, "cnProbes.rda"))
	cnProbes <- get("cnProbes")
	cnps <- rownames(cnProbes)
	cnps <- cnps[cnps %in% featureNames(object)]
	index <- match(cnps, featureNames(object), nomatch=0)
	index <- index[index != 0]	
	featureNames(object)[index]	
  })


setMethod("snpIndex", "eSet", function(object){
	index <- match(snpNames(object), featureNames(object), nomatch=0)
	index[index != 0]
})

setMethod("snpNames", "eSet", function(object){
	path <- system.file("extdata", package=paste(annotation(object), "Crlmm", sep=""))
	load(file.path(path, "snpProbes.rda"))
	snpProbes <- get("snpProbes")
	snps <- rownames(snpProbes)
	snps <- snps[snps %in% featureNames(object)]
	index <- match(snps, featureNames(object), nomatch=0)
	index <- index[index != 0]
	featureNames(object)[index]
})

setMethod("chromosome", "eSet", function(object) fData(object)$chromosome)
setMethod("position", "eSet", function(object) fData(object)$position)

##setMethod("getParam", signature(object="eSet",
##				name="character",
##				batch="ANY"),
##	  function(object, name, batch){
##		  if(length(batch) > 1){
##			  warning("batch argument to getParam should have length 1.  Only using the first")
##			  batch <- batch[1]
##		  }
##		  getParam.SnpSuperSet(object, name, batch)
##})

##setMethod("batch", "eSet", function(object){
##	if("batch" %in% varLabels(object)){
##		return(object$batch)
##	} else {
##		return(protocolData(object)$batch)
##	}
##})

##setMethod("combine", signature=signature(x="eSet", y="eSet"),
##	  function(x, y, ...){
##		  ##Check that both x and y are valid objects
##		  if(!validObject(x)) stop("x is not a valid object")
##		  if(!validObject(y)) stop("y is not a valid object")
##
##		  if(annotation(x) != annotation(y)) stop("must have same value for annotation slot")
##		  if(class(x) != class(y)) stop("objects must have the same class")
##		  if(storageMode(assayData(x)) != storageMode(assayData(y))){
##			  stop("objects must have same storage mode for assayData")
##		  }
##		  fd <- combine(featureData(x), featureData(y))
##		  pd <- combine(phenoData(x), phenoData(y))            
##		  ad.x <- as.list(assayData(x))
##		  ad.y <- as.list(assayData(y))
##		  ad.xy <- mapply(rbind, ad.x, ad.y, SIMPLIFY=FALSE)
##		  id.x <- match(rownames(ad.xy[[1]]), featureNames(fd))
##		  ee <- combine(experimentData(x), experimentData(y))
##		  assayData(x) <- ad.xy
##		  storageMode(assayData(x)) <- storageMode(assayData(y))            
##		  experimentData(x) <- ee
##		  featureData(x) <- fd
##		  phenoData(x) <- pd
##		  x
##          })



##setMethod("pr", signature(object="eSet",
##			  name="character",
##			  batch="ANY",
##			  value="numeric"), 
##	  function(object, name, batch, value){
##		  label <- paste(name, batch, sep="_")
##		  colindex <- match(label, fvarLabels(object))
##		  if(length(colindex) == 1){
##			  fData(object)[, colindex] <- value
##		  } 
##		  if(is.na(colindex)){
##			  stop(paste(label, " not found in object"))
##		  }
##		  if(length(colindex) > 1){
##			  stop(paste(label, " not unique"))
##		  }
##		  object
##	  })
