setMethod("getOptions", "IlluminaAlleleSet", function(object){
	copynumber <- TRUE
	nonpolymorphic <- TRUE
	verbose <- TRUE
	seed <- 1
	readOpts <- readOptions.IlluminaAlleleSet()
	snprmaOpts <- snprmaOptions.AffymetrixAlleleSet()
	crlmmOpts <- crlmmOptions.AffymetrixAlleleSet()
	cnOpts <- cnOptions()
	list(verbose=verbose,
	     copynumber=copynumber,
	     nonpolymorphic=nonpolymorphic,
	     seed=seed,
	     readOpts=readOpts,
	     snprmaOpts=snprmaOpts,
	     crlmmOpts=crlmmOpts,
	     cnOpts=cnOpts)
})
readOptions.IlluminaAlleleSet <- function(sampleSheet=NULL,
					  fileExt=list(green="Grn.idat", red="Red.idat", sep="_"),
					  ids=NULL,
					  arrayInfoColNames=list(barcode="SentrixBarcode_A", position="SentrixPosition_A"),
					  highDensity=FALSE,
					  path=".",
					  stripNorm=TRUE,
					  useTarget=TRUE){
	list(sampleSheet=sampleSheet,
	     fileExt=fileExt,
	     ids=ids,
	     arrayInfoColNames=arrayInfoColNames,
	     highDensity=highDensity,
	     path=path,
	     stripNorm=stripNorm,
	     useTarget=useTarget)
}

construct.IlluminaAlleleSet <- function(object, filenames){
	message("reading first idat file to extract feature data")
	M <- .illuminaFeatureInfo(object, filenames)
	nr <- nrow(M)
	dns <- list(featureNames(M), basename(filenames))
	RG <- new("IlluminaRGSet",
		  R=initializeBigMatrix(dns),
		  G=initializeBigMatrix(dns),
		  zero=initializeBigMatrix(dns),
		  featureData=M,
		  annotation=annotation(object),
		  options=crlmmOptions(object))
	sampleNames(RG) <- basename(filenames)
	storageMode(RG) <- "environment"
	RG##featureData=ops$illuminaOpts[["featureData"]])
}
setMethod("construct", "IlluminaAlleleSet", construct.IlluminaAlleleSet)
.illuminaFeatureInfo <- function(object, filenames){
        fileExt <- crlmmOptions(object)$readOpts[["fileExt"]]
        grnfile <- paste(filenames[1], fileExt$green, sep=fileExt$sep)
        if(!file.exists(grnfile)){
                stop(paste(grnfile, " does not exist. Check fileExt argument"))
        }
        G <- readIDAT(grnfile)
        idsG = rownames(G$Quants)
        nr <- length(idsG)
	fD <- new("AnnotatedDataFrame", data=data.frame(row.names=idsG))##, varMetadata=data.frame(labelDescript
        return(fD)
}

