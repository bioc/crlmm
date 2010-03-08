setAs("AffymetrixAlleleSet", "CallSet", function(from, to){
	new("CallSet",
	    alleleA=A(from),
	    alleleB=B(from),		  
	    call=initializeBigMatrix(list(featureNames(from), sampleNames(from))),
	    callProbability=initializeBigMatrix(list(featureNames(from), sampleNames(from))),
	    phenoData=phenoData(from),
	    experimentData=experimentData(from),
	    annotation=annotation(from),
	    protocolData=protocolData(from),
	    featureData=featureData(from),
	    genomeAnnotation=genomeAnnotation(from),
	    options=crlmmOptions(from))
})

setAs("IlluminaAlleleSet", "CallSet", function(from, to){
	new("CallSet",
	    alleleA=A(from),
	    alleleB=B(from),		  
	    call=initializeBigMatrix(list(featureNames(from), sampleNames(from))),
	    callProbability=initializeBigMatrix(list(featureNames(from), sampleNames(from))),
	    phenoData=phenoData(from),
	    experimentData=experimentData(from),
	    annotation=annotation(from),
	    protocolData=protocolData(from),
	    featureData=featureData(from),
	    genomeAnnotation=genomeAnnotation(from),
	    options=crlmmOptions(from))
})



setMethod("allele", "AlleleSet",
          function(object, allele, strand){
            stopifnot(!missing(allele))
            allele <- match.arg(allele, c("A", "B"))
            both <- bothStrands(object)
            if (!both){
              what <- paste("allele", allele, sep="")
            }else{
              stopifnot(!missing(strand))
              strand <- match.arg(strand, c("sense", "antisense"))
              what <- paste(strand, "Allele", allele, sep="")
            }
            assayDataElement(object, what)
          })
setMethod("A", "AlleleSet", function(object) allele(object, "A"))
setMethod("B", "AlleleSet", function(object) allele(object, "B"))

setReplaceMethod("A", c("AlleleSet", "matrix_or_ff"),
		 function(object, value){
			 assayDataElementReplace(object, "alleleA", value)
		 })
setReplaceMethod("B", c("AlleleSet", "matrix_or_ff"),
		 function(object, value){
			 assayDataElementReplace(object, "alleleB", value)
		 })

##setMethod("rma", "AffymetrixAlleleSet", function(object){
##	object <- snprma(object)
##	ops <- crlmmOptions(object)
##	if(ops$copynumber & ops$nonpolymorphic)
##		object <- cnrma(object)
##	return(object)
##})

##setMethod("rma", "SmallDataSet", function(ops){
##	RG <- new("NChannelSet",
##		  R=getR(ops),
##		  G=getG(ops),
##		  X=getR(ops),##doesn't really matter
##		  Y=getR(ops),##doesn't really matter
##		  zero=getZero(ops),
##		  phenoData=getPhenoData(ops),
##		  featureData=ops$illuminaOpts[["featureData"]])
####	XY <- new("NChannelSet",
####		  X=tmpmat,
####		  Y=tmpmat,
####		  zero=tmpmat, # Xnb=tmpmat, Ynb=tmpmat, Xse=tmpmat, Yse=tmpmat, zero=tmpmat,
####		  annotation=chipType,
####		  phenoData=phenoData(RG),
####		  protocolData=protocolData(RG),
####		  storage.mode="environment")
##	storageMode(RG) <- "environment"
##	RG <- readIdatFiles(RG, ops)
##
##	XY <- RGtoXY(RG, chipType=ops$cdfName)
##	object <- preprocessInfinium2(XY)
##	return(object)
##})




