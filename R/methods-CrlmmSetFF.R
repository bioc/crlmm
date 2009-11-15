##setAs("SnpCallSetPlus", "CrlmmSetFF",
##	  function(from, to){
##		  CA <- CB <- matrix(NA, nrow(from), ncol(from))
##		  dimnames(CA) <- dimnames(CB) <- list(featureNames(from), sampleNames(from))		  
##		  new("CrlmmSetFF",
##		      calls=calls(from),
##		      callsConfidence=callsConfidence(from),
##		      senseThetaA=A(from),
##		      senseThetaB=B(from),
##		      CA=CA,
##		      CB=CB,
##		      phenoData=phenoData(from),
##		      experimentData=experimentData(from),
##		      annotation=annotation(from),
##		      protocolData=protocolData(from),
##		      featureData=featureData(from))
##	  })

##setMethod("[", "CrlmmSetFF", function(x, i, j, ..., drop = FALSE) {
	## does NOT subset any of the assayDataElements
	## does subset the featureData, phenoData, protocolData

	## accessors to the assayData work as follows
	## object2 <- object[1:10, ]
	##
	## A(object2)[1:10, ]
	##  - i.  A(object2) gets the file
	##  - ii. A is subset by the dimNames information as follows
	##        dimNames[[1]] %in% featureNames(object2)
	##        dimNames[[2]] %in% sampleNames(object2)
	##    iii. first 10 elements are retrieved
##})

##setMethod("A", "CrlmmSetFF", function(object) assayData(object)[["senseThetaA"]])
##setMethod("B", "CrlmmSetFF", function(object) assayData(object)[["senseThetaB"]])
##setMethod("CA", "CrlmmSetFF", function(object) assayData(object)[["CA"]]/100)
##setMethod("CB", "CrlmmSetFF", function(object) assayData(object)[["CB"]]/100)


