##setValidity("SnpQSet", function(object) {
##	##msg <- validMsg(NULL, Biobase:::isValidVersion(object, "CopyNumberSet"))
##	msg <- validMsg(NULL, assayDataValidMembers(assayData(object), c("A", "B")))
##	if (is.null(msg)) TRUE else msg
##})

setMethod("initialize", "SnpQSet",
          function(.Object,
                   assayData = assayDataNew(senseThetaA=senseThetaA, senseThetaB=senseThetaB),
                   senseThetaA=new("matrix"),
                   senseThetaB=new("matrix"),
                   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   experimentData=new("MIAME"),
                   annotation=new("character")){
            .Object <- callNextMethod(.Object,
				      assayData = assayDataNew(senseThetaA=senseThetaA, senseThetaB=senseThetaB),
				      phenoData=phenoData,
				      experimentData=experimentData,
				      annotation=annotation)
            .Object
    })

setValidity("SnpQSet",
            function(object)
            assayDataValidMembers(assayData(object),
                                  c("senseThetaA",
                                    "senseThetaB")
            ))
setMethod("A", "SnpQSet", function(object) senseThetaA(object))
setMethod("B", "SnpQSet", function(object) senseThetaB(object))
setReplaceMethod("A", signature(object="SnpQSet", value="matrix"),
		 function(object, value){
			 assayDataElementReplace(object, "senseThetaA", value)			
		 })
setReplaceMethod("B", signature(object="SnpQSet", value="matrix"),
		 function(object, value){
			 assayDataElementReplace(object, "senseThetaB", value)			
		 })
setMethod("isSnp", "SnpQSet", function(object) {
	isSnp.SnpQSet(object)
})


