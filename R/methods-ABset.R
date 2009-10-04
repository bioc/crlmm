setValidity("ABset", function(object) {
	##msg <- validMsg(NULL, Biobase:::isValidVersion(object, "CopyNumberSet"))
	msg <- validMsg(NULL, assayDataValidMembers(assayData(object), c("A", "B")))
	if (is.null(msg)) TRUE else msg
})
setMethod("A", "ABset", function(object) assayData(object)[["A"]])
setMethod("B", "ABset", function(object) assayData(object)[["B"]])
setReplaceMethod("A", signature(object="ABset", value="matrix"),
		 function(object, value){
			 assayDataElementReplace(object, "A", value)			
		 })
setReplaceMethod("B", signature(object="ABset", value="matrix"),
		 function(object, value){
			 assayDataElementReplace(object, "B", value)			
		 })



