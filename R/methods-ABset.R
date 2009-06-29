setValidity("ABset", function(object) {
	##msg <- validMsg(NULL, Biobase:::isValidVersion(object, "CopyNumberSet"))
	msg <- validMsg(NULL, assayDataValidMembers(assayData(object), c("A", "B")))
	if (is.null(msg)) TRUE else msg
})
setMethod("A", "ABset", function(object) assayData(object)[["A"]])
setMethod("B", "ABset", function(object) assayData(object)[["B"]])




