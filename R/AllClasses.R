setOldClass("ffdf")
setOldClass("ff_matrix")
setClassUnion("matrix_or_ff", c("matrix", "ff_matrix"))
setClassUnion("list_or_ffdf", c("list", "ffdf"))
setClass("CrlmmContainer", contains="eSet",
	 representation(options="list",
			genomeAnnotation="ANY",
			"VIRTUAL"))

setMethod("show", "CrlmmContainer", function(object){
	callNextMethod(object)
	cat("options: \n")
	print(names(crlmmOptions(object)))
	cat("\n")
	cat("genomeAnnotation:", nrow(genomeAnnotation(object)), " rows, ", ncol(genomeAnnotation(object)), " columns\n")
	print(genomeAnnotation(object)[1:5, ])
	cat("\n")
})
setClass("AlleleSet", contains="CrlmmContainer")
setClass("CallSet", contains="AlleleSet")
setClass("CNSet", contains="CallSet",
	 representation(lM="list_or_ffdf"))

setClass("IlluminaRGSet", contains="CrlmmContainer")
setClass("IlluminaXYSet", contains="CrlmmContainer")

setClass("AffymetrixAlleleSet", contains="AlleleSet")  ##AffymetrixAlleleSet
setClass("IlluminaAlleleSet", contains="AlleleSet")
##setClass("AffymetrixBigData", contains="AffymetrixAlleleSet")
##setClass("AffymetrixSmallData", contains="AffymetrixAlleleSet")
##setClass("IlluminaSmallData", contains="IlluminaAlleleSet")
##setClass("IlluminaBigData", contains="IlluminaAlleleSet")
##setMethod("initialize", "AffymetrixBigData", function(.Object, annotation){
##	.Object <- callNextMethod(.Object)
##	if(!missing(annotation)) annotation(.Object) <- annotation
##	.Object
##})
##setClass("AffymetrixCallSet", contains="CallSet")
##setClass("IlluminaCallSet", contains="CallSet")
setMethod("initialize", "AlleleSet", function(.Object, alleleA=new("matrix"), alleleB=new("matrix"), ...){
	.Object <- callNextMethod(.Object, alleleA=alleleA, alleleB=alleleB, ...)
	storageMode(.Object) <- "environment"
	.Object
})
setMethod("initialize", "CallSet", function(.Object, call=new("matrix"), callProbability=new("matrix"), ...){
	.Object <- callNextMethod(.Object, call=call, callProbability=callProbability, ...)
	storageMode(.Object) <- "environment"
	.Object
})
setMethod("initialize", "CNSet", function(.Object, CA=new("matrix"), CB=new("matrix"), lM=new("list"), ...){
	.Object <- callNextMethod(.Object, CA=CA, CB=CB, lM=lM,...)
	storageMode(.Object) <- "environment"
	.Object
})
setValidity("AlleleSet", function(object) {
	assayDataValidMembers(assayData(object), c("alleleA", "alleleB"))
})
setValidity("IlluminaRGSet", function(object) {
	assayDataValidMembers(assayData(object), c("R", "G", "zero"))
})
setValidity("IlluminaXYSet", function(object) {
	assayDataValidMembers(assayData(object), c("X", "Y", "zero"))
})

setValidity("CallSet", function(object) {
	assayDataValidMembers(assayData(object), c("alleleA", "alleleB", "call", "callProbability"))
})
setValidity("CNSet", function(object) {
	assayDataValidMembers(assayData(object), c("alleleA", "alleleB", "call", "callProbability", "CA", "CB"))
})
	 



