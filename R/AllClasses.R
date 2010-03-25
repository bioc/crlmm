setOldClass("ffdf")
setOldClass("ff_matrix")
##setClassUnion("matrix_or_ff", c("matrix", "ff_matrix"))
setClassUnion("list_or_ffdf", c("list", "ffdf"))
setClassUnion("ff_or_matrix", c("ff_matrix", "matrix", "ffdf"))
setClass("CNSetLM", contains="CNSet", representation(lM="list_or_ffdf"))
setMethod("initialize", "CNSetLM", function(.Object, lM=new("list"), ...){
	.Object@lM <- lM
	.Object <- callNextMethod(.Object, ...)
})
