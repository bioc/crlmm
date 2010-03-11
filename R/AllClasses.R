setOldClass("ffdf")
setOldClass("ff_matrix")
##setClassUnion("matrix_or_ff", c("matrix", "ff_matrix"))
setClassUnion("list_or_ffdf", c("list", "ffdf"))
setClass("CNSetLM", contains="CNSet", representation(lM="list_or_ffdf"))
setMethod("initialize", "CNSetLM", function(.Object, CA=new("matrix"), CB=new("matrix"), lM=new("list"), ...){
	.Object <- callNextMethod(.Object, CA=CA, CB=CB, lM=lM, ...)
	.Object
})
