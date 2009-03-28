## Methods for crlmm

setGeneric("calls", function(x) standardGeneric("calls"))
setMethod("calls", "SnpSet", function(x) assayData(x)$call)

setGeneric("confs", function(x) standardGeneric("confs"))
setMethod("confs", "SnpSet", function(x) assayData(x)$callProbability)
