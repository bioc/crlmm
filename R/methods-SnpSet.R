## Methods for crlmm
setMethod("calls", "SnpSet", function(object) assayData(object)$call)
setReplaceMethod("calls", signature(object="SnpSet", value="matrix"), function(object, value) assayDataElementReplace(object, "call", value))
setMethod("confs", "SnpSet", function(object) assayData(object)$callProbability)
setReplaceMethod("confs", signature(object="SnpSet", value="matrix"), function(object, value) assayDataElementReplace(object, "callProbability", value))
