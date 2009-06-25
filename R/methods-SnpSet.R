## Methods for crlmm


setMethod("calls", "SnpSet", function(object) assayData(object)$call)
setMethod("confs", "SnpSet", function(x) assayData(x)$callProbability)
