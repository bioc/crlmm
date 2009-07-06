## Methods for crlmm


setMethod("calls", "SnpSet", function(object) assayData(object)$call)
setMethod("confs", "SnpSet", function(object) assayData(object)$callProbability)
