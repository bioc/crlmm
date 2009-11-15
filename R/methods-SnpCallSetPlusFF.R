setMethod("featureNames", "SnpCallSetPlusFF", function(object) sampleNames(featureData(object)))
setMethod("A", "SnpCallSetPlusFF", function(object){})

