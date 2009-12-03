setMethod("A", "AlleleSet", function(object) allele(object, "A"))
setMethod("B", "AlleleSet", function(object) allele(object, "B"))
##setReplaceMethod("A", signature(object="AlleleSet", value="matrix"),
##		 function(object, value){
##			 assayDataElementReplace(object, "senseThetaA", value)			
##		 })
##setReplaceMethod("B", signature(object="AlleleSet", value="matrix"),
##		 function(object, value){
##			 assayDataElementReplace(object, "senseThetaB", value)			
##		 })
setMethod("isSnp", "AlleleSet", function(object) {
	isSnp.AlleleSet(object)
})


