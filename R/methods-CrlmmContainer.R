setMethod("crlmmOptions", "CrlmmContainer", function(object) object@options)
setReplaceMethod("crlmmOptions", "CrlmmContainer", function(object, value){
	object@options <- value
	return(object)
})
setMethod("[", "CrlmmContainer", function(x, i, j, ..., drop=FALSE){
	if(!missing(i))
		if(is.logical(i)) i <- which(i)
	if(!missing(j))
		if(is.logical(j)) j <- which(j)
	x <- callNextMethod(x, i, j, ..., drop=FALSE)
	if(!missing(i)){
		genomeAnnotation(x) <- genomeAnnotation(x)[i, ]
	}
	x
})
setMethod("getGenomeAnnotation", signature(object="CrlmmContainer"),
	  function(object, featurenames){
		  getGenomeAnnotation(crlmmOptions(object), featureNames(object))
	  })
setMethod("genomeAnnotation", signature(object="CrlmmContainer"), function(object) object@genomeAnnotation)
setReplaceMethod("genomeAnnotation", signature(object="CrlmmContainer", "matrix_or_ff"),
		 function(object, value){
			 object@genomeAnnotation <- value
			 object
		 })
setMethod("chromosome", "CrlmmContainer", function(object)
	  genomeAnnotation(object)[, "chromosome"]
	  )
setMethod("position", "CrlmmContainer", function(object)
	  genomeAnnotation(object)[, "position"]
	  )
setMethod("isSnp", "CrlmmContainer", function(object)
	  genomeAnnotation(object)[, "isSnp"]
	  )

setMethod("close", "CrlmmContainer", function(con, ...){
	object <- con
	names <- ls(assayData(object))
	L <- length(names)
	for(i in 1:L) close(eval(substitute(assayData(object)[[NAME]], list(NAME=names[i]))))
	close(genomeAnnotation(object))
	return()
})
setMethod("open", "CrlmmContainer", function(con, ...){
	object <- con
	names <- ls(assayData(object))
	L <- length(names)
	for(i in 1:L) open(eval(substitute(assayData(object)[[NAME]], list(NAME=names[i]))))
	open(genomeAnnotation(object))
	return()
})


