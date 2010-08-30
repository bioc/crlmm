setClassUnion("integerOrMissing", c("integer", "missing", "numeric"))
setGeneric("totalCopyNumber", function(object, i, j, ...) standardGeneric("totalCopyNumber"))
setMethod("totalCopyNumber",
	  signature=signature(object="CNSet", i="integerOrMissing", j="integerOrMissing"),
	  function(object, i, j, ...){
	if(missing(i) & missing(j)){
		if(inherits(CA(object), "ff") | inherits(CA(object), "ffdf")) stop("Must specify i and/or j for ff objects")
	}
	if(missing(i) & !missing(j)){
		snp.index <- which(isSnp(object))	
		cn.total <- as.matrix(CA(object)[, j])
		if(length(snp.index) > 0){
			cb <- as.matrix(CB(object)[snp.index, j])
			snps <- (1:nrow(cn.total))[i %in% snp.index]
			cn.total[snps, ] <- cn.total[snps, j] + cb				
		}
	}
	if(!missing(i) & missing(j)){
		snp.index <- intersect(which(isSnp(object)), i)
		cn.total <- as.matrix(CA(object)[i, ])
		if(length(snp.index) > 0){
			cb <- as.matrix(CB(object)[snp.index, ])
			snps <- (1:nrow(cn.total))[i %in% snp.index]
			cn.total[snps, ] <- cn.total[snps, ] + cb				
		}
	}
	if(!missing(i) & !missing(j)){
		snp.index <- intersect(which(isSnp(object)), i)		
		cn.total <- as.matrix(CA(object)[i, j])
		if(length(snp.index) > 0){
			cb <- as.matrix(CB(object)[snp.index, j])
			snps <- (1:nrow(cn.total))[i %in% snp.index]
			cn.total[snps, ] <- cn.total[snps, ] + cb
		}
	}
	cn.total <- cn.total/100
	dimnames(cn.total) <- NULL
	return(cn.total)
})



