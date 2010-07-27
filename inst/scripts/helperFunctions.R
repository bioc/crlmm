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


constructFeatureData <- function(gns, cdfName){
	pkgname <- paste(cdfName, "Crlmm", sep="")	
	path <- system.file("extdata", package=pkgname)
	load(file.path(path, "cnProbes.rda"))
	load(file.path(path, "snpProbes.rda"))
	cnProbes$chr <- chromosome2integer(cnProbes$chr)
	cnProbes <- as.matrix(cnProbes)
	snpProbes$chr <- chromosome2integer(snpProbes$chr)
	snpProbes <- as.matrix(snpProbes)
	mapping <- rbind(snpProbes, cnProbes, deparse.level=0)
	mapping <- mapping[match(gns, rownames(mapping)), ]
	isSnp <- 1L-as.integer(gns %in% rownames(cnProbes))
	mapping <- cbind(mapping, isSnp, deparse.level=0)
	stopifnot(identical(rownames(mapping), gns))
	colnames(mapping) <- c("chromosome", "position", "isSnp")
	new("AnnotatedDataFrame",
	    data=data.frame(mapping),
	    varMetadata=data.frame(labelDescription=colnames(mapping)))	
}
constructAssayData <- function(np, snp, object, storage.mode="environment", order.index){
	stopifnot(identical(snp$gns, featureNames(object)))
	gns <- c(featureNames(object), np$gns)
	sns <- np$sns
	np <- np[1:2]
	snp <- snp[1:2]
	stripnames <- function(x) {
		dimnames(x) <- NULL
		x
	}
	np <- lapply(np, stripnames)
	snp <- lapply(snp, stripnames)
	A <- rbind(snp[[1]], np[[1]], deparse.level=0)[order.index, ]
	B <- rbind(snp[[2]], np[[2]], deparse.level=0)[order.index, ]
	gt <- stripnames(calls(object))
	emptyMatrix <- matrix(integer(), nrow(np[[1]]), ncol(A))
	gt <- rbind(gt, emptyMatrix, deparse.level=0)[order.index,]
	pr <- stripnames(snpCallProbability(object))
	pr <- rbind(pr, emptyMatrix, deparse.level=0)[order.index, ]
	emptyMatrix <- matrix(integer(), nrow(A), ncol(A))	
	aD <- assayDataNew(storage.mode,
			   alleleA=A,
			   alleleB=B,
			   call=gt,
			   callProbability=pr,
			   CA=emptyMatrix,
			   CB=emptyMatrix)
}
