setMethod("[", "CrlmmSetList", function(x, i, j, ..., drop = FALSE){
            if (missing(drop)) drop <- FALSE
            if (missing(i) && missing(j))
              {
                 if (length(list(...))!=0)
                   stop("specify genes or samples to subset; use '",
                        substitute(x), "$", names(list(...))[[1]],
                        "' to access phenoData variables")
                 return(x)
               }
            if (!missing(j)){
              f1 <- function(x, j){
                x <- x[, j]
              }
              x <- lapply(x, f1, j)
            }
            if(!missing(i)){
              f2 <- function(x, i){
                x <- x[i, ]
              }
              x <- lapply(x, f2, i)
            }
	as(x, "CrlmmSetList")	
})

##setReplaceMethod("[[", "CrlmmSetList", function(x, i, j, ..., value) {
##	browser()
##	##otherwise infinite recursion
##	x <- as(x, "list")
##	x[[i]] <- value
##	x <- as(x, "CrlmmSetList")
##	x <- .harmonizeDimnames(x)
##	stopifnot(identical(featureNames(x[[i]]), featureNames(x[[1]])))
##	return(x)
##})

setMethod("A", "CrlmmSetList", function(object) A(object[[1]]))

setMethod("annotation", "CrlmmSetList", function(object) annotation(object[[1]]))

setMethod("B", "CrlmmSetList", function(object) B(object[[1]]))

setMethod("CA", "CrlmmSetList", function(object, ...) CA(object[[3]], ...))
setMethod("CB", "CrlmmSetList", function(object, ...) CB(object[[3]], ...))

setMethod("calls", "CrlmmSetList", function(object) calls(object[[2]]))
setMethod("cnIndex", "CrlmmSetList", function(object, ...) {
	match(cnNames(object[[1]], annotation(object)), featureNames(object))
})

setMethod("copyNumber", "CrlmmSetList", function(object) copyNumber(object[[3]]))

setMethod("combine", signature=signature(x="CrlmmSetList", y="CrlmmSetList"),
          function(x, y, ...){
		  x.abset <- x[[1]]
		  y.abset <- y[[1]]

		  x.snpset <- x[[2]]
		  y.snpset <- y[[2]]

		  abset <- combine(x.abset, y.abset)

		  ##we have hijacked the featureData slot to store parameters.  Biobase will not allow combining our 'feature' data.
		  warning("the featureData is not easily combined...  removing the featureData")
		  ##fd1 <- featureData(x.snpset)
		  ##fd2 <- featureData(y.snpset)
		  featureData(x.snpset) <- annotatedDataFrameFrom(calls(x.snpset), byrow=TRUE)
		  featureData(y.snpset) <- annotatedDataFrameFrom(calls(y.snpset), byrow=TRUE)
		  snpset <- combine(x.snpset, y.snpset)
		  merged <- list(abset, snpset)
		  merged <- as(merged, "CrlmmSetList")
		  merged
	  })
		  


##setMethod("fData", "CrlmmSetList", function(object) featureNames(object[[1]]))

setMethod("featureNames", "CrlmmSetList", function(object) featureNames(object[[1]]))

setMethod("ncol", signature(x="CrlmmSetList"), function(x) ncol(x[[1]]))
setMethod("nrow", signature(x="CrlmmSetList"), function(x) nrow(x[[1]]))
setMethod("plot", signature(x="CrlmmSetList"),
	  function(x, y, ...){
		  A <- log2(A(x))
		  B <- log2(B(x))
		  plot(A, B, ...)
	  })

setMethod("points", signature(x="CrlmmSetList"),
	  function(x, y, ...){
		  A <- log2(A(x))
		  B <- log2(B(x))
		  points(A, B, ...)
	  })

setMethod("sampleNames", "CrlmmSetList", function(object) sampleNames(object[[1]]))

setMethod("show", "CrlmmSetList", function(object){
	##for(i in seq(along=object)) show(object[[i]])
	cat("\n Elements in CrlmmSetList object: \n")
	cat("\n")
	for(i in 1:length(object)){
		cat("class: ", class(object[[i]]), "\n")
		cat("assayData elements: ", ls(assayData(object[[i]])), "\n")
		cat("Dimensions: ", dim(object[[i]]))		
		cat("\n \n")
	}
##	cat("\n")
##	cat("Dimensions:\n")
##	print(dims(object))
})

setMethod("chromosome", "CrlmmSetList", function(object) chromosome(object[[3]]))
setMethod("position", "CrlmmSetList", function(object) position(object[[3]]))
setMethod("confs", "CrlmmSetList", function(object) confs(object[[2]]))


setMethod(".harmonizeDimnames", "CrlmmSetList", function(object){
	i <- length(object)
	while(i > 1){
		object[[i-1]] <- harmonizeDimnamesTo(object[[i-1]], object[[i]])
		i <- i-1
	}
	object
})
setMethod("dims", "CrlmmSetList", function(object) sapply(object, dim))
setMethod("batch", "CrlmmSetList", function(object) batch(object[[3]]))
setMethod("$", "CrlmmSetList", function(x, name) {
	##if(!(name %in% .parameterNames()[output(x) != 0])){
	if(length(x) != 3){
		stop("'$' operature reserved for accessing parameter names in CopyNumberSet object.  CrlmmSetList must be of length 3")
	}
	j <- grep(name, fvarLabels(x[[3]]))	
	if(length(j) < 1)
		stop(name, " not in fvarLabels of CopyNumberSet object")
	if(length(j) > 1){
		warning("Multiple instances of ", name, " in fvarLabels.  Using the first instance")
		j <- j[1]
	}
	param <- fData(x[[3]])[, j]
	param
})

	
setMethod("snpIndex", "CrlmmSetList", function(object, ...){
	match(snpNames(object[[1]], annotation(object)), featureNames(object))
})
setMethod("splitByChromosome", "CrlmmSetList", function(object, cdfName, outdir){
	path <- system.file("extdata", package=paste(cdfName, "Crlmm", sep=""))	
	load(file.path(path, "snpProbes.rda"))
	load(file.path(path, "cnProbes.rda"))
	k <- grep("chr", colnames(snpProbes))
	if(length(k) < 1) stop("chr or chromosome not in colnames(snpProbes)")
	for(CHR in 1:24){
		cat("Chromosome ", CHR, "\n")
		snps <- rownames(snpProbes)[snpProbes[, k] == CHR]
		cnps <- rownames(cnProbes)[cnProbes[, k] == CHR]
		index <- c(match(snps, featureNames(object)),
			   match(cnps, featureNames(object)))
		index <- index[!is.na(index)]
		crlmmSetList <- object[index, ]
		save(crlmmSetList, file=file.path(outdir, paste("crlmmSetList_", CHR, ".rda", sep="")))
	}
})


