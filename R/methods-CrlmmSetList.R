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

setMethod("A", "CrlmmSetList", function(object) A(object[[1]]))
setMethod("B", "CrlmmSetList", function(object) B(object[[1]]))
setMethod("calls", "CrlmmSetList", function(object) calls(object[[2]]))
setMethod("cnIndex", "CrlmmSetList", function(object) match(cnNames(object[[1]]), featureNames(object)))

setMethod("combine", signature=signature(x="CrlmmSetList", y="CrlmmSetList"),
          function(x, y, ...){
		  x.abset <- x[[1]]
		  y.abset <- y[[1]]

		  x.snpset <- x[[2]]
		  y.snpset <- y[[2]]

		  abset <- combine(x.abset, y.abset)

		  ##we have hijacked the featureData slot to store parameters.  Biobase will not allow combining our 'feature' data.
		  warning("removing featureData")
		  ##fd1 <- featureData(x.snpset)
		  ##fd2 <- featureData(y.snpset)
		  featureData(x.snpset) <- annotatedDataFrameFrom(calls(x.snpset), byrow=TRUE)
		  featureData(y.snpset) <- annotatedDataFrameFrom(calls(y.snpset), byrow=TRUE)
		  snpset <- combine(x.snpset, y.snpset)
		  merged <- list(abset, snpset)
		  merged <- as(merged, "CrlmmSetList")
		  merged
	  })
		  


setMethod("featureNames", "CrlmmSetList", function(object) featureNames(object[[1]]))
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
setMethod("scanDates", "CrlmmSetList", function(object) scanDates(object[[1]]))
setMethod("show", "CrlmmSetList", function(object){
	show(object[[1]])
	show(object[[2]])
})
setMethod("snpIndex", "CrlmmSetList", function(object) match(snpNames(object[[1]]), featureNames(object)))
setMethod("splitByChromosome", "CrlmmSetList", function(object, cdfName, outdir){
	path <- system.file("extdata", package=paste(cdfName, "Crlmm", sep=""))	
	load(file.path(path, "snpProbes.rda"))
	load(file.path(path, "cnProbes.rda"))				
	for(CHR in 1:24){
		cat("Chromosome ", CHR, "\n")
		snps <- rownames(snpProbes)[snpProbes[, "chrom"] == CHR]
		cnps <- rownames(cnProbes)[cnProbes[, "chrom"] == CHR]
		index <- c(match(snps, featureNames(object)),
			   match(cnps, featureNames(object)))
		crlmmResults <- object[index, ]
		save(crlmmResults, file=file.path(outdir, paste("crlmmResults_", CHR, ".rda", sep="")))
	}
})


