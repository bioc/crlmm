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
setMethod(".harmonizeDimnames", "CrlmmSetList", function(object){
	i <- length(object)
	while(i > 1){
		object[[i-1]] <- harmonizeDimnamesTo(object[[i-1]], object[[i]])
		i <- i-1
	}
	object
})

setMethod("A", "CrlmmSetList", function(object) A(object[[1]]))
setMethod("annotation", "CrlmmSetList", function(object) annotation(object[[1]]))
setMethod("B", "CrlmmSetList", function(object) B(object[[1]]))
setMethod("batch", "CrlmmSetList", function(object) batch(object[[3]]))
setMethod("CA", "CrlmmSetList", function(object, ...) CA(object[[3]], ...))
setMethod("CB", "CrlmmSetList", function(object, ...) CB(object[[3]], ...))
setMethod("calls", "CrlmmSetList", function(object) calls(object[[2]]))
setMethod("chromosome", "CrlmmSetList", function(object) chromosome(object[[3]]))
setMethod("cnIndex", "CrlmmSetList", function(object, ...) {
	match(cnNames(object[[1]], annotation(object)), featureNames(object))
})
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
setMethod("confs", "CrlmmSetList", function(object) confs(object[[2]]))
setMethod("copyNumber", "CrlmmSetList", function(object) copyNumber(object[[3]]))
setMethod("dims", "CrlmmSetList", function(object) sapply(object, dim))
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
setMethod("position", "CrlmmSetList", function(object) position(object[[3]]))
setMethod("sampleNames", "CrlmmSetList", function(object) sampleNames(object[[1]]))
setMethod("show", "CrlmmSetList", function(object){
	cat("\n Elements in CrlmmSetList object: \n")
	cat("\n")
	for(i in 1:length(object)){
		cat("class: ", class(object[[i]]), "\n")
		cat("assayData elements: ", ls(assayData(object[[i]])), "\n")
		cat("Dimensions: ", dim(object[[i]]))		
		cat("\n \n")
	}
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
setMethod("update", "CrlmmSetList", function(object, ...){
	computeCopynumber(object, ...)
})


setMethod("boxplot", "CrlmmSetList", function(x, ...){
##boxplot.CrlmmSetList <- function(x, ...){
	if(length(x) != 3) stop("elements of list should be of class ABset, SnpSet, and CopyNumberSet, respectively.")
	genotypes <- calls(x)-1
	A1 <- A(x)
	B1 <- B(x)
	Alist <- split(A1, genotypes)
	Alist <- rev(Alist)
	Blist <- split(B1, genotypes)
	ylim <- range(unlist(Alist))
	boxplot(Alist, xaxt="n", ylab=expression(I[A]), 
		cex.axis=0.6,
		xlab="",
		ylim=range(unlist(Alist), na.rm=TRUE),	
		border="grey50", xaxs="i", at=0:2, xlim=c(-0.5, 2.5),
		cex.main=0.9, xaxt="n",
		yaxt="n",
		col=cols)
	axis(2, at=pretty(ylim), cex.axis=0.8)
	axis(1, at=0:2, labels=rev(c("2A (AA genotype)", "1A (AB genotype)", "0A (BB genotype)")), cex.axis=0.7)
	##extracts nuA for first batch
	suppressWarnings(nuA <- x$nuA)
	segments(-1, nuA, 
		 0, nuA, lty=2, col="blue")
	suppressWarnings(phiA <- x$phiA)
	##phiA <- fData(x)[["phiA_A"]]
	segments(0, nuA,
		 2.5, nuA+2.5*phiA, lwd=2, col="blue")
	axis(2, at=nuA, labels=expression(hat(nu[A])), cex.axis=0.9)
	text(0, ylim[1], labels=paste("n =", length(Alist[[1]])), cex=0.8)
	text(1, ylim[1], labels=paste("n =", length(Alist[[2]])), cex=0.8)
	text(2, ylim[1], labels=paste("n =", length(Alist[[3]])), cex=0.8)
##	segments(0.5, nuA+0.5*phiA,
##		 1, pretty(ylim, n=5)[2], lty=2, col="blue")
##	segments(1, pretty(ylim, n=5)[2],
##		 1.2, pretty(ylim)[2], lty=2, col="blue")
##	text(1.25, pretty(ylim)[2], pos=4, 
##	     labels=expression(hat(nu[A]) + c[A]*hat(phi[A])))
	legend("topleft", fill=rev(cols), legend=c("AA", "AB", "BB"))##, title="diallelic genotypes")	

	ylim <- range(unlist(Blist))
	boxplot(Blist, xaxt="n", ylab=expression(I[B]), 
		##xlab="diallelic genotypes", 
		cex.axis=0.6, 
		ylim=range(unlist(Blist), na.rm=TRUE),
		border="grey50", xaxs="i", 
		at=0:2, xlim=c(-0.5, 2.5),
		cex.main=0.9, yaxt="n", 
		col=rev(cols))
	axis(2, at=pretty(ylim), cex.axis=0.8)
	axis(1, at=0:2, labels=c("0B (AA genotype)", "1B (AB genotype)", "2B (BB genotype)"), cex.axis=0.7)
	##nuB <- fData(x)[["nuB_A"]]
	suppressWarnings(nuB <- x$nuB)
	segments(-1, nuB, 
		 0, nuB, lty=2, col="blue")
	##phiB <- fData(x[i,])[["phiB_A"]]
	suppressWarnings(phiB <- x$phiB)
	segments(0, nuB,
		 2.5, nuB+2.5*phiB, lwd=2, col="blue")
	axis(2, at=nuB, labels=expression(hat(nu[B])), cex.axis=0.9)
	text(0, ylim[1], labels=paste("n =", length(Blist[[1]])), cex=0.8)
	text(1, ylim[1], labels=paste("n =", length(Blist[[2]])), cex=0.8)
	text(2, ylim[1], labels=paste("n =", length(Blist[[3]])), cex=0.8)
##	segments(0.5, nuB+0.5*phiB,
##		 1, pretty(ylim,n=5)[2], lty=2, col="blue")
##	segments(1, pretty(ylim, n=5)[2],
##		 1.2, pretty(ylim,n=5)[2], lty=2, col="blue")
##	text(1.25, pretty(ylim,n=5)[2], pos=4, labels=expression(hat(nu[B]) + c[B]*hat(phi[B])))	
	mtext(featureNames(x)[i], 3, outer=TRUE, line=1)	
})
