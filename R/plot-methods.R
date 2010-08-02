setOldClass("ellipse")
setMethod("plot", c("numeric", "CNSetLM"), function(x, y, copynumber=2, batch, plot.it=TRUE, ...){
	plotCNSetLM(x, y, copynumber, batch, plot.it, ...)
	})
plotCNSetLM <- function(x, y, copynumber, batch, plot.it, ...){
	require(ellipse)
	object <- y
	stopifnot(length(batch) == 1)
	stopifnot(batch %in% batch(object))
	sample.index <- which(batch(object) %in% batch)
	I <- x
	if(missing(batch)) batch <- seq(along=unique(protocolData(object)$batch))
	logA <- log2(as.matrix(A(object)[I, sample.index]))
	dimnames(logA) <- NULL
	logB <- log2(as.matrix(B(object)[I, sample.index]))
	dimnames(logB) <- NULL
	if(!"ylim" %in% names(list(...))){
		ylim <- range(c(as.numeric(logA), as.numeric(logB)), na.rm=TRUE)
	}
	hasxlab <- function(...) !all(is.na(pmatch(names(list(...)), "lab")))
	hascol <- function(...) !all(is.na(pmatch(names(list(...)), "col")))
	hasbg <- function(...) !all(is.na(pmatch(names(list(...)), "bg")))	
##	if(color.by.genotype){
##		gt <- as.matrix(calls(object)[I, ])
##		if(!hascol(...)) col <- gt
##	}
	calls.class <- class(calls(object)[[1]])
	ffIsLoaded <- calls.class[1] == "ff_matrix" | calls.class[1] == "ffdf" | calls.class[1]=="ff"
	column <- grep(batch, unique(batch(object)))
	stopifnot(length(column) == 1)
	if(ffIsLoaded){
		nuA <- (physical(lM(object))[["nuA"]])[I, column , drop=FALSE]
		nuB <- (physical(lM(object))[["nuB"]])[I, column , drop=FALSE]
		phiA <- (physical(lM(object))[["phiA"]])[I, column ,drop=FALSE]
		phiB <- (physical(lM(object))[["phiB"]])[I, column ,drop=FALSE]
		tau2A <- (physical(lM(object))[["tau2A"]])[I, column ,drop=FALSE]
		tau2B <- (physical(lM(object))[["tau2B"]])[I, column ,drop=FALSE]
		sigma2A <- (physical(lM(object))[["sig2A"]])[I, column ,drop=FALSE]
		sigma2B <- (physical(lM(object))[["sig2B"]])[I, column ,drop=FALSE]
		corrAB <- (physical(lM(object))[["corrAB"]])[I, column ,drop=FALSE]
		corrAA <- (physical(lM(object))[["corrAA"]])[I, column ,drop=FALSE]
		corrBB <- (physical(lM(object))[["corrBB"]])[I, column ,drop=FALSE]
	} else {
		nuA <- lM(object)[["nuA"]][I, column , drop=FALSE]
		nuB <- lM(object)[["nuB"]][I, column , drop=FALSE]
		phiA <- lM(object)[["phiA"]][I, column ,drop=FALSE]
		phiB <- lM(object)[["phiB"]][I, column ,drop=FALSE]
		tau2A <- lM(object)[["tau2A"]][I, column ,drop=FALSE]
		tau2B <- lM(object)[["tau2B"]][I, column ,drop=FALSE]
		sigma2A <- lM(object)[["sig2A"]][I, column ,drop=FALSE]
		sigma2B <- lM(object)[["sig2B"]][I, column ,drop=FALSE]
		corrAB <- lM(object)[["corrAB"]][I, column ,drop=FALSE]
		corrAA <- lM(object)[["corrAA"]][I, column ,drop=FALSE]
		corrBB <- lM(object)[["corrBB"]][I, column ,drop=FALSE]
	}
	for(i in seq(along=I)){
		if(plot.it){
			if(hasxlab(...)){
				plot(logA[i, ], logB[i, ], ...)
			} else{
				plot(logA[i, ], logB[i, ], xlab="log(A)", ylab="log(B)", ...)
			}
		}
		if(all(is.na(nuA))) {
			message("Parameter estimates for batch ", b, " not available")
			next()
		}
		for(CN in copynumber){
			for(CA in 0:CN){
				CB <- CN-CA
				A.scale <- sqrt(tau2A[i, ]*(CA==0) + sigma2A[i, ]*(CA > 0))
				B.scale <- sqrt(tau2B[i, ]*(CB==0) + sigma2B[i, ]*(CB > 0))
				scale <- c(A.scale, B.scale)
				if(CA == 0 & CB > 0) rho <- corrBB[i, ]
				if(CA > 0 & CB == 0) rho <- corrAA[i, ]
				if(CA > 0 & CB > 0) rho <- corrAB[i, ]
				if(CA == 0 & CB == 0) rho <- 0
				lines(ellipse(x=rho, centre=c(log2(nuA[i, ]+CA*phiA[i, ]),
						     log2(nuB[i, ]+CB*phiB[i,])),
					      scale=scale), ...)
			}
		}	
	}
}
