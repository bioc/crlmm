setOldClass("ellipse")
setMethod("plot", c("numeric", "CNSetLM"), function(x, y, copynumber=2, batch, plot.it=TRUE, ...){
	plotCNSetLM(x, y, copynumber, batch, plot.it, ...)
	})
plotCNSetLM <- function(x, y, copynumber, batch, plot.it, ...){
	require(ellipse)
	object <- y
	I <- x
	if(missing(batch)) batch <- seq(along=unique(protocolData(object)$batch))
	logA <- log2(as.matrix(A(object)[I, ]))
	logB <- log2(as.matrix(B(object)[I, ]))
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
	ffIsLoaded <- class(calls(object))[[1]] == "ff"
	if(ffIsLoaded){
		nuA <- (physical(lM(object))[["nuA"]])[I, , drop=FALSE]
		nuB <- (physical(lM(object))[["nuB"]])[I, , drop=FALSE]
		phiA <- (physical(lM(object))[["phiA"]])[I, ,drop=FALSE]
		phiB <- (physical(lM(object))[["phiB"]])[I, ,drop=FALSE]
		tau2A <- (physical(lM(object))[["tau2A"]])[I, ,drop=FALSE]
		tau2B <- (physical(lM(object))[["tau2B"]])[I, ,drop=FALSE]
		sigma2A <- (physical(lM(object))[["sig2A"]])[I, ,drop=FALSE]
		sigma2B <- (physical(lM(object))[["sig2B"]])[I, ,drop=FALSE]
		corrAB <- (physical(lM(object))[["corrAB"]])[I, ,drop=FALSE]
		corrAA <- (physical(lM(object))[["corrAA"]])[I, ,drop=FALSE]
		corrBB <- (physical(lM(object))[["corrBB"]])[I, ,drop=FALSE]
	} else {
		nuA <- lM(object)[["nuA"]][I, , drop=FALSE]
		nuB <- lM(object)[["nuB"]][I, , drop=FALSE]
		phiA <- lM(object)[["phiA"]][I, ,drop=FALSE]
		phiB <- lM(object)[["phiB"]][I, ,drop=FALSE]
		tau2A <- lM(object)[["tau2A"]][I, ,drop=FALSE]
		tau2B <- lM(object)[["tau2B"]][I, ,drop=FALSE]
		sigma2A <- lM(object)[["sig2A"]][I, ,drop=FALSE]
		sigma2B <- lM(object)[["sig2B"]][I, ,drop=FALSE]
		corrAB <- lM(object)[["corrAB"]][I, ,drop=FALSE]
		corrAA <- lM(object)[["corrAA"]][I, ,drop=FALSE]
		corrBB <- lM(object)[["corrBB"]][I, ,drop=FALSE]
	}
	for(i in seq(along=I)){
		for(b in batch){
			jj <- which(batch(object) == batch)
			if(plot.it){
				if(hasxlab(...)){
					plot(logA[i, jj], logB[i, jj], ...)
				} else{
					plot(logA[i, jj], logB[i, jj], xlab="log(A)", ylab="log(B)", ...)
				}
			}
			if(all(is.na(nuA[, b]))) {
				message("Parameter estimates for batch ", b, " not available")
				next()
			}
			for(CN in copynumber){
				for(CA in 0:CN){
					CB <- CN-CA
					A.scale <- sqrt(tau2A[i, b]*(CA==0) + sigma2A[i, b]*(CA > 0))
					B.scale <- sqrt(tau2B[i, b]*(CB==0) + sigma2B[i, b]*(CB > 0))
					scale <- c(A.scale, B.scale)
					if(CA == 0 & CB > 0) rho <- corrBB[i, b]
					if(CA > 0 & CB == 0) rho <- corrAA[i, b]
					if(CA > 0 & CB > 0) rho <- corrAB[i, b]
					if(CA == 0 & CB == 0) rho <- 0
					lines(ellipse(x=rho, centre=c(log2(nuA[i, b]+CA*phiA[i, b]),
							     log2(nuB[i, b]+CB*phiB[i,b])),
						      scale=scale), ...)
				}
			}	
		}
	}
}
