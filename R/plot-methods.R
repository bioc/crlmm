setOldClass("ellipse")
setMethod("lines", c("CNSetLM"), function(x, y, batch, copynumber, ...){
	linesCNSetLM(x, y, batch, copynumber, ...)
})
linesCNSetLM <- function(x, y, batch, copynumber, x.axis="A", ...){
	require(ellipse)
	object <- x
	I <- y
	stopifnot(length(I) == 1) 
	calls.class <- class(calls(object)[[1]])
	ffIsLoaded <- calls.class[1] == "ff_matrix" | calls.class[1] == "ffdf" | calls.class[1]=="ff"
	column <- grep(batch, unique(batch(object)))
	stopifnot(length(column) == 1)
	if(ffIsLoaded){
		nuA <- (physical(lM(object))[["nuA"]])[I, column , drop=TRUE]
		nuB <- (physical(lM(object))[["nuB"]])[I, column , drop=TRUE]
		phiA <- (physical(lM(object))[["phiA"]])[I, column ,drop=TRUE]
		phiB <- (physical(lM(object))[["phiB"]])[I, column ,drop=TRUE]
		tau2A <- (physical(lM(object))[["tau2A"]])[I, column ,drop=TRUE]
		tau2B <- (physical(lM(object))[["tau2B"]])[I, column ,drop=TRUE]
		sigma2A <- (physical(lM(object))[["sig2A"]])[I, column ,drop=TRUE]
		sigma2B <- (physical(lM(object))[["sig2B"]])[I, column ,drop=TRUE]
		corrAB <- (physical(lM(object))[["corrAB"]])[I, column ,drop=TRUE]
		corrAA <- (physical(lM(object))[["corrAA"]])[I, column ,drop=TRUE]
		corrBB <- (physical(lM(object))[["corrBB"]])[I, column ,drop=TRUE]
	} else {
		nuA <- lM(object)[["nuA"]][I, column , drop=TRUE]
		nuB <- lM(object)[["nuB"]][I, column , drop=TRUE]
		phiA <- lM(object)[["phiA"]][I, column ,drop=TRUE]
		phiB <- lM(object)[["phiB"]][I, column ,drop=TRUE]
		tau2A <- lM(object)[["tau2A"]][I, column ,drop=TRUE]
		tau2B <- lM(object)[["tau2B"]][I, column ,drop=TRUE]
		sigma2A <- lM(object)[["sig2A"]][I, column ,drop=TRUE]
		sigma2B <- lM(object)[["sig2B"]][I, column ,drop=TRUE]
		corrAB <- lM(object)[["corrAB"]][I, column ,drop=TRUE]
		corrAA <- lM(object)[["corrAA"]][I, column ,drop=TRUE]
		corrBB <- lM(object)[["corrBB"]][I, column ,drop=TRUE]
	}
	if(all(is.na(nuA))) {
		message("Parameter estimates for batch ", b, " not available")
		next()
	}
	for(CN in copynumber){
		for(CA in 0:CN){
			CB <- CN-CA
			A.scale <- sqrt(tau2A*(CA==0) + sigma2A*(CA > 0))
			B.scale <- sqrt(tau2B*(CB==0) + sigma2B*(CB > 0))
			scale <- c(A.scale, B.scale)
			if(CA == 0 & CB > 0) rho <- corrBB
			if(CA > 0 & CB == 0) rho <- corrAA
			if(CA > 0 & CB > 0) rho <- corrAB
			if(CA == 0 & CB == 0) rho <- 0
			if(x.axis=="A"){
				lines(ellipse(x=rho, centre=c(log2(nuA+CA*phiA),
						     log2(nuB+CB*phiB)),
					      scale=scale), ...)
			} else {
				lines(ellipse(x=rho, centre=c(log2(nuB+CB*phiB),
						     log2(nuA+CA*phiA)),
					      scale=rev(scale)), ...)
			}
		}
	}	
}
