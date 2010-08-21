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
	nuA <- nu(object, "A")[I, column]
	nuB <- nu(object, "B")[I, column]
	phiA <- phi(object, "A")[I, column]
	phiB <- phi(object, "B")[I, column]		
	tau2A <- tau2(object, "A")[I, column]
	tau2B <- tau2(object, "B")[I, column]
	sigma2A <- sigma2(object, "A")[I, column]
	sigma2B <- sigma2(object, "B")[I, column]
	corrAB <- corr(object, "AB")[I, column]
	corrAA <- corr(object, "AA")[I, column]
	corrBB <- corr(object, "BB")[I, column]
	if(all(is.na(nuA))) {
		message("Parameter estimates for batch ", batch, " not available")
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
