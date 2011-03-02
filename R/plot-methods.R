setMethod("lines", signature=signature(x="CNSet"),
	  function(x, y, batch, copynumber, grid=FALSE, ...){
		  linesCNSet(x, y, batch, copynumber, grid=grid, ...)
	  })

linesCNSet <- function(x, y, batch, copynumber, x.axis="A", grid=FALSE, ...){
	require(ellipse)
	object <- x
	marker.index <- y
	stopifnot(length(marker.index) == 1)
	batch.index <- match(batch, batchNames(object))
	stopifnot(length(batch.index) == 1)
	nuA <- nu(object, "A")[marker.index, batch.index]
	nuB <- nu(object, "B")[marker.index, batch.index]
	phiA <- phi(object, "A")[marker.index, batch.index]
	phiB <- phi(object, "B")[marker.index, batch.index]

	taus <- tau2(object, i=marker.index, j=batch.index)
	tau2A <- taus[, "A", "BB", ]
	tau2B <- taus[, "B", "AA", ]
	sigma2A <- taus[, "A", "AA", ]
	sigma2B <- taus[, "B", "BB", ]
	cors <- corr(object, i=marker.index, j=batch.index)[, , ]
	corrAB <- cors[["AB"]]
	corrAA <- cors[["AA"]]
	corrBB <- cors[["BB"]]
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
				dat.ellipse <- ellipse(x=rho, centre=c(log2(nuA+CA*phiA),
							      log2(nuB+CB*phiB)),
						       scale=scale)
				if(!grid){
					lines(dat.ellipse, ...)
				} else {
					llines(dat.ellipse[, 1], dat.ellipse[, 2], ...)
				}
			} else {
				dat.ellipse <- ellipse(x=rho, centre=c(log2(nuB+CB*phiB),
							      log2(nuA+CA*phiA)), scale=rev(scale))
				if(!grid){
					lines(dat.ellipse, ...)
				} else {
					llines(dat.ellipse[, 1], dat.ellipse[, 2], ...)
				}
			}
		}
	}
}
