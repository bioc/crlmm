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


cnPanel <- function(x, y, ..., pch.cols, gt, cbs.segs, hmm.segs=NULL, shades, subscripts, add.ideogram=TRUE){
	##if(panel.number() == 2) browser()
	add.ideogram <- add.ideogram[[panel.number()]]
	##cbs.segs <- cbs.segs[[panel.number()]]
	draw.hmm.states <- ifelse(panel.number() <= length(hmm.segs), TRUE, FALSE)
	panel.grid(h=6, v=10)
	which.hom <- which(gt[subscripts] == 1 | gt[subscripts]==3)
	which.het <- which(gt[subscripts] == 2)
	panel.xyplot(x, y, col="grey60", ...)
	lpoints(x[which.hom], y[which.hom], col=pch.cols[1], ...)
	lpoints(x[which.het], y[which.het], col=pch.cols[2], ...)
	lsegments(x0=start(cbs.segs)/1e6, x1=end(cbs.segs)/1e6,
		  y0=cbs.segs$seg.mean,
		  y1=cbs.segs$seg.mean, ...)
	if(draw.hmm.states){
		hmm.segs <- hmm.segs[order(width(hmm.segs), decreasing=TRUE), ]
		lrect(xleft=start(hmm.segs)/1e6,
		      xright=end(hmm.segs)/1e6,
		      ybottom=-0.4,
		      ytop=0,
		      border=shades[hmm.segs$state],
		      col=shades[hmm.segs$state])
	}
	ltext(130, 5, paste("MAD:", round(mad(y, na.rm=TRUE), 2)))
	if(add.ideogram){
		pathto <- system.file("hg18", package="SNPchip")
		cytoband <- read.table(file.path(pathto, "cytoBand.txt"), as.is=TRUE)
		cytoband$V2 <- cytoband$V2/1e6
		cytoband$V3 <- cytoband$V3/1e6
		colnames(cytoband) <- c("chrom", "start", "end", "name", "gieStain")
		plotCytoband(unique(hmm.segs$chrom),
			     cytoband=cytoband,
			     cytoband.ycoords=c(5.6, 5.9),
			     new=FALSE,
			     label.cytoband=FALSE,
			     build="hg18",
			     use.lattice=TRUE)
	}
}
