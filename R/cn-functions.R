copynumber <- function(filenames,
		       batch,
		       summaries=c("lrr", "baf", "ca", "cb", "gt", "gtconf"),
		       write=FALSE,
		       outdir=".",
		       onefile=FALSE,
		       rda=TRUE){
	fnamelist <- split(filenames, batch)
	batchlist <- split(batch, batchlist)
	## for each element in list.  Avoid
	foreach(i=fnamelist) %do% simpleusage(filenames=fnameslist[[i]], batch=batchlist[[i]])
}

simpleusage <- function(filenames, batch, ...){
	object <- genotype(fnamelist, batch, ...)
	## run genotype on each element in fnamelist
	object <- genotypeSummary(object)
	## different than currently implemented
	##    - shrink to saved batch medians, etc.
	##    - shinkage across markers (?)
	object <- shrinkSummary(object)
	results <- array(NA, nrow(object), ncol(object), length(summaries))
	if(c("ca", "cb") %in% summaries){
		## estimateCnParameters
		##results <- estimatecnParameters(...)
	} else results <- NULL
	is.baf <- "baf" %in% summaries || "lrr" %in% summaries
	if(is.baf){
		results2 <- calculateRtheta(object) ## returns list
		## keep as list, or coerce to array
	}
	if(write){
		##write2(, onefile=onefile)
	}
	return(results)
}


imputeTheta <- function(ia, ib, theta, S=100){
##	y <- rbind(y1, y2)
	y <- theta
	N <- nrow(y); p <- ncol(y)

	mu0 <- apply(y, 2, mean, na.rm=TRUE)
	sd0 <- mu0/5
	L0 <- matrix(.1, p, p); diag(L0) <- 1; L0 <- L0*outer(sd0, sd0)
	nu0 <- p+2; S0 <- L0

	## starting values
	Sigma <- S0
	Y.full <- y
	O <- 1*(!is.na(y))
	if(all(O[1, ] == 0)){
		return(rep(NA, length(ia)))
	}
	for(j in seq_len(p)) Y.full[is.na(Y.full[, j]), j] <- mu0[j]
	## Gibbs sampler
	##THETA <- SIGMA <- Y.MISS <- NULL
	## problems:  approx. 90 observations for the means in the other batches
	##   -- only 3 observations for the mean in the current batch.
##	THETA <- matrix(NA, iter, p)
##	SIGMA <- matrix(NA, iter, p^2)
	Y.MISS <- matrix(NA, S, sum(O[1,]==0))
	bafs <- matrix(NA, S, length(a))
	for(s in seq_len(S)){
		## update lambda
		ybar <- apply(Y.full, 2, mean)
		Ln <- solve(solve(L0) + N*solve(Sigma))
		mun <- Ln%*%(solve(L0)%*%mu0 + N*solve(Sigma)%*%ybar)
		lambda <- rmvnorm(1, mun, Ln)

		##update sigma
		Sn <- S0+(t(Y.full)-c(lambda))%*%t(t(Y.full)-c(lambda))
		Sigma <- solve(rwish(nu0+N, solve(Sn)))

		## update missing data (only care about the first row.)
		a <- O[1, ] == 1
		b <- O[1, ] == 0
		iSa <- solve(Sigma[a,a])
		beta.j <- Sigma[b,a]%*%iSa
		Sigma.j <- Sigma[b,b]-Sigma[b,a]%*%iSa%*%Sigma[a,b]
		lambda.j <- lambda[b]+beta.j%*%(t(Y.full[1,a, drop=FALSE]-lambda[a]))
		Y.full[1,b] <- rmvnorm(1, lambda.j, Sigma.j)
		##SIGMA[s, ] <- c(Sigma)
		Y.MISS[s, ] <- Y.full[1, O[1, ]==0]
	}
	na.cols <- which(is.na(y[1, ]))
	THETA <- matrix(NA, S, 3)
	THETA[, na.cols] <- Y.MISS
	if(length(na.cols) < length(ia))
		THETA[, -na.cols] <- y[1, -na.cols]

	obs.theta <- atan2(ib, ia)*2/pi
	theta.aa <- matrix(THETA[, 1], S, 3)
	theta.ab <- matrix(THETA[, 2], S, 3)
	theta.bb <- matrix(THETA[, 3], S, 3)
	lessAA <- obs.theta < theta.aa
	lessAB <- obs.theta < theta.ab
	lessBB <- obs.theta < theta.bb
	grAA <- !lessAA ## >= theta.aa
	grAB <- !lessAB ## >= theta.ab
	grBB <- !lessBB ## >= theta.bb
	##not.na <- !is.na(theta.aa)
	I1 <- grAA & lessAB
	I2 <- grAB & lessBB
	##mu <- apply(Y.MISS, 2, mean)
	##bf <- matrix(NA, S, length(ia))
	bf <- rep(NA, S*length(ia))
	I1 <- as.logical(I1); I2 <- as.logical(I2)
	bf[I1] <- 0.5 * as.numeric(((obs.theta-theta.aa)/(theta.ab-theta.aa)))[I1]
	bf[I2] <- as.numeric((.5 * (obs.theta - theta.ab) / (theta.bb - theta.ab)))[I2] + 0.5
	bf[as.logical(lessAA)] <- 0
	bf[as.logical(grBB)] <- 1
	bf <- matrix(bf, S, 3, byrow=TRUE)
	pm.bf <- apply(bf, 2, mean)
	return(pm.bf)
}
