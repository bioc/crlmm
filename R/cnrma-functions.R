rowCovs <- function(x, y, ...){
	notna <- !is.na(x)
	N <- rowSums(notna)
	x <- suppressWarnings(log2(x))
	if(!missing(y)) y <- suppressWarnings(log2(y))
	return(rowSums((x - rowMeans(x, ...)) * (y - rowMeans(y, ...)), ...)/(N-1))
}

rowMAD <- function(x, y, ...){
	notna <- !is.na(x)
	mad <- 1.4826*rowMedians(abs(x-rowMedians(x, ...)), ...)
	return(mad)
}

rowCors <- function(x, y, ...){
	N <- rowSums(!is.na(x))
	x <- suppressWarnings(log2(x))
	y <- suppressWarnings(log2(y))
	sd.x <- rowSds(x, ...)
	sd.y <- rowSds(y, ...)
	covar <- rowSums((x - rowMeans(x, ...)) * (y - rowMeans(y, ...)), ...)/(N-1)
	return(covar/(sd.x*sd.y))
}

nuphiAllele <- function(allele, Ystar, W, envir, p){
	Ns <- get("Ns", envir)
	NOHET <- mean(Ns[, p, 2], na.rm=TRUE) < 0.05
	chrom <- get("chrom", envir)
	if(missing(allele)) stop("must specify allele")
	if(chrom == 23){
		gender <- get("gender", envir)
		if(allele == "A" & gender == "female")
			X <- cbind(1, 2:0)
		if(allele == "B" & gender == "female")
			X <- cbind(1, 0:2)
		if(allele == "A" & gender == "male")
			X <- cbind(1, 1:0)
		if(allele == "B" & gender == "male")
			X <- cbind(1, 0:1)
	} else {##autosome
		if(allele == "A") X <- cbind(1, 2:0) else X <- cbind(1, 0:2)
		if(NOHET) X <- X[-2, ] ##more than 1 X chromosome, but all homozygous		
	}
	if(any(!is.finite(W))){## | any(!is.finite(V))){
		i <- which(rowSums(!is.finite(W)) > 0)
		browser()
		stop("Inf values in W or V")
	}
	##How to quickly generate Xstar, Xstar = diag(W) %*% X
	generateX <- function(w, X) as.numeric(diag(w) %*% X)
	##Not instant
	Xstar <- apply(W, 1, generateX, X)
	generateIXTX <- function(x, nrow=3) {
		X <- matrix(x, nrow=nrow)
		XTX <- crossprod(X)
		solve(XTX)
	}
	##a little slow
	IXTX <- apply(Xstar, 2, generateIXTX, nrow=nrow(X))
	betahat <- matrix(NA, 2, nrow(Ystar))
	ses <- matrix(NA, 2, nrow(Ystar))
	for(i in 1:nrow(Ystar)){
		betahat[, i] <- crossprod(matrix(IXTX[, i], 2, 2), crossprod(matrix(Xstar[, i], nrow=nrow(X)), Ystar[i, ]))
		ssr <- sum((Ystar[i, ] - matrix(Xstar[, i], nrow(X), 2) %*% matrix(betahat[, i], 2, 1))^2)
		ses[, i] <- sqrt(diag(matrix(IXTX[, i], 2, 2) * ssr))
	}
	nu <- betahat[1, ]
	phi <- betahat[2, ]
	nu.se <- ses[1,]
	phi.se <- ses[2,]
	##dimnames(nu) <- dimnames(phi) <- dimnames(nu.ses) <- dimnames(phi.ses)
	assign("nu", nu, envir=parent.frame(1))
	assign("phi", phi, envir=parent.frame(1))
	assign("nu.se", nu.se, envir=parent.frame(1))
	assign("phi.se", phi.se, envir=parent.frame(1))	
}



celDatesFrom <- function(celfiles, path=""){
	require(affyio)
	celdates <- vector("character", length(celfiles))
	celtimes <- vector("character", length(celfiles))
	for(i in seq(along=celfiles)){
		if(i %% 100 == 0) cat(".")
		tmp <- read.celfile.header(file.path(path, celfiles[i]), info="full")$DatHeader
		tmp <- strsplit(tmp, "\ +")
		celdates[i] <- tmp[[1]][6]
		celtimes[i] <- tmp[[1]][7]
	}
	tmp <- paste(celdates, celtimes)
	celdts <- strptime(tmp, "%m/%d/%y %H:%M:%S")
	return(celdts)
}


cnrma <- function(filenames, sns, cdfName, seed=1, verbose=FALSE){
##	require(preprocessCore) || stop("Package preprocessCore not available")
	if (missing(sns)) sns <- basename(filenames)
	if (missing(cdfName))
		cdfName <- read.celfile.header(filenames[1])$cdfName
	##  stuffDir <- changeToCrlmmAnnotationName(cdfName)
	pkgname <- crlmm:::getCrlmmAnnotationName(cdfName)
	if(!require(pkgname, character.only=TRUE, quietly=!verbose)){
		suggCall <- paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep="")
		msg <- paste("If", pkgname, "is installed on an alternative location, please load it manually by using", suggCall)
		message(strwrap(msg))
		stop("Package ", pkgname, " could not be found.")
		rm(suggCall, msg)
	}
	
	##I'm using the pd.mapping package
	data(npProbesFid, package=pkgname)
	fid <- npProbesFid
	gc()
	set.seed(seed)
	idx2 <- sample(length(fid), 10^5) ##for skewness. no need to do everything
	SKW <- vector("numeric", length(filenames))
	NP <- matrix(NA, nrow(map), length(filenames))
	verbose <- TRUE
	if(verbose){
		message("Processing ", length(filenames), " files.")
		if (getRversion() > '2.7.0') pb <- txtProgressBar(min=0, max=length(filenames), style=3)
	}
	##load reference distribution obtained from hapmap
	data(list="1m_reference_cn", pkgname)
	for(i in seq(along=filenames)){
		y <- as.matrix(read.celfile(filenames[i], intensity.means.only=TRUE)[["INTENSITY"]][["MEAN"]][fid])
		x <- log2(y[idx2])
		SKW[i] <- mean((x-mean(x))^3)/(sd(x)^3)
		rm(x)
		NP[, i] <- as.integer(normalize.quantiles.use.target(y, target=reference))
		##A[, i] <- intMedianSummaries(y[aIndex, 1, drop=FALSE], pnsa)
		##B[, i] <- intMedianSummaries(y[bIndex, 1, drop=FALSE], pnsb)
		##Now to fit the EM
		if (verbose)
			if (getRversion() > '2.7.0') setTxtProgressBar(pb, i)
			else cat(".")
	}
	dimnames(NP) <- list(names(fid), sns)
	##dimnames(NP) <- list(map[, "man_fsetid"], sns)
	res3 <- list(NP=NP, SKW=SKW)
	return(res3)
}

getFlags <- function(phi.thr, envir){
	nuA <- get("nuA", envir)
	nuB <- get("nuB", envir)
	phiA <- get("phiA", envir)
	phiB <- get("phiB", envir)
	
	negativeNus <- nuA < 1 | nuB < 1
	negativeNus <- negativeNus > 0

	negativePhis <- phiA < phi.thr | phiB < phi.thr
	negativePhis <- negativePhis > 0
	negativeCoef <- negativeNus | negativePhis

	notfinitePhi <- !is.finite(phiA) | !is.finite(phiB)
	notfinitePhi <- notfinitePhi > 0
	flags <- negativeCoef | notfinitePhi  
}

goodSnps <- function(phi.thr, envir, fewAA=20, fewBB=20){
	Ns <- get("Ns", envir)
	flags <- getFlags(phi.thr=phi.thr, envir)
	fewAA <- Ns[, , 1] < fewAA
	fewBB <- Ns[, , 3] < fewBB
	flagsA <- flags | fewAA
	flagsB <- flags | fewBB
	flags <- list(A=flagsA, B=flagsB)
	return(flags)
}

instantiateObjects <- function(calls, NP, plate, envir, chrom){
	if(missing(chrom)) stop("must specify chrom")
	assign("chrom", chrom, envir)
	snps <- rownames(calls)
	cnvs <- rownames(NP)
	sns <- basename(colnames(calls))
	stopifnot(identical(colnames(calls), colnames(NP)))
	assign("sns", sns, envir)	
	assign("snps", snps, envir)
	assign("cnvs", cnvs, envir)	

	CA <- CB <- matrix(NA, nrow(calls), ncol(calls))
	assign("plate", plate, envir)
	uplate <- unique(plate)	
	Ns <- array(NA, dim=c(nrow(calls), length(uplate), 3))
	dimnames(Ns)[[3]] <- c("AA", "AB", "BB")
	muA.AA <- matrix(NA, nrow=nrow(calls), length(uplate))
	##Sigma.AA <- array(NA, dim=c(nrow(calls), 3, length(uplate)))
	##dimnames(Sigma.AA)[2:3] <- list(c("varA", "varB", "cov"), uplate)

	NP.CT <- matrix(NA, nrow(NP), ncol(NP))
	NP.sds <- matrix(NA, nrow(NP), ncol(NP))
	assign("NP.CT", NP.CT, envir)
	assign("NP.sds", NP.sds, envir)
	nus <- matrix(NA, nrow(NP), length(uplate))
	assign("nus", nus, envir=envir)  
	assign("phis", nus, envir=envir)

	plates.completed <- rep(FALSE, length(uplate))
	assign("plates.completed", plates.completed, envir)
	
	fit.variance <- NULL
	snpflags <- NULL
	assign("fit.variance", fit.variance, envir=envir)
	assign("snpflags", snpflags, envir=envir)
	
	assign("Ns", Ns, envir=envir)		
	assign("uplate", uplate, envir=envir)	
	assign("muA.AA", muA.AA, envir=envir)
	assign("muA.AB", muA.AA, envir=envir)
	assign("muA.BB", muA.AA, envir=envir)
	assign("muB.AA", muA.AA, envir=envir)
	assign("muB.AB", muA.AA, envir=envir)
	assign("muB.BB", muA.AA, envir=envir)
	
	assign("tau2A", muA.AA, envir=envir)
	assign("tau2B", muA.AA, envir=envir)
	assign("sig2A", muA.AA, envir=envir)
	assign("sig2B", muA.AA, envir=envir)
	assign("corr", muA.AA, envir=envir)
	assign("corrA.BB", muA.AA, envir=envir)
	assign("corrB.AA", muA.AA, envir=envir)			
	
	assign("nuA", muA.AA, envir=envir)
	assign("nuB", muA.AA, envir=envir)
	assign("phiA", muA.AA, envir=envir)
	assign("phiB", muA.AA, envir=envir)
	assign("nuA.se", muA.AA, envir=envir)
	assign("nuB.se", muA.AA, envir=envir)
	assign("phiA.se", muA.AA, envir=envir)
	assign("phiB.se", muA.AA, envir=envir)
	assign("sd.CT", CA, envir=envir)
	assign("CA", CA, envir=envir)
	assign("CB", CB, envir=envir)
}

computeCopynumber <- function(A,
			      B,
			      calls,
			      conf,
			      NP,
			      plate,
			      fit.variance=NULL,
			      MIN.OBS=3,
			      envir,
			      chrom, P, DF.PRIOR=50, CONF.THR=0.99, trim=0, upperTail=TRUE, ...){
	if(length(ls(envir)) == 0) instantiateObjects(calls, NP, plate, envir, chrom)
	##require(genefilter)
	require(splines)
	require(Biobase)##needs rowMedians
	##will be updating these objects
	uplate <- get("uplate", envir)
	##AAA <- A
	##BBB <- B
	message("Sufficient statistics")
	if(missing(P)) P <- seq(along=uplate)
	for(p in P){
		cat(".")
		if(sum(plate == uplate[p]) < 10) next()
		oneBatch(plateIndex=p,
			 G=calls[, plate==uplate[p]],
			 A=A[, plate==uplate[p]],
			 B=B[, plate==uplate[p]],
			 envir=envir,
			 MIN.OBS=MIN.OBS,
			 DF.PRIOR=DF.PRIOR,
			 trim=trim, upperTail=upperTail)
	}
	message("\nEstimating coefficients")	
	for(p in P){
		cat(".")
		coefs(plateIndex=p, conf=conf[, plate==uplate[p]],
		      envir=envir, CONF.THR=CONF.THR)
	}
	message("\nAllele specific copy number")	
	for(p in P){
		cat(".")
		polymorphic(plateIndex=p,
			   A=A[, plate==uplate[p]],
			   B=B[, plate==uplate[p]],
			   envir=envir)
	}
	message("\nCopy number for nonpolymorphic probes...")	
	for(p in P){
		cat(".")
		nonpolymorphic(plateIndex=p,
			       NP=NP[, plate==uplate[p]],
			       envir=envir)
	}
}

nonpolymorphic <- function(plateIndex, NP, envir){
	require(genefilter) 
	p <- plateIndex
	plate <- get("plate", envir)
	uplate <- get("uplate", envir)	
	plates.completed <- get("plates.completed", envir)
	if(!plates.completed[p]) return()
	##snpflags <- get("snpflags", envir)
	##if(is.null(snpflags)){
	snpflags <- goodSnps(phi.thr=10, envir=envir, fewAA=10, fewBB=10)
	##assign("snpflags", snpflags, envir=envir)
	flagsA <- snpflags$A[, p]
	flagsB <- snpflags$B[, p]
	if(all(flagsA) | all(flagsB)) stop("all snps are flagged")
	nuA <- get("nuA", envir)
	nuB <- get("nuB", envir)
	phiA <- get("phiA", envir)
	phiB <- get("phiB", envir)
	uplate <- get("uplate", envir)
	sns <- get("sns", envir)
	muA.AA <- get("muA.AA", envir)
	muA.AB <- get("muA.AB", envir)
	muA.BB <- get("muA.BB", envir)
	muB.AA <- get("muB.AA", envir)
	muB.AB <- get("muB.AB", envir)
	muB.BB <- get("muB.BB", envir)
	NP.CT <- get("NP.CT", envir)
	nus <- get("nus", envir)
	phis <- get("phis", envir)
	NP.sds <- get("NP.sds", envir)
	##fit.variance <- get("fit.variance", envir)
	NP.CT <- get("NP.CT", envir)
	##---------------------------------------------------------------------------
	## Train on unflagged SNPs
	##---------------------------------------------------------------------------		
	plateInd <- plate == uplate[p]
	muA.AAp <- muA.AA[!flagsA, p]
	muB.BBp <- muB.BB[!flagsB, p]
		
	##From SNP data
	X <- cbind(1, log(c(muA.AAp, muB.BBp)))
	Y <- log(c(phiA[!flagsA, p], phiB[!flagsB, p]))
	##Y.nu <- log(c(nuA[!snpflags$A, p], nuB[!snpflags$B, p]))
	##Y <- cbind(Y.nu, Y.phi)
	betahat <- solve(crossprod(X), crossprod(X, Y))
	##Prediction
	mus <- rowMedians(NP, na.rm=TRUE)
	X <- cbind(1, log(mus))
	Yhat <- as.numeric(X %*% betahat)
	##NP.nus[, p] <- exp(Yhat[, 1])
	##NP.phis[, p] <- exp(Yhat[, 2])
	phi <- exp(Yhat)
	nu <- mus - 2*phi
	
	phis[, p] <- as.integer(phi)
	nus[, p] <- as.integer(nu)
	##plot(NP.phis[, p], NP.nus[, p], pch=".")
	CT <- 1/phi*(NP-nu)
	tmp <- matrix(as.integer(CT*100), nrow(CT), ncol(CT))
	NP.CT[, plateInd] <- tmp
	if(FALSE){
		##should be centered at 2
		tmp <- rowMeans(NP.CT[, plateInd])
		hist(tmp/100, breaks=200)
	}

	##---------------------------------------------------------------------------
	## For NA SNPs, treat as nonpolymorphic
	##---------------------------------------------------------------------------
	CA <- get("CA", envir)
	CB <- get("CB", envir)	
	tmpA <- CA[, plate==uplate[p]]
	tmpB <- CB[, plate==uplate[p]]	
	indexA <- which(rowSums(is.na(tmpA)) > 0)
	indexB <- which(rowSums(is.na(tmpB)) > 0)
	index <- union(indexA, indexB)
	if(length(index) > 0) browser()
	
	##---------------------------------------------------------------------------
	## this part could stand improvement
	##  - estimate the uncertainty
	##---------------------------------------------------------------------------
	
	##---------------------------------------------------------------------------
	## Estimate variance for copy number probes
	## VAR(CT) = var(1/phi * (NP - nu))
	##         = 1/phi^2 * VAR(NP)
	##         = 1/phi^2 * f(mu) * sigma
	## ** Phi is predicted from a regression using the row medians as predictors
	##   -- this may underestimate the uncertainty
	##---------------------------------------------------------------------------		
##	log.mus <- log2(mus)
##	log.var <- log2(rowVars(NP))
##	f <- predict(fit.variance, newdata=data.frame(log.mus=log.mus))
##	resid <- log.var - f
##	sds <- 1/phi*sqrt(2^(predict(fit.variance, newdata=data.frame(log.mus=log.mus+resid))))
##
##	robustSD <- function(X) diff(quantile(X, probs=c(0.16, (1-0.16)), na.rm=TRUE))/2
##	sample.sds <- apply(CT, 2, robustSD)
##	sample.sds <- sample.sds/median(sample.sds, na.rm=TRUE)
##	NP.sds[, plate==uplate[p]] <- sds %*% t(sample.sds)
	assign("phis", phis, envir)
	assign("nus", nus, envir)
	assign("NP.CT", NP.CT, envir)
##	assign("NP.sds", NP.sds, envir)
##	firstPass.NP <- list(phis=phis, nus=nus, CT=NP.CT, sds=sds,
##			     gns=gns, sns=sns, bns=bns)
##	firstPass.NP
}

oneBatch <- function(plateIndex, G, A, B, MIN.OBS=3, DF.PRIOR, envir, trim, upperTail, ...){
	p <- plateIndex
	plate <- get("plate", envir)
	AA <- G == 1
	AB <- G == 2
	BB <- G == 3
	Ns <- get("Ns", envir)
	Ns[, p, "AA"] <- rowSums(AA)
	Ns[, p, "AB"] <- rowSums(AB)
	Ns[, p, "BB"] <- rowSums(BB)
	AA[AA == FALSE] <- NA
	AB[AB == FALSE] <- NA
	BB[BB == FALSE] <- NA
	Ip <- array(NA, dim=c(nrow(G), ncol(G), 3, 2))
	dimnames(Ip) <- list(rownames(G), colnames(G), c("AA", "AB", "BB"), c("A", "B"))
	Ip[, , "AA", "A"] <- AA*A
	Ip[, , "AB", "A"] <- AB*A
	Ip[, , "BB", "A"] <- BB*A
	Ip[, , "AA", "B"] <- AA*B
	Ip[, , "AB", "B"] <- AB*B
	Ip[, , "BB", "B"] <- BB*B
	##---------------------------------------------------------------------------
	## Sufficient statistics (plate-specific)
	##---------------------------------------------------------------------------
	##These matrices instantiated in the calling function are updated
	muA.AA <- get("muA.AA", envir)
	muA.AB <- get("muA.AB", envir)
	muA.BB <- get("muA.BB", envir)
	muB.AA <- get("muB.AA", envir)
	muB.AB <- get("muB.AB", envir)
	muB.BB <- get("muB.BB", envir)
	tau2A <- get("tau2A", envir)
	tau2B <- get("tau2B", envir)
	sig2A <- get("sig2A", envir)
	sig2B <- get("sig2B", envir)
	corr <- get("corr", envir)
	corrA.BB <- get("corrA.BB", envir)  ## B allele
	corrB.AA <- get("corrB.AA", envir)

	AA.A <- AA*A
	AB.A <- AB*A
	BB.A <- BB*A
	AA.B <- AA*B
	AB.B <- AB*B
	BB.B <- BB*B
	if(trim > 0){
		##rowMedians is not robust enough when a variant is common
		##Try a trimmed rowMedian, trimming only one tail (must specify which)
		##  - for genotypes with 3 or more observations, exclude the upper X% 
		##        - one way to do this is to replace these observations with NAs
		##         e.g, exclude round(Ns * X%, 0) observations
		##  - for genotypes with fewer than 10 observations,
		##  - recalculate rowMedians
		replaceWithNAs <- function(x, trim, upperTail){
			##put NA's last if trimming the upperTail
			if(upperTail) decreasing <- TRUE else decreasing <- FALSE
			NN <- round(sum(!is.na(x)) * trim, 0)
			ix <- order(x, decreasing=decreasing, na.last=TRUE)[1:NN]
			x[ix] <- NA
			return(x)
		}
		##which rows should be trimmed
		rowsToTrim <- which(round(Ns[, p, "AA"] * trim, 0) > 0)
		##replace values in the tail of A with NAs
		AA.A[rowsToTrim, ] <- t(apply(AA.A[rowsToTrim, ], 1, replaceWithNAs, trim=trim, upperTail=upperTail))
		AA.B[rowsToTrim, ] <- t(apply(AA.B[rowsToTrim, ], 1, replaceWithNAs, trim=trim, upperTail=upperTail))
		rowsToTrim <- which(round(Ns[, p, "AB"] * trim, 0) > 0)
		AB.A[rowsToTrim, ] <- t(apply(AB.A[rowsToTrim, ], 1, replaceWithNAs, trim=trim, upperTail=upperTail))
		AB.B[rowsToTrim, ] <- t(apply(AB.B[rowsToTrim, ], 1, replaceWithNAs, trim=trim, upperTail=upperTail))
		rowsToTrim <- which(round(Ns[, p, "BB"] * trim, 0) > 0)
		BB.A[rowsToTrim, ] <- t(apply(BB.A[rowsToTrim, ], 1, replaceWithNAs, trim=trim, upperTail=upperTail))
		BB.B[rowsToTrim, ] <- t(apply(BB.B[rowsToTrim, ], 1, replaceWithNAs, trim=trim, upperTail=upperTail))

		##Should probably recompute the Ns -- change the
		##degrees of freedom
		Ns[, p, "AA"] <- rowSums(!is.na(AA.A))
		Ns[, p, "AB"] <- rowSums(!is.na(AB.A))
		Ns[, p, "BB"] <- rowSums(!is.na(BB.A))		
	}
	muA.AA[, p] <- rowMedians(AA.A, na.rm=TRUE)
	muA.AB[, p] <- rowMedians(AB.A, na.rm=TRUE)
	muA.BB[, p] <- rowMedians(BB.A, na.rm=TRUE)
	muB.AA[, p] <- rowMedians(AA.B, na.rm=TRUE)
	muB.AB[, p] <- rowMedians(AB.B, na.rm=TRUE)
	muB.BB[, p] <- rowMedians(BB.B, na.rm=TRUE)
	sigmaA <- matrix(NA, nrow(AA), 3)
	sigmaB <- matrix(NA, nrow(AA), 3)	
	colnames(sigmaB) <- colnames(sigmaA) <- c("AA", "AB", "BB")


	sigmaA[, "AA"] <- rowMAD(AA.A, na.rm=TRUE)
	sigmaA[, "AB"] <- rowMAD(AB.A, na.rm=TRUE)
	sigmaA[, "BB"] <- rowMAD(BB.A, na.rm=TRUE)
	sigmaB[, "AA"] <- rowMAD(AA.B, na.rm=TRUE)
	sigmaB[, "AB"] <- rowMAD(AB.B, na.rm=TRUE)
	sigmaB[, "BB"] <- rowMAD(BB.B, na.rm=TRUE)

##	sigmaA[, "AA"] <- rowSds(AA.A, na.rm=TRUE)
##	sigmaA[, "AB"] <- rowSds(AB.A, na.rm=TRUE)
##	sigmaA[, "BB"] <- rowSds(BB*A, na.rm=TRUE)
##	sigmaB[, "AA"] <- rowSds(AA*B, na.rm=TRUE)
##	sigmaB[, "AB"] <- rowSds(AB*B, na.rm=TRUE)
##	sigmaB[, "BB"] <- rowSds(BB*B, na.rm=TRUE)
	##shrink
	DF <- Ns[, p, ]-1
	DF[DF < 1] <- 1
	medsA <- apply(sigmaA, 2, "median", na.rm=TRUE)
	medsB <- apply(sigmaB, 2, "median", na.rm=TRUE)			
	for(m in 1:3){
		sigmaA[, m] <- (sigmaA[, m]*DF[, m] + medsA[m]*DF.PRIOR)/(DF.PRIOR+DF[, m])
		sigmaA[is.na(sigmaA[, m]), m] <- medsA[m]
		sigmaB[, m] <- (sigmaB[, m]*DF[, m] + medsB[m]*DF.PRIOR)/(DF.PRIOR+DF[, m])
		sigmaB[is.na(sigmaB[, m]), m] <- medsB[m]		
	}
	index.AA <- which(Ns[, p, "AA"] >= 2)
	index.AB <- which(Ns[, p, "AB"] >= 2)
	index.BB <- which(Ns[, p, "BB"] >= 2)
	
	x <- Ip[index.BB, , "BB", "A"]
	##tau2A[index.BB, p] <- rowVars(x, na.rm=TRUE)##var(log(IA)| BB)
	tau2A[index.BB, p] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2
	DF <- Ns[, p, "BB"]-1
	DF[DF < 1] <- 1
	med <- median(tau2A[, p], na.rm=TRUE)
	tau2A[, p] <- (tau2A[, p] * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
	tau2A[is.na(tau2A[, p]), p] <- med
	
	x <- Ip[index.BB, , "BB", "B"]
##	sig2B[index.BB, p] <- rowVars(log2(x), na.rm=TRUE) ##var(log(IB)|BB)
	sig2B[index.BB, p] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2	
	##system.time(tmp <- rowCovs(x, x, na.rm=TRUE)) ##var(log(IB)|BB)
	##system.time(tmp2 <- rowVars(log2(x), na.rm=TRUE))
	med <- median(sig2B[, p], na.rm=TRUE)
	sig2B[, p] <- (sig2B[, p] * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
	sig2B[is.na(sig2B[, p]), p] <- med
	
	x <- Ip[index.AA, , "AA", "B"]
	##tau2B[index.AA, p] <- rowCovs(x, x, na.rm=TRUE)##var(log(IB)| AA)
	tau2B[index.AA, p] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2		
	DF <- Ns[, p, "AA"]-1
	DF[DF < 1] <- 1
	med <- median(tau2B[, p], na.rm=TRUE)
	tau2B[, p] <- (tau2B[, p] * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
	tau2B[is.na(tau2B[, p]), p] <- med
	
	x <- Ip[index.AA, , "AA", "A"]
	##sig2A[index.AA, p] <- rowVars(log2(x), na.rm=TRUE)##var(log(IA)|AA)
	sig2A[index.AA, p] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2##var(log(IA)|AA)	
	##sig2A[index.AA, p] <- rowCovs(x, x, na.rm=TRUE) ##var(log(IA)|AA)
	med <- median(sig2A[, p], na.rm=TRUE)
	sig2A[, p] <- (sig2A[, p]*DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
	sig2A[is.na(sig2A[, p]), p] <- med	

	##estimate the correlation where there is at least 1 or more copies (AB genotypes)
	if(length(index.AB) > 0){ ##all homozygous is possible	
		x <- Ip[index.AB, , "AB", "A"]
		y <- Ip[index.AB, , "AB", "B"]
		corr[index.AB, p] <- rowCors(x, y, na.rm=TRUE)
		corr[corr < 0] <- 0
		DF <- Ns[, p, "AB"]-1
		DF[DF<1] <- 1
		med <- median(corr[, p], na.rm=TRUE)
		corr[, p] <- (corr[, p]*DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
		corr[is.na(corr[, p]), p] <- med
	}
	##estimate the correlation of the background errors and AA, BB genotypes
	backgroundB <- Ip[index.AA, , "AA", "B"]
	signalA <- Ip[index.AA, , "AA", "A"]
	corrB.AA[index.AA, p] <- rowCors(backgroundB, signalA, na.rm=TRUE)
	DF <- Ns[, p, "AA"]-1
	DF[DF < 1] <- 1
	med <- median(corrB.AA[, p], na.rm=TRUE)
	corrB.AA[, p] <- (corrB.AA[, p]*DF + med*DF.PRIOR)/(DF.PRIOR + DF)
	corrB.AA[is.na(corrB.AA[, p]), p] <- med

	backgroundA <- Ip[index.BB, , "BB", "A"]
	signalB <- Ip[index.BB, , "BB", "B"]
	corrA.BB[index.BB, p] <- rowCors(backgroundA, signalB, na.rm=TRUE)
	DF <- Ns[, p, "BB"]-1
	DF[DF < 1] <- 1
	med <- median(corrA.BB[, p], na.rm=TRUE)
	corrA.BB[, p] <- (corrA.BB[, p]*DF + med*DF.PRIOR)/(DF.PRIOR + DF)
	corrA.BB[is.na(corrA.BB[, p]), p] <- med	
	
	##---------------------------------------------------------------------------
	## Predict sufficient statistics for unobserved genotypes (plate-specific)
	##---------------------------------------------------------------------------
	correct.orderA <- muA.AA[, p] > muA.BB[, p]
	correct.orderB <- muB.BB[, p] > muB.AA[, p]
	if(length(index.AB) > 0){
		nobs <- rowSums(Ns[, p, ] >= MIN.OBS) == 3
	} else nobs <- rowSums(Ns[, p, c(1,3)] >= MIN.OBS) == 2
	index.complete <- which(correct.orderA & correct.orderB & nobs) ##be selective here
	size <- min(5000, length(index.complete))
	if(size == 5000) index.complete <- sample(index.complete, 5000)
	if(length(index.complete) < 200){
		warning("too few snps pass my criteria for predicting the sufficient statistics")
		##message("Plate has too few samples to call many genotypes...can we use prior information to predict")
		browser()
	}
	index.AA <- which(Ns[, p, "AA"] == 0 & Ns[, p, "AB"] >= MIN.OBS)
	index.AB <- which(Ns[, p, "AB"] == 0 & (Ns[, p, "AA"] >= MIN.OBS & Ns[, p, "BB"] >= MIN.OBS))
	index.BB <- which(Ns[, p, "BB"] == 0 & Ns[, p, "AB"] >= MIN.OBS)
	if(FALSE){
		png("figures/predictedMusBB%02d.png", width=800, height=500)
		par(mfrow=c(1, 2), las=1, mar=c(3, 1, 1, 1), oma=c(1, 2, 1, 1), pty="s")				
		plot(log(muA.BB[index.complete, 1]), log(muB.BB[index.complete, 1]), pch=".", ylim=c(5,12), xlim=c(5,12),
		     xlab="A", ylab="B", main="good")
		points(log(muA.AB[index.complete, 1]), log(muB.AB[index.complete, 1]), pch=".")
		points(log(muA.AA[index.complete, 1]), log(muB.AA[index.complete, 1]), pch=".")
		plot(log(muA.BB[index.BB, 1]), log(muB.BB[index.BB, 1]), pch=".", ylim=c(5,12), xlim=c(5,12),
		     xlab="A", ylab="B", main="missing BB")
		points(log(muA.AB[index.BB, 1]), log(muB.AB[index.BB, 1]), pch=".")
		points(log(muA.AA[index.BB, 1]), log(muB.AA[index.BB, 1]), pch=".")				
	}
	if(length(index.AA) > 0){
		X.AA <- cbind(1, muA.AB[index.complete, p], muA.BB[index.complete, p],
			      muB.AB[index.complete, p], muB.BB[index.complete, p])
		Y <- cbind(muA.AA[index.complete, p], muB.AA[index.complete, p])
		XtY <- crossprod(X.AA, Y)
		betahat <- solve(crossprod(X.AA), XtY)
		X.AA <- cbind(1, muA.AB[index.AA, p], muA.BB[index.AA, p],
			      muB.AB[index.AA, p], muB.BB[index.AA, p])
		musAA <- X.AA %*% betahat
		musAA[musAA <= 100] <- 100
		muA.AA[index.AA, p] <- musAA[, 1]
		muB.AA[index.AA, p] <- musAA[, 2]
	}
	if(length(index.AB) > 0){
		X.AB <- cbind(1, muA.AA[index.complete, p], muA.BB[index.complete, p],
			      muB.AA[index.complete, p], muB.BB[index.complete, p])
		Y <- cbind(muA.AB[index.complete, p], muB.AB[index.complete, p])
		XtY <- crossprod(X.AB, Y)
		betahat <- solve(crossprod(X.AB), XtY)
		X.AB <- cbind(1, muA.AA[index.AB, p], muA.BB[index.AB, p],
			      muB.AA[index.AB, p], muB.BB[index.AB, p])
		musAB <- X.AB %*% betahat
		musAB[musAB <= 100] <- 100					
		muA.AB[index.AB, p] <- musAB[, 1]
		muB.AB[index.AB, p] <- musAB[, 2]				
	}
	if(length(index.BB) > 0){
		X.BB <- cbind(1, muA.AA[index.complete, p], muA.AB[index.complete, p],
			      muB.AA[index.complete, p], muB.AB[index.complete, p])
		Y <- cbind(muA.BB[index.complete, p], muB.BB[index.complete, p])
		XtY <- crossprod(X.BB, Y)
		betahat <- solve(crossprod(X.BB), XtY)
		X.BB <- cbind(1, muA.AA[index.BB, p], muA.AB[index.BB, p],
			      muB.AA[index.BB, p], muB.AB[index.BB, p])
		musBB <- X.BB %*% betahat
		musBB[musBB <= 100] <- 100				
		muA.BB[index.BB, p] <- musBB[, 1]
		muB.BB[index.BB, p] <- musBB[, 2]								
	}
	if(FALSE){
		points(log(muA.BB[index.BB, 1]), log(muB.BB[index.BB, 1]), pch=".", col="blue")
		dev.off()
	}
	##missing two genotypes
	noAA <- Ns[, p, "AA"] <= 2
	noAB <- Ns[, p, "AB"] <= 2
	noBB <- Ns[, p, "BB"] <= 2
	##---------------------------------------------------------------------------
	## Two genotype clusters not observed -- would sequence help? (didn't seem that helpful)
	## 1 extract index of complete data
	## 2 Regress  mu1,mu3 ~ sequence + mu2
	## 3 Predict mu1*, mu3* for missing genotypes
	##---------------------------------------------------------------------------			
	if(sum(noAA & noAB) > 0){
		##predict other cluster centers using only the observed genotypes
		##predict AA:
		X <- cbind(1, muA.BB[index.complete, p], muB.BB[index.complete, p])
		Y <- cbind(muA.AA[index.complete, p],
			   muB.AA[index.complete, p],
			   muA.AB[index.complete, p],
			   muB.AB[index.complete, p])
		XtY <- crossprod(X, Y)
		betahat <- solve(crossprod(X), XtY)
		X <- cbind(1, muA.BB[noAA & noAB, p], muB.BB[noAA & noAB, p])
		mus <- X %*% betahat
		mus[mus <= 100] <- 100				
		muA.AA[noAA & noAB, p] <- mus[, 1]
		muB.AA[noAA & noAB, p] <- mus[, 2]
		muA.AB[noAA & noAB, p] <- mus[, 3]
		muB.AB[noAA & noAB, p] <- mus[, 4]				
	}
	if(FALSE){
		i <- which(noAA & noAB)
		par(mfrow=c(1, 1), las=1, mar=c(3, 3, 1, 1), pty="s")
		plot(log(muA.BB[i, 1]), log(muB.BB[i, 1]), pch=".", ylim=c(5,12), xlim=c(5,12),
		     xlab="A", ylab="B", main="missing BB")
		points(log(muA.AB[i, 1]), log(muB.AB[i, 1]), pch=".", col="blue")
		points(log(muA.AA[i, 1]), log(muB.AA[i, 1]), pch=".", col="red")
	}
	if(sum(noBB & noAB) > 0){
		##predict other cluster centers using only the observed genotypes
		X <- cbind(1, muA.AA[index.complete, p], muB.AA[index.complete, p])
		Y <- cbind(muA.AB[index.complete, p],
			   muB.AB[index.complete, p],
			   muA.BB[index.complete, p], #muA.AB[index.complete, p],
			   muB.BB[index.complete, p])#, muB.AB[index.complete, p])
		XtY <- crossprod(X, Y)
		betahat <- solve(crossprod(X), XtY)
		X <- cbind(1, muA.AA[noBB & noAB, p], muB.AA[noBB & noAB, p])
		mus <- X %*% betahat
		mus[mus <= 100] <- 100
		muA.AB[noBB & noAB, p] <- mus[, 1]
		muB.AB[noBB & noAB, p] <- mus[, 2]								
		muA.BB[noBB & noAB, p] <- mus[, 3]
		muB.BB[noBB & noAB, p] <- mus[, 4]				
		if(FALSE){
			i <- which(noBB & noAB)
			par(mfrow=c(1, 1), las=1, mar=c(3, 3, 1, 1), pty="s")
			plot(log(muA.BB[i, 1]), log(muB.BB[i, 1]), pch=".", ylim=c(5,12), xlim=c(5,12),
			     xlab="A", ylab="B", main="missing BB")
			points(log(muA.AB[i, 1]), log(muB.AB[i, 1]), pch=".", col="blue")
			points(log(muA.AA[i, 1]), log(muB.AA[i, 1]), pch=".", col="red")
		}				
	}
	dn.Ns <- dimnames(Ns)
	Ns <- array(as.integer(Ns), dim=dim(Ns))
	dimnames(Ns)[[3]] <- dn.Ns[[3]]

	assign("Ns", Ns, envir)	
	assign("muA.AA", muA.AA, envir)
	assign("muA.AB", muA.AB, envir)
	assign("muA.BB", muA.BB, envir)
	assign("muB.AA", muB.AA, envir)
	assign("muB.AB", muB.AB, envir)
	assign("muB.BB", muB.BB, envir)
	assign("sigmaA", sigmaA, envir)
	assign("sigmaB", sigmaB, envir)	
	assign("tau2A", tau2A, envir)
	assign("tau2B", tau2B, envir)
	assign("sig2A", sig2A, envir)
	assign("sig2B", sig2B, envir)
	assign("corr", corr, envir)
	assign("corrB.AA", corrB.AA, envir)
	assign("corrA.BB", corrA.BB, envir)		
	plates.completed <- get("plates.completed", envir)
	plates.completed[p] <- TRUE
	assign("plates.completed", plates.completed, envir)
}

coefs <- function(plateIndex, conf, MIN.OBS=3, envir, CONF.THR=0.99){
	p <- plateIndex
	##if(p == 5) browser()
	plates.completed <- get("plates.completed", envir)
	if(!plates.completed[p]) return()
	plate <- get("plate", envir)
	nuA <- get("nuA", envir)
	nuB <- get("nuB", envir)
	nuA.se <- get("nuA.se", envir)
	nuB.se <- get("nuB.se", envir)
	phiA <- get("phiA", envir)
	phiB <- get("phiB", envir)
	phiA.se <- get("phiA.se", envir)
	phiB.se <- get("phiB.se", envir)
	muA.AA <- get("muA.AA", envir)
	muA.AB <- get("muA.AB", envir)
	muA.BB <- get("muA.BB", envir)
	muB.AA <- get("muB.AA", envir)
	muB.AB <- get("muB.AB", envir)
	muB.BB <- get("muB.BB", envir)
	Ns <- get("Ns", envir)
	uplate <- get("uplate", envir)
	sigmaA <- get("sigmaA", envir)
	sigmaB <- get("sigmaB", envir)	
	##fit.variance <- get("fit.variance", envir)
	IA <- cbind(muA.AA[, p], muA.AB[, p], muA.BB[, p])
	IB <- cbind(muB.AA[, p], muB.AB[, p], muB.BB[, p])
	NOHET <- mean(Ns[, p, 2], na.rm=TRUE) < 0.05
	##predict missing variance (do only once)
	is.complete <- rowSums(Ns[, p, ] >= 1) == 3
	if(NOHET) is.complete <- rowSums(Ns[, p, c(1,3)] >= MIN.OBS) == 2		
	correct.orderA <- muA.AA[, p] > muA.BB[, p]
	correct.orderB <- muB.BB[, p] > muB.AA[, p]
	highConf <- 1-exp(-conf/1000)
	confInd <- rowMeans(highConf) > CONF.THR
	keep <- confInd & is.complete & correct.orderA & correct.orderB
##	if(is.null(fit.variance)){
##		log.sigma2 <- as.numeric(log2(sigma2A[keep, ]))
##		log.mus <- as.numeric(log2(IA[keep, ]))
##		##log.mus <- as.numeric(log(mus[, p, , "A"]))
##		nodups <- which(!duplicated(log.sigma2))
##		i <- sample(nodups, 5000)
##		include <- intersect(which(!is.na(log.sigma2) & !is.na(log.mus)), nodups)
##		fit.variance <- lm(log.sigma2 ~ ns(log.mus, 3), subset=sample(include, 5000))
##		assign("fit.variance", fit.variance, envir)
##	} 
##	i <- which(rowSums(is.na(sigma2A)) > 0)
##	log.mus <- as.numeric(log2(IA[i, ]))
##	predictS2 <- predict(fit.variance, newdata=data.frame(log.mus=log.mus))
##	predictS2 <- matrix(2^(predictS2), ncol=3)
##	sigma2A[i, ] <- predictS2
	##---------------------------------------------------------------------------
	## Estimate nu and phi
	##---------------------------------------------------------------------------
	Np <- Ns[, p, ]
	sigma2A <- sigmaA^2
	sigma2B <- sigmaB^2	
	if(NOHET){
		##only homozygous
		Np <- Np[, -2]
		Np[Np < 1] <- 1
		IA <- IA[, c(1, 3)]
		sigma2A <- sigma2A[, c(1,3)]
	}else 	Np[Np < 1] <- 1
	##Using MAD instead of 
	V <- sigma2A/Np
	W <- sqrt(1/V)
	Ystar <- IA*W
	is.complete <- rowSums(is.na(W)) == 0
	##is.complete <- is.complete & correct.orderA & correct.orderB & confInd & notmissing
	nuphiAllele(allele="A", Ystar=Ystar[is.complete, ], W=W[is.complete, ], envir=envir, p=p)
	nuA[is.complete, p] <- nu
	phiA[is.complete, p] <- phi
	nuA.se[is.complete, p] <- nu.se
	phiA.se[is.complete, p] <- phi.se
	
	if(NOHET){
		IB <- IB[, c(1, 3)]
		sigma2B <- sigma2B[, c(1,3)]
	}
##	i <- which(rowSums(is.na(sigma2B)) > 0)
##	log.mus <- as.numeric(log2(IB[i, ]))
##	predictS2 <- predict(fit.variance, newdata=data.frame(log.mus=log.mus))
##	predictS2 <- matrix(2^(predictS2), ncol=ncol(sigma2B))
##	sigma2B[i, ] <- predictS2		
	V <- sigma2B/Np
	W <- sqrt(1/V)
	Ystar <- IB*W
	is.complete <- rowSums(is.na(W)) == 0
	nuphiAllele(allele="B", Ystar=Ystar[is.complete, ], W=W[is.complete, ], envir=envir, p=p)
	nuB[is.complete, p] <- nu
	phiB[is.complete, p] <- phi
	nuB.se[is.complete, p] <- nu.se
	phiB.se[is.complete, p] <- phi.se
	phiA <- matrix(as.integer(phiA), nrow(phiA), ncol(phiA))
	phiB <- matrix(as.integer(phiB), nrow(phiA), ncol(phiA))

	assign("nuA", nuA, envir)
	assign("nuB", nuB, envir)
	assign("phiA", phiA, envir)
	assign("phiB", phiB, envir)
	assign("nuA.se", nuA.se, envir)
	assign("nuB.se", nuB.se, envir)
	assign("phiA.se", phiA.se, envir)
	assign("phiB.se", phiB.se, envir)
}

polymorphic <- function(plateIndex, A, B, envir){
	p <- plateIndex
	plates.completed <- get("plates.completed", envir)
	if(!plates.completed[p]) return()
	plate <- get("plate", envir)	
	nuA <- get("nuA", envir)
	nuB <- get("nuB", envir)
	phiA <- get("phiA", envir)
	phiB <- get("phiB", envir)
	uplate <- get("uplate", envir)
	muA.AA <- get("muA.AA", envir)
	muA.AB <- get("muA.AB", envir)
	muA.BB <- get("muA.BB", envir)
	muB.AA <- get("muB.AA", envir)
	muB.AB <- get("muB.AB", envir)
	muB.BB <- get("muB.BB", envir)
	sd.CT <- get("sd.CT", envir)
	sigmaA <- get("sigmaA", envir)
	sigmaB <- get("sigmaB", envir)	
	Ns <- get("Ns", envir)
	CA <- get("CA", envir)
	CB <- get("CB", envir)
	
	NOHET <- mean(Ns[, p, 2], na.rm=TRUE) < 0.05
	if(NOHET){
		IA <- cbind(muA.AA[, p], muA.BB[, p])
		IB <- cbind(muB.AA[, p], muB.BB[, p])		
		sigma2A <- cbind(sigmaA[, c(1, 3)]^2)
		sigma2B <- cbind(sigmaB[, c(1, 3)]^2)		
##		sigma2A <- cbind(Sigma.AA[, "varA", p], Sigma.BB[, "varA", p])
##		colnames(sigma2A) <- c("AA", "BB")
##		sigma2B <- cbind(Sigma.AA[, "varB", p], Sigma.BB[, "varB", p])
##		colnames(sigma2B) <- c("AA", "BB")
	} else{
		IA <- cbind(muA.AA[, p], muA.AB[, p], muA.BB[, p])
		IB <- cbind(muB.AA[, p], muB.AB[, p], muB.BB[, p])		
		sigma2A <- sigmaA^2
		sigma2B <- sigmaB^2		
##		sigma2A <- cbind(Sigma.AA[, "varA", p], Sigma.AB[, "varA", p], Sigma.BB[, "varA", p])
##		colnames(sigma2A) <- c("AA", "AB", "BB")
##		sigma2B <- cbind(Sigma.AA[, "varB", p], Sigma.AB[, "varB", p], Sigma.BB[, "varB", p])
##		colnames(sigma2B) <- c("AA", "AB", "BB")
	}
	##sigma2A <- get("sigma2A", envir)
	##sigma2B <- get("sigma2B", envir)
##	fit.variance <- get("fit.variance", envir)

	##---------------------------------------------------------------------------
	## Estimate CA, CB
	##---------------------------------------------------------------------------
	phiA[phiA < 1] <- 1	
	tmp <- (1/phiA[, p])*(A - nuA[, p])
	tmp <- matrix(as.integer(tmp*100), nrow(tmp), ncol(tmp))
	if(FALSE) hist(rowMeans(tmp/100), breaks=1000, xlim=c(-2, 8))
	CA[, plate==uplate[p]] <- tmp
	phiB[phiB < 1] <- 1
	tmp <- (1/phiB[, p])*(B - nuB[, p])
	tmp <- matrix(as.integer(tmp*100), nrow(tmp), ncol(tmp))
	if(FALSE) hist(rowMeans(tmp/100), breaks=1000, xlim=c(-2, 8))		
	CB[, plate==uplate[p]] <- tmp

	nuA[nuA < 1] <- 1
	nuB[nuB < 1] <- 1

	tmpA <- CA[, plate==uplate[p]]
	tmpB <- CB[, plate==uplate[p]]	
	indexA <- which(rowSums(is.na(tmpA)) > 0)
	indexB <- which(rowSums(is.na(tmpB)) > 0)
	index <- union(indexA, indexB)
	if(length(index) > 0) browser()	
	##---------------------------------------------------------------------------
	## Estimate var(CA), var(CB), var(CA+CB)
	##---------------------------------------------------------------------------
	##var(CA) = 1/phiA^2*var(IA - nuA)
	##        = 1/phiA^2*var(IA)
	##        = 1/phiA^2*f(mu)
##	log.musA <- as.numeric(log2(IA))		
##	log.sigma2A <- as.numeric(log2(sigma2A))
##	log.musB <- as.numeric(log2(IB))
##	log.sigma2B <- as.numeric(log2(sigma2B))
##	f <- predict(fit.variance, newdata=data.frame(log.mus=log.musA))
##	resid <- matrix(log.sigma2A - f, ncol=ncol(sigma2A))
##	ls2 <- rowMedians(resid, na.rm=TRUE)		
##	log.musA <- matrix(log.musA, ncol=ncol(sigma2A))
##	tmp <- apply(log.musA + ls2, 2, function(log.mus, fit) sqrt(2^(predict(fit, newdata=data.frame(log.mus=log.mus)))), fit=fit.variance)
##	sdA <- 1/phiA[, p]*rowMedians(tmp)
##	##If all three points are below the spline, the 3 residuals
##	##will be negative.  adding the median residual to the log mus
##	##protects against over/under estimating the variance
##	##sdA[tmpIndex, p, ] <- apply(log.musA[tmpIndex, ] + ls2, 2,
##	##function(log.mus, fit) sqrt(exp(predict(fit.logVariance,
##	##newdata=data.frame(log.mus=log.mus)))))
##	robustSD <- function(X) diff(quantile(X, probs=c(0.16, (1-0.16)), na.rm=TRUE))/2	
##	sample.sd <- apply(CA[, plate==uplate[p]]/100, 2, robustSD)
##	sd.ca <- sdA %*% t(sample.sd)
##	##where missing, use the residual from the genotype with the most observations
##	f <- predict(fit.variance, newdata=data.frame(log.mus=log.musB))
##	resid <- matrix(log.sigma2B - f, ncol=ncol(sigma2B))
##	ls2 <- rowMedians(resid, na.rm=TRUE)		
##	log.musB <- matrix(log.musB, ncol=ncol(sigma2B))
##	tmp <- apply(log.musB + ls2, 2, function(log.mus, fit) sqrt(2^(predict(fit, newdata=data.frame(log.mus=log.mus)))), fit=fit.variance)
##	sdB <- 1/phiB[, p]*rowMedians(tmp)
##	sample.sd <- apply(CB[, plate==uplate[p]]/100, 2, robustSD)
##	sd.cb <- sdB %*% t(sample.sd)
##	covar <- rowCovs(1/phiA[, p]*log.musA, 1/phiB[, p]*log.musB)
##	##tmp1 <- sqrt(sdA^2 + sdB^2 - 2*covar) ##probe-specific variance
##	tmp <- sqrt(sd.ca^2 + sd.cb^2 - 2*covar) ##probe-specific variance
##	sd.CT[, plate==uplate[p]] <- matrix(as.integer(100*tmp), nrow(tmp), ncol(tmp))
	##another approach:
	##tmp <- apply(log(A), 2, function(log.mus, fit) sqrt(exp(predict(fit, newdata=data.frame(log.mus=log.mus)))), fit=fit.variance)	
	##tmp2 <- 1/phiA[, p]*tmp
##	CA <- matrix(as.integer(CA*100), nrow(CA), ncol(CA))
##	CB <- matrix(as.integer(CB*100), nrow(CB), ncol(CB))
##	sd.CT <- matrix(as.integer(sd.CT*100), nrow(sd.CT), ncol(sd.CT))	
	assign("CA", CA, envir)
	assign("CB", CB, envir)
##	assign("sd.CT", sd.CT, envir)	
	##---------------------------------------------------------------------------
	## nonpolymorphic probes
	##---------------------------------------------------------------------------
}
