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

generateX <- function(w, X) as.numeric(diag(w) %*% X)
generateIXTX <- function(x, nrow=3) {
	X <- matrix(x, nrow=nrow)
	XTX <- crossprod(X)
	solve(XTX)
}

nuphiAllele <- function(p, allele, Ystar, W, envir){
	Ns <- envir[["Ns"]]
	CHR <- envir[["chrom"]]
	nuA <- envir[["nuA"]]
	nuB <- envir[["nuB"]]
	nuA.se <- envir[["nuA.se"]]
	nuB.se <- envir[["nuB.se"]]
	phiA <- envir[["phiA"]]
	phiB <- envir[["phiB"]]
	phiA.se <- envir[["phiA.se"]]
	phiB.se <- envir[["phiB.se"]]
	if(CHR==23){
		phiAx <- envir[["phiAx"]]
		phiBx <- envir[["phiBx"]]
	}
	complete <- rowSums(is.na(W)) == 0
	NOHET <- mean(Ns[, p, "AB"], na.rm=TRUE) < 0.05
	if(missing(allele)) stop("must specify allele")
	if(CHR == 23){
		gender <- envir[["gender"]]
		if(allele == "A") X <- cbind(1, c(1, 0, 2, 1, 0), c(0, 1, 0, 1, 2))
		if(allele == "B") X <- cbind(1, c(0, 1, 0, 1, 2), c(1, 0, 2, 1, 0))			
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
	Xstar <- apply(W, 1, generateX, X)
	IXTX <- apply(Xstar, 2, generateIXTX, nrow=nrow(X))
	if(CHR == 23){
		betahat <- matrix(NA, 3, nrow(Ystar))
		ses <- matrix(NA, 3, nrow(Ystar))		
	} else{
		betahat <- matrix(NA, 2, nrow(Ystar))
		ses <- matrix(NA, 2, nrow(Ystar))
	}
	for(i in 1:nrow(Ystar)){
		betahat[, i] <- crossprod(matrix(IXTX[, i], ncol(X), ncol(X)), crossprod(matrix(Xstar[, i], nrow=nrow(X)), Ystar[i, ]))
		ssr <- sum((Ystar[i, ] - matrix(Xstar[, i], nrow(X), ncol(X)) %*% matrix(betahat[, i], ncol(X), 1))^2)
		ses[, i] <- sqrt(diag(matrix(IXTX[, i], ncol(X), ncol(X)) * ssr))
	}
	if(allele == "A"){
		nuA[complete, p] <- betahat[1, ]
		phiA[complete, p] <- betahat[2, ]
		nuA.se[complete, p] <- ses[1, ]
		phiA.se[complete, p] <- ses[2, ]
		envir[["nuA"]] <- nuA
		envir[["phiA"]] <- phiA
		envir[["nuA.se"]] <- nuA.se
		envir[["phiA.se"]] <- phiA.se
		if(CHR == 23){
			phiAx[complete, p] <- betahat[3, ]
			envir[["phiAx"]] <- phiAx
		}
	}
	if(allele=="B"){
		nuB[complete, p] <- betahat[1, ]
		phiB[complete, p] <- betahat[2, ]
		nuB.se[complete, p] <- ses[1, ]
		phiB.se[complete, p] <- ses[2, ]
		envir[["nuB"]] <- nuB
		envir[["phiB"]] <- phiB
		envir[["nuB.se"]] <- nuB.se
		envir[["phiB.se"]] <- phiB.se
		if(CHR == 23){
			phiBx[complete, p] <- betahat[3, ]
			envir[["phiBx"]] <- phiBx
		}
	}
}

celDates <- function(celfiles){
	if(!all(file.exists(celfiles))) stop("1 or more cel file does not exist")
	celdates <- vector("character", length(celfiles))
	celtimes <- vector("character", length(celfiles))
	for(i in seq(along=celfiles)){
		if(i %% 100 == 0) cat(".")
		tmp <- read.celfile.header(celfiles[i], info="full")$DatHeader
		tmp <- strsplit(tmp, "\ +")
		celdates[i] <- tmp[[1]][6]
		celtimes[i] <- tmp[[1]][7]
	}
	tmp <- paste(celdates, celtimes)
	celdts <- strptime(tmp, "%m/%d/%y %H:%M:%S")
	return(celdts)
}

cnrma <- function(filenames, cdfName="genomewidesnp6", sns, seed=1, verbose=FALSE){
	if(cdfName != "genomewidesnp6") stop("Only genomewidesnp6 supported at this time")
	## BC: 03/14/09
	## getting pkgname from cdfName, in the future this might be useful
	## as the method might be implemented for other platforms
	pkgname <- getCrlmmAnnotationName(cdfName)
	require(pkgname, character.only=TRUE) || stop("Package ", pkgname, " not available")
	if (missing(sns)) sns <- basename(filenames)
	## Loading data in .crlmmPkgEnv and extracting from there
        loader("npProbesFid.rda", .crlmmPkgEnv, pkgname)
##	data("npProbesFid", package=pkgname, envir=.crlmmPkgEnv)
	fid <- getVarInEnv("npProbesFid")
	set.seed(seed)
	idx2 <- sample(length(fid), 10^5) ##for skewness. no need to do everything
	SKW <- vector("numeric", length(filenames))
	NP <- matrix(NA, length(fid), length(filenames))
	verbose <- TRUE
	if(verbose){
		message("Processing ", length(filenames), " files.")
		if (getRversion() > '2.7.0') pb <- txtProgressBar(min=0, max=length(filenames), style=3)
	}
	##load reference distribution obtained from hapmap
        loader("1m_reference_cn.rda", .crlmmPkgEnv, pkgname)
##	data(list="1m_reference_cn", package="genomewidesnp6Crlmm", envir=.crlmmPkgEnv)
	reference <- getVarInEnv("reference")
	for(i in seq(along=filenames)){
		y <- as.matrix(read.celfile(filenames[i], intensity.means.only=TRUE)[["INTENSITY"]][["MEAN"]][fid])
		x <- log2(y[idx2])
		SKW[i] <- mean((x-mean(x))^3)/(sd(x)^3)
		rm(x)
		NP[, i] <- as.integer(normalize.quantiles.use.target(y, target=reference))
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
	negativePhis <- phiA < phi.thr | phiB < phi.thr
	negativeCoef <- negativeNus | negativePhis

	notfinitePhi <- !is.finite(phiA) | !is.finite(phiB)
	flags <- negativeCoef | notfinitePhi
	return(flags)
}

goodSnps <- function(phi.thr, envir, fewAA=20, fewBB=20){
	Ns <- get("Ns", envir)
	flags <- getFlags(phi.thr=phi.thr, envir)
	fewAA <- Ns[, , "AA"] < fewAA
	fewBB <- Ns[, , "BB"] < fewBB
	flagsA <- flags | fewAA
	flagsB <- flags | fewBB
	flags <- list(A=flagsA, B=flagsB)
	return(flags)
}

instantiateObjects <- function(calls, NP, plate, envir, chrom, A, B,
			       gender, SNRmin=5, SNR,
                               pkgname="genomewidesnp6Crlmm"){
	envir[["chrom"]] <- chrom
	CHR_INDEX <- paste(chrom, "index", sep="")
        fname <- paste(CHR_INDEX, ".rda", sep="")
        loader(fname, .crlmmPkgEnv, pkgname)
        index <- get("index", envir=.crlmmPkgEnv)
##	data(list=CHR_INDEX, package="genomewidesnp6Crlmm")
	A <- A[index[[1]], SNR > SNRmin]
	B <- B[index[[1]], SNR > SNRmin]
	calls <- calls[index[[1]], SNR > SNRmin]
	conf <- conf[index[[1]], SNR > SNRmin]
	NP <- NP[index[[2]], SNR > SNRmin]
	plate <- plate[SNR > SNRmin]
	uplate <- unique(plate)
	SNR <- SNR[SNR > SNRmin]
	
	envir[["uplate"]] <- uplate
	envir[["plate"]] <- plate	
	envir[["NP"]] <- NP
	envir[["A"]] <- A
	envir[["B"]] <- B
	envir[["calls"]] <- calls
	envir[["conf"]] <- conf
	snps <- rownames(calls)
	cnvs <- rownames(NP)
	sns <- basename(colnames(calls))
	stopifnot(identical(colnames(calls), colnames(NP)))
	envir[["sns"]] <- sns
	envir[["snps"]] <- snps
	envir[["cnvs"]] <- cnvs

	if(chrom == 23){
		if(is.null(gender)){
			message("Estimating gender")
			XMedian <- apply(log2(A[, , drop=FALSE]) + log2(B[, , drop=FALSE]), 2, median)/2
			gender <- kmeans(XMedian, c(min(XMedian[SNR>SNRmin]), max(XMedian[SNR>SNRmin])))[["cluster"]]			
			##gender <- getGender(res)
			gender[gender==2] <- "female"
			gender[gender=="1"] <- "male"
			envir[["gender"]] <- gender
		} else envir[["gender"]] <- gender
		phiAx <- matrix(NA, nrow(calls), length(uplate))
		envir[["phiAx"]] <- phiAx
		envir[["phiBx"]] <- phiAx
	}
	CA <- CB <- matrix(NA, nrow(calls), ncol(calls))
	envir[["CA"]] <- CA
	envir[["CB"]] <- CB
	
	Ns <- array(NA, dim=c(nrow(calls), length(uplate), 5))
	dimnames(Ns)[[3]] <- c("A", "B", "AA", "AB", "BB")

	envir[["Ns"]] <- envir[["muB"]] <- envir[["muA"]] <- Ns
	envir[["vB"]] <- envir[["vA"]] <- Ns

	CT.sds <- CT <- matrix(NA, nrow(NP), length(sns))
	##NP.CT <- matrix(NA, nrow(NP), ncol(NP))
	##NP.sds <- matrix(NA, nrow(NP), ncol(NP))
	envir[["CT"]] <- CT
	envir[["CT.sds"]] <- CT.sds

	##assign("NP.CT", NP.CT, envir)
	##assign("NP.sds", NP.sds, envir)
	nuT <- matrix(NA, nrow(NP), length(uplate))
	phiT <- nuT
	envir[["nuT"]] <- nuT
	envir[["phiT"]] <- phiT
	##assign("nus", nus, envir=envir)  
	##assign("phis", nus, envir=envir)

	plates.completed <- rep(FALSE, length(uplate))
	envir[["plates.completed"]] <- plates.completed

	steps <- rep(FALSE, 4)
	names(steps) <- c("suffStats", "coef", "snp-cn", "np-cn")
	envir[["steps"]] <- steps

	snpflags <- matrix(FALSE, length(snps), length(uplate))
	npflags <- matrix(FALSE, length(cnvs), length(uplate))
	##assign("snpflags", snpflags, envir=envir)
	##assign("npflags", npflags, envir=envir)
	envir[["snpflags"]] <- snpflags
	envir[["npflags"]] <- npflags

	tau2A <- matrix(NA, nrow(calls), length(uplate))
	envir[["tau2A"]] <- tau2A
	envir[["tau2B"]] <- tau2A
	envir[["sig2A"]] <- tau2A
	envir[["sig2B"]] <- tau2A
	sig2T <- matrix(NA, nrow(NP), ncol(NP))
	envir[["sig2T"]] <- sig2T
	envir[["corr"]] <- tau2A
	envir[["corrA.BB"]] <- tau2A
	envir[["corrB.AA"]] <- tau2A
	envir[["nuB"]] <- envir[["nuA"]] <- tau2A
	envir[["phiB"]] <- envir[["phiA"]] <- tau2A
	envir[["nuB.se"]] <- envir[["nuA.se"]] <- tau2A
	envir[["phiB.se"]] <- envir[["phiA.se"]] <- tau2A
	normal <- matrix(TRUE, nrow(A), ncol(A))
	normalNP <- matrix(TRUE, nrow(NP), ncol(NP))	
	envir[["normal"]] <- normal
	envir[["normalNP"]] <- normalNP
}



computeCopynumber <- function(chrom,
			      A,
			      B,
			      calls,
			      conf,
			      NP,
			      plate,
			      MIN.OBS=1,
			      envir,
			      P,
			      DF.PRIOR=50,
			      CONF.THR=0.99,
			      bias.adj=FALSE,
			      priorProb,
			      gender=NULL,
			      SNR,
			      SNRmin=5, seed=123, verbose=TRUE, ...){
	set.seed(seed)
	if(missing(chrom)) stop("must specify chromosome")
	if(length(ls(envir)) == 0) {
		instantiateObjects(calls=calls, NP=NP, plate=plate,
				   envir=envir, chrom=chrom, A=A, B=B,
				   gender=gender, SNR=SNR, SNRmin=SNRmin)
	}
	plate <- envir[["plate"]]
	uplate <- envir[["uplate"]]	
	calls <- envir[["calls"]]
	A <- envir[["A"]]
	B <- envir[["B"]]
	conf <- envir[["conf"]]
	NP <- envir[["NP"]]
	if(bias.adj){
		##assign uniform priors for total copy number states
		if(missing(priorProb)) priorProb <- rep(1/4, 4)
		envir[["steps"]] <- rep(FALSE, 4)
	}
	##will be updating these objects
	message("Sufficient statistics")
	if(missing(P)) P <- seq(along=uplate)
	steps <- envir[["steps"]]
	if(!steps[1]){
		for(p in P){
			cat(".")
			if(sum(plate == uplate[p]) < 10) next()
			J <- plate==uplate[p]
			oneBatch(plateIndex=p,
				 A=A[, J],
				 B=B[, J],
				 calls=calls[, J],
				 conf=conf[, J],
				 gender=NULL,
				 NP[, J],
				 plate[J],
				 MIN.OBS=1,
				 envir=envir,
				 DF.PRIOR=DF.PRIOR,
				 CONF.THR=CONF.THR,
				 bias.adj=bias.adj,
				 priorProb=priorProb,...)
		}
		steps[1] <- TRUE
		envir[["steps"]] <- steps
	}
	if(!steps[2]){
		message("\nEstimating coefficients")	
		for(p in P){
			cat(".")
			coefs(plateIndex=p, conf=conf[, plate==uplate[p]],
			      envir=envir, CONF.THR=CONF.THR, MIN.OBS=MIN.OBS)
		}
		steps[2] <- TRUE
		envir[["steps"]] <- steps		
	}
	if(!steps[3]){
		message("\nAllele specific copy number")	
		for(p in P){
			cat(".")
			polymorphic(plateIndex=p,
				    A=A[, plate==uplate[p]],
				    B=B[, plate==uplate[p]],
				    envir=envir)
		}
		steps[3] <- TRUE
		envir[["steps"]] <- steps						
	}
	if(!steps[4]){
		message("\nCopy number for nonpolymorphic probes...")	
		for(p in P){
			cat(".")
			nonpolymorphic(plateIndex=p,
				       NP=NP[, plate==uplate[p]],
				       envir=envir)
		}
		steps[4] <- TRUE
		envir[["step"]] <- steps
	}
}

nonpolymorphic <- function(plateIndex, NP, envir, CONF.THR=0.99, DF.PRIOR=50, pkgname="genomewidesnp6Crlmm"){
	p <- plateIndex
	CHR <- envir[["chrom"]]
	plate <- envir[["plate"]]
	uplate <- envir[["uplate"]]
	plates.completed <- envir[["plates.completed"]]
	if(!plates.completed[p]) return()
	snpflags <- goodSnps(phi.thr=2^6, envir=envir, fewAA=10, fewBB=10)
	flagsA <- snpflags$A[, p]
	flagsB <- snpflags$B[, p]
	if(all(flagsA) | all(flagsB)) stop("all snps are flagged")
	nuA <- envir[["nuA"]][, p]
	nuB <- envir[["nuB"]][, p]
	phiA <- envir[["phiA"]][, p]
	phiB <- envir[["phiB"]][, p]
	uplate <- envir[["uplate"]]
	sns <- envir[["sns"]]
	muA <- envir[["muA"]][, p, ]
	muB <- envir[["muB"]][, p, ]
	nuT <- envir[["nuT"]]
	phiT <- envir[["phiT"]]
	sig2T <- envir[["sig2T"]]
	A <- envir[["A"]][, plate==uplate[p]]
	B <- envir[["B"]][, plate==uplate[p]]
	CA <- envir[["A"]][, plate==uplate[p]]
	CB <- envir[["B"]][, plate==uplate[p]]
	if(CHR == 23){
		phiAx <- envir[["phiAx"]][, p]
		phiBx <- envir[["phiBx"]][, p]
	}
	##---------------------------------------------------------------------------
	## Train on unflagged SNPs
	##---------------------------------------------------------------------------
	##Might be best to train using the X chromosome, since for the
	##X phi and nu have been adjusted for cross-hybridization
	plateInd <- plate == uplate[p]
	##muA <- muA[!flagsA, p, c("A", "AA")]
	##muB <- muB[!flagsB, p, c("B", "BB")]
	muA <- muA[!flagsA, "AA"]
	muB <- muB[!flagsB, "BB"]
	X <- cbind(1, log2(c(muA, muB)))
	Y <- log2(c(phiA[!flagsA], phiB[!flagsB]))
	if(nrow(X) > 5000) ix <- sample(1:nrow(X), 5000)
	##X <- X[ix, ]
	##Y <- Y[ix]
	betahat <- solve(crossprod(X[ix, ]), crossprod(X[ix, ], Y[ix]))
##	Yhat <- X%*%betahat
##	phihat <- 2^Yhat
##	nuhat <- 2^X[, 2] - 2*phihat
##	nuAB <- c(nuA[!flagsA], nuB[!flagsB])
##	plot(log2(nuhat), log2(nuAB), pch=".")
##	plot(log2(nuhat)-log2(nuAB), pch=".")
##	hist(log2(nuAB))
	##plot(Y-Yhat, pch=".")
	##plot(Y, Yhat, pch=".")
	##calculate R2
	if(CHR == 23){
		cnvs <- envir[["cnvs"]]
                loader("cnProbes", pkgname, .crlmmPkgEnv)
                cnProbes <- get("cnProbes", envir=.crlmmPkgEnv)
##		data(cnProbes, package="genomewidesnp6Crlmm")
		cnProbes <- cnProbes[match(cnvs, rownames(cnProbes)), ]
		par <- cnProbes[, "position"] < 2709520 | (cnProbes[, "position"] > 154584237 & cnProbes[, "position"] < 154913754)
		gender <- envir[["gender"]]
		mu1 <- rowMedians(NP[, gender=="male"], na.rm=TRUE)
		mu2 <- rowMedians(NP[, gender=="female"], na.rm=TRUE)
		mus <- log(cbind(mu1, mu2))
		X <- cbind(1, mus[, 1])
		##For build Hg18
		##http://genome.ucsc.edu/cgi-bin/hgGateway
		##pseudo-autosomal regions on X
		##chrX:1-2,709,520 and chrX:154584237-154913754, respectively
		Yhat1 <- as.numeric(X %*% betahat)
		X <- cbind(1, mus[, 2])		
		Yhat2 <- as.numeric(X %*% betahat)
		phi1 <- exp(Yhat1)
		phi2 <- exp(Yhat2)
		nu1 <- exp(mus[, 1]) - phi1
		nu1[par] <- exp(mus[par, 1]) - 2*phi1[par]
		nu2 <- exp(mus[, 2]) - 2*phi2
		CT1 <- 1/phi1*(NP[, gender=="male"]-nu1)
		CT2 <- 1/phi2*(NP[, gender=="female"]-nu2)
		CT1 <- matrix(as.integer(100*CT1), nrow(CT1), ncol(CT1))
		CT2 <- matrix(as.integer(100*CT2), nrow(CT2), ncol(CT2))
		CT <- envir[["CT"]]
		CT[, plate==uplate[p] & gender=="male"] <- CT1
		CT[, plate==uplate[p] & gender=="female"] <- CT2
		envir[["CT"]] <- CT
	} else {
		normalNP <- envir[["normalNP"]]
		normalNP <- normalNP[, plate==uplate[p]]
		mus <- rowMedians(NP * normalNP, na.rm=TRUE)
		crosshyb <- median(muA) - median(mus)
		X <- cbind(1, log2(mus+crosshyb))
		##X <- cbind(1, log2(mus))
		logPhiT <- X %*% betahat
		phiT[, p] <- 2^(logPhiT)
		nuT[, p] <- mus - 2*phiT[, p]
		T <- 1/phiT[, p]*(NP - nuT[, p])
		CT <- envir[["CT"]]
		CT[, plate==uplate[p]] <- matrix(as.integer(100*T), nrow(T), ncol(T))

		##Variance for prediction region
		sig2T[, plate==uplate[[p]]] <- rowMAD(log2(NP*normalNP), na.rm=TRUE)^2	
		envir[["sig2T"]] <- sig2T
		envir[["CT"]] <- CT
		envir[["phiT"]] <- phiT
		envir[["nuT"]] <- nuT
	}
	##---------------------------------------------------------------------------
	## For NA SNPs, treat as nonpolymorphic
	##---------------------------------------------------------------------------
}

withinGenotypeMoments <- function(p, A, B, calls, conf, CONF.THR, DF.PRIOR, envir){
	CHR <- envir[["chrom"]]
	Ns <- envir[["Ns"]]
	muA <- envir[["muA"]]
	muB <- envir[["muB"]]
	vA <- envir[["vA"]]
	vB <- envir[["vB"]]
	plate <- envir[["plate"]]
	uplate <- envir[["uplate"]]
	normal <- envir[["normal"]][, plate==uplate[p]]
	G <- calls; rm(calls); gc()	
	if(is.null(normal)) normal <- matrix(TRUE, nrow(G), ncol(G))

	highConf <- 1-exp(-conf/1000)
	highConf <- highConf > CONF.THR
	if(CHR == 23){
		gender <- envir[["gender"]]
		IX <- matrix(gender, nrow(G), ncol(G), byrow=TRUE)
		IX <- IX == "female"
	} else IX <- matrix(TRUE, nrow(G), ncol(G))
	
	index <- GT.B <- GT.A <- vector("list", 3)
	names(index) <- names(GT.B) <- names(GT.A) <- c("AA", "AB", "BB")
	##--------------------------------------------------
	##within-genotype sufficient statistics
	##--------------------------------------------------
	GT.B <- GT.A <- list()
	for(j in 1:3){
		GT <- G==j & highConf & IX & normal
		Ns[, p, j+2] <- rowSums(GT, na.rm=TRUE)		
		GT[GT == FALSE] <- NA
		GT.A[[j]] <- GT*A
		GT.B[[j]] <- GT*B
		index[[j]] <- which(Ns[, p, j+2] > 0)
		
		muA[, p, j+2] <- rowMedians(GT.A[[j]], na.rm=TRUE)
		muB[, p, j+2] <- rowMedians(GT.B[[j]], na.rm=TRUE)
		vA[, p, j+2] <- rowMAD(GT.A[[j]], na.rm=TRUE)
		vB[, p, j+2] <- rowMAD(GT.B[[j]], na.rm=TRUE)

		##Shrink towards the typical variance
		DF <- Ns[, p, j+2]-1
		DF[DF < 1] <- 1
		v0A <- median(vA[, p, j+2], na.rm=TRUE)
		v0B <- median(vB[, p, j+2], na.rm=TRUE)
		vA[, p, j+2] <- (vA[, p, j+2]*DF + v0A*DF.PRIOR)/(DF.PRIOR+DF)
		vA[is.na(vA[, p, j+2]), p, j+2] <- v0A
		vB[, p, j+2] <- (vB[, p, j+2]*DF + v0B*DF.PRIOR)/(DF.PRIOR+DF)
		vB[is.na(vB[, p, j+2]), p, j+2] <- v0B
	}
	if(CHR == 23){
		k <- 1
		for(j in c(1,3)){
			GT <- G==j & highConf & !IX 
			Ns[, p, k] <- rowSums(GT)
			GT[GT == FALSE] <- NA
			muA[, p, k] <- rowMedians(GT*A, na.rm=TRUE)
			muB[, p, k] <- rowMedians(GT*B, na.rm=TRUE)
			vA[, p, k] <- rowMAD(GT*A, na.rm=TRUE)
			vB[, p, k] <- rowMAD(GT*B, na.rm=TRUE)
			
			DF <- Ns[, p, k]-1
			DF[DF < 1] <- 1
			v0A <- median(vA[, p, k], na.rm=TRUE)
			v0B <- median(vB[, p, k], na.rm=TRUE)
			vA[, p, k] <- (vA[, p, k]*DF + v0A*DF.PRIOR)/(DF.PRIOR+DF)
			vA[is.na(vA[, p, k]), p, k] <- v0A
			vB[, p, k] <- (vB[, p, k]*DF + v0B*DF.PRIOR)/(DF.PRIOR+DF)
			vB[is.na(vB[, p, k]), p, k] <- v0B			
			k <- k+1
		}
	}
	envir[["GT.A"]] <- GT.A
	envir[["GT.B"]] <- GT.B
	envir[["Ns"]] <- Ns
	envir[["index"]] <- index
	envir[["muA"]] <- muA
	envir[["muB"]] <- muB
	envir[["vA"]] <- vA
	envir[["vB"]] <- vB
}
	

oneBatch <- function(plateIndex,
		     A,
		     B,
		     calls,
		     conf,
		     gender,
		     NP,
		     plate,
		     MIN.OBS=1,
		     envir,
		     DF.PRIOR=50,
		     CONF.THR=0.99,
		     trim=0,
		     bias.adj=FALSE, priorProb, ...){
	p <- plateIndex
	CHR <- envir[["chrom"]]
	if(bias.adj){
		nuA <- envir[["nuA"]]
		if(all(is.na(nuA))) {
			message("Background and signal coefficients have not yet been estimated -- can not do bias correction yet")
			stop("Must run computeCopynumber a second time with bias.adj=TRUE to do the adjustment")
		}
		message("running bias adjustment")		
		##adjustment for nonpolymorphic probes
		biasAdjNP(plateIndex=p, envir=envir, priorProb=priorProb)
		##adjustment for SNPs
		biasAdj(plateIndex=p, envir=envir, priorProb=priorProb)
		message("Recomputing location and scale parameters")		
	}
	withinGenotypeMoments(p=p, A=A, B=B, calls=calls, conf=conf, CONF.THR=CONF.THR,
			      DF.PRIOR=DF.PRIOR, envir=envir)
	GT.A <- envir[["GT.A"]]
	GT.B <- envir[["GT.B"]]
	index <- envir[["index"]]
	locationAndScale(p=p, GT.A=GT.A, GT.B=GT.B, index=index, envir=envir, DF.PRIOR=DF.PRIOR)
	muA <- envir[["muA"]]
	muB <- envir[["muB"]]
	Ns <- envir[["Ns"]]

	##---------------------------------------------------------------------------
	## Predict sufficient statistics for unobserved genotypes (plate-specific)
	##---------------------------------------------------------------------------
	index.AA <- which(Ns[, p, "AA"] >= 3)
	index.AB <- which(Ns[, p, "AB"] >= 3)
	index.BB <- which(Ns[, p, "BB"] >= 3)
	correct.orderA <- muA[, p, "AA"] > muA[, p, "BB"]
	correct.orderB <- muB[, p, "BB"] > muB[, p, "AA"]
	##For chr X, this will ignore the males 
	nobs <- rowSums(Ns[, p, 3:5] >= MIN.OBS, na.rm=TRUE) == 3
	index.complete <- which(correct.orderA & correct.orderB & nobs) ##be selective here
	size <- min(5000, length(index.complete))
	if(size == 5000) index.complete <- sample(index.complete, 5000)
	if(length(index.complete) < 200){
		warning("fewer than 200 snps pass criteria for predicting the sufficient statistics")
		stop()
	}
	index[[1]] <- which(Ns[, p, "AA"] == 0 & (Ns[, p, "AB"] >= MIN.OBS & Ns[, p, "BB"] >= MIN.OBS))
	index[[2]] <- which(Ns[, p, "AB"] == 0 & (Ns[, p, "AA"] >= MIN.OBS & Ns[, p, "BB"] >= MIN.OBS))
	index[[3]] <- which(Ns[, p, "BB"] == 0 & (Ns[, p, "AB"] >= MIN.OBS & Ns[, p, "AA"] >= MIN.OBS))
	mnA <- muA[, p, 3:5]
	mnB <- muB[, p, 3:5]
	for(j in 1:3){
		if(length(index[[j]]) == 0) next()
		X <- cbind(1, mnA[index.complete,  -j], mnB[index.complete,  -j])
		Y <- cbind(mnA[index.complete, j], mnB[index.complete,  j])
		betahat <- solve(crossprod(X), crossprod(X,Y))
		X <- cbind(1, mnA[index[[j]],  -j],  mnB[index[[j]],  -j])
		mus <- X %*% betahat
		muA[index[[j]], p, j+2] <- mus[, 1]
		muB[index[[j]], p, j+2] <- mus[, 2]
	}
	nobsA <- Ns[, p, "A"] > 10
	nobsB <- Ns[, p, "B"] > 10
	complete <- list()
	complete[[1]] <- which(correct.orderA & correct.orderB & nobsA) ##be selective here
	complete[[2]] <- which(correct.orderA & correct.orderB & nobsB) ##be selective here	
	size <- min(5000, length(complete[[1]]))
	if(size == 5000) complete <- lapply(complete, function(x) sample(x, size))
	if(CHR == 23){
		index[[1]] <- which(Ns[, p, "A"] == 0)
		index[[2]] <- which(Ns[, p, "B"] == 0)
		cols <- 2:1
		for(j in 1:2){
			if(length(index[[j]]) == 0) next()
			X <- cbind(1, muA[complete[[j]], p, cols[j]], muB[complete[[j]], p, cols[j]])
			Y <- cbind(muA[complete[[j]], p, j], muB[complete[[j]], p, j])
			betahat <- solve(crossprod(X), crossprod(X,Y))
			X <- cbind(1, muA[index[[j]], p, cols[j]],  muB[index[[j]], p, cols[j]])
			mus <- X %*% betahat
			muA[index[[j]], p, j] <- mus[, 1]
			muB[index[[j]], p, j] <- mus[, 2]
		}
	}
	##missing two genotypes
	noAA <- Ns[, p, "AA"] < MIN.OBS
	noAB <- Ns[, p, "AB"] < MIN.OBS
	noBB <- Ns[, p, "BB"] < MIN.OBS
	index[[1]] <- noAA & noAB
	index[[2]] <- noBB & noAB
	index[[3]] <- noAA & noBB
	snpflags <- envir[["snpflags"]]
	snpflags[, p] <- index[[1]] | index[[2]] | index[[3]]

	##---------------------------------------------------------------------------
	## Two genotype clusters not observed -- would sequence help? (didn't seem that helpful)
	## 1 extract index of complete data
	## 2 Regress  mu1,mu3 ~ sequence + mu2
	## 3 Predict mu1*, mu3* for missing genotypes
	##---------------------------------------------------------------------------
	cols <- c(3, 1, 2)
	for(j in 1:3){
		if(sum(index[[1]]) == 0) next()
		k <- cols[j]
		X <- cbind(1, mnA[index.complete, k], mnB[index.complete, k])
		Y <- cbind(mnA[index.complete,  -k],
			   mnB[index.complete,  -k])
		betahat <- solve(crossprod(X), crossprod(X,Y))
		X <- cbind(1, mnA[index[[j]],  k], mnB[index[[j]],  k])
		mus <- X %*% betahat
		muA[index[[j]], p, -c(1, 2, k+2)] <- mus[, 1:2]
		muB[index[[j]], p, -c(1, 2, k+2)] <- mus[, 3:4]
	}
	
	negA <- rowSums(muA[, p, ] < 0) > 0
	negB <- rowSums(muB[, p, ] < 0) > 0	
	snpflags[, p] <- snpflags[, p] | negA | negB | rowSums(is.na(muA[, p, 3:5]), na.rm=TRUE) > 0
	envir[["snpflags"]] <- snpflags
	dn.Ns <- dimnames(Ns)
	Ns <- array(as.integer(Ns), dim=dim(Ns))
	dimnames(Ns)[[3]] <- dn.Ns[[3]]
	envir[["Ns"]] <- Ns
	envir[["muA"]] <- muA
	envir[["muB"]] <- muB
	plates.completed <- envir[["plates.completed"]]
	plates.completed[p] <- TRUE
	envir[["plates.completed"]] <- plates.completed
}

locationAndScale <- function(p, GT.A, GT.B, index, envir, DF.PRIOR){
	tau2A <- envir[["tau2A"]]
	tau2B <- envir[["tau2B"]]
	sig2A <- envir[["sig2A"]]
	sig2B <- envir[["sig2B"]]

	corr <- envir[["corr"]]
	corrA.BB <- envir[["corrA.BB"]]
	corrB.AA <- envir[["corrB.AA"]]
	Ns <- get("Ns", envir)
	
	index.AA <- index[[1]]
	index.AB <- index[[2]]
	index.BB <- index[[3]]
	rm(index); gc()

	AA.A <- GT.A[[1]]
	AB.A <- GT.A[[2]]
	BB.A <- GT.A[[3]]
	
	AA.B <- GT.B[[1]]
	AB.B <- GT.B[[2]]
	BB.B <- GT.B[[3]]	
	x <- BB.A[index.BB, ]
	tau2A[index.BB, p] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2
	DF <- Ns[, p, "BB"]-1
	DF[DF < 1] <- 1
	med <- median(tau2A[, p], na.rm=TRUE)
	tau2A[, p] <- (tau2A[, p] * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
	tau2A[is.na(tau2A[, p]), p] <- med
	
	x <- BB.B[index.BB, ]
	sig2B[index.BB, p] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2	
	med <- median(sig2B[, p], na.rm=TRUE)
	sig2B[, p] <- (sig2B[, p] * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
	sig2B[is.na(sig2B[, p]), p] <- med
	
	x <- AA.B[index.AA, ]
	tau2B[index.AA, p] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2		
	DF <- Ns[, p, "AA"]-1
	DF[DF < 1] <- 1
	med <- median(tau2B[, p], na.rm=TRUE)
	tau2B[, p] <- (tau2B[, p] * DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
	tau2B[is.na(tau2B[, p]), p] <- med
	
	x <- AA.A[index.AA, ]
	sig2A[index.AA, p] <- rowMAD(log2(x), log2(x), na.rm=TRUE)^2##var(log(IA)|AA)	
	med <- median(sig2A[, p], na.rm=TRUE)
	sig2A[, p] <- (sig2A[, p]*DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
	sig2A[is.na(sig2A[, p]), p] <- med	

	if(length(index.AB) > 0){ ##all homozygous is possible
		x <- AB.A[index.AB, ]
		y <- AB.B[index.AB, ]
		corr[index.AB, p] <- rowCors(x, y, na.rm=TRUE)
		corr[corr < 0] <- 0
		DF <- Ns[, p, "AB"]-1
		DF[DF<1] <- 1
		med <- median(corr[, p], na.rm=TRUE)
		corr[, p] <- (corr[, p]*DF  +  med * DF.PRIOR)/(DF.PRIOR + DF)
		corr[is.na(corr[, p]), p] <- med
	}
	backgroundB <- AA.B[index.AA, ]
	signalA <- AA.A[index.AA, ]
	corrB.AA[index.AA, p] <- rowCors(backgroundB, signalA, na.rm=TRUE)
	DF <- Ns[, p, "AA"]-1
	DF[DF < 1] <- 1
	med <- median(corrB.AA[, p], na.rm=TRUE)
	corrB.AA[, p] <- (corrB.AA[, p]*DF + med*DF.PRIOR)/(DF.PRIOR + DF)
	corrB.AA[is.na(corrB.AA[, p]), p] <- med

	backgroundA <- BB.A[index.BB, ]
	signalB <- BB.B[index.BB, ]
	corrA.BB[index.BB, p] <- rowCors(backgroundA, signalB, na.rm=TRUE)
	DF <- Ns[, p, "BB"]-1
	DF[DF < 1] <- 1
	med <- median(corrA.BB[, p], na.rm=TRUE)
	corrA.BB[, p] <- (corrA.BB[, p]*DF + med*DF.PRIOR)/(DF.PRIOR + DF)
	corrA.BB[is.na(corrA.BB[, p]), p] <- med

	envir[["tau2A"]] <- tau2A
	envir[["tau2B"]] <- tau2B
	envir[["sig2A"]] <- sig2A
	envir[["sig2B"]] <- sig2B
	envir[["corr"]] <- corr
	envir[["corrB.AA"]] <- corrB.AA
	envir[["corrA.BB"]] <- corrA.BB	
}

coefs <- function(plateIndex, conf, MIN.OBS=3, envir, CONF.THR=0.99){
	p <- plateIndex
	plates.completed <- envir[["plates.completed"]]
	if(!plates.completed[p]) return()
	CHR <- envir[["chrom"]]
	plate <- envir[["plate"]]
	muA <- envir[["muA"]]
	muB <- envir[["muB"]]
	vA <- envir[["vA"]]
	vB <- envir[["vB"]]
	Ns <- envir[["Ns"]]
	uplate <- envir[["uplate"]]
	if(CHR != 23){
		IA <- muA[, p, 3:5]
		IB <- muB[, p, 3:5]
		vA <- vA[, p, 3:5]
		vB <- vB[, p, 3:5]
		Np <- Ns[, p, 3:5]
	} else {
		IA <- muA[, p, ]
		IB <- muB[, p, ]
		vA <- vA[, p, ]
		vB <- vB[, p, ]
		Np <- Ns[, p, ]
	}
	NOHET <- mean(Ns[, p, "AB"], na.rm=TRUE) < 0.05
	##---------------------------------------------------------------------------
	## Estimate nu and phi
	##---------------------------------------------------------------------------
	if(NOHET){
		##only homozygous
		Np <- Np[, -2]
		Np[Np < 1] <- 1
		IA <- IA[, c(1, 3)]
		IB <- IB[, c(1, 3)]		
		vA <- vA[, c(3,5)]
		vB <- vB[, c(3,5)]
	}else 	Np[Np < 1] <- 1
	vA2 <- vA^2/Np
	vB2 <- vB^2/Np
	wA <- sqrt(1/vA2)
	wB <- sqrt(1/vB2)
	YA <- IA*wA
	YB <- IB*wB
	nuphiAllele(p=p, allele="A", Ystar=YA, W=wA, envir=envir)
	nuphiAllele(p=p, allele="B", Ystar=YB, W=wB, envir=envir)

	##---------------------------------------------------------------------------
	##Estimate crosshyb using X chromosome and sequence information
	##---------------------------------------------------------------------------
	##browser()
	####data(sequences, package="genomewidesnp6Crlmm")
	##snpflags <- envir[["snpflags"]]
	##muA <- envir[["muA"]][, p, 3:5]
	##muB <- envir[["muB"]][, p, 3:5]
	##Y <- envir[["phiAx"]]
	##load("sequences.rda")
	##seqA <- sequences[, "A", ][, 1]
	##seqA <- seqA[match(snps, names(seqA))]
	##X <- cbind(1, sequenceDesignMatrix(seqA))
	##X <- cbind(X, nuA[, p], phiA[, p], nuB[, p], phiB[, p])
	##missing <- rowSums(is.na(X)) > 0
	##betahat <- solve(crossprod(X[!missing, ]), crossprod(X[!missing, ], Y[!missing]))
}

polymorphic <- function(plateIndex, A, B, envir){
	p <- plateIndex
	plates.completed <- envir[["plates.completed"]]
	if(!plates.completed[p]) return()
	CHR <- envir[["chrom"]]
	plate <- envir[["plate"]]
	uplate <- envir[["uplate"]]
	vA <- envir[["vA"]]
	vB <- envir[["vB"]]
	nuA <- envir[["nuA"]][, p]
	nuB <- envir[["nuB"]][, p]
	nuA.se <- envir[["nuA.se"]]
	nuB.se <- envir[["nuB.se"]]
	phiA <- envir[["phiA"]][, p]
	phiB <- envir[["phiB"]][, p]
	phiA.se <- envir[["phiA.se"]]
	phiB.se <- envir[["phiB.se"]]
	Ns <- get("Ns", envir)
	CA <- get("CA", envir)
	CB <- get("CB", envir)
	NOHET <- mean(Ns[, p, "AB"], na.rm=TRUE) < 0.05
	##---------------------------------------------------------------------------
	## Estimate CA, CB
	##---------------------------------------------------------------------------
	if(CHR == 23){
		phiAx <- as.matrix(envir[["phiAx"]])
		phiBx <- as.matrix(envir[["phiBx"]])
		phiAx <- phiAx[, p]
		phiBx <- phiBx[, p]
		phistar <- phiBx/phiA
		tmp <- (B-nuB - phistar*A + phistar*nuA)/phiB
		copyB <- tmp/(1-phistar*phiAx/phiB)
		copyA <- (A-nuA-phiAx*copyB)/phiA
		CB[, plate==uplate[p]] <- matrix(as.integer(100*copyB), nrow(copyB), ncol(copyB))
		CA[, plate==uplate[p]] <- matrix(as.integer(100*copyA), nrow(copyA), ncol(copyA))
	} else{
		CA[, plate==uplate[p]] <- matrix(as.integer(100*1/phiA*(A-nuA)), nrow(A), ncol(A))
		CB[, plate==uplate[p]] <- matrix(as.integer(100*1/phiB*(B-nuB)), nrow(A), ncol(A))
	}
	assign("CA", CA, envir)
	assign("CB", CB, envir)
}




biasAdj <- function(plateIndex, envir, priorProb){
	CHR <- envir[["chrom"]]
	if(CHR == 23){
		phiAx <- envir[["phiAx"]]
		phiBx <- envir[["phiBx"]]
	}
	A <- envir[["A"]]
	B <- envir[["B"]]
	sig2A <- envir[["sig2A"]]
	sig2B <- envir[["sig2B"]]
	tau2A <- envir[["tau2A"]]
	tau2B <- envir[["tau2B"]]
	corrA.BB <- envir[["corrA.BB"]]
	corrB.AA <- envir[["corrB.AA"]]
	corr <- envir[["corr"]]
	nuA <- envir[["nuA"]]
	nuB <- envir[["nuB"]]
	phiA <- envir[["phiA"]]
	phiB <- envir[["phiB"]]
	normal <- envir[["normal"]]
	p <- plateIndex
	plate <- envir[["plate"]]
	if(missing(priorProb)) priorProb <- rep(1/4, 4) ##uniform
	emit <- array(NA, dim=c(nrow(A), ncol(A), 10))##SNPs x sample x 'truth'	

	lA <- log2(A)
	lB <- log2(B)	
	X <- cbind(lA, lB)
	counter <- 1##state counter								
	for(CT in 0:3){
		for(CA in 0:CT){
			cat(".")
			CB <- CT-CA
			A.scale <- sqrt(tau2A[, p]*(CA==0) + sig2A[, p]*(CA > 0))
			B.scale <- sqrt(tau2B[, p]*(CA==0) + sig2B[, p]*(CA > 0))
			if(CA == 0 & CB == 0) rho <- 0
			if(CA == 0 & CB > 0) rho <- corrA.BB[, p]
			if(CA > 0 & CB == 0) rho <- corrB.AA[, p]
			if(CA > 0 & CB > 0) rho <- corr[, p]
			if(CHR == 23){
				means <- cbind(suppressWarnings(log2(nuA[, p]+CA*phiA[, p] + CB*phiAx[, p])), suppressWarnings(log2(nuB[, p]+CB*phiB[, p] + CA*phiBx[, p])))
			} else{
				means <- cbind(suppressWarnings(log2(nuA[, p]+CA*phiA[, p])), suppressWarnings(log2(nuB[, p]+CB*phiB[, p])))
			}
			covs <- rho*A.scale*B.scale
			A.scale2 <- A.scale^2
			B.scale2 <- B.scale^2
			m <- 1			
			for(i in 1:nrow(A)){
				Sigma <- matrix(c(A.scale2[i], covs[i], covs[i], B.scale2[i]), 2,2)
				xx <- matrix(X[i, ], ncol=2)
				tmp <- dmvnorm(xx, mean=means[i, ], sigma=Sigma) 
				emit[m, , counter] <- tmp
				m <- m+1				
			}
			counter <- counter+1
		}
	}
	homDel <- priorProb[1]*emit[, , 1]
	hemDel <- priorProb[2]*emit[, , c(2, 3)] # + priorProb[3]*emit[, c(4, 5, 6)] + priorProb[4]*emit[, c(7:10)]
	norm <- priorProb[3]*emit[, , 4:6]
	amp <- priorProb[4]*emit[, , 7:10]
	##sum over the different combinations within each copy number state
	hemDel <- apply(hemDel, c(1,2), sum)
	norm <- apply(norm, c(1, 2), sum)
	amp <- apply(amp, c(1,2), sum)

	tmp <- array(NA, dim=c(nrow(A), ncol(A), 4))
	tmp[, , 1] <- homDel
	tmp[, , 2] <- hemDel
	tmp[, , 3] <- norm
	tmp[, , 4] <- amp
	tmp2 <- apply(tmp, c(1, 2), function(x) order(x, decreasing=TRUE)[1])
	##Adjust for SNPs that have less than 80% of the samples in an altered state
	##flag the remainder?
	if(CHR != 23){
		tmp3 <- tmp2 != 3
	}else{
		##should also consider pseud
		gender <- envir[["gender"]]
		tmp3 <- tmp2
		tmp3[, gender=="male"] <- tmp2[, gender=="male"] != 2
		tmp3[, gender=="female"] <- tmp2[, gender=="female"] != 3
	}
	##Those near 1 have NaNs for nu and phi.  this occurs by NaNs in the muA[,, "A"] or muA[, , "B"] for X chromosome
	propAlt <- rowMeans(tmp3)##prop normal
	ii <- propAlt < 0.75
	##only exclude observations from one tail, depending on
	##whether more are up or down
	if(CHR != 23){
		moreup <- rowSums(tmp2 > 3) > rowSums(tmp2 < 3)
		notUp <-  tmp2[ii & moreup, ] <= 3
		notDown <- tmp2[ii & !moreup, ] >= 3
		NORM <- matrix(NA, nrow(A), ncol(A))
		NORM[ii & moreup, ] <- notUp
		NORM[ii & !moreup, ] <- notDown
		##Define NORM so that we can iterate this step
		##NA's in the previous iteration (normal) will be propogated		
		normal <- NORM*normal
	} else{
		fem <- tmp2[, gender=="female"]
		mal <- tmp2[, gender=="male"]
		moreupF <- rowSums(fem > 3) > rowSums(fem < 3)
		moreupM <- rowSums(mal > 2) > rowSums(mal < 2)
		notUpF <-  fem[ii & moreupF, ] <= 3
		notUpM <-  fem[ii & moreupM, ] <= 2	       
		notDownF <- fem[ii & !moreupF, ] >= 3
		notDownM <- mal[ii & !moreupM, ] >= 2
		normalF <- matrix(TRUE, nrow(fem), ncol(fem))
		normalF[ii & moreupF, ] <- notUpF
		normalF[ii & !moreupF, ] <- notDownF
		normalM <- matrix(TRUE, nrow(mal), ncol(mal))
		normalM[ii & moreupM, ] <- notUpM
		normalM[ii & !moreupM, ] <- notDownM
		normal <- matrix(TRUE, nrow(A), ncol(A))
		normal[, gender=="female"] <- normalF
		normal[, gender=="male"] <- normalM
	}
	flagAltered <- which(propAlt > 0.5)
	envir[["flagAltered"]] <- flagAltered
	normal[normal == FALSE] <- NA
	envir[["normal"]] <- normal
}


posteriorNonpolymorphic <- function(plateIndex, envir, priorProb, cnStates=0:6){
	p <- plateIndex
	CHR <- envir[["chrom"]]
	if(missing(priorProb)) priorProb <- rep(1/length(cnStates), length(cnStates)) ##uniform	
	plate <- envir[["plate"]]
	uplate <- envir[["plate"]]
	NP <- envir[["NP"]][, plate==uplate[p]]
	nuT <- envir[["nuT"]][, p]
	phiT <- envir[["phiT"]][, p]
	sig2T <- envir[["sig2T"]][, p]
	##Assuming background variance for np probes is the same on the log-scale
	emit <- array(NA, dim=c(nrow(NP), ncol(NP), length(cnStates)))##SNPs x sample x 'truth'
	lT <- log2(NP)
	sds <- sqrt(sig2T)
	counter <- 1##state counter	
	for(CT in cnStates){
		cat(".")
		if(CHR == 23) browser()
		means <- suppressWarnings(log2(nuT + CT*phiT))
		emit[, , counter] <- dnorm(lT, mean=means, sd=sds)
		counter <- counter+1
	}
	for(j in seq(along=cnStates)){
		emit[, , j] <- priorProb[j]*emit[, , j]
	}
	homDel <- emit[, , 1]
	hemDel <- emit[, , 2]
	norm <- emit[, , 3]
	amp <- emit[, , 4]
	amp4 <- emit[, , 5]
	amp5 <- emit[, , 6]
	amp6 <- emit[, , 7]
	total <- homDel+hemDel+norm+amp+amp4+amp5+amp6
	weights <- array(NA, dim=c(nrow(NP), ncol(NP), length(cnStates)))
	weights[, , 1] <- homDel/total
	weights[, , 2] <- hemDel/total
	weights[, , 3] <- norm/total
	weights[, , 4] <- amp/total
	weights[, , 5] <- amp4/total
	weights[, , 6] <- amp5/total
	weights[, , 7] <- amp6/total
	##posterior mode
	posteriorMode <- apply(weights, c(1, 2), function(x) order(x, decreasing=TRUE)[1])
	posteriorMode <- posteriorMode-1
	##sns <- envir[["sns"]]
	##colnames(posteriorMode) <- sns
	##envir[["np.posteriorMode"]] <- posteriorMode
	##envir[["np.weights"]] <- weights
	posteriorMeans <- 0*homDel/total + 1*hemDel/total + 2*norm/total + 3*amp/total + 4*amp4/total + 5*amp5/total + 6*amp6/total
	##colnames(posteriorMeans) <- sns
	##envir[["np.posteriorMeans"]] <- posteriorMeans
	return(posteriorMode)
}

posteriorWrapper <- function(envir){
	snp.PM <- matrix(NA, length(envir[["snps"]]), length(envir[["sns"]]))
	np.PM <- matrix(NA, length(envir[["cnvs"]]), length(envir[["sns"]]))
	plate <- envir[["plate"]]
	uplate <- envir[["uplate"]]
	for(p in seq(along=uplate)){
		tmp <- expectedC(plateIndex=p, envir=envir)
		snp.PM[, plate==uplate[p]] <- tmp
		##snp.pm <- env[["posteriorMode"]]
		##trace(posteriorNonpolymorphic, browser)
		tmp <- posteriorNonpolymorphic(plateIndex=p, envir=envir)
		np.PM[, plate==uplate[p]] <- tmp##env[["np.posteriorMode"]]
		##pMode <- rbind(snp.pm, np.pm)
		##rownames(pMode) <- c(env[["snps"]], env[["cnvs"]])
		##dn <- dimnames(pMode)
		##pMode <- matrix(as.integer(pMode), nrow(pMode), ncol(pMode))
	}
	PM <- rbind(snp.PM, np.PM)
	PM <- matrix(as.integer(PM), nrow(PM), ncol(PM))
	dns <- list(c(envir[["snps"]], envir[["cnvs"]]), envir[["sns"]])
	dimnames(PM) <- dns
	return(PM)
}





##for polymorphic probes
expectedC <- function(plateIndex, envir, priorProb, cnStates=0:6){
	p <- plateIndex
	CHR <- envir[["chrom"]]
	if(missing(priorProb)) priorProb <- rep(1/length(cnStates), length(cnStates)) ##uniform	
	plate <- envir[["plate"]]
	uplate <- envir[["uplate"]]
	A <- envir[["A"]]
	B <- envir[["B"]]
	A <- A[, plate==uplate[p]]
	B <- B[, plate==uplate[p]]
	calls <- envir[["calls"]]	
	calls <- calls[, plate==unique(plate)[p]]
	probA <- sqrt(rowMeans(calls == 1, na.rm=TRUE))
	probB <- sqrt(rowMeans(calls == 3, na.rm=TRUE))
	sig2A <- envir[["sig2A"]]
	sig2B <- envir[["sig2B"]]
	tau2A <- envir[["tau2A"]]
	tau2B <- envir[["tau2B"]]
	corrA.BB <- envir[["corrA.BB"]]
	corrB.AA <- envir[["corrB.AA"]]
	corr <- envir[["corr"]]
	nuA <- envir[["nuA"]]
	nuB <- envir[["nuB"]]
	phiA <- envir[["phiA"]]
	phiB <- envir[["phiB"]]
	emit <- array(NA, dim=c(nrow(A), ncol(A), 28))##SNPs x sample x 'truth'
	##AAAA, AAAB, AABB, ABBB, BBBB
	##AAAAA, AAAAB, AAABB, AABBB, ABBBB, BBBBB
	##AAAAAA, AAAAAB, AAAABB, AAABBB, AABBBB, ABBBBB, BBBBBB
	lA <- log2(A)
	lB <- log2(B)	
	X <- cbind(lA, lB)	
	counter <- 1##state counter
	for(CT in cnStates){
		cat(".")
		for(CA in 0:CT){
			CB <- CT-CA
			A.scale <- sqrt(tau2A[, p]*(CA==0) + sig2A[, p]*(CA > 0))
			B.scale <- sqrt(tau2B[, p]*(CB==0) + sig2B[, p]*(CB > 0))
			scale <- c(A.scale, B.scale)
			if(CA == 0 & CB == 0) rho <- 0
			if(CA == 0 & CB > 0) rho <- corrA.BB[, p]
			if(CA > 0 & CB == 0) rho <- corrB.AA[, p]
			if(CA > 0 & CB > 0) rho <- corr[, p]
			if(CHR == 23) browser()
			means <- cbind(suppressWarnings(log2(nuA[, p]+CA*phiA[, p])), suppressWarnings(log2(nuB[, p]+CB*phiB[, p])))
			covs <- rho*A.scale*B.scale
			A.scale2 <- A.scale^2
			B.scale2 <- B.scale^2			
			##ensure positive definite			
			##Sigma <- as.matrix(nearPD(matrix(c(A.scale^2, covs,
			##covs, B.scale^2), 2, 2))[[1]])
			m <- 1##snp counter				
			for(i in 1:nrow(A)){
				Sigma <- matrix(c(A.scale2[i], covs[i], covs[i], B.scale2[i]), 2,2)
				xx <- matrix(X[i, ], ncol=2)
				tmp <- dmvnorm(xx, mean=means[i, ], sigma=Sigma) 				
				##Using HWE: P(CA=ca, CB=cb|CT=c)				
				ptmp <- (probA[i]^CA)*(probB[i]^CB)*tmp
				emit[m, , counter] <- ptmp
				m <- m+1				
			}
			counter <- counter+1			
		}
	}
	##priorProb=P(CT=c)
	homDel <- priorProb[1]*emit[, , 1]
	hemDel <- priorProb[2]*emit[, , c(2, 3)] # + priorProb[3]*emit[, c(4, 5, 6)] + priorProb[4]*emit[, c(7:10)]
	norm <- priorProb[3]*emit[, , 4:6]
	amp <- priorProb[4]*emit[, , 7:10]
	amp4 <- priorProb[5]*emit[, , 11:15]
	amp5 <- priorProb[6]*emit[, , 16:21]
	amp6 <- priorProb[7]*emit[, , 22:28]	
	##sum over the different combinations within each copy number state
	hemDel <- apply(hemDel, c(1,2), sum)
	norm <- apply(norm, c(1, 2), sum)
	amp <- apply(amp, c(1,2), sum)
	amp4 <- apply(amp4, c(1,2), sum)
	amp5 <- apply(amp5, c(1,2), sum)
	amp6 <- apply(amp6, c(1,2), sum)
	total <- homDel+hemDel+norm+amp+amp4+amp5+amp6
	weights <- array(NA, dim=c(nrow(homDel), ncol(A), 7))
	weights[, , 1] <- homDel/total
	weights[, , 2] <- hemDel/total
	weights[, , 3] <- norm/total
	weights[, , 4] <- amp/total
	weights[, , 5] <- amp4/total
	weights[, , 6] <- amp5/total
	weights[, , 7] <- amp6/total
	##posterior mode
	posteriorMode <- apply(weights, c(1, 2), function(x) order(x, decreasing=TRUE)[1])
	posteriorMode <- posteriorMode-1
	##This is for one plate.  Need to instantiate a much bigger
	##object in the environment
	
	##envir[["posteriorMode"]] <- posteriorMode
	##weights <- list(homDel/total, hemDel/total, norm/total, amp/total, amp4/total, amp5/total, amp6/total)
	##names(weights) <- c(cnStates)
	##envir[["weights"]] <- weights
	posteriorMeans <- 0*homDel/total + 1*hemDel/total + 2*norm/total + 3*amp/total + 4*amp4/total + 5*amp5/total + 6*amp6/total
	##sns <- envir[["sns"]]
	##colnames(posteriorMeans) <- sns
	##envir[["posteriorMeans"]] <- posteriorMeans
	return(posteriorMode)
}





biasAdjNP <- function(plateIndex, envir, priorProb){
	p <- plateIndex
	normalNP <- envir[["normalNP"]]
	CHR <- envir[["chrom"]]
	NP <- envir[["NP"]]
	plate <- envir[["plate"]]
	uplate <- envir[["uplate"]]
	sig2T <- envir[["sig2T"]]
	NP <- NP[, plate==uplate[p]]
	sig2T <- sig2T[, plate==uplate[p]]

	##Assume that on the log-scale, that the background variance is the same...
	tau2T <- sig2T	
	nuT <- envir[["nuT"]]
	phiT <- envir[["phiT"]]
	p <- plateIndex
	plate <- envir[["plate"]]
	if(missing(priorProb)) priorProb <- rep(1/4, 4) ##uniform
	emit <- array(NA, dim=c(nrow(NP), ncol(NP), 4))##SNPs x sample x 'truth'	
	lT <- log2(NP)
	counter <- 1 ##state counter
	for(CT in 0:3){
		sds <- sqrt(tau2T*(CT==0) + sig2T*(CT > 0))
		means <- suppressWarnings(log2(nuT[, p]+CT*phiT[, p]))
		tmp <- dnorm(lT, mean=means, sd=sds)
		emit[, , counter] <- tmp
		counter <- counter+1
	}
	tmp2 <- apply(emit, c(1, 2), function(x) order(x, decreasing=TRUE)[1])
	##Adjust for SNPs that have less than 80% of the samples in an altered state
	##flag the remainder?
	if(CHR != 23){
		tmp3 <- tmp2 != 3
	}else{
		browser()
		##should also consider pseudoautosomal
		gender <- envir[["gender"]]
		tmp3 <- tmp2
		tmp3[, gender=="male"] <- tmp2[, gender=="male"] != 2
		tmp3[, gender=="female"] <- tmp2[, gender=="female"] != 3
	}
	##Those near 1 have NaNs for nu and phi.  this occurs by NaNs in the muA[,, "A"] or muA[, , "B"] for X chromosome
	propAlt <- rowMeans(tmp3)##prop normal
	ii <- propAlt < 0.75
	##only exclude observations from one tail, depending on
	##whether more are up or down
	if(CHR != 23){
		moreup <- rowSums(tmp2 > 3) > rowSums(tmp2 < 3)
		notUp <-  tmp2[ii & moreup, ] <= 3
		notDown <- tmp2[ii & !moreup, ] >= 3
		NORM <- matrix(TRUE, nrow(NP), ncol(NP))
		NORM[ii & moreup, ] <- notUp
		NORM[ii & !moreup, ] <- notDown
		normalNP <- normalNP*NORM
	}
	flagAltered <- which(propAlt > 0.5)
	envir[["flagAlteredNP"]] <- flagAltered
	normalNP[normalNP == FALSE] <- NA
	envir[["normalNP"]] <- normalNP
}
