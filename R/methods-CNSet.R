setMethod("posteriorMean", signature(object="CNSet"), function(object) assayDataElement(object, "posteriorMean"))
setReplaceMethod("posteriorMean", signature(object="CNSet", value="matrix"), function(object, value) assayDataElementReplace(object, "posteriorMean", value))

linearParamElementReplace <- function(obj, elt, value) {
    storage.mode <- storageMode(batchStatistics(obj))
    switch(storage.mode,
           "lockedEnvironment" = {
               aData <- copyEnv(batchStatistics(obj))
               if (is.null(value)) rm(list=elt, envir=aData)
               else aData[[elt]] <- value
               Biobase:::assayDataEnvLock(aData)
               batchStatistics(obj) <- aData
           },
           "environment" = {
               if (is.null(value)) rm(list=elt, envir=batchStatistics(obj))
               else batchStatistics(obj)[[elt]] <- value
           },
           list = batchStatistics(obj)[[elt]] <- value)
    obj
}

## parameters

## allele A
##   autosome SNPs
##   autosome NPs
##   chromosome X NPs for women
C1 <- function(object, marker.index, batch.index, sample.index){
	acn <- matrix(NA, nrow=length(marker.index), ncol=length(sample.index))
	for(k in seq_along(batch.index)){
		l <- batch.index[k]
		## calculate cn for all the samples that are in this batch
		jj <- sample.index[as.character(batch(object))[sample.index] == batchNames(object)[l]]
		bg <- nuA(object)[marker.index, l]
		slope <- phiA(object)[marker.index, l]
		I <- as.matrix(A(object)[marker.index, jj])
		acn[, match(jj, sample.index)] <- 1/slope * (I - bg)
	}
	return(as.matrix(acn))
}

## allele B  (treated allele 'A' for chromosome X NPs)
##   autosome SNPs
##   chromosome X for male nonpolymorphic markers
C2 <- function(object, marker.index, batch.index, sample.index, NP.X=FALSE){
	acn <- matrix(NA, nrow=length(marker.index), ncol=length(sample.index))
	for(k in seq_along(batch.index)){
		l <- batch.index[k]
		jj <- sample.index[as.character(batch(object))[sample.index] == batchNames(object)[l]]
		bg <- nuB(object)[marker.index, l]
		slope <- phiB(object)[marker.index, l]
		if(!NP.X){
			I <- as.matrix(B(object)[marker.index, jj])
		} else I <- as.matrix(A(object)[marker.index, jj])
		acn[, match(jj, sample.index)] <- 1/slope * (I - bg)
	}
##	if(length(acn) > 1){
##		acn <- do.call("cbind", acn)
##	} else acn <- acn[[1]]
	return(as.matrix(acn))
}

## Chromosome X SNPs
C3 <- function(object, allele, marker.index, batch.index, sample.index){
##	acn <- vector("list", length(batch.index))
	acn <- matrix(NA, nrow=length(marker.index), ncol=length(sample.index))
	for(k in seq_along(batch.index)){
		l <- batch.index[k]
		##j <- which(as.character(batch(object))[sample.index] == batchNames(object)[l])
		jj <- sample.index[as.character(batch(object))[sample.index] == batchNames(object)[l]]
		##phiA2 and phiB2 are not always estimable  -- need both men and women
		phiA2 <- phiPrimeA(object)[marker.index, l]
		phiB2 <- phiPrimeB(object)[marker.index, l]
		phiA <- phiA(object)[marker.index, l]
		phiB <- phiB(object)[marker.index, l]
		nuA <- nuA(object)[marker.index, l]
		nuB <- nuB(object)[marker.index, l]
		phiA <- phiA(object)[marker.index, l]
		IA <- as.matrix(A(object)[marker.index, jj])
		IB <- as.matrix(B(object)[marker.index, jj])
		## I = nu + acn * phi
		## acn = 1/phi* (I-nu)
		phistar <- phiB2/phiA
		tmp <- (IB - nuB - phistar*IA + phistar*nuA)/phiB
		CB <- tmp/(1-phistar*phiA2/phiB)
		index <- which(is.na(phiB2) | is.na(phiA2))
		if(length(index) > 0){
			cb <- 1/phiB[index] * (IB[index, ] - nuB[index])
			CB[index, ] <- cb
		}
		if(allele == "B"){
			acn[, match(jj, sample.index)] <- CB
			##acn[[k]] <- CB
		}
		if(allele == "A"){
			CA <- (IA-nuA-phiA2*CB)/phiA
			if(length(index) > 0){
				ca <- 1/phiA[index] * (IA[index, ] - nuA[index])
 				CA[index, ] <- ca
			}
			acn[, match(jj, sample.index)] <- CA
		}
	}
##	if(length(acn) > 1){
##		acn <- do.call("cbind", acn)
##	} else acn <- acn[[1]]
	return(as.matrix(acn))
}




ACN <- function(object, allele, i , j){
	if(missing(i) & missing(j)) stop("must specify rows (i) or columns (j)")
	is.ff <- is(calls(object), "ff") | is(calls(object), "ffdf")
	missing.i <- missing(i)
	missing.j <- missing(j)
	if(!missing.i){
		is.ann <- !is.na(chromosome(object)[i])
		is.X <- chromosome(object)[i]==23 & is.ann
		is.auto <- chromosome(object)[i] < 23 & is.ann
		is.snp <- isSnp(object)[i] & is.ann
	} else{
		is.ann <- !is.na(chromosome(object))
		is.X <- chromosome(object)==23 & is.ann
		is.auto <- chromosome(object) < 23 & is.ann
		is.snp <- isSnp(object) & is.ann
		i <- 1:nrow(object)
	}
	## Define batch.index and sample.index
	if(!missing.j) {
		batches <- unique(as.character(batch(object))[j])
		##batches <- as.character(batch(object)[j])
		batch.index <- match(batches, batchNames(object))
	} else {
		batch.index <- seq_along(batchNames(object))
		j <- 1:ncol(object)
	}
	nr <- length(i)
	nc <- length(j)
	acn <- matrix(NA, nr, nc)
	dimnames(acn) <- list(featureNames(object)[i],
			      sampleNames(object)[j])
	if(allele == "A"){
		if(is.ff){
			open(nuA(object))
			open(phiA(object))
			open(A(object))
		}
		## --
		## 4 types of markers for allele A
		##--
		## 1. autosomal SNPs or autosomal NPs
		if(any(is.auto)){
			auto.index <- which(is.auto)
			marker.index <- i[is.auto]
			acn[auto.index, ] <- C1(object, marker.index, batch.index, j)
		}
		if(any(is.X)){
			##2. CHR X SNPs (men and women)
			if(any(is.snp)){
				if(is.ff) {
					open(phiPrimeA(object))
					open(phiPrimeB(object))
					open(phiB(object))
					open(nuB(object))
					open(B(object))
				}
				marker.index <- i[is.X & is.snp]
				acn.index <- which(is.X & is.snp)
				acn[acn.index, ] <- C3(object, allele="A", marker.index, batch.index, j)
				if(is.ff) {
					close(phiPrimeA(object))
					close(phiPrimeB(object))
					close(phiB(object))
					close(nuB(object))
					close(B(object))
				}
			}
			if(any(!is.snp)){
				marker.index <- i[is.X & !is.snp]
				acn.index <- which(is.X & !is.snp)
				acn[acn.index, ] <- NA
				female.index <- j[object$gender[j] == 2]
				## 3. CHR X NPs: women
				if(length(female.index) > 0){
					female.batch.index <- match(unique(as.character(batch(object))[female.index]), batchNames(object))
					jj <- which(object$gender[j] == 2)
					acn[acn.index, jj] <- C1(object, marker.index, female.batch.index, female.index)
				}
				male.index <- j[object$gender[j] == 1]
				if(length(male.index) > 0){
					if(is.ff){
						open(nuB(object))
						open(phiB(object))
					}
					male.batch.index <- match(unique(as.character(batch(object))[male.index]), batchNames(object))
					jj <- which(object$gender[j] == 1)
					acn[acn.index, jj] <- C2(object, marker.index, male.batch.index, male.index, NP.X=TRUE)
					if(is.ff){
						close(nuB(object))
						close(phiB(object))
					}
				}
			}
		}
		if(is.ff){
			close(nuA(object))
			close(phiA(object))
			close(A(object))
		}
	}
	if(allele == "B"){
		if(is.ff){
			open(nuB(object))
			open(phiB(object))
			open(B(object))
		}
		if(any(!is.snp)){
			acn.index <- which(!is.snp)
			acn[acn.index, ] <- 0
		}
		if(any(is.auto)){
			auto.index <- which(is.auto & is.snp)
			if(length(auto.index) > 0){
				marker.index <- i[auto.index]
				acn[auto.index, ] <- C2(object, marker.index, batch.index, j)
			}
		}
		if(any(is.X)){
			if(is.ff){
				open(phiPrimeA(object))
				open(phiPrimeB(object))
				open(phiA(object))
				open(nuA(object))
				open(A(object))
			}
			marker.index <- i[is.X & is.snp]
			acn.index <- which(is.X & is.snp)
			acn[acn.index, ] <- C3(object, allele="B", marker.index, batch.index, j)
			if(is.ff){
				close(phiPrimeA(object))
				close(phiPrimeB(object))
				close(phiA(object))
				close(nuA(object))
				close(A(object))
			}
			if(any(!is.snp)){
				acn.index <- which(!is.snp)
				marker.index <- i[!is.snp]
				acn[acn.index, ] <- 0
			}
		}
	}
	if(is.ff){
		close(nuB(object))
		close(phiB(object))
		close(B(object))
	}
	return(acn)
}

setMethod("CA", signature=signature(object="CNSet"),
	  function(object, ...){
		  ca <- ACN(object, allele="A", ...)
		  ca[ca < 0] <- 0
		  ca[ca > 5] <- 5
		  return(ca)
	  })
setMethod("CB", signature=signature(object="CNSet"),
	  function(object, ...) {
		  cb <- ACN(object, allele="B", ...)
		  cb[cb < 0] <- 0
		  cb[cb > 5] <- 5
		  return(cb)
	  })

setMethod("totalCopynumber", signature=signature(object="CNSet"),
	  function(object, ...){
		  ca <- CA(object, ...)
		  cb <- CB(object, ...)
		  return(ca+cb)
	  })
rawCopynumber <- totalCopynumber

setMethod("posteriorProbability", signature(object="CNSet"),
	  function(object, predictRegion, copyNumber=0:4){
		  logI <- array(NA, dim=c(nrow(object), ncol(object), 2), dimnames=list(NULL, NULL, LETTERS[1:2]))
		  getIntensity <- function(object){
			  logI[, , 1] <- log2(A(object))
			  logI[, , 2] <- log2(B(object))
			  return(logI)
		  }
		  logI <- getIntensity(object)
		  ##gts <- lapply(as.list(0:4), genotypes)
		  prob <- array(NA, dim=c(nrow(object), ncol(object), length(copyNumber)),
				dimnames=list(NULL,NULL, paste("copynumber",copyNumber, sep="")))
		  for(i in seq_along(copyNumber)){
			  G <- genotypes(copyNumber[i])
			  P <- array(NA, dim=c(nrow(object), ncol(object), length(G)),
				     dimnames=list(NULL,NULL,G))
			  for(g in seq_along(G)){
				  gt <- G[g]
				  P[, , gt] <- dbvn(x=logI, mu=predictRegion[[gt]]$mu, Sigma=predictRegion[[gt]]$cov)
			  }
			  ## the marginal probability for the total copy number
			  ##  -- integrate out the probability for the different genotypes
			  for(j in 1:ncol(P)){
				  prob[, j, i] <- rowSums(as.matrix(P[, j, ]), na.rm=TRUE)
			  }
		  }
		  ## divide by normalizing constant.  Probs across copy number states must sum to one.
		  ##nc <- apply(prob, c(1,3), sum, na.rm=TRUE)
		  nc <- matrix(NA, nrow(object), ncol(object))
		  for(j in seq_len(ncol(object))){
			  nc <- rowSums(as.matrix(prob[,j, ]), na.rm=TRUE)
			  prob[, j, ] <- prob[, j, ]/nc
		  }
		  return(prob)
	  })

setMethod("calculatePosteriorMean", signature(object="CNSet"),
	  function(object, posteriorProb, copyNumber=0:4, w, ...){
		  if(missing(w)) w <- rep(1/length(copyNumber),length(copyNumber))
		  stopifnot(sum(w)==1)
		  stopifnot(dim(posteriorProb)[[3]] == length(copyNumber))
		  pm <- matrix(0, nrow(object), ncol(object))
		  for(i in seq_along(copyNumber)){
			  pm <- pm + posteriorProb[, , i] * copyNumber[i] * w[i]
		  }
		  return(pm)
	  })


## for a given copy number, return a named list of bivariate normal prediction regions
##   - elements of list are named by genotype
##   'AAA'  list should have 'mu' and 'Sigma'
##                mu is a R x 2 matrix, R is number of features
##                Sigma is a R x 4 matrix.  Each row can be coerced to a 2 x 2 matrix
##          For nonpolymorphic markers, 2nd column of mu is NA and elements 2-4 of Sigma are NA
##
setMethod("predictionRegion", signature(object="CNSet", copyNumber="integer"),
	  function(object, copyNumber=0:4){
		  ## would it be more efficient to store as an array?
		  ## mu: features x genotype x allele
		  ## Sigma: features x genotype x covariance
		  stopifnot(all(copyNumber %in% 0:4))
		  gts <- lapply(as.list(copyNumber), genotypes)
		  nms <- unlist(gts)
		  res <- vector("list", length(nms))
		  ##names(res) <- paste("copyNumber", copyNumber, sep="")
		  names(res) <- nms
		  nu <- getNu(object)
		  phi <- getPhi(object)
		  mus <- matrix(NA, nrow(nu), 2, dimnames=list(NULL, LETTERS[1:2]))
		  Sigma <- matrix(NA, nrow(nu), 3)
		  bivariateCenter <- function(nu, phi){
			  ##  lexical scope for mus, CA, CB
			  mus[,1] <- log2(nu[1] + CA * phi[1])
			  mus[,2] <- log2(nu[2] + CB * phi[2])
			  mus
		  }
		  for(i in seq_along(copyNumber)){
			  G <- genotypes(copyNumber[i])
			  tmp <- vector("list", length(G))
			  names(tmp) <- G
			  CN <- copyNumber[i]
			  for(g in seq_along(G)){
				  gt <- G[g]
				  CB <- g-1
				  CA <- CN-CB
				  gt.corr <- genotypeCorrelation(gt)
				  nma <- ifelse(CA == 0, "tau2A.BB", "tau2A.AA")
				  nmb <- ifelse(CB == 0, "tau2B.AA", "tau2B.BB")
				  Sigma[,1] <- getVar(object, nma)
				  Sigma[,3] <- getVar(object, nmb)
				  Sigma[,2] <- getCor(object, gt.corr)
				  res[[gt]]$mu <- bivariateCenter(nu, phi)
				  res[[gt]]$cov <- Sigma
			  }
		  }
		  return(res)
	  })
