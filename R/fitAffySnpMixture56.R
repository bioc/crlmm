fitAffySnpMixture56 <- function(S, M, knots, probs=rep(1/3, 3), eps=.01, maxit=10, verbose=FALSE){
  ##56 stands for 5 and 6 arrays but will also work for Illumina
  ##Note the unfortunate choice of numbering:
  ##1 is BB, 2 AB, and 3 AA. Opposite to everything else!
  ##this is legacy code I decided not to change.
  ## this why at the end we report -coefs: F1 is the negative f
  mus <- append(quantile(M, c(1, 5)/6, names=FALSE), 0, 1)
  sigmas <- rep(mad(c(M[M<mus[1]]-mus[1], M[M>mus[3]]-mus[3])), 3)
  sigmas[2] <- sigmas[2]/2
 
  weights <- apply(cbind(mus, sigmas), 1, function(p) dnorm(M, p[1], p[2]))
  previousF1 <- -Inf
  change <- eps+1
  it <- 0
 
  if(verbose) message("Max change must be under ", eps, ".")
  matS <- stupidSplineBasis(S, knots)
  while (change > eps & it < maxit){
    it <- it+1
    ## E
    z <- sweep(weights, 2, probs, "*")
    LogLik <- rowSums(z)
    z <- sweep(z, 1, LogLik, "/")
    probs <- colMeans(z)
 
    ## M
    fit1 <- crossprod(chol2inv(chol(crossprod(sweep(matS, 1, z[, 1], FUN="*"), matS))), crossprod(matS, z[, 1]*M))
 
    fit2 <- sum(z[, 2]*M)/sum(z[, 2])
    F1 <- matS%*%fit1
    sigmas[c(1, 3)] <- sqrt(sum(z[, 1]*(M-F1)^2)/sum(z[, 1]))
    sigmas[2] <- sqrt(sum(z[, 2]*(M-fit2)^2)/sum(z[, 2]))
 
    weights[, 1] <- dnorm(M, F1, sigmas[1])
    weights[, 2] <- dnorm(M, fit2, sigmas[2])
    weights[, 3] <- dnorm(M, -F1, sigmas[3])
    
    change <- max(abs(F1-previousF1))
    previousF1 <- F1
    if(verbose) message("Iter ", it, ": ", change, ".")
  }
  medF1 <- median(-F1)
 return(list(coef= -fit1, medF1=medF1, sigma1=sigmas[1], sigma2=sigmas[2]))
}

