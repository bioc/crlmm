\name{crlmmCopynumber}
\alias{crlmmCopynumber}
\alias{crlmmCopynumber2}
\alias{crlmmCopynumberLD}
\title{Locus- and allele-specific estimation of copy number}
\description{
  Locus- and allele-specific estimation of copy number.
}
\usage{
crlmmCopynumber(object, chromosome = 1:23, which.batches, MIN.SAMPLES
= 10, SNRMin = 5, MIN.OBS = 3, DF.PRIOR = 50, bias.adj = FALSE,
prior.prob = rep(1/4, 4), seed = 1, verbose = TRUE, GT.CONF.THR =
0.99, PHI.THR = 2^6, nHOM.THR = 5, MIN.NU = 2^3, MIN.PHI = 2^3,
THR.NU.PHI = TRUE, thresholdCopynumber = TRUE)

crlmmCopynumberLD(object, which.batches, MIN.SAMPLES = 10, SNRMin = 5,
MIN.OBS = 1, DF.PRIOR = 50, bias.adj = FALSE, prior.prob = rep(1/4,
4), seed = 1, verbose = TRUE, GT.CONF.THR = 0.99, PHI.THR = 2^6,
nHOM.THR = 5, MIN.NU = 2^3, MIN.PHI = 2^3, THR.NU.PHI = TRUE,
thresholdCopynumber = TRUE)

}
\arguments{
  \item{object}{object of class \code{SnpSuperSet}.
}

  \item{which.batches}{ Character vector with length equal to the number of
  samples.  Used to adjust for batch effects.  Chemistry plate or
  date often work well.  See examples.

  Ignored in crlmmCopynumberLD.

}

 \item{chromosome}{Numeric vector indicating which chromosomes to
 process (length <= 23). For chromosome X, use 23. A copy number
 method for chromosome Y is not yet available.
}

  \item{MIN.SAMPLES}{ 'Integer'.  The minimum number of samples in a
  batch.  Bathes with fewer than MIN.SAMPLES are skipped.  Therefore,
  samples in batches with fewer than MIN.SAMPLES have NA's for the
  allele-specific copy number and NA's for the linear model
  parameters.

}

  \item{SNRMin}{ Samples with low signal to noise ratios are
  excluded.  
}

  \item{MIN.OBS}{ 

  For genotypes with fewer than \code{MIN.OBS}, the within-genotype
  median is imputed from the observed genotypes.  For example, assume
  at at a given SNP genotypes AA and AB were observed and BB is an
  unobserved genotype.  For SNPs in which all 3 genotypes were
  observed, we fit the model E(mean_BB) = beta0 + beta1*mean_AA +
  beta2*mean_AB, obtaining estimates; of beta0, beta1, and beta2.  The
  imputed mean at the SNP with unobserved BB is then beta0hat +
  beta1hat * mean_AA of beta2hat * mean_AB.

}
  \item{DF.PRIOR}{

  The 2 x 2 covariance matrix of the background and signal variances
  is estimated from the data at each locus.  This matrix is then
  smoothed towards a common matrix estimated from all of the loci.
  DF.PRIOR controls the amount of smoothing towards the common matrix,
  with higher values corresponding to greater smoothing.  Currently,
  DF.PRIOR is not estimated from the data.  Future versions may
  estimate DF.PRIOR empirically.

}
  \item{bias.adj}{ 

  If \code{TRUE}, initial estimates of the linear model are updated
  after excluding samples that have a low posterior probability of
  normal copy number.  Excluding samples that have a low posterior
  probability can be helpful at loci in which a substantial fraction
  of the samples have a copy number alteration.  For additional
  information, see Scharpf et al., 2009.

  This argument is ignored in crlmmCopynumberLD.

}
  \item{prior.prob}{

  A numerical vector providing prior probabilities for copy number
  states corresponding to homozygous deletion, hemizygous deletion,
  normal copy number, and amplification, respectively.

}
  \item{seed}{ Seed for random number generation.}

  \item{verbose}{ Logical. }

  \item{GT.CONF.THR}{ 

    Confidence threshold for genotype calls (0, 1).  Calls with
    confidence scores below this theshold are not used to estimate the
    within-genotype medians.

}

  \item{PHI.THR}{ 
    SNPs with slopes (phi values) below this value are flagged.
    Flagged SNPs are not used in a regression to impute background and
    slope coefficients at nonpolymorphic loci.
}

  \item{nHOM.THR}{ 

  If fewer than \code{nHOM.THR} homozygous genotypes (AA or BB) are
    observed, the SNPs is flagged.  Flagged SNPs are not used in a
    regression to impute background and slope coefficients at
    nonpolymorphic loci.

}

  \item{MIN.NU}{ numeric. Minimum threshold for background. Ignored if \code{THR.NU.PHI} is \code{FALSE}.
}
  \item{MIN.PHI}{numeric. Minimum threshold for slope. Ignored if \code{THR.NU.PHI} is \code{FALSE}.
}
  \item{THR.NU.PHI}{
  If \code{THR.NU.PHI} is \code{FALSE}, \code{MIN.NU} and
  \code{MIN.PHI} are ignored.
}
  \item{thresholdCopynumber}{

  If \code{TRUE}, allele-specific number estimates are truncated.
  Values less than 0.05 are assigned the value 0.05; values exceeding
  5 are assigned the value 5.  

  Ignored in crlmmCopynumberLD.  Extreme values are automatically
  truncated.

}
}

\details{

	The function crlmmCopynumber uses matrices instead of ff
	objects if the ff library is not loaded.

	The function crlmmCopynumberLD allows parallel processing via
	and requires large data support via the ff package.

	We plan to phase out crlmmCopynumber and replace this function
	by crlmmCopynumberLD.

}
\author{R. Scharpf}
\examples{
## data(example.callSet)
## cnSet <- crlmmCopynumberLD(example.callSet)
## total copy number
## cn <- copyNumber(cnSet)
## allele-specific copy number
## ca <- CA(cnSet)/100 ## A dosage
## cb <- CB(cnSet)/100 ## B dosage
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
