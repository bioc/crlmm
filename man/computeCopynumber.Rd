\name{computeCopynumber}
\alias{computeCopynumber}
\title{Computes copy number}
\description{
  
  A function that transforms the quantile-normalized fluorescence
  intensities of the polymorphic and nonpolymorphic probes to a copy
  number scale.  
  
}
\usage{
computeCopynumber(chrom, A, B, calls, conf, NP, plate, MIN.OBS=1, 
envir, P, DF.PRIOR = 50, CONF.THR = 0.99, bias.adj=FALSE,
priorProb, gender=NULL, SNR, SNRmin, seed=123, verbose=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chrom}{Chromosome (an integer).  Use 23 for X and 24 for Y.}
  \item{A}{The A allele intensities from \code{snprma}}
  \item{B}{The B allele intensities from \code{snprma}}
  \item{calls}{The genotype calls from \code{crlmm}}
  \item{conf}{The genotype confidence scores from \code{crlmm}}
  \item{NP}{The quantile normalized intensities of the nonpolymorphic probes}
  \item{plate}{The bach variable.  Should be the same length as the
    number of columns in A}
  \item{MIN.OBS}{Integer: The minimum number of observations in a genotype
    cluster for which a SNP is deemed complete.}
  \item{envir}{An environment to save intermediate objects}
  \item{P}{Mainly for debugging a particular plate/batch.}
  \item{DF.PRIOR}{The degrees of freedom for the prior.  Higher numbers
    will shrink the variance and correlation more.}
  \item{CONF.THR}{A threshold for the genotype confidence scores.
    Genotypes with scores below the threshold are ignored when computing
  SNP-specific  within-genotype estimates of location and scale.}
  \item{bias.adj}{Logical: whether to adjust the location and scale
    parameters to account for biases due to common copy number
    variants.  This is a SNP-specific adjustment.  Parameters for
    background and slope must have already been estimated and available
    from the environment variable.}
  \item{priorProb}{Numerical vector of length 4.  The prior probability
    of each copy number state (0, 1, 2, 3, and 4). The default is a
    uniform prior.  Ignored if bias.adj=FALSE}
  \item{gender}{Gender of subjects. If not specified, we predict the
    gender from the X chromosome.}
  \item{SNR}{Signal to noise ratio from crlmm.}
  \item{SNRmin}{The minimum value for the SNR -- we suggest 5. Samples
    with SNR below SNRmin are excluded.}
  \item{seed}{Seed used for random samples}
  \item{verbose}{Logical: verbose output}
  \item{\dots}{Currently ignored}
}

\details{
  
  This function essentially transforms intensities from the raw scale to
  a copy number scale.  We assume that the median within-SNP intensity
  across samples is 2.  We make no assumption about the chromosomal copy
  number. This function is mostly useful for detecting rare variants
  (e.g., variants that would not affect the SNP-specific quantile-based
  estimators of location and scale).  A correction for more common
  variants is coming, as well as improved estimates of the uncertainty.
  
}

\value{

  For now, nothing is returned. All objects created by this function are
  stored in the environment passed to this function.  See the copynumber
  vignette for details on accessing these objects.
  
}
\references{Nothing yet.}
\author{Rob Scharpf}
\keyword{manip}

