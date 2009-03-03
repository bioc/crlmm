\name{computeCopynumber}
\alias{computeCopynumber}
\title{Computes copy number}
\description{
  
  A function that transforms the quantile-normalized fluorescence
  intensities of the polymorphic and nonpolymorphic probes to a copy
  number scale.  
  
}
\usage{
computeCopynumber(A, B, calls, conf, NP, plate, fit.variance = NULL, MIN.OBS = 3, envir, chrom, P, DF.PRIOR = 50, CONF.THR = 0.99, trim = 0, upperTail = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{The A allele intensities from \Rfunction{snprma}}
  \item{B}{The B allele intensities from \Rfunction{snprma}}
  \item{calls}{The genotype calls from \Rfunction{crlmm}}
  \item{conf}{The genotype confidence scores from \Rfunction{crlmm}}
  \item{NP}{The quantile normalized intensities of the nonpolymorphic probes}
  \item{plate}{The bach variable.  Should be the same length as the
    number of columns in A}
  \item{fit.variance}{Temporary option to ignore.}
  \item{MIN.OBS}{The minimum number of observations in a genotype
    cluster for which a SNP is deemed complete.}
  \item{envir}{An environment to save intermediate objects}
  \item{chrom}{Is mainly useful for indicating sex chromosomes.}
  \item{P}{Mainly for debugging a particular plate/batch.}
  \item{DF.PRIOR}{The degrees of freedom for the prior.  Higher numbers
    will shrink the variance and correlation more.}
  \item{CONF.THR}{A threshold for the genotype confidence scores.
    Genotypes with scores below the threshold are ignored when computing
  SNP-specific  within-genotype estimates of location and scale.}
  \item{trim}{Temporary option: a number between 0 and 1.  Lops off a percentage of one tail of the
    intensities for the A and B alleles}
  \item{upperTail}{Logical: if TRUE and trim = p, quantiles greater than
  1-p are replaced with NAs}
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

