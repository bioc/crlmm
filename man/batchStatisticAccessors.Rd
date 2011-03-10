\name{batchStatisticAccessors}
\alias{Ns}
\alias{corr}
\alias{tau2}
\alias{mads}
\alias{medians}

\title{

	Accessors for batch-specific summary statistics.

}

\description{

	The summary statistics stored here are used by the tools for
	copy number estimation.

}

\usage{
corr(object, ...)
tau2(object, ...)
mads(object,...)
medians(object,...)
Ns(object,...)
}

\arguments{
  \item{object}{  An object of class \code{CNSet}.}
  \item{\dots}{
  An additional argument named 'i' can be passed to subset the markers
  and an argument 'j' can be passed to subset the batches.  Other
  arguments are ignored.
  }
}

\value{
	
An array with dimension R x A x G x C, or R x G x C.  

R: number of markers
A: number of alleles (2)
G: number of biallelic genotypes (3)
C: number of batches

\code{Ns} returns an array of genotype frequencies stratified by
batch.  Dimension R x G x C.

\code{corr} returns an array of within-genotype correlations
(log2-scale) stratified by batch. Dimension R x G x C.

\code{medians} returns an array of the within-genotype medians
(intensity-scale) stratified by batch and allele. Dimension R x A x G
x C.

\code{mads} returns an array of the within-genotype median absolute
deviations (intensity-scale) stratified by batch and allele. Dimension
is the same as for \code{medians}.

\code{tau2} returns an array of the squared within-genotype median
absolute deviation on the log-scale.  Only the mads for AA and BB
genotypes are stored.  Dimension is R x A x G x C, where G is AA or
BB.  Note that the mad for allele A/B for subjects with genotype BB/AA
is a robust estimate of the background variance, whereas the the mad
for allele A/B for subjects with genotype AA/BB is a robust estimate
of the variance for copy number greater than 0 (we assume that on the
log-scale the variance is rougly constant for CA, CB > 0).

}


\seealso{

	\code{\link{batchStatistics}}

}

\examples{
data(sample.CNSet)
## All NAs. Need to replace sample.CNSet with a HapMap example
Ns(cnSet, i=1:5, j=1:2)
corr(cnSet, i=1:5, j=1:2)
medians(cnSet, i=1:5, j=1:2)
mads(cnSet, i=1:5, j=1:2)
tau2(cnSet, i=1:5, j=1:2)
}
\keyword{manip}
