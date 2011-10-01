\name{sample.CNSet}
\alias{sample.CNSet}
\alias{cnSet}
\docType{data}
\title{
	Object of class 'CNSet'
}
\description{

	The data for the first 16 polymorphic markers in the HapMap analysis.

}
\details{
  This object was created from the copynumber vignette in inst/scripts.
}
\usage{data(sample.CNSet)}
\format{
  The data illustrates the \code{CNSet-class}, with
	\code{assayData} containing the quantile-normalized
	intensities for the A and B alleles, genotype calls and
	confidence scores.  New slots that specific to copy number
	estimation are \code{batch} and \code{batchStatistics}.

}
\examples{
\dontrun{
      ## hapmap phase 3 data
      data(hapmapSet, package="CnvScripts")
      marker.index <- which(chromosome(hapmapSet) == 8)
      marker.index <- marker.index[1:60e3]
      sample.CNSet <- hapmapSet[marker.index, c(1168:1169)]
      ## 2 samples, many markers
      save(sample.CNSet, file="~/Software/crlmm/data/sample.CNSet.rda")

      ## all samples, a few markers
      snp.index <- which(isSnp(hapmapSet))[1:100]
      np.index <- which(!isSnp(hapmapSet))[1:100]
      marker.index <- c(snp.index, np.index)
      sample.CNSet2 <- hapmapSet[marker.index, ]
      save(sample.CNSet2, file="~/Software/crlmm/data/sample.CNSet2.rda")
}
data(sample.CNSet)
## --------------------------------------------------
## accessors for the feature-level info
## --------------------------------------------------
chromosome(sample.CNSet)[1:5]
position(sample.CNSet)[1:5]
isSnp(sample.CNSet)[1:5]
table(isSnp(sample.CNSet))
## --------------------------------------------------
## sample-level statistics computed by crlmm
## --------------------------------------------------
varLabels(sample.CNSet)
## accessors for sample-level statistics
## The signal to noise ratio (SNR)
sample.CNSet$SNR[1:5]
## the skew
sample.CNSet$SKW[1:5]
## the gender (gender is imputed unless specified in the call to crlmm)
table(sample.CNSet$gender)  ## 1=male, 2=female
## --------------------------------------------------
## batchStatistics
## -------------------------------------------------- estimate of
## intercept from linear model
dim(nu(sample.CNSet, "A"))
## background for the A allele in the 2 batches for the
## first 5 markers
nu(sample.CNSet, "A")[1:5, ]
## background for the B allele in the 2 batches for the
## first 5 markers
nu(sample.CNSet, "B")[1:5, ]
## the slope
phi(sample.CNSet, "A")[1:5, ]
## correlation within genotype cluster AA
##corr(sample.CNSet, "AA")[1:5, ]
#### correlation within genotype cluster AB
##corr(sample.CNSet, "AB")[1:5, ]
#### correlation within genotype cluster BB
##corr(sample.CNSet, "BB")[1:5, ]
## --------------------------------------------------

## --------------------------------------------------
## calculating allele-specific copy number
## --------------------------------------------------
## copy number for allele A, first 5 markers, first 2 samples
(ca <- CA(sample.CNSet, i=1:5, j=1:2))
## copy number for allele B, first 5 markers, first 2 samples
(cb <- CB(sample.CNSet, i=1:5, j=1:2))
## total copy number for first 5 markers, first 2 samples
(cn1 <- ca+cb)

## total copy number at first 5 nonpolymorphic loci
index <- which(!isSnp(sample.CNSet))[1:5]
cn2 <- CA(sample.CNSet, i=index, j=1:2)
## note, cb is NA at nonpolymorphic loci
(cb <- CB(sample.CNSet, i=index, j=1:2))
## note, ca+cb will give NAs at nonpolymorphic loci
CA(sample.CNSet, i=index, j=1:2) + cb
## A shortcut for total copy number
cn3 <- totalCopynumber(sample.CNSet, i=1:5, j=1:2)
all.equal(cn3, cn1)
cn4 <- totalCopynumber(sample.CNSet, i=index, j=1:2)
all.equal(cn4, cn2)

## markers 1-5, all samples
cn5 <- totalCopynumber(sample.CNSet, i=1:5)
## all markers, samples 1-5
cn6 <- totalCopynumber(sample.CNSet, j=1:5)
}
\keyword{datasets}
