\name{sample.CNSet}
\alias{sample.CNSet}
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
data(sample.CNSet)
## --------------------------------------------------
## accessors for the feature-level info
## --------------------------------------------------
chromosome(cnSet)[1:5]
position(cnSet)[1:5]
isSnp(cnSet)[1:5]
table(isSnp(cnSet))
## --------------------------------------------------
## sample-level statistics computed by crlmm
## --------------------------------------------------
varLabels(cnSet)
## accessors for sample-level statistics
## The signal to noise ratio (SNR)
cnSet$SNR[1:5]
## the skew
cnSet$SKW[1:5]
## the gender (gender is imputed unless specified in the call to crlmm)
table(cnSet$gender)  ## 1=male, 2=female
## --------------------------------------------------
## batchStatistics
## -------------------------------------------------- estimate of
## intercept from linear model
dim(nu(cnSet, "A"))
## background for the A allele in the 2 batches for the
## first 5 markers
nu(cnSet, "A")[1:5, ]
## background for the B allele in the 2 batches for the
## first 5 markers
nu(cnSet, "B")[1:5, ]
## the slope
phi(cnSet, "A")[1:5, ]
## correlation within genotype cluster AA
##corr(cnSet, "AA")[1:5, ]
#### correlation within genotype cluster AB
##corr(cnSet, "AB")[1:5, ]
#### correlation within genotype cluster BB
##corr(cnSet, "BB")[1:5, ]
## --------------------------------------------------

## --------------------------------------------------
## calculating allele-specific copy number
## --------------------------------------------------
## copy number for allele A, first 5 markers, first 2 samples
(ca <- CA(cnSet, i=1:5, j=1:2))
## copy number for allele B, first 5 markers, first 2 samples
(cb <- CB(cnSet, i=1:5, j=1:2))
## total copy number for first 5 markers, first 2 samples
(cn1 <- ca+cb)

## total copy number at first 5 nonpolymorphic loci
index <- which(!isSnp(cnSet))[1:5]
cn2 <- CA(cnSet, i=index, j=1:2)
## note, cb is NA at nonpolymorphic loci
(cb <- CB(cnSet, i=index, j=1:2))
## note, ca+cb will give NAs at nonpolymorphic loci
CA(cnSet, i=index, j=1:2) + cb
## A shortcut for total copy number
cn3 <- totalCopynumber(cnSet, i=1:5, j=1:2)
all.equal(cn3, cn1)
cn4 <- totalCopynumber(cnSet, i=index, j=1:2)
all.equal(cn4, cn2)

## markers 1-5, all samples
cn5 <- totalCopynumber(cnSet, i=1:5)
## all markers, samples 1-5
cn6 <- totalCopynumber(cnSet, j=1:5)
}
\keyword{datasets}
