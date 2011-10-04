\name{sample.CNSet}
\alias{sample.CNSet}
\alias{sample.CNSet2}
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
\usage{

data(sample.CNSet)
data(sample.CNSet2)

}
\format{
  The data illustrates the \code{CNSet-class}, with
	\code{assayData} containing the quantile-normalized
	intensities for the A and B alleles, genotype calls and
	confidence scores.  New slots that specific to copy number
	estimation are \code{batch} and \code{batchStatistics}.

}
\examples{
\dontshow{
\dontrun{

      ## hapmap phase 3 data
      data(hapmapSet, package="CnvScripts")
      marker.index <- which(chromosome(hapmapSet) == 8)
      marker.index <- marker.index[1:60e3]
      sample.CNSet <- hapmapSet[marker.index, c(1168:1169)]
      ## 2 samples, many markers

      ## all samples, a few markers
      snp.index <- which(isSnp(hapmapSet))[1:100]
      np.index <- which(!isSnp(hapmapSet))[1:100]
      marker.index <- c(snp.index, np.index)
      sample.CNSet2 <- hapmapSet[marker.index, ]
}
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
sample.CNSet$SNR[]
## the skew
sample.CNSet$SKW[]
## the gender (gender is imputed unless specified in the call to crlmm)
table(sample.CNSet$gender[])  ## 1=male, 2=female
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

## --------------------------------------------------
## calculating allele-specific copy number
## --------------------------------------------------
(ca <- CA(sample.CNSet, i=1:5, j=1:2))
## copy number for allele B, first 5 markers, first 2 samples
(cb <- CB(sample.CNSet, i=1:5, j=1:2))
index <- which(!isSnp(sample.CNSet))[1:5]
cn2 <- CA(sample.CNSet, i=index, j=1:2)
## note, cb is 0 at nonpolymorphic loci
CB(sample.CNSet, i=index, j=1:2)
## A shortcut for total copy number
cn3 <- totalCopynumber(sample.CNSet, i=1:5, j=1:2)
all.equal(cn3, ca+cb)
cn4 <- totalCopynumber(sample.CNSet, i=index, j=1:2)
all.equal(cn4, cn2)

## markers 1-5, all samples
cn5 <- totalCopynumber(sample.CNSet, i=1:5)
## all markers, samples 1-2
cn6 <- totalCopynumber(sample.CNSet, j=1:2)
}
\keyword{datasets}
