\name{cnSetExample}
\alias{cnSetExample}
\alias{cnSetExample2}
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

data(cnSetExample)
data(cnSetExample2)

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
      cnSetExample <- sample.CNSet
      save(cnSetExample, file="~/Software/crlmm/data/cnSetExample.rda")
      ## 2 samples, many markers

      ## all samples, a few markers
      snp.index <- which(isSnp(hapmapSet))[1:100]
      np.index <- which(!isSnp(hapmapSet))[1:100]
      marker.index <- c(snp.index, np.index)
      sample.CNSet2 <- hapmapSet[marker.index, ]
      sample.CNSet2$gender <- sample.CNSet2$gender[]
      sample.CNSet2$SKW <- sample.CNSet2$SKW[]
      sample.CNSet2$SNR <- sample.CNSet2$SNR[]
      cnSetExample2 <- sample.CNSet2
      save(cnSetExample2, file="~/Software/crlmm/data/cnSetExample2.rda")
}
}
data(cnSetExample)
## --------------------------------------------------
## accessors for the feature-level info
## --------------------------------------------------
chromosome(cnSetExample)[1:5]
position(cnSetExample)[1:5]
isSnp(cnSetExample)[1:5]
table(isSnp(cnSetExample))
## --------------------------------------------------
## sample-level statistics computed by crlmm
## --------------------------------------------------
varLabels(cnSetExample)
## accessors for sample-level statistics
## The signal to noise ratio (SNR)
cnSetExample$SNR[]
## the skew
cnSetExample$SKW[]
## the gender (gender is imputed unless specified in the call to crlmm)
table(cnSetExample$gender[])  ## 1=male, 2=female
## --------------------------------------------------
## batchStatistics
## -------------------------------------------------- estimate of
## intercept from linear model
dim(nu(cnSetExample, "A"))
## background for the A allele in the 2 batches for the
## first 5 markers
nu(cnSetExample, "A")[1:5, ]
## background for the B allele in the 2 batches for the
## first 5 markers
nu(cnSetExample, "B")[1:5, ]
## the slope
phi(cnSetExample, "A")[1:5, ]

## --------------------------------------------------
## calculating allele-specific copy number
## --------------------------------------------------
(ca <- CA(cnSetExample, i=1:5, j=1:2))
## copy number for allele B, first 5 markers, first 2 samples
(cb <- CB(cnSetExample, i=1:5, j=1:2))
index <- which(!isSnp(cnSetExample))[1:5]
cn2 <- CA(cnSetExample, i=index, j=1:2)
## note, cb is 0 at nonpolymorphic loci
CB(cnSetExample, i=index, j=1:2)
## A shortcut for total copy number
cn3 <- totalCopynumber(cnSetExample, i=1:5, j=1:2)
all.equal(cn3, ca+cb)
cn4 <- totalCopynumber(cnSetExample, i=index, j=1:2)
all.equal(cn4, cn2)

## markers 1-5, all samples
cn5 <- totalCopynumber(cnSetExample, i=1:5)
## all markers, samples 1-2
cn6 <- totalCopynumber(cnSetExample, j=1:2)
}
\keyword{datasets}
