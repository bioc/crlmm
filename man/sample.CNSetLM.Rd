\name{sample.CNSetLM}
\alias{sample.CNSetLM}
\docType{data}
\title{
	Dataset of class 'CNSetLM'
}
\description{

	The data for 2119 polymorphic and nonpolymorphic markers on
	chromosome 1 for the CEPH and Yoruban HapMap samples.

}
\usage{data(sample.CNSetLM)}
\format{

	This class has been deprecated.  See example below for how to
	update an existing 'CNSetLM' object to class 'CNSet'.

	The data illustrates the \code{CNSetLM-class}, with
	\code{assayData} containing the quantile-normalized
	intensities for the A and B alleles, genotype calls and
	confidence scores (call and callProbability), and
	allele-specific copy number (CA and CB). The parameters from
	the linear model are stored in the lM slot.

}
\examples{
## class CNSetLM has been deprecated
data(sample.CNSetLM)
## update to class CNSet
cnSet <- as(sample.CNSetLM, "CNSet")
all(isCurrent(cnSet)) ## is the cnSet object current?

\dontrun{
	## updating class CNSetLM using ff objects
	## a bigger object with multiple batches
	if(require(ff)){
		outdir <- "/amber1/scratch/rscharpf/jss/hapmap2"
		load(file.path(outdir, "container.rda"))
		container <- object; rm(object); gc()
		container2 <- as(container, "CNSet")
		all(isCurrent(container2))
		## test replacement methods, subset methods
		table(batch(container2))
		##generates warning ... would need open, close in the '[' method
		invisible(open(nuA(container2)))
		xx <- nu(container2, "A")[1:5, ]
		nuA(container2)[1:5, ] <- xx
		invisible(close(nuA(container2)))
	}
}
## --------------------------------------------------
## accessors for the feature-level info
## --------------------------------------------------
chromosome(cnSet)[1:5]
position(cnSet)[1:5]
isSnp(cnSet)[1:5]
## 980 nonpolymorphic markers and 1139 polymoprhic markers
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
## --------------------------------------------------
##
## accessors for parameters estimated from the linear model for copy
## number (note that the parameters have dimension R x C, where R
## corresponds to the number of features and C corresponds to the
## number of batches)
## -------------------------------------------------- estimate of
## background
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

## NOTE: subsetting the object before extracting copy number
##       can be very inefficient when the data set is very large,
##       particularly if using ff objects.  IN particular, subsetting
##       the CNSet object will subset all of the assay data elements
##       and all of the elements in the LinearModelParameter slot
\dontrun{
	cnsubset <- cnSet[1:5, ]
}
}
\keyword{datasets}
