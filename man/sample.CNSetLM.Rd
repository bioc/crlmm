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
all(isCurrent(cnSet))
}
\keyword{datasets}
