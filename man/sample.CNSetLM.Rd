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

	The data illustrates the \code{\link{CNSetLM-class}}, with
	\code{assayData} containing the quantile-normalized
	intensities for the A and B alleles, genotype calls and
	confidence scores (call and callProbability), and
	allele-specific copy number (CA and CB). The parameters from
	the linear model are stored in the lM slot.
	
}
\examples{
data(sample.CNSetLM)
}
\keyword{datasets}
