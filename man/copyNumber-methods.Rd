\name{copyNumber-methods}
\docType{methods}
\alias{copyNumber-methods}
\alias{copyNumber,CNSet-method}
\title{Computes the total copy number}
\description{

	This function calculates the total copy number by adding the
	copy number for the two alleles at polymorphic loci, and
	returning the copy number for one allele at nonpolymorphic
	loci.

}
\section{Methods}{
\describe{
\item{\code{signature(object = "CNSet")}}{Calculates CA + CB, where CA is the copy number for allele A and CB is the copy number for allele B}
}}
\seealso{\code{\linkS4class{oligoClasses}{CNSet}}}
\keyword{methods}

