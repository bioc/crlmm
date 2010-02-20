\name{AlleleSet-methods}
\docType{methods}

\title{Indicator for polymorphic probes}
\description{
  
  This functions uses the annotation slot of the AlleleSet object to load
  the corresponding annotation package and determine whether each probe
  in the object interrogates a polymorphic or nonpolymorphic allele. For
  instance, in the Affy 6.0 platform roughly 900,000 of the 1.8 million
  markers are for polymorphic alleles.
  
}
\usage{
A(object)
B(object)
}
\arguments{
  \item{object}{AlleleSet object}
}
\value{
	matrix of normalized intensities
}
\keyword{manip}
