\name{SnpQSet-methods}
\alias{A}
\alias{A,SnpQSet-method}
\alias{A<-}
\alias{A<-,SnpQSet,matrix-method}
\alias{B}
\alias{B,SnpQSet-method}
\alias{B<-}
\alias{B<-,SnpQSet,matrix-method}
\alias{isSnp}
\alias{isSnp,SnpQSet-method}
\title{Indicator for polymorphic probes}
\description{
  
  This functions uses the annotation slot of the SnpQSet object to load
  the corresponding annotation package and determine whether each probe
  in the object interrogates a polymorphic or nonpolymorphic allele. For
  instance, in the Affy 6.0 platform roughly 900,000 of the 1.8 million
  markers are for polymorphic alleles.
  
}
\usage{
isSnp(object)
}
\arguments{
  \item{object}{SnpQSet object}
}
\value{
  an indicator: 0=nonpolymorphic, 1=polymorphic
}
\keyword{manip}
