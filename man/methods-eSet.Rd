\name{methods-eSet}
\alias{chromosome}
\alias{chromosome,eSet-method}
\alias{cnIndex}
\alias{cnIndex,eSet-method}
\alias{position}
\alias{position,eSet-method}
\alias{snpIndex}
\alias{snpIndex,eSet-method}

\title{Methods for eSet derivatives}

\description{
  \code{chromosome} return the chromosome number for SNP and NP probes.
    
  \code{cnIndex, ...} return the row indices of copy number probes.

  \code{position} return the physical position of the SNP or the
  physical position for the first index of the NP probe.
  
  \code{snpIndex, ...} return the row indices of polymorphic loci.

}

\usage{

  chromosome(object)
  cnIndex(object, ...)
  position(object)
  snpIndex(object, ...)

}

\arguments{
  
  \item{object}{Any object extending eSet (ABset, CopyNumberSet,
    CrlmmSetList, SnpSet}

  \item{\dots}{Usually the character string indicating annotation package.}

}

\value{

  Vector of indices.

}

\seealso{

  \code{"\linkS4class{ABset}"}
  \code{"\linkS4class{CopyNumberSet}"}
  \code{"\linkS4class{CrlmmSetList}"}
  
}
  

\keyword{manip}