\name{methods-eSet}
\alias{cnIndex}
\alias{cnIndex,eSet-method}
\alias{snpIndex}
\alias{snpIndex,eSet-method}


\title{Methods for eSet derivatives}

\description{

  \code{snpIndex} return the row indices of polymorphic loci.

  \code{cnIndex} return the row indices of copy number probes.

}

\usage{

  cnIndex(object)

  snpIndex(object)

}

\arguments{
  
  \item{object}{Any object extending eSet (ABset, CopyNumberSet,
    CrlmmSetList, SnpSet}

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