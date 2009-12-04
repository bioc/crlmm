\name{methods-eSet}
%%\alias{chromosome}
\alias{chromosome,eSet-method}
%%\alias{position}
\alias{position,eSet-method}
\title{Methods for eSet derivatives}

\description{
  \code{chromosome} return the chromosome number for SNP and NP probes.
  \code{position} return the physical position of the SNP or the
  physical position for the first index of the NP probe.
}

\section{Methods}{
  \describe{
    \item{\code{chromosome(object)}:}{ accessor for chromosome}
    \item{\code{position(object)}:}{accessor for genomic position}
  }
}
\keyword{manip}