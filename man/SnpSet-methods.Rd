\name{SnpSet-methods}
\docType{methods}
%\alias{calls}
\alias{calls,SnpSet-method}
\alias{calls<-,SnpSet,matrix-method}
%\alias{confs}
\alias{confs,SnpSet-method}
%\alias{confs<-}
\alias{confs<-,SnpSet,matrix-method}
\title{Accessors for Calls and Confidences on a SnpSet object}

\description{
  \code{calls} returns the genotype calls. CRLMM stores genotype calls
  as integers (1 - AA; 2 - AB; 3 - BB).

  \code{confs} returns the confidences associated to the genotype
  calls. The current implementation of CRLMM stores the confidences as
  integers by using the transformation:

  conf = round(-1000*log2(1-p)),

  where 'p' is the posterior probability of the call.
}

\section{Methods}{
\describe{
  \item{\code{initialize(SnpSet)}:}{Object instantiation, used by
    \code{new}; not be be called directly by the user.}  
  \item{\code{calls(object)}}{accessor for genotype calls}
  \item{\code{confs(object)}}{accessor for crlmm genotype confidence scores}  
}}

\examples{
  theCalls <- matrix(sample(1:3, 20, rep=TRUE), nc=2)
  p <- matrix(runif(20), nc=2)
  theConfs <- round(-1000*log2(1-p))
  obj <- new("SnpSet", call=theCalls, callProbability=theConfs)
  calls(obj)
  confs(obj)
}
\keyword{manip}
