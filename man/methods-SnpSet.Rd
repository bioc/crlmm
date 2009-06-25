\name{calls}
\alias{calls}
\alias{calls,SnpSet-method}
\alias{confs}
\alias{confs,SnpSet-method}
\title{Accessors for Calls and Confidences on a SnpSet object}
\description{
  \code{calls} returns the genotype calls. CRLMM stores genotype calls
  as integers (1 - AA; 2 - AB; 3 - BB).

  \code{confs} returns the confidences associated to the genotype
  calls. THe current implementation of CRLMM stores the confidences as
  integers by using the transformation:

  conf = round(-1000*log2(1-p)),

  where 'p' is the posterior probability of the call.
}
\usage{
calls(object)
confs(object)
}
\arguments{
  \item{object}{SnpSet object}
}
\value{
  Matrix of genotype calls or confidences.
  }
\examples{
  set.seed(1)
  theCalls <- matrix(sample(1:3, 20, rep=TRUE), nc=2)
  p <- matrix(runif(20), nc=2)
  theConfs <- round(-1000*log2(1-p))
  obj <- new("SnpSet", call=theCalls, callProbability=theConfs)
  calls(obj)
  confs(obj)
}
\keyword{manip}
