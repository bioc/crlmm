\name{CNSetLM-methods}
\alias{lM}

\title{Methods for CNSetLM class}
\description{ Accessors for CNSetLM class}
\usage{
	lM(object)
}

\arguments{
  \item{object}{\code{CNSetLM}}
}
\details{
	\code{lM} returns a list (or an ffdf object if large data support is enabled) of the parameters estimated from a linear model fit for each SNP.  The parameters are batch and locus-specific.
}
\value{object of class \code{ffdf} or \code{list}}

\author{R. Scharpf}

\seealso{\code{\link{crlmmCopynumber}}, \code{\link{crlmmCopynumber2}}, \code{\link{CNSetLM-class}}}
\keyword{manip}
