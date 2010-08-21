\name{constructIlluminaCNSet}
\alias{constructIlluminaCNSet}
\title{
	Construct an instance of CNSetLM after preprocessing Illumina files
}
\description{

	Assemble the preprocessed data and genotype calls from
	\code{crlmmIllumina} to initialize a \code{CNSetLM} object.

}
\usage{
constructIlluminaCNSet(crlmmResult, path, snpFile, cnFile)
}
\arguments{
  \item{crlmmResult}{

  A \code{SnpSet} object returned by function \code{crlmmIllumina} or
  \code{crlmmIllumina2}.

}

  \item{path}{path to files created by \code{crlmmIllumina}}
	
  \item{snpFile}{
  The \code{snpFile} filename specified in \code{crlmmIllumina}.
}
  \item{cnFile}{
  The \code{cnFile} filename specified in \code{crlmmIllumina}.
}
}

\value{
	An object of class \code{CNSetLM}.
}
\author{
R. Scharpf
}
\seealso{
	\code{\link{CNSetLM-class}}, \code{\link{crlmmIllumina}}
}

\keyword{manip}

