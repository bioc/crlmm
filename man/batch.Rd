\name{batch}
\alias{batch}
\alias{batch,eSet-method}
\title{
	Function to extract batch information.
}
\description{
	Checks the phenoData and protocolData for a variable named
	batch and, if present, returns the vector.
}
\usage{
batch(object)
}
\arguments{
  \item{object}{
  An object extending the \code{eSet} class.
}
}
\details{

	For copy number estimation, a batch variable must be
	specified.  Currently, we suggest storing this variable in the
	protocolData.

	Batch represents groups of samples that were processed (DNA
	preparation and collection, PCR amplification, scan date) at
	similar times. Often, the 96 well chemistry plate or scan date
	is a useful surrogate for batch.

}
\value{

	Vector indicating batch.

}

\author{
R. Scharpf
}

\seealso{
	\code{\link{genotype}}, \code{\link{genotype2}}
}
\keyword{manip}

