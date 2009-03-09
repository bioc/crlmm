\name{cnrma}
\alias{cnrma}
\title{Quantile normalizes the intensities for the nonpolymorphic probe
  to a HapMap reference distribution.}
\description{Quantile normalizes the intensities for the nonpolymorphic probe
  to a HapMap reference distribution.
}
\usage{
cnrma(filenames, sns, cdfName, seed = 1, verbose=FALSE)
}
\arguments{
  \item{filenames}{filenames with complete path}
  \item{sns}{the sample names. If missing, the basename of the filenames is used.}
  \item{cdfName}{Should be 'genomewidesnp6Crlmm'}
  \item{seed}{seed for sampling the intensities to calculate skewness}
  \item{verbose}{logical}
}
\value{
  A list.  First element is the matrix of quantile-normalized
  intensities.  The second element is the skew.
}
\author{Rob Scharpf}
\note{Not tested}
\examples{
	\dontrun{
		library(genomewidesnp6Crlmm)
		library(hapmapsnp6)
		path <- system.file("celFiles", package="hapmapsnp6")
		celFiles <- list.celfiles(path, full.names=TRUE)
		cnrmaResult <- cnrma(celFiles)
	}
}
\keyword{robust}

