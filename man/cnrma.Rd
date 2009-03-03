\name{cnrma}
\alias{cnrma}
\title{Quantile normalizes the intensities for the nonpolymorphic probe
  to a HapMap reference distribution.}
\description{
}
\usage{
cnrma(filenames, sns, cdfName, seed = 1, verbose=FALSE)
}
\arguments{
  \item{filenames}{filenames with complete path}
  \item{sns}{the sample names. If missing, the basename of the filenames is used.}
  \item{cdfName}{Should be 'genomewidesnp6Crlmm'}
  \item{seed}{seed for sampling the intensities to calculate skewness}
}
\details{
}
\value{
  A list.  First element is the matrix of quantile-normalized
  intensities.  The second element is the skew.
}
\author{Rob Scharpf}
\note{Not tested}
\seealso{\code{\link{snprma}}}
\examples{
	\dontrun{
		library(genomewidesnp6Crlmm)
		library(hapmapsnp6)
		path <- system.file("celFiles", package="hapmapsnp5")
		celFiles <- list.celfiles(path, full.names=TRUE)
		cnrmaResult <- cnrma(celFiles)
	}
}
\keyword{robust}

