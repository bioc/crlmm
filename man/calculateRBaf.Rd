\name{calculateRBaf}
\alias{calculateRBaf}
\alias{calculateRBaf,CNSet-method}
\title{Calculate log R ratios and B allele frequencies.}
\description{Calculate log R ratios and B allele frequencies from
  a \code{CNSet} object}
\usage{
calculateRBaf(object, batch.name)
}

\arguments{
  \item{object}{A \code{CNSet} object.}
  \item{batch.name}{A character string.  See details.}
}
\details{\code{batch.name} must be a value in
  \code{batch(object)}. Currently, one must specify a single
  \code{batch.name}. If a character vector for \code{batch.name} is
  supplied, only the first is evaluated.

  TODO: A description of how these values are calculated.

}

\value{
  A list.

  \code{baf}: A matrix of B allele frequencies.

  \code{lrr}: A matrix of log R ratios.

}

\references{

  Reference for BAFs, LRRs.

}
\author{Lynn Mireless}

\examples{
data(cnSetExample)
baf.lrr <- calculateRBaf(cnSetExample, "SHELF")
hist(baf.lrr[["baf"]], breaks=100)
hist(baf.lrr[["lrr"]], breaks=100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{list}
\keyword{methods}