\name{calculatePosteriorMean}
\alias{calculatePosteriorMean}
\alias{calculatePosteriorMean,CNSet-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the posterior mean of the integer copy numbers}
\description{

  Computes the posterior mean copy number from the bivariate normal
  prediction regions.

}

\usage{
calculatePosteriorMean(object, posteriorProb, copyNumber = 0:4, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A \code{CNSet} object.
}
  \item{posteriorProb}{
    The posterior probability for copy numbers 0-4.
  }

  \item{copyNumber}{
    Integer vector.
}
  \item{\dots}{
    Ignored
}
}
\details{
}
\value{
  An array  (features x samples x copy number)
}

\author{
R. Scharpf
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{predictionRegion}}, \code{\link{posteriorProbability}}
}
\examples{
data(sample.CNSet)
pr <- predictionRegion(sample.CNSet, copyNumber=0:4)
pp <- posteriorProbability(sample.CNSet, predictRegion=pr)
pm <- calculatePosteriorMean(sample.CNSet, posteriorProb=pp)

## multiple batches
data(sample.CNSet2)
pr <- predictionRegion(sample.CNSet2, copyNumber=0:4)
pp <- posteriorProbability(sample.CNSet2, predictRegion=pr)
## prediction regions not available for some batches
pm <- calculatePosteriorMean(sample.CNSet2, posteriorProb=pp)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{array}
\keyword{models}
