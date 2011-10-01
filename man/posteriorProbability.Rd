\name{posteriorProbability}
\alias{posteriorProbability}
\alias{posteriorProbability,CNSet-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the posterior probability for integer copy numbers.}

\description{
  Calculate the posterior probability for integer copy numbers using the
  bivariate normal prediction regions.
}

\usage{
posteriorProbability(object, predictRegion, copyNumber = 0:4, w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A \code{CNSet} object.
}
  \item{predictRegion}{
    A list containing the bivariate normal prediction region for each of
    the possible genotypes.
}
  \item{copyNumber}{
    Integer vector.
  }

  \item{w}{ numeric vector of prior probabilities for each of the copy
    number states.  Must be the same length as \code{copyNumber} and sum
    to 1.}

}
\details{

}
\value{
An array (features x samples x copy number)
}

\author{
R. Scharpf
}
\note{
  %Remark about nonpolymorphic probes
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{predictionRegion}}, \code{\link{genotypes}}, \code{\link{calculatePosteriorMean}}
}
\examples{
data(sample.CNSet)
pr <- predictionRegion(sample.CNSet, copyNumber=0:4)
pp <- posteriorProbability(sample.CNSet, predictRegion=pr)
dim(pp)

## multiple batches
data(sample.CNSet2)
pr <- predictionRegion(sample.CNSet2, copyNumber=0:4)
pp <- posteriorProbability(sample.CNSet2, predictRegion=pr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{array}