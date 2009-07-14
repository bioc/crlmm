\name{crlmmWrapper}
\Rdversion{1.1}
\alias{crlmmWrapper}
\title{ Function to quantile normalize the cell files and genotype the
  polymorphic loci.
}
\description{
  Quantile normalizes the cel files to a target reference distribution
  (both polymorphic and nonpolymorphic loci).  Genotype calls and
  confidence scores are assigned using the crlmm algorithm.
}
\usage{
crlmmWrapper(filenames, outdir = "./", cdfName = "genomewidesnp6", save.it=FALSE, splitByChr=TRUE, ...)
}
\arguments{
  \item{filenames}{
    Complete path to cel files
  }
  \item{outdir}{
    Where to store results.
  }
  \item{cdfName}{
    Annotation package.
  }
  \item{save.it}{
    Whether to save intermediate results. Sometimes this can be helpful.
  }
  \item{splitByChr}{
    Logical. Whether to split results by chromosome (default is TRUE)
  }
  \item{\dots}{
    additional arguments to function \code{crlmm}
  }
}
\value{
  An object of class \code{CrlmmSetList}
}

\author{
  R. Scharpf
}
\seealso{
  \code{\linkS4class{CrlmmSetList}}
  \code{\link{crlmm}}
}
\examples{
\dontrun{
library(crlmm)
library(genomewidesnp6Crlmm)
celFiles <- list.files("PATH_TO_CELS", full.name=TRUE)
crlmmList <- crlmmWrapper(celFiles, outdir="./", save.it=TRUE, intensityFile="intensities.rda")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{robust}
\keyword{methods}
\keyword{list}
