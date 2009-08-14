\name{crlmmWrapper}
\Rdversion{1.1}
\alias{crlmmWrapper}
\title{ Wrapper to quantile normalize and genotype Affymetrix 6.0 cel files.
}
\description{
  Quantile normalizes the cel files to a target reference distribution
  (both polymorphic and nonpolymorphic loci).  Genotype calls and
  confidence scores are assigned using the crlmm algorithm.
}
\usage{
crlmmWrapper(filenames, cdfName = "genomewidesnp6", load.it=FALSE,
save.it=FALSE, splitByChr=TRUE, crlmmFile, intensityFile, rgFile,
platform=c("affymetrix", "illumina")[1], ...)
}
\arguments{
  \item{filenames}{
    Complete path to cel files
  }
  \item{cdfName}{
    Annotation package (currently, only 'genomewidesnp6' is supported).
  }
  \item{load.it}{ Logical.  If TRUE, the wrapper function will try to
    load the 'intensityFile' and 'crlmmFile'. When FALSE, the quantile
    normalization and genotype calling steps are recomputed regardless
    of whether these objects exist.}
  \item{save.it}{ Logical. Whether to save intermediate files.  When
    TRUE, the SnpSet object containing genotype calls is saved to the
    file indicated by 'crlmmFile'. This can be useful if an error occurs
    in the crlmmWrapper function after the genotyping has been completed.
  }
  \item{splitByChr}{
    Logical. Whether to split results by chromosome (default is TRUE).
    Helpful for datasets containing 50+ samples.  Files are saved as
    \code{paste(dirname(crlmmFile), "crlmmSetList_", CHR, ".rda", sep="")}.
  }
  \item{crlmmFile}{Character string.  Name of file (including the
    complete path) containing the genotype calls.  }
  \item{intensityFile}{Character string.  Name of file (with
    complete path) containing the quantile-normalized intensities.  }  
  \item{rgFile}{Character string.  Name of file (including the
    complete path) containing the RG intensities (only applicable for
    illumina platform).  }
  \item{platform}{Character string.  Currently only \code{affymetrix}
    and \code{illumina} are supported.}
  \item{\dots}{
    additional arguments to function \code{readIdatFiles} (illumina
    platform only)
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
crlmmSetList <- crlmmWrapper(celFiles, intensityFile="intensities.rda")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{robust}
\keyword{methods}
\keyword{list}
