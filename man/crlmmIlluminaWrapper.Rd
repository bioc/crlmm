\name{crlmmIlluminaWrapper}
\Rdversion{1.1}
\alias{crlmmIlluminaWrapper}
\title{
  Reads idat files, quantile normalizes the raw intensities, and
  provides diallelic genotype calls at polymorphic loci.
}
\description{
  This function is a wrapper for utilities used to preprocess and
  genotype files from the illumina platform.   See Section see also.
}
\usage{
crlmmIlluminaWrapper(outdir = "./", cdfName, save.intermediate = FALSE, splitByChr = TRUE, intensityFile, ...)
}
\arguments{
  \item{outdir}{Character string.  The complete path to the directory
    where output will be saved.}
  \item{cdfName}{Character string.  Annotation package.  See
    crlmm:::validCdfNames() for supported packages.
  }
  \item{save.intermediate}{Logical.  Whether to save intermediate
    files. This is used primarily for debugging.
  }
  \item{splitByChr}{Logical.  Whether to split output files by
    chromosome.  This is recommended for datasets with 50 or more samples.
  }
  \item{intensityFile}{Character string.  Name of intermediate file for storing
    quantile normalized intensities.  This is always saved regardless of
    the value of save.intermediate.
  }
  \item{\dots}{Additional arguments to readIdatFiles.
  }
}
\details{
}
\value{
  An object of class \code{CrlmmSetList}.  This object will have the
  following two elements:

  - an object of class \code{ABset}
  - an object of class \code{SnpSet}
  
}
\author{
  Rob Scharpf
}
\note{
}
\seealso{
  \code{\linkS4class{CrlmmSetList}}
  \code{\linkS4class{ABset}}
  \code{\link[Biobase]{SnpSet}}
}
\examples{
}
\keyword{manip}

