\name{crlmmSet-class}
\Rdversion{1.1}
\docType{class}
\alias{crlmmSet-class}
\alias{calls}
\alias{calls,crlmmSet-method}
\alias{confs}
\alias{confs,crlmmSet-method}
\alias{list2crlmmSet}
\alias{list2crlmmSet,list-method}

\title{Class "crlmmSet"}
\description{An eSet extension to store genotype results.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("crlmmSet",
    assayData, phenoData, featureData, experimentData, annotation,
    ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{assayData}:}{Object of class \code{"AssayData"}. For an
      object to be valid, the assayData slot must have two matrices
      named 'calls' and 'confs'.}
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}}
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}}
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}}
    \item{\code{annotation}:}{Object of class \code{"character"}}
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}}
  }
}
\section{Extends}{
Class \code{"\link[Biobase]{eSet-class}"}, directly.
Class \code{"\link[Biobase]{VersionedBiobase-class}"}, by class "eSet", distance 2.
Class \code{"\link[Biobase]{Versioned-class}"}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{
    \item{calls}{\code{signature(x = "crlmmSet")}: extractor for
      genotype calls}
    \item{confs}{\code{signature(x = "crlmmSet")}: extractor for confidences}
	 }
}
\examples{
if (require(genomewidesnp5Crlmm) & require(hapmapsnp5)){
  path <- system.file("celFiles", package="hapmapsnp5")

  ## the filenames with full path...
  ## very useful when genotyping samples not in the working directory
  cels <- list.celfiles(path, full.names=TRUE)
  crlmmOutput <- crlmm(cels)
  crlmmSetObject <- list2crlmmSet(crlmmOutput)
  crlmmSetObject
  calls(crlmmSetObject)[1:10, 1:2]
  confs(crlmmSetObject)[1:10, 1:2]
}
}
\keyword{classes}
\keyword{methods}
