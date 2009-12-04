\name{CNSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CNSet-class}
\alias{CA}
\alias{CA<-}
\alias{CB}
\alias{CB<-}
\alias{CA,CNSet-method}
\alias{CA<-,CNSet,matrix-method}
\alias{CB,CNSet-method}
\alias{CB<-,CNSet,matrix-method}
\alias{computeHmm}
\alias{computeHmm,CNSet-method}
\alias{copyNumber,CNSet-method}
\alias{ellipse,CNSet-method}
\alias{emissionPr}
\alias{emissionPr<-}
\alias{emissionPr,CNSet-method}
\alias{emissionPr<-,CNSet,array-method}
\alias{segmentData}
\alias{segmentData<-}
\alias{segmentData,CNSet-method}
\alias{segmentData<-,CNSet,RangedData-method}

\title{Container for allele-specific copy number and lower level SNP /
  NP summaries}
\description{Container for allele-specific copy number and lower level SNP /
  NP summaries}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CNSet", phenoData, featureData, annotation, experimentData, protocolData, calls, callsConfidence, senseThetaA, senseThetaB, CA, CB, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SnpCallSetPlus}"}, directly.
Class \code{"\linkS4class{SnpQSet}"}, by class "SnpCallSetPlus", distance 2.
Class \code{"\linkS4class{SnpCallSet}"}, by class "SnpCallSetPlus", distance 2.
Class \code{"\linkS4class{QuantificationSet}"}, by class "SnpCallSetPlus", distance 3.
Class \code{"\linkS4class{SnpLevelSet}"}, by class "SnpCallSetPlus", distance 3.
Class \code{"\linkS4class{eSet}"}, by class "SnpCallSetPlus", distance 4.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "SnpCallSetPlus", distance 5.
Class \code{"\linkS4class{Versioned}"}, by class "SnpCallSetPlus", distance 6.
}
\section{Methods}{
  \describe{
    \item{CA}{\code{signature(object = "CNSet")}: ... }
    \item{CA<-}{\code{signature(object = "CNSet", value = "matrix")}: ... }
    \item{CB}{\code{signature(object = "CNSet")}: ... }
    \item{CB<-}{\code{signature(object = "CNSet", value = "matrix")}: ... }
    \item{computeHmm}{\code{signature(object = "CNSet")}: ... }
    \item{copyNumber}{\code{signature(object = "CNSet")}: ... }
    \item{ellipse}{\code{signature(x = "CNSet")}: ... }
	 }
}
\author{R. Scharpf}
\examples{
showClass("CNSet")
}
\keyword{classes}
