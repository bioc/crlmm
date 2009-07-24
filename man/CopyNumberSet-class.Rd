\name{CopyNumberSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CopyNumberSet-class}
\alias{batch,CopyNumberSet-method}
\alias{CA}
\alias{CA<-}
\alias{CA,CopyNumberSet-method}
\alias{CA<-,CopyNumberSet,matrix-method}
\alias{CB}
\alias{CB<-}
\alias{CB,CopyNumberSet-method}
\alias{CB<-,CopyNumberSet,matrix-method}
\alias{chromosome,CopyNumberSet-method}
\alias{copyNumber,CopyNumberSet-method}
\alias{ellipse,CopyNumberSet-method}
\alias{plot,ABset,CopyNumberSet-method}
\alias{position,CopyNumberSet-method}
\title{Class "CopyNumberSet"}
\description{Container for allele-specific estimates of copy number and
  the corresponding uncertainty}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CopyNumberSet", assayData, phenoData, featureData, experimentData, annotation, protocolData, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{assayData}:}{Object of class \code{"AssayData"}}
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}}
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}  }
    \item{\code{annotation}:}{Object of class \code{"character"}  }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}  }
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{

    \item{batch}{\code{signature(object="CopyNumberSet")}: Extracts the
    batch information used to estimate copy number.}

    \item{CA}{\code{signature(object = "CopyNumberSet", ...)}: Extracts the
      copy number for allele 'A' at polymorphic loci. At nonpolymorphic
      loci, the total copy number is returned.}

    \item{CA<-}{\code{signature(object = "CopyNumberSet",
    value="matrix")}: Replaces CA matrix with supplied value.}    

    \item{CB}{\code{signature(object = "CopyNumberSet", ...)}:
      Extracts copy number for allele 'B' at polymorphic loci. NAs are
      returned for nonpolymorphic loci.}    

    \item{CB<-}{\code{signature(object = "CopyNumberSet",
	value="matrix")}: Replaces CB matrix with supplied value .}

    
    \item{chromosome}{\code{signature(object = "CopyNumberSet")}:
      Extract chromosome.}

    \item{copyNumber}{\code{signature(object = "CopyNumberSet")}:
      Returns CA + CB.}

    \item{ellipse}{\code{signature(x = "CopyNumberSet", ...)}:
      Extracts parameters from \code{featureData} slot and draws
      prediction regions for the supplied copy number.}

    \item{plot}{\code{signature(x = "ABset", y = "CopyNumberSet")}:
      Physical position}
    \item{position}{\code{signature(object = "CopyNumberSet")}: physical position}
  }
}

\author{R. Scharpf}
\examples{
showClass("CopyNumberSet")

\dontrun{
##returns the copy number for allele A at polymorphic loci
CA(object[snpIndex(object), ])
##returns the total copy number at nonpolymorphic loci
CA(object[cnIndex(object), ])
}
}
\keyword{classes}
