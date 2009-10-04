\name{ABset-class}
\Rdversion{1.1}
\docType{class}
\alias{ABset-class}
\alias{A}
\alias{B}
\alias{A,ABset-method}
\alias{B,ABset-method}
\alias{A<-,ABset,matrix-method}
\alias{B<-,ABset,matrix-method}
\title{Class "ABset"}
\description{Containter for quantile-normalized intensities}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ABset", assayData, phenoData, featureData, experimentData, annotation, protocolData, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{assayData}:}{Object of class \code{"AssayData"} }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} }
    \item{\code{annotation}:}{Object of class \code{"character"}  }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}  }
  }
}
\section{Extends}{
Class \code{\link[Biobase:class.eSet]{eSet}}, directly.
Class \code{\link[Biobase:class.VersionedBiobase]{VersionedBiobase}}, by class "eSet", distance 2.
Class \code{\link[Biobase:class.Versioned]{Versioned}}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{
    \item{A}{\code{signature(object="ABset")}: accessor for the
      quantile-normalized intensities of allele A for polymorphic probes
      and the quantile normalized intensities for the copy number
      probes.}

    \item{"A<-"}{\code{signature(object="ABset", value="matrix")}:
      replacement method for the A allele intensities.}

    \item{"B<-"}{\code{signature(object="ABset", value="matrix")}:
    replacement method for the B allele intensities.}    
    
    \item{B}{\code{signature(object="ABset")}: accessor for the
      quantile-normalized intensities of allele B for polymorphic probes.}    
    \item{plot}{\code{signature(x = "ABset", y = "CopyNumberSet")}:
... }}

}
\author{R. Scharpf}
\examples{
showClass("ABset")
}
\keyword{classes}
