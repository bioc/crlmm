\name{CNSetLM-class}
\Rdversion{1.1}
\docType{class}
\alias{CNSetLM}
\alias{CNSetLM-class}
\alias{[,CNSetLM-method}
\alias{corr,CNSetLM,character-method}
\alias{lines,CNSetLM-method}
\alias{lM,CNSetLM-method}
\alias{lM<-,CNSetLM,list_or_ffdf-method}
\alias{open,CNSetLM-method}
\alias{nu,CNSetLM,character-method}
\alias{phi,CNSetLM,character-method}
\alias{show,CNSetLM-method}
\alias{sigma2,CNSetLM,character-method}
\alias{tau2,CNSetLM,character-method}

\title{CNSetLM class}
\description{Container for allele-specific copy number and linear model
parameters}
\section{Objects from the Class}{
		 Objects from the class can be created by calls of the
		 form \code{ new("CNSetLM", CA=matrix(), CB=matrix(), alleleA=matrix(), alleleB=matrix(), call=matrix(), callProbability=matrix())}
		 
}
\section{Slots}{
  \describe{
    \item{\code{lM}:}{Object of class \code{"list_or_ffdf"}}
    \item{\code{assayData}:}{Object of class \code{"AssayData"} }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} }
    \item{\code{annotation}:}{Object of class \code{"character"}  }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} }
  }
}
\section{Extends}{
Class \code{"\linkS4class{CNSet}"}, directly.
Class \code{"\linkS4class{SnpSuperSet}"}, by class "CNSet", distance 2.
Class \code{"\linkS4class{AlleleSet}"}, by class "CNSet", distance 3.
Class \code{"\linkS4class{SnpSet}"}, by class "CNSet", distance 3.
Class \code{"\linkS4class{eSet}"}, by class "CNSet", distance 4.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "CNSet", distance 5.
Class \code{"\linkS4class{Versioned}"}, by class "CNSet", distance 6.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "CNSetLM")}: subset \code{CNSetLM} objects}
    \item{lines}{\code{signature(x="CNSetLM")}: for drawing prediction regions on A vs B scatterplots}
    \item{lM}{\code{signature(object = "CNSetLM")}: Extract list or
    ffdf object containing linear model parameters}
    \item{nu}{\code{signature(object = "CNSetLM"), \code{signature(allele="character")}}: 
            intercept for linear model. See \code{\link{nu}}.}
  \item{open}{\code{signature(con = "CNSetLM")}: opens file connects
    to ff objects for assayData elements and linear model parameters}
    \item{phi}{\code{signature(object = "CNSetLM"), \code{signature(allele="character")}}: 
            slope for linear model.  See \code{\link{phi}}.}
    \item{show}{\code{signature(object = "CNSetLM")}: print method
    for the class }
    \item{sigma2}{\code{signature(object = "CNSetLM"), \code{signature(allele="character")}}: 
    The within genotype         }   
	 }
}
\author{ R. Scharpf}
\seealso{
	\code{\linkS4class{SnpSuperSet}}, \code{\linkS4class{CNSet}}
}

\examples{
showClass("CNSetLM")
}
\keyword{classes}
