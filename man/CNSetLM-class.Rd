\name{CNSetLM-class}
\Rdversion{1.1}
\docType{class}
\alias{CNSetLM}
\alias{CNSetLM-class}
\alias{[,CNSetLM-method}
\alias{lM,CNSetLM-method}
\alias{lM<-,CNSetLM,list_or_ffdf-method}
\alias{open,CNSetLM-method}
\alias{show,CNSetLM-method}

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
    \item{lM}{\code{signature(object = "CNSetLM")}: Extract list or
    ffdf object containing linear model parameters}
##    \item{lM<-}{\code{signature(object = "CNSetLM", value = "list_or_ffdf")}: ... }
    \item{open}{\code{signature(con = "CNSetLM")}: opens file connects
    to ff objects for assayData elements and linear model parameters}
    \item{show}{\code{signature(object = "CNSetLM")}: print method
    for the class }
	 }
}
\author{ R. Scharpf}
\seealso{
	\code{\linkS4class{SnpSuperSet}}, \code{\linkS4class{CNSet}}
}
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("CNSetLM")
}
\keyword{classes}
