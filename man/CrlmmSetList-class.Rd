\name{CrlmmSetList-class}
\Rdversion{1.1}
\docType{class}
\alias{CrlmmSetList-class}
\alias{[,CrlmmSetList-method}
\alias{A,CrlmmSetList-method}
\alias{B,CrlmmSetList-method}
\alias{calls,CrlmmSetList-method}
\alias{cnIndex,CrlmmSetList-method}
\alias{combine,CrlmmSetList,CrlmmSetList-method}
\alias{plot,CrlmmSetList,ANY-method}
\alias{points,CrlmmSetList-method}
\alias{sampleNames,CrlmmSetList-method}
\alias{scanDates,CrlmmSetList-method}
\alias{show,CrlmmSetList-method}
\alias{snpIndex,CrlmmSetList-method}


\title{Class "CrlmmSetList"}
\description{Container for quantile normalized intensities and genotype calls.}
\section{Objects from the Class}{
  Objects from the class are created by calls to the \code{crlmmWrapper} function.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
  }    
%	   \item{\code{ABset}:}{Object of class \code{"ABset"}.
%	     Contains quantile normalized A and B intensities.  Quantile
%	   normalized intensities for nonpolymorphic probes are stored
%	   in the A slot.}
%	 \item{\code{crlmmResult}:}{Object of class
%	   \code{"SnpSet"}. Contains genotype calls and the
%	   corresponding confidence estimates.}
}

\section{Details}{

  Instances of \code{CrlmmSetList} are a list with two elements:

  1.  an object of class \code{ABset}

  2.  an object of class \code{SnpSet}

  The \code{featureNames} and \code{sampleNames} of both objects are
  required to be identical.

  Quantile-normalized intensities for the copy number probes are stored
  in the assay data element A of the ABset object.  Corresponding rows
  for the B allele are recorded as NAs.  Genotype calls for copy number
  probes in the SnpSet object are recorded as NAs.

}
  
\section{Extends}{
  Class \code{"\linkS4class{list}"}, from data part.
  Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
  Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{"["}{\code{signature(x = "CrlmmSetList")}: subsets both the }
  
    \item{"A"}{\code{signature(object = "CrlmmSetList")}: extracts the
      quantile normalized intensities for the A allele in the
      \code{ABset} element.}
    
  
    \item{"B"}{\code{signature(object = "CrlmmSetList")}: extracts the
      quantile normalized intensities for the B allele in the
      \code{ABset} element.}

    \item{"calls"}{\code{signature(object = "CrlmmSetList")}: extracts
      the genotype calls from the \code{SnpSet} element.}

    \item{"cnIndex"}{\code{signature(object = "CrlmmSetList")}: returns
      the row indices of the copy number probes.}

    \item{"combine"}{\code{signature(x = "CrlmmSetList", y = "CrlmmSetList")}: combines
      objects of the class.  }

    \item{"plot"}{\code{signature(x = "CrlmmSetList")}: For A versus B
      scatterplots.  }

    \item{"points"}{\code{signature(x = "CrlmmSetList")}: For A versus B
      scatterplots.  }

    \item{"sampleNames"}{\code{signature(object = "CrlmmSetList")}:
      Accessor for the column identifiers of the assay data elements.  }

    \item{"show"}{\code{signature(object = "CrlmmSetList")}:
      Shows the \code{ABset} and \code{SnpSet} elements of the list.  }


    \item{"snpIndex"}{\code{signature(object = "CrlmmSetList")}: returns
      the row indices of the polymorphic probes.}
    
    
  }
}


\section{Warnings}{

   \code{combine} will produce a warning as the \code{featureData} are
   not easily combined and therefore removed.  This method needs
   improvement.

}
\author{R. Scharpf}
\seealso{
  See also \code{"\link{crlmmWrapper}"},
  \code{"\linkS4class{ABset}"}
}
\examples{
showClass("CrlmmSetList")
}
\keyword{classes}
