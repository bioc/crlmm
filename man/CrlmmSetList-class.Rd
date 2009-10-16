\name{CrlmmSetList-class}
\Rdversion{1.1}
\docType{class}
\alias{CrlmmSetList-class}
\alias{[,CrlmmSetList-method}
\alias{$,CrlmmSetList-method}
\alias{addFeatureAnnotation}
\alias{addFeatureAnnotation,CrlmmSetList-method}
\alias{A<-}
\alias{B<-}
\alias{A,CrlmmSetList-method}
\alias{A<-,CrlmmSetList,matrix-method}
\alias{B,CrlmmSetList-method}
\alias{B<-,CrlmmSetList,matrix-method}
\alias{batch,CrlmmSetList-method}
\alias{CA,CrlmmSetList-method}
\alias{CB,CrlmmSetList-method}
\alias{CA<-,CrlmmSetList,matrix-method}
\alias{CB<-,CrlmmSetList,matrix-method}
\alias{calls,CrlmmSetList-method}
\alias{chromosome,CrlmmSetList-method}
\alias{cnIndex,CrlmmSetList-method}
\alias{combine,CrlmmSetList,CrlmmSetList-method}
\alias{confs,CrlmmSetList-method}
\alias{copyNumber,CrlmmSetList-method}
\alias{nrow,CrlmmSetList-method}
\alias{ncol,CrlmmSetList-method}
\alias{plot,CrlmmSetList,ANY-method}
\alias{points,CrlmmSetList-method}
\alias{position,CrlmmSetList-method}
\alias{sampleNames,CrlmmSetList-method}
\alias{show,CrlmmSetList-method}
\alias{snpIndex,CrlmmSetList-method}
\alias{update,CrlmmSetList-method}
\alias{update,character-method}


\title{Class "CrlmmSetList"}
\description{Container for quantile normalized intensities and genotype calls.}
\section{Objects from the Class}{
  Objects from the class are created by calls to the \code{crlmmWrapper}
  function. 
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
  Class \code{\link[Biobase:assayData]{assayData}}, by class
  "list", distance 2. 
}
\section{Methods}{
  \describe{
    \item{"["}{\code{signature(x = "CrlmmSetList")}: subsets both the
      features and samples of each element in th elist}
  
    \item{"$"}{\code{signature(x = "CrlmmSetList")}: used to access
    element from the \code{featureData} of the \code{CopyNumberSet}
    element.  In particular, useful for extract the SNP- and
    batch-specific parameters used to estimate copy number that are
    currently stored in the \code{featureData} of \code{CopyNumberSet}
    objects. }

  \item{"addFeatureAnnotation"}{\code{signature(object =
  "CrlmmSetList")}: Creates an object of class AnnotatedDataFrame from
  the CrlmmSetList object.  The new annotation object contains
  information on chromosome and physical position and has the same
  ordering of rows as the CrlmmSetList object. }

    \item{"A"}{\code{signature(object = "CrlmmSetList")}: extracts the
      quantile normalized intensities for the A allele in the
      \code{ABset} element.}

    \item{"A<-"}{\code{signature(object = "CrlmmSetList", value =
    "matrix")}: Replacement method for the A intensities.}    
  
    \item{"B"}{\code{signature(object = "CrlmmSetList")}: extracts the
      quantile normalized intensities for the B allele in the
      \code{ABset} element.}
    
    \item{"B<-"}{\code{signature(object = "CrlmmSetList", value =
    "matrix")}: Replacement method for the B intensities.}        
    
    \item{"batch"}{\code{signature(object = "CrlmmSetList")}: extracts the
      batch information used to estimate copy number.}    

    \item{"CA"}{\code{signature(object = "CrlmmSetList")}: extracts the
      copy number for allele A at polymorphic loci. For nonpolymorphic
      probes, CA returns the total copy number. }
  
    \item{"CB"}{\code{signature(object = "CrlmmSetList")}: extracts the
      copy number for allele B at polymorphic loci. For nonpolymorphic
      probes, CB returns 'NA'.}

    \item{"CA<-"}{\code{signature(object = "CrlmmSetList", value =
	"matrix")}: Replacement method for the CA estimates.}

    \item{"CB<-"}{\code{signature(object = "CrlmmSetList", value =
    "matrix")}: Replacement method for the CB estimates.}                

    \item{"calls"}{\code{signature(object = "CrlmmSetList")}: extracts
      the genotype calls from the \code{SnpSet} element.}

    \item{"cnIndex"}{\code{signature(object = "CrlmmSetList")}: returns
      the row indices of the copy number probes.}

    \item{"combine"}{\code{signature(x = "CrlmmSetList", y = "CrlmmSetList")}: combines
      objects of the class.  }

    \item{"confs"}{\code{signature(x = "CrlmmSetList")}: confidence
    scores for crlmm genotype calls. }    
    
    \item{"nrow"}{\code{signature(x = "CrlmmSetList")}: Number of rows
    (features) of each element in the list.  }
    
    \item{"ncol"}{\code{signature(x = "CrlmmSetList")}: Number of columns
    (samples) of each element in the list.  }    
	
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

    \item{"update"}{\code{signature(object = "CrlmmSetList")}: This
    method calls \code{computeCopynumber} to compute copy number for the
    \code{crlmmSetList} object.  The object returned by this method is
    itself an instance of class \code{CrlmmSetList}.  The third element
    of the list is an instance of \code{CopyNumberSet} containing the
    locus-level, allele-specific estimates of copynumber.  }

  \item{"update"}{\code{signature(object = "character")}: This method
    loads the CrlmmSetList file given as the argument to object --
    the object is a character string or a vector of character strings
    that provides the complete path to the CrlmmSetList file.
    The CrlmmSetList object is loaded and updated with copy number
    estimates -- an object of class CopyNumberSet becomes the third
    element in the CrlmmSetList object.  After updating the CrlmmSetList
    object, the object is saved to the same path.}
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
