\name{CNSet-methods}
\Rdversion{1.1}
\docType{methods}
\alias{CA,CNSet-method}
\alias{CB,CNSet-method}
\alias{lines,CNSet-method}
\alias{totalCopynumber,CNSet-method}

\title{crlmm methods for class "CNSet"}
\description{

	CNSet is a container defined in the oligoClasses package for
	 storing normalized intensities for genotyping platforms,
	 genotype calls, and parameters estimated for copy
	 number. Accessors for data that an object of this class
	 contains are largely defined in the package
	 oligoClasses. CNSet methods that involve more complex
	 calculations that are specific to the crlmm package, such as
	 computing allele-specific copy number, are included in crlmm
	 and described here.

}

\section{Methods}{
  \describe{
    \item{CA}{\code{signature(object="CNSet")}: ...}
    \item{CB}{\code{signature(object="CNSet")}: ...}
    \item{lines}{\code{signature(x="CNSet")}: ...}
    \item{totalCopynumber}{\code{signature(object="CNSet")}: ...}
	 }
}

\seealso{
  \code{\link{CNSet-class}}, \code{\link{CA}}, \code{\link{CB}}, \code{\link{totalCopynumber}}
}
\keyword{methods}
