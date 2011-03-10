\name{CNSet-methods}
\Rdversion{1.1}
\docType{methods}
\alias{CA,CNSet-method}
\alias{CB,CNSet-method}
\alias{lines,CNSet-method}
\alias{totalCopynumber,CNSet-method}
\alias{rawCopynumber,CNSet-method}
\alias{nuA,CNSet-method}
\alias{nuB,CNSet-method}
\alias{Ns,CNSet-method}
\alias{corr,CNSet-method}
\alias{mads,CNSet-method}
\alias{medians,CNSet-method}
\alias{phiA,CNSet-method}
\alias{phiB,CNSet-method}
\alias{tau2,CNSet-method}

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
    \item{rawCopynumber}{\code{signature(object="CNSet")}: ...}
    \item{nuA}{\code{signature(object="CNSet")}: ...}
    \item{nuB}{\code{signature(object="CNSet")}: ...}
    \item{phiA}{\code{signature(object="CNSet")}: ...}
    \item{phiB}{\code{signature(object="CNSet")}: ...}
    \item{Ns}{\code{signature(object="CNSet")}: ...}
    \item{corr}{\code{signature(object="CNSet")}: ...}
    \item{mads}{\code{signature(x="CNSet")}: ...}
    \item{medians}{\code{signature(object="CNSet")}: ...}
    \item{tau2}{\code{signature(object="CNSet")}: ...}
	 }
}

\seealso{
  \code{\link{CNSet-class}}, \code{\link{CA}}, \code{\link{CB}},
  \code{\link{totalCopynumber}}, \code{\link{rawCopynumber}}
}
\keyword{methods}
