\name{computeCopynumber}
\alias{computeCopynumber}
\title{Computes copy number}
\description{
  A function that transforms the quantile-normalized fluorescence
  intensities of the polymorphic and nonpolymorphic probes to a copy
  number scale.  
  
}
\usage{
computeCopynumber(object, CHR, bias.adj=FALSE, batch, SNRmin=5, cdfName,
...)
}
\arguments{
  \item{object}{object of class \code{CrlmmSetList}.}
  \item{CHR}{Chromosome (an integer).  Use 23 for X and 24 for Y.}
  \item{bias.adj}{Logical.  When TRUE, an extra iteration is performed
    (see details) .}
  \item{batch}{Factor.  Must be the same length as the number of samples
    in the object.}
  \item{SNRmin}{The minimum value for the SNR -- we suggest 5. Samples
    with SNR below SNRmin are dropped.}
  \item{cdfName}{Annotation package }
  \item{\dots}{arguments to \code{.computeCopynumber}.}
}

\details{

  One may altenatively use the \code{update} method (a wrapper for
  \code{computeCopynumber}) to add a \code{CopyNumberSet} object to the
  \code{CrlmmSetList} object.

  This function requires 10 or more samples to estimate model
  parameters.  Preferably, 70+ samples would be processed together in a
  batch.
  
  This function translates the quantile-normalized fluorescence
  intensities to the scale of copy number.  We assume that for any given
  locus the median copy number is two for each batch. When
  \code{bias.adj=TRUE}, an extra iteration is performed whereby samples
  with a high posterior probability of having a non-normal copy number
  are excluded.  The resulting within-genotype estimators of location
  and scale are more robust to a large number of samples having a copy
  number variant.

  For details, see the technical report:

Scharpf RB, Ruczinski I, Carvalho B, Doan B, Chakravarti A, and Irizarry
  R. A multilevel model to address batch effects in copy number
  estimation using SNP arrays
  (http://www.bepress.com/cgi/viewcontent.cgi?article=1193&context=jhubiostat).


}

\value{
  
  An object of class \code{CrlmmSetList}.

}

\seealso{
  \code{\link{update}},  \code{\linkS4class{CrlmmSetList}}, \code{\linkS4class{CopyNumberSet}},   \code{.computeCopynumber}
}
\author{Rob Scharpf}
\keyword{manip}
\keyword{models}
\keyword{regression}
\keyword{robust}

