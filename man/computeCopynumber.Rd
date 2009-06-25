\name{computeCopynumber}
\alias{computeCopynumber}
\title{Computes copy number}
\description{
  A function that transforms the quantile-normalized fluorescence
  intensities of the polymorphic and nonpolymorphic probes to a copy
  number scale.  
  
}
\usage{
computeCopynumber(object, CHR, bias.adj=FALSE, batch, SNRmin=5, cdfName="genomewidesnp6", ...)
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
  \item{\dots}{arguments to \code{.computeCopynumber} -- not called by
    the user.}
}

\details{
  
  This function transforms the intensities to a copy number scale.  We
  assume that the median within-SNP intensity across samples is 2.  We
  make no assumption about the chromosomal copy number. This function is
  useful for detecting rare variants (e.g., variants that would not
  effect the SNP-specific quantile-based estimators of location and
  scale).  When \code{bias.adj=TRUE}, an extra iteration is performed
  whereby samples with a high posterior probability of having a
  non-normal copy number are excluded.  The resulting within-genotype
  estimators of location and scale are more robust to a large number of
  samples having a copy number variant.

}

\value{
  
  An object of class \code{CopyNumberSet}.

}

\seealso{
  \code{\linkS4class{CopyNumberSet-class}}
}
\references{Nothing yet.}
\author{Rob Scharpf}
\keyword{manip}
\keyword{models}
\keyword{regression}
\keyword{robust}

