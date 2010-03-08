\name{crlmmCopynumber}
\alias{crlmmCopynumber}
\title{
	Estimates copy number for many Affymetrix and Illumina Platforms
}

\description{ 

	      This function consists of the following steps:

	1. quantile normalizes the intensities at
	   polymorphic and nonpolymorphic loci

	2. summarizes the intensities for multiple probes at a single locus by
	   taking the median

	3. calls biallelic genotypes at all polymorphic loci using the
	   crlmm algorithm.  Confidence scores of the genotype calls are
	   also computed

	4. Develops prediction regions for integer copy number at both
	polymorphic and nonpolymorphic loci. 

	- For polymorphic loci, the prediction regions are bivariate
	normal ellipses on the log A versus log B intensities.

	- For nonpolymorphic loci, the prediction regions are
          univariate intervals based on the log 'A' intensities

	  The prediction regions can be used directly by downstream
	  hidden Markov models (e.g., the VanillaICE package)

	5.  Transforms the normalized intensities to allele-specific
	estimates of copy number at polymorphic loci (see the
	accessors CA and CB).  At nonpolymorphic loci, estimates of
	total copy number are stored in the CA slot.

	The total copy number can be used for downstream segmentation
	algorithms (e.g., circular binary segmentation).

}	

\usage{
crlmmCopynumber(cnOptions, ...)
}

\arguments{
  \item{cnOptions}{
  All the options for genotyping and copy number estimation.
}
  \item{\dots}{
  Not implemented  
}
}
\details{

The cnOptions is a list containing all of the options for this
function.  In particular, the following options can be useful:

\code{outdir}: where to store output files

\code{cdfName}: package to use for extracting annotation information
(see the 'see also' section below)

 \code{load.it}: when TRUE, an attempt is made to load previously
saved normalized intensities and genotype calls.  In particular, one
could redo the copy number estimation steps without having to
re-normalize or genotype all the SNPs.

\code{chromosome}: An integer 1-23, where 23 corresponds to chromosome
'X'.  This is only relevant for \code{crlmmCopynumber} and is ignored
during preprocessing and genotyping.

}
\value{

	Nothing is returned.  Output is saved to .rda files.

}

\author{
	R. Scharpf
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

	See \code{\link{cnOptions}} for all of the options available
	for the estimation procedure.  See \code{\link{crlmm}} for the
	genotyping algorithm.  See \code{\link{snprma}} for the
	normalization.

	For valid \code{cdfName}, see  \code{\link{validCdfNames}}

	For preprocessing and genotyping, see \code{\link{preprocessAndGenotype}}.

	For copy number estimation, see \code{\link{crlmmCopynumber}}.


}

\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

