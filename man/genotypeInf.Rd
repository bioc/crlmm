\name{genotypeInf}
\alias{genotypeInf}

\title{

	Genotyping of Illumina Infinium II arrays.

}
\description{

	Genotyping of Illumina Infinium II arrays. This function
	provides CRLMM genotypes and confidence scores for the the
	polymorphic markers and is a required step prior to copy
	number estimation.

}
\usage{
genotypeInf(cnSet, mixtureParams, probs = rep(1/3, 3), SNRMin = 5, recallMin = 10, recallRegMin = 1000, verbose = TRUE, returnParams = TRUE, badSNP = 0.7, gender = NULL, DF = 6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cnSet}{An object of class \code{CNSet}}
  \item{mixtureParams}{
    \code{data.frame} containing mixture model parameters needed for
    genotyping.  The mixture model parameters are estimated from the
    \code{preprocessInf} function.
  }
  \item{probs}{'numeric' vector with priors for AA, AB and BB.}
  \item{SNRMin}{'numeric' scalar defining the minimum SNR used to filter
    out samples.}
  \item{recallMin}{Minimum number of samples for recalibration. }
  \item{recallRegMin}{Minimum number of SNP's for regression.}
  \item{verbose}{  'logical.'  Whether to print descriptive messages during processing.}
  \item{returnParams}{'logical'. Return recalibrated parameters from crlmm.}
  \item{badSNP}{'numeric'. Threshold to flag as bad SNP (affects
  batchQC)}
  \item{gender}{  integer vector (  male = 1, female =2 ) or missing,
  with same length as filenames.  If missing, the gender is
  predicted.}
  \item{DF}{'integer' with number of degrees of freedom to use with t-distribution.}
}

\details{

	The CRLMM genotype calls and confidence scores are written to
	file using \code{ff} protocols for I/O.  For the most part,
	the calls and confidence scores can be accessed as though the
	data is in memory through the methods \code{snpCall} and
	\code{snpCallProbability}, respectively.

	The genotype calls are stored using an integer representation: 1
	- AA, 2 - AB, 3 - BB.  Similarly, the call probabilities are
	stored using an integer representation to reduce file size using
	the transformation 'round(-1000*log2(1-p))', where p is the
	probability.  The function \code{i2P} can be used to convert the
	integers back to the scale of probabilities.
      }

\value{

Logical.  If the genotyping is completed, the value 'TRUE' is
returned.  Note that \code{assayData} elements 'call' and
'callProbability' are updated on disk. Therefore, the genotypes and
confidence scores can be retrieved using accessors for the
\code{CNSet}  class.

}
\author{
R. Scharpf
}
\seealso{
	\code{\link{crlmm}}, \code{\link{snpCall}}, \code{\link{snpCallProbability}}
}
\examples{
	## See the 'illumina_copynumber' vignette in inst/scripts of
	## the source package
}
\keyword{classif}