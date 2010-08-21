\name{linearModelParams}
\alias{lM}
\alias{corr}
\alias{nu}
\alias{phi}
\alias{sigma2}
\alias{tau2}
\title{
	Accessors for linear model parameters.
}
\description{
	
	Accessors for linear model parameters.

}

\usage{
lM(object)
corr(object, allele)
nu(object, allele)
phi(object, allele)
sigma2(object, allele)
tau2(object, allele)
}

\arguments{
  \item{object}{  An object of class \code{CNSetLM}}
  \item{allele}{
   Character string. Valid entries for most accessors are 'A' or 'B'.
   For the \code{corr} accessor, valid entries are 'AA', 'AB', or 'BB'.
   }
}

\details{

\code{lM}:  Extracts entire list of linear model parameters.

\code{corr}: The within-genotype correlation of log2(A) and log2(B) intensities.

\code{nu}: The intercept for the linear model.  The linear model is
fit to the A and B alleles independently.

\code{phi}: The slope for the linear model.  The linear model is fit
independently to the A and B alleles.

\code{sigma2}: For allele A, sigma2 is calculated as the squared MAD
of the log2(intensity) for allele 'A' among subjects with genotype AA.
For allele B, sigma2 is calculated as the squared MAD of the
log2(intensity) for allele 'B' among subjects with genotype BB.
sigma2 can be interpreted as a robust estimate of the signal variance.

\code{tau2}: For allele A, tau2 is calculated as the squared MAD of
the log2(intensity) for allele 'A' among subjects with genotype BB.
For allele B, tau2 is calculated as the squared MAD of the
log2(intensity) for allele 'B' among subjects with genotype AA.  tau2
can be interpeted as a robust estimate of the background variance.

}

\value{

	A matrix or \code{ff} object.

}

\author{
R. Scharpf
}

\seealso{
	\code{\link{CNSetLM-class}}
}
\examples{
## object with ff class
if(require("ff")){
	data(sample.CNSetLMff)
	invisible(open(sample.CNSetLMff))
	class(lM(sample.CNSetLMff))
	params <- lM(sample.CNSetLMff)
	nuA <- nu(sample.CNSetLMff, "A")
	nuB <- nu(sample.CNSetLMff, "B")
	phA <- phi(sample.CNSetLMff, "A")
	phB <- phi(sample.CNSetLMff, "B")
	sig2A <- sigma2(sample.CNSetLMff, "A") 
	sig2B <- sigma2(sample.CNSetLMff, "B")
	tau2A <- tau2(sample.CNSetLMff, "A")
	tau2B <- tau2(sample.CNSetLMff, "B")
	corrAA <- corr(sample.CNSetLMff, "AA")
	corrBB <- corr(sample.CNSetLMff, "BB")
	corrAB <- corr(sample.CNSetLMff, "AB")
	invisible(close(sample.CNSetLMff))
}
## object with matrix class
data(sample.CNSetLM)
class(lM(sample.CNSetLM))
nuA <- nu(sample.CNSetLM, "A")
nuB <- nu(sample.CNSetLM, "B")
phA <- phi(sample.CNSetLM, "A")
phB <- phi(sample.CNSetLM, "B")
sig2A <- sigma2(sample.CNSetLM, "A") 
sig2B <- sigma2(sample.CNSetLM, "B")
tau2A <- tau2(sample.CNSetLM, "A")
tau2B <- tau2(sample.CNSetLM, "B")
corrAA <- corr(sample.CNSetLM, "AA")
corrBB <- corr(sample.CNSetLM, "BB")
corrAB <- corr(sample.CNSetLM, "AB")
}
\keyword{manip}

