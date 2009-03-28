\name{crlmm-package}
\alias{crlmm-package}
\docType{package}
\title{
Genotype Calling via CRLMM Algorithm
}
\description{
Faster implementation of CRLMM specific to SNP 5.0 and 6.0 arrays.
}
\details{
Index:
\preformatted{
crlmm-package           New implementation of the CRLMM Algorithm.
crlmm                   Genotype SNP 5.0 or 6.0 samples.
calls                   Accessor for genotype calls.
confs                   Accessor for confidences.
}
The 'crlmm' package reimplements the CRLMM algorithm present in the
'oligo' package. This implementation primes for efficient genotyping of
samples on SNP 5.0 and SNP 6.0 Affymetrix arrays.

To use this package, the user must have additional data packages:
'genomewidesnp5Crlmm' - SNP 5.0 arrays
'genomewidesnp6Crlmm' - SNP 6.0 arrays
}
\author{
Rafael A Irizarry
Maintainer: Benilton S Carvalho <bcarvalh@jhsph.edu>
}
\references{
  Carvalho B, Louis TA, Irizarry RA. Describing Uncertainty in
  Genome-wide Genotype Calling. (in prep)
}
\keyword{ package }
