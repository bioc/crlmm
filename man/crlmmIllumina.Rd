\name{crlmmIllumina}
\alias{crlmmIllumina}
\title{Genotype Illumina Infinium II BeadChip data with CRLMM}
\description{
  This implementation of the CRLMM is especially designed for
  data from Illumina Infinium II BeadChips.
}
\usage{

crlmmIllumina(RG, XY, stripNorm=TRUE, useTarget=TRUE,
      row.names=TRUE, col.names=TRUE,
      probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5,
      gender=NULL, seed=1, mixtureSampleSize=10^5,
      eps=0.1, verbose=TRUE, cdfName, sns, recallMin=10,
      recallRegMin=1000, returnParams=FALSE, badSNP=0.7)
}

\arguments{
  \item{RG}{\code{NChannelSet} containing R and G bead intensities}
  \item{XY}{\code{NChannelSet} containing X and Y bead intensities}
  \item{stripNorm}{'logical'.  Should the data be strip-level normalized?}
  \item{useTarget}{'logical' (only used when \code{stripNorm=TRUE}).
    Should the reference HapMap intensities be used in strip-level normalization?}
  \item{row.names}{'logical'. Use rownames - SNP names?}
  \item{col.names}{'logical'. Use colnames - Sample names?}
  \item{probs}{'numeric' vector with priors for AA, AB and BB.}
  \item{DF}{'integer' with number of degrees of freedom to use with t-distribution.}
  \item{SNRMin}{'numeric' scalar defining the minimum SNR used to filter
  out samples.}
  \item{gender}{'integer' vector, with same length as 'filenames',
    defining sex. (1 - male; 2 - female)}
  \item{seed}{'integer' scalar for random number generator (used to
    sample \code{mixtureSampleSize} SNPs for mixture model.}
  \item{mixtureSampleSize}{'integer'. The number of SNP's to be used
    when fitting the mixture model.}
  \item{eps}{Minimum change for mixture model.}
  \item{verbose}{'logical'.}
  \item{cdfName}{'character' defining the chip annotation (manifest) to use
    ('human370v1c', human550v3b', 'human650v3a', 'human1mv1c',
    'human370quadv3c', 'human610quadv1b', 'human660quadv1a' 'human1mduov3b')}
  \item{sns}{'character' vector with sample names to be used.}
  \item{recallMin}{'integer'. Minimum number of samples for recalibration.}
  \item{recallRegMin}{'integer'. Minimum number of SNP's for regression.}
  \item{returnParams}{'logical'. Return recalibrated parameters.}
  \item{badSNP}{'numeric'. Threshold to flag as bad SNP (affects batchQC)}
}
\value{
  \item{A}{normalized A allele intensity}
  \item{B}{normalized B allele intensity}
  \item{calls}{Genotype calls (1 - AA, 2 - AB, 3 - BB)}
  \item{confs}{Confidence scores 'round(-1000*log2(1-p))'}
  \item{SNPQC}{SNP Quality Scores}
  \item{batchQC}{Batch Quality Score}
  \item{params}{Recalibrated parameters}
}

\details{
  Note: The user should specify either the \code{RG} or \code{XY}
  intensities, not both.
}

\references{
  Carvalho B, Bengtsson H, Speed TP, Irizarry RA. Exploration,
  normalization, and genotype calls of high-density oligonucleotide SNP
  array data. Biostatistics. 2007 Apr;8(2):485-99. Epub 2006 Dec
  22. PMID: 17189563.

  Carvalho B, Louis TA, Irizarry RA. Describing Uncertainty in
  Genome-wide Genotype Calling. (in prep)
}

\author{Matt Ritchie}

\examples{
## crlmmOut = crlmmIllumina(RG)
}
\keyword{classif}
