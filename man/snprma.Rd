\name{snprma}
\Rdversion{1.1}
\alias{snprma}

\title{
  Preprocessing tool for SNP arrays.
}
\description{
  SNPRMA will preprocess SNP chips. The preprocessing consists of
  quantile normalization to a known target distribution and
  summarization to the SNP-Allele level.
}
\usage{
snprma(filenames, mixtureSampleSize = 10^5, fitMixture = FALSE, eps = 0.1, verbose = TRUE, seed = 1, cdfName, sns)
}
\arguments{
  \item{filenames}{
    'character' vector with file names.
}
  \item{mixtureSampleSize}{
    Sample size to be use when fitting the mixture model.
}
  \item{fitMixture}{
    'logical'. Fit the mixture model?
}
  \item{eps}{
    Stop criteria.
}
  \item{verbose}{
    'logical'.
}
  \item{seed}{
    Seed to be used when sampling.
}
  \item{cdfName}{
    cdfName: 'GenomeWideSnp_5', 'GenomeWideSnp_6'
}
  \item{sns}{
    Sample names.
}
}
\value{
  \item{A}{Summarized intensities for Allele A}
  \item{B}{Summarized intensities for Allele B}
  \item{sns}{Sample names}
  \item{gns}{SNP names}
  \item{SNR}{Signal-to-noise ratio}
  \item{SKW}{Skewness}
  \item{mixtureParams}{Parameters from mixture model}
  \item{cdfName}{Name of the CDF}
}
\examples{
if (require(genomewidesnp5Crlmm) & require(hapmapsnp5)){
  path <- system.file("celFiles", package="hapmapsnp5")

  ## the filenames with full path...
  ## very useful when genotyping samples not in the working directory
  cels <- list.celfiles(path, full.names=TRUE)
  snprmaOutput <- snprma(cels)
  snprmaOutput[["A"]][1:10,]
  snprmaOutput[["B"]][1:10,]
}
}
\keyword{manip}
\keyword{classif}