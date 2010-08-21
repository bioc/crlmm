\name{genotype}
\alias{genotype}
\alias{genotype2}
\alias{genotypeLD}
\title{
	Preprocessing and genotyping of Affymetrix arrays.
}
\description{
	Preprocessing and genotyping of Affymetrix arrays.	
}
\usage{
genotype(filenames, cdfName, batch, mixtureSampleSize = 10^5, eps =0.1,
         verbose = TRUE, seed = 1, sns, probs = rep(1/3, 3),
         DF = 6, SNRMin = 5, recallMin = 10, recallRegMin = 1000,
         gender = NULL, returnParams = TRUE, badSNP = 0.7)


}
\arguments{
  \item{filenames}{ complete path to CEL files}
  \item{cdfName}{  annotation package  (see also
  \code{validCdfNames})}
  \item{batch}{ batch variable. See details. }
  \item{mixtureSampleSize}{    Sample size to be use when fitting the mixture model.}
 \item{eps}{   Stop criteria.}
  \item{verbose}{  Logical.  Whether to print descriptive messages during processing.}
  \item{seed}{ Seed to be used when sampling. Useful for reproducibility}
  \item{sns}{The sample identifiers.  If missing, the default sample names are \code{basename(filenames)}}
  \item{probs}{'numeric' vector with priors for AA, AB and BB.}
  \item{DF}{'integer' with number of degrees of freedom to use with t-distribution.}
  \item{SNRMin}{'numeric' scalar defining the minimum SNR used to filter
  out samples.}
  \item{recallMin}{Minimum number of samples for recalibration. }
  \item{recallRegMin}{Minimum number of SNP's for regression.}
  \item{gender}{  integer vector (  male = 1, female =2 ) or missing,
  with same length as filenames.  If missing, the gender is predicted.}
  \item{returnParams}{'logical'. Return recalibrated parameters from crlmm.}
  \item{badSNP}{'numeric'. Threshold to flag as bad SNP (affects batchQC)}
}

\details{

	For large datasets it is important to utilize the large data
	support by installing and loading the ff package before calling
	the \code{genotype} function. In previous versions of the
	\code{crlmm} package, we useed different functions for
	genotyping depending on whether the ff package is loaded, namely
	\code{genotype} and \code{genotype2}.  The \code{genotype}
	function now handles both instances.

}

\value{	A \code{SnpSuperSet} instance.}
\references{

  Carvalho B, Bengtsson H, Speed TP, Irizarry RA. Exploration,
  normalization, and genotype calls of high-density oligonucleotide SNP
  array data. Biostatistics. 2007 Apr;8(2):485-99. Epub 2006 Dec
  22. PMID: 17189563.

  Carvalho BS, Louis TA, Irizarry RA. 
  Quantifying uncertainty in genotype calls.
  Bioinformatics. 2010 Jan 15;26(2):242-9.

}
\author{R. Scharpf}
\note{For large datasets, load the 'ff' package prior to genotyping --
this will greatly reduce the RAM required for big jobs.  See
\code{ldPath} and \code{ocSamples}.}

\seealso{
	\code{\link{snprma}}, \code{\link{crlmm}},
	\code{\link[oligoClasses]{ocSamples}},
	\code{\link[oligoClasses]{ldOpts}},
	\code{\link{batch}},
	\code{\link{crlmmCopynumber}}
}
\examples{
if (require(ff) & require(genomewidesnp6Crlmm) & require(hapmapsnp6)){

  path <- system.file("celFiles", package="hapmapsnp6")
  ## the filenames with full path...
  ## very useful when genotyping samples not in the working directory
  cels <- list.celfiles(path, full.names=TRUE)

  ## Note: one would need at least 10 CEL files for copy number estimation
  ## To use less RAM, specify a smaller argument to ocProbesets
  ocProbesets(50e3)
  batch <- as.factor(rep("A", length(cels)))
  (cnSet <- genotype(cels, cdfName="genomewidesnp6"))

  dim(cnSet)
  table(isSnp(cnSet))

## The above is a trivial example.  Typically you may have a large
##    number of cel files, many of which were processed at different
##    times. For such datasets, it is important to specify a batch
##    variable.  
      batch(cnSet)
      protocolData(cnSet)$ScanDate
}
}
\keyword{ classif }

