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
genotype(filenames, cdfName, batch, mixtureSampleSize = 10^5, eps =
0.1, verbose = TRUE, seed = 1, sns, copynumber =TRUE, probs =
rep(1/3, 3), DF = 6, SNRMin = 5, recallMin = 10, recallRegMin = 1000,
gender = NULL, returnParams = TRUE, badSNP = 0.7)

genotypeLD(filenames, cdfName, batch, mixtureSampleSize = 10^5, eps = 0.1, verbose = TRUE, seed = 1, sns, copynumber = TRUE, probs = rep(1/3, 3), DF = 6, SNRMin = 5, recallMin = 10, recallRegMin = 1000, gender = NULL, returnParams = TRUE, badSNP = 0.7)
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
  \item{copynumber}{ Whether to quantile normalize the nonpolymorphic probes.  If TRUE, the quantile normalized intensities for nonpolymorphic markers are included in the 'A' matrix.}
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
	support by installing and loading the ff package before
	calling the genotype or genotypeLD function.

	Currently, two functions are provided for preprocessing and
	genotyping Affymetrix platforms: genotype and genotypeLD. For
	small datasets, genotype and genotypeLD are identical.  For
	large datasets, genotypeLD provides large data support (via ff)
	and permits the use of clusters or multiple cores (via snow
	package) to speed up genotyping (similar to \code{crlmm2}).
	The \code{genotype} function will be phased out in the future
	and replaced by \code{genotypeLD}.

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

  ## To use less RAM, specify a smaller argument to ocProbesets
  ocProbesets(50e3)
  (cnSet <- genotypeLD(cels, cdfName="genomewidesnp6",
  copynumber=TRUE))

  dim(cnSet)
  table(isSnp(cnSet))

## The above is a trivial example.  Typically you may have a large
##    number of cel files, many of which were processed at different
##    times. For such datasets, it is important to set a batch
##    variable.  If not specified, the scan date of the file is used
##    as the batch variable.  
      batch(cnSet)
      protocolData(cnSet)$ScanDate
}
}
\keyword{ classif }

