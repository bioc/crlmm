\name{crlmm}
\alias{crlmm}

\title{Genotype oligonucleotide arrays with CRLMM}
\description{
  This is a faster and more efficient implementation of the CRLMM
  algorithm, especially designed for Affymetrix SNP 5 and 6 arrays (to
  be soon extended to other platforms).
}
\usage{
crlmm(filenames, row.names=TRUE, col.names=TRUE,
      probs=c(1/3, 1/3, 1/3), DF=6, SNRMin=5,
      gender=NULL, save.it=FALSE, load.it=FALSE, snprmaFile, callsFile,
      confsFile, AFile, BFile,
      mixtureSampleSize=10^5,
      eps=0.1, verbose=TRUE, cdfName, sns, recallMin=10,
      recallRegMin=1000, returnParams=FALSE, badSNP=0.7, crlmmBatchSize=1000)
}

\arguments{
  \item{filenames}{'character' vector with CEL files to be genotyped.}
  \item{row.names}{'logical'. Use rownames - SNP names?}
  \item{col.names}{'logical'. Use colnames - Sample names?}
  \item{probs}{'numeric' vector with priors for AA, AB and BB.}
  \item{DF}{'integer' with number of degrees of freedom to use with t-distribution.}
  \item{SNRMin}{'numeric' scalar defining the minimum SNR used to filter
  out samples.}
  \item{gender}{'integer' vector, with same length as 'filenames',
    defining sex. (1 - male; 2 - female)}
  \item{save.it}{'logical'. Save preprocessed data?}
  \item{load.it}{'logical'. Load preprocessed data to speed up analysis?}
  \item{snprmaFile}{'character' with filename to be saved/loaded.  Statistics saved in this file include the skew (SKW), signal to noise ratio (SNR), mixture model parameters (mixtureParams), and gender.  }
  \item{callsFile}{'character' string -- filename to be saved/loaded for biallelic genotype calls.}
  \item{confsFile}{'character' string -- filename to be saved/loaded for confidence scores of the biallelic genotype calls.}
  \item{AFile}{'character' string -- filename to be saved/loaded for quantile normalized intensities of the A allele}
  \item{BFile}{'character' string -- filename to be saved/loaded for quantile normalized intensities of the B allele}
  \item{mixtureSampleSize}{Number of SNP's to be used with the mixture model.}
  \item{eps}{Minimum change for mixture model.}
  \item{verbose}{'logical'.}
  \item{cdfName}{'character' defining the CDF name to use
    ('GenomeWideSnp5', 'GenomeWideSnp6')}
  \item{sns}{'character' vector with sample names to be used.}
  \item{recallMin}{Minimum number of samples for recalibration.}
  \item{recallRegMin}{Minimum number of SNP's for regression.}
  \item{returnParams}{'logical'. Return recalibrated parameters.}
  \item{badSNP}{'numeric'. Threshold to flag as bad SNP (affects batchQC)}
  \item{crlmmBatchSize}{'integer', indicating the number of samples to be genotyped in a batch by crlmm.  The purpose of this argument is to keep the RAM at manageable levels.}
}
\value{
	
	Nothing is returned.  The following files are created as a result of this function:

	- AFile
	- BFile
	- CAFile
	- CBFile
	- callsFile
	- confsFile
	- snprmaFile

	See the \code{cnOptions} function for the default filenames.
	
}

\references{
  Carvalho B, Bengtsson H, Speed TP, Irizarry RA. Exploration,
  normalization, and genotype calls of high-density oligonucleotide SNP
  array data. Biostatistics. 2007 Apr;8(2):485-99. Epub 2006 Dec
  22. PMID: 17189563.

  Carvalho BS, Louis TA, Irizarry RA. 
  Quantifying uncertainty in genotype calls.
  Bioinformatics. 2010 Jan 15;26(2):242-9.
}
\examples{
## this can be slow
if (require(genomewidesnp5Crlmm) & require(hapmapsnp5)){
   ## For small datasets:
  path <- system.file("celFiles", package="hapmapsnp5")
  ## the filenames with full path...
  ## very useful when genotyping samples not in the working directory
  cels <- list.celfiles(path, full.names=TRUE)
  obj <- constructClass("genomewidesnp5")
  class(obj)
  str(crlmmOptions(obj))  ## default settings can be modified by changing the elements of this list
  callSet <- crlmm(obj, cels)
  ##Accessors
  A(callSet)[1:5,]
  B(callSet)[1:5,]
  calls(callSet)[1:5,]
  confs(callSet)[1:5,]
  callProbability(callSet)[1:5,]
  ## For big datasets
  if(require(ff)){
	  options(fffinalizer="close", fftempdir=".", fffinonexit="TRUE")
	  obj <- constructClass("genomewidesnp5")
	  class(obj)
	  callSet <- crlmm(obj, cels)
	  A(callSet)
	  B(callSet)
	  calls(callSet)
	  confs(callSet)
  }
}
}
\seealso{\code{\link{cnOptions}}, \code{\link{preprocessAndGenotype}}}
\keyword{classif}

