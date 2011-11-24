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
  \item{returnParams}{'logical'. Return recalibrated parameters from
    crlmm.}
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

	\code{genotype} is essentially a wrapper of the \code{crlmm}
	function for genotyping.  Differences include (1) that the copy
	number probes (if present) are also quantile-normalized and (2)
	the class of object returned by this function, \code{CNSet}, is
	needed for subsequent copy number estimation.  Note that the
	batch variable that must be passed to this function has no
	effect on the normalization or genotyping steps.  Rather,
	\code{batch} is required in order to initialize a \code{CNSet}
	container with the appropriate dimensions.
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
  (cnSet <- genotype(cels, cdfName="genomewidesnp6", batch=batch))

Segment faults that occur with the above step can often be traced to a
corrupt cel file. To check if any of the files are corrupt, try reading
the files in one at a time:

\dontrun{
require(affyio)
for(i in seq_along(cels)) {
	cat("Reading ", i, "\n")
	invisible(read.celfile(cels[i], intensity.means.only=TRUE))
}
}

  ## when gender is not specified (as in the above example), crlmm tries
  ## to predict the gender from SNPs on chromosome X
  cnSet$gender

  ## If gender is known, one should check that the assigned gender is
  ## correct. Alternatively, one can pass gender as an argument to the
  ## genotype function.
  gender <- c("female", "female", "male")
  gender[gender == "female"] <- 2
  gender[gender == "male"] <- 1
  dim(cnSet)
  table(isSnp(cnSet))
  \dontrun{
          ##---------------------------------------------------------------------------
          ##
	  ## Genotype Imputation
          ##
          ##---------------------------------------------------------------------------
     if(require("ff")){
	  path <- system.file("celFiles", package="hapmapsnp6")
	  cels <- list.celfiles(path, full.names=TRUE)
	  gtSet <- genotype(cels, batch=rep(1L,3))

          gt1 <- calls(crlmmOutput)
          gt2 <- calls(gtSet)[isSnp(gtSet), ]
          gt2 <- gt2[match(rownames(gt2), rownames(gt1)), ]
	  stopifnot(all.equal(gt1, gt2))

	  XIndex <- which(chromosome(gtSet)==23 & isSnp(gtSet))
	  YIndex <- which(chromosome(gtSet)==24 & isSnp(gtSet))
	  A.X <- log2(A(gtSet)[XIndex,,drop=FALSE])
	  B.X <- log2(B(gtSet)[XIndex,,drop=FALSE])
	  meds.X <- apply(A.X+B.X, 2, median)/2
	  A.Y <- log2(A(gtSet)[YIndex,,drop=FALSE])
	  B.Y <- log2(B(gtSet)[YIndex,,drop=FALSE])
	  meds.Y <- apply(A.Y+B.Y, 2, median)/2
	  R <- meds.X - meds.Y
	  SNR <- gtSet$SNR[]
	  SNRMin <- 5
	  gender <- kmeans(R, c(min(R[SNR[]>SNRMin]), max(R[SNR[]>SNRMin])))[["cluster"]]
	  plot(seq_along(R), R, pch=21, bg=c("royalblue", "red")[gender],
	       xlim=c(0, 4), xaxt='n', xlab="sample index")
	  legend("topright", fill=c("royalblue", "red"), legend=c("male", "female"))
	  delete(A(gtSet))
	  delete(B(gtSet))
	  delete(calls(gtSet))
	  delete(assayDataElement(gtSet, "snpCallProbability"))
	  ff.filenames <- list.files(".", pattern=".ff")
	  unlink(ff.filename)
     } ## end if(require("ff"))
  }
}
}
\keyword{ classif }

