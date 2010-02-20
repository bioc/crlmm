\name{cnOptions}
\alias{cnOptions}
\title{
	Options for copy number estimation
}
\description{
	This function returns all the user-modifiable arguments to the crlmm copy number function.
}
\usage{
    cnOptions(outdir = "./", cdfName, 
              AFile="A",
	      BFile="B",
	      callsFile="genotypes",
	      rgFile="rgFile.rda",
	      cnFile="cnSet",
	      cnrmaFile="cn_rmaResult.rda",
	      snprmaFile="snp_rmaResult.rda",
	      confsFile="confsFile",
	      save.it = TRUE, 
	      save.cnset = TRUE, 
	      load.it = TRUE, 
	      MIN.OBS = 3, MIN.SAMPLES = 10, batch = NULL, DF.PRIOR = 50, 
	      bias.adj = FALSE, prior.prob = rep(1/4, 4), 
	      SNRmin = 4, chromosome = 1:24, seed = 123, verbose = TRUE, 
	      GT.CONF.THR = 0.99, PHI.THR = 2^6, nHOM.THR = 5, 
	      MIN.NU = 2^3, MIN.PHI = 2^3, 
	      THR.NU.PHI = TRUE, thresholdCopynumber = TRUE, unlink = TRUE, 
xo	      use.poe=FALSE, crlmmBatchSize=1100, ...)
}
\arguments{
  \item{outdir}{
	Path to store output from genotyping / copy number algorithms
}
  \item{cdfName}{
        Character string indicating array-type.  See crlmm:::validCdfNames().
}

  \item{AFile}{

  File name for quantile normalized intensities for the A allele is <outdir>/AFile<.rda>

}

  \item{BFile}{

  File name for quantile normalized intensities for the B allele is <outdir>/BFile<.rda>

}

  \item{callsFile}{

  File name for storing biallelic genotypes from the crlmm algorithm is <outdir>/genotypes<.rda>

}

  \item{rgFile}{

  Applicable for the Illumina platform only.

  File name for storing the normalized R and G summaries is <outdir>/rgFile<.rda>

}

  \item{cnFile}{

  Filename for CNSet object is <outdir>/cnSet_<batch>_<CHR><.rda>
       
}

\item{cnrmaFile}{

  Filename for summary statistics from the quantile normalization of
  the nonpolymorphic probes is <outdir>/cnrmaFile<.rda>.  Currently,
  this file contains only the skew and is not used elswhere. This file
  is not meant to be loaded by the user.

}

\item{snprmaFile}{

  Filename for summary statistics from the quantile normalization of
  the A and B alleles is <outdir>/snprmaFile<.rda>.  Currently, this
  file contains skew, signal to noise ratio, mixture model parameter,
  the platform ids for the SNPs, and the sample names.  This file is
  not meant to be loaded by the user.

}

\item{confsFile}{

	The crlmm confidence scores for the biallelic genotype calls
	are stored in <outdir>/consFile<.rda>.  
	
	This file is not meant to be loaded by the user.

}


  \item{save.it}{
  
  Whether to save intermediate files.  Not currently implemented.
  
}

  \item{save.cnset}{

  This argument is ignored for the processing of more than one
  chromosome (e.g., \code{length(chromosome) > 1}).  If \code{TRUE},
  results from the copy number estimation for each chromosome is saved
  to a file.  The format of the saved file is:

  \code{<outdir>/cnSet_<chromosome>.rda}

}

  \item{load.it}{

  If \code{TRUE}, intermediate files are loaded from \code{outdir}.
  See \code{rgFile}, \code{intensityFile}, \code{crlmmFile}.

}

  \item{MIN.OBS}{

  For genotypes with fewer than \code{MIN.OBS}, the within-genotype
  median is imputed from the observed within-genotype mediants at that
  loci.  The parameters used in the regression are estimated from
  polymorphic loci where the genotypes AA, AB, and BB each have a
  frequency greater than \code{MIN.OBS}.
  
}

  \item{MIN.SAMPLES}{

  The minimum number of samples required for each batch. Batches
  wither fewer than \code{MIN.SAMPLES} are skipped.

}

  \item{batch}{

  Character string or factor denoting the batch for each file to be
  processed.  The length of this argument should be the same as the
  number of files. The batch covariate is a surrogate for experimental
  conditions that change over calendar time. Typically, batch can be
  denoted by the 96 well chemistry plate or the month / year.  

}

  \item{DF.PRIOR}{

  The 2 x 2 covariance matrix of the background and signal variances
  is estimated from the data at each locus.  This matrix is then
  smoothed towards a common matrix estimated from all of the loci.
  DF.PRIOR controls the amount of smoothing towards the common matrix,
  with higher values corresponding to greater smoothing.  Currently,
  DF.PRIOR is not estimated from the data.  Future versions may
  estimate DF.PRIOR empirically.

}

  \item{bias.adj}{

  If \code{TRUE}, initial estimates of the linear model are updated
  after excluding samples that have a low posterior probability of
  normal copy number.  Excluding samples that have a low posterior
  probability can be helpful at loci in which a substantial fraction
  of the samples have a copy number alteration.  For additional
  information, see Scharpf et al., 2009.

}

  \item{prior.prob}{

  A numerical vector providing prior probabilities for copy number
  states corresponding to homozygous deletion, hemizygous deletion,
  normal copy number, and amplification, respectively.

}

  \item{SNRmin}{

  The signal to noise ratio (SNR) estimated during the CRLMM
  genotyping is a summary measure of sample quality based on the
  separation of the genotype clusters.  Smaller values of the SNR
  correspond to samples of lower quality.  Samples are excluded from
  the copy number estimation step if SNR values are less than this
  value. A SNR less than 5 for the Affymetrix platform generally
  corresponds to low quality.  For the Illumina platform, we have
  observed samples with poor quality with SNR < 32.  More specific
  recommendations for Illumina are being evaluated.

}

  \item{chromosome}{

  The chromosome(s) to estimate copy number.  Valid entries are 1-23,
  where 23 corresponds to chromosome X.  Quantile normalization and
  genotyping are performed for all SNPs and nonpolymorphic features,
  irrespective of the \code{chromosome} argument.

}

  \item{seed}{

  Seed for random number generation (integer).  Used only for reproducibility.

}

  \item{verbose}{

  If \code{TRUE}, verbose output.

}

  \item{GT.CONF.THR}{

    Confidence threshold for genotype calls (0, 1).  Calls with
    confidence scores below this theshold are not used to estimate the
    within-genotype medians.

}

  \item{PHI.THR}{

    SNPs with slopes (phi values) below this value are flagged.
    Flagged SNPs are not used in a regression to impute background and
    slope coefficients at nonpolymorphic loci.

}

  \item{nHOM.THR}{

  If fewer than \code{nHOM.THR} homozygous genotypes (AA or BB) are
    observed, the SNPs is flagged.  Flagged SNPs are not used in a
    regression to impute background and slope coefficients at
    nonpolymorphic loci.

}

  \item{MIN.NU}{

  Minimum value for the estimate of background in the linear model.
  Negative values are permissible from the estimation, but not
  plausible.  Ignored if \code{THR.NU.PHI} is \code{FALSE}.

}

  \item{MIN.PHI}{

  Minimum value for the estimate of signal (phi) in the linear model.
  Negative values are permissible from the estimation, but not
  plausible.   Ignored if \code{THR.NU.PHI} is \code{FALSE}.

}

  \item{THR.NU.PHI}{

  If \code{THR.NU.PHI} is \code{FALSE}, \code{MIN.NU} and
  \code{MIN.PHI} are ignored.

}

  \item{thresholdCopynumber}{

  If \code{TRUE}, allele-specific number estimates are truncated.
  Values less than 0.05 are assigned the value 0.05; values exceeding
  5 are assigned the value 5.  

}

  \item{unlink}{

  Whether to remove intermediate files storing the quantile normalized
  intensities.  Ignored if \code{save.it} is \code{FALSE}.

}

\item{use.poe}{
	
	Under development.  Whether to use the EM implementation of
	POE (unsupervised) in the R package metaArray to estimate copy
	number at nonpolymorphic loci.  A supervised POE could be used
	in, for instance, a case-control study.  Otherwise, we assume
	the marginal mode is copy number 2 (normal).

}

\item{crlmmBatchSize}{

	'integer'.  The maximum number of samples that can be process
	in a single batch.  When using the ff package, the maximum
	number of elements in a matrix is 2^31.  For Affy 6.0, the
	limit is therefore approximately 1100 samples (1.8e6 * 1100 < 2^31).
}

  \item{\dots}{

  Additional arguments are passed to \code{readIdatFiles} (Illumina
  platforms only).

}
}
\details{

	The minimum required arguments when calling this function are
	\code{cdfName} and \code{batch}.  The user will generally want
	to assign a valid path to \code{outdir} that specifies where
	the intermediate files and processed data are saved.
	
}

\value{

	In general, nothing is returned and all results are saved to
	\code{outdir}.

}


\references{

	RB Scharpf, I Ruczinski, B Carvalho, B Doan, A Chakravarti,
	and R Irizarry (2009), A multilevel model to address batch
	effects in copy number estimation using SNP arrays (Technical
	Report).

}

\author{  R. Scharpf}


\seealso{
}
\examples{
require(hapmapsnp6)
path <- system.file("celFiles", package="hapmapsnp6")
celfiles <- list.celfiles(path)
## the different populations were run in different batches.  For the
##  files in this package, the batch is indicated by the 13th character
##  in the string
batch <- substr(celfiles, 13, 13)
}
\keyword{manip}

