\name{crlmmCopynumber}
\alias{crlmmCopynumber}
\alias{crlmmCopynumber2}
\alias{crlmmCopynumberLD}
\title{Locus- and allele-specific estimation of copy number}
\description{
  Locus- and allele-specific estimation of copy number.
}
\usage{
crlmmCopynumber(object, MIN.SAMPLES=10, SNRMin = 5, MIN.OBS = 1, 
	        DF.PRIOR = 50, bias.adj = FALSE,
                prior.prob = rep(1/4, 4), seed = 1, verbose = TRUE,
                GT.CONF.THR = 0.95, MIN.NU = 2^3, MIN.PHI = 2^3,
                THR.NU.PHI = TRUE, type=c("SNP", "NP", "X.SNP", "X.NP"))
}

\arguments{
  
  \item{object}{object of class \code{CNSet}.}
  
 \item{MIN.SAMPLES}{ 'Integer'.  The minimum number of samples in a
  batch.  Bathes with fewer than MIN.SAMPLES are skipped.  Therefore,
  samples in batches with fewer than MIN.SAMPLES have NA's for the
  allele-specific copy number and NA's for the linear model
  parameters.

}

  \item{SNRMin}{ Samples with low signal to noise ratios are  excluded.  }

  \item{MIN.OBS}{ 

  For a SNP with with fewer than \code{MIN.OBS} of a genotype in a given
  batch, the within-genotype median is imputed.  The imputation is based
  on a regression using SNPs for which all three biallelic genotypes are
  observed.  For example, assume at at a given SNP genotypes AA and AB
  were observed and BB is an unobserved genotype.  For SNPs in which all
  3 genotypes were observed, we fit the model E(mean_BB) = beta0 +
  beta1*mean_AA + beta2*mean_AB, obtaining estimates; of beta0, beta1,
  and beta2.  The imputed mean at the SNP with unobserved BB is then
  beta0hat + beta1hat * mean_AA of beta2hat * mean_AB.

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

  \code{bias.adj} is currently ignored (as well as the prior.prob
  argument).  We plan to add this feature back to the crlmm package in
  the near future. This feature, when \code{TRUE}, updated initial
  estimates from the linear model after excluding samples with a low
  posterior probability of normal copy number.  Excluding samples that
  have a low posterior probability can be helpful at loci in which a
  substantial fraction of the samples have a copy number alteration.
  For additional information, see Scharpf et al., 2010.

}
  \item{prior.prob}{

    This argument is currently ignored.  A numerical vector providing
  prior probabilities for copy number states corresponding to homozygous
  deletion, hemizygous deletion, normal copy number, and amplification,
  respectively.

}
  \item{seed}{ Seed for random number generation.}

  \item{verbose}{ Logical. }

  \item{GT.CONF.THR}{ 

    Confidence threshold for genotype calls (0, 1).  Calls with
    confidence scores below this theshold are not used to estimate the
    within-genotype medians. See Carvalho et al., 2007 for information
    regarding confidence scores of biallelic genotypes.

}


  \item{MIN.NU}{ numeric. Minimum value for background
    intensity. Ignored if \code{THR.NU.PHI} is \code{FALSE}. }
  
  \item{MIN.PHI}{numeric. Minimum value for slope. Ignored if
    \code{THR.NU.PHI} is \code{FALSE}.}
  
  \item{THR.NU.PHI}{  If \code{THR.NU.PHI} is \code{FALSE},
    \code{MIN.NU} and \code{MIN.PHI} are ignored. When TRUE, background
    (nu) and slope (phi) coefficients below MIN.NU and MIN.PHI are set
    to MIN.NU and MIN.PHI, respectively.}

  \item{type}{ Character string vector that must be one or more of
    "SNP", "NP", "X.SNP", or "X.NP". Type refers to a set of
    markers. See details below}
    
 }
  
 \references{
   
  Carvalho B, Bengtsson H, Speed TP, Irizarry RA. Exploration,
  normalization, and genotype calls of high-density oligonucleotide SNP
  array data. Biostatistics. 2007 Apr;8(2):485-99. Epub 2006 Dec
  22. PMID: 17189563.

  Carvalho BS, Louis TA, Irizarry RA. 
  Quantifying uncertainty in genotype calls.
  Bioinformatics. 2010 Jan 15;26(2):242-9.

  Scharpf RB, Ruczinski I, Carvalho B, Doan B, Chakravarti A, and
  Irizarry RA, Biostatistics.  Biostatistics, Epub July 2010.
  
}

\details{

	The function crlmmCopynumber uses matrices instead of ff objects
	if the ff library is not loaded. When the ff package is loaded,
	large data support is enabled.  Normalized intensities
	(\code{alleleA} and \code{alleleB}), genotype calls and
	confidence scores (\code{snpCall} and \code{snpCallProbability})
	are stored in \code{assayData} slot.  Summary statistics for
	each batch, including the linear model paramters for copy
	number, are stored in the \code{batchStatistics} slot.  Both the
	\code{assayData} and \code{batchStatistics} slot are of class
	\code{AssayData} with elements that are ff objects (if ff
	package is loaded) or matrices.

	The functions \code{crlmmCopynumberLD} and
	\code{crlmmCopynumber2} have been deprecated.

	The argument \code{type} can be used to specify a subset of
	markers for which the copy number estimation algorithm is run.
	One or more of the following possible entries are valid: 'SNP',
	'NP', 'X.SNP', and 'X.NP'.

	'SNP' referers to autosomal SNPs.

	'NP' refers to autosomal nonpolymorphic markers.

	'X.SNP' refers to SNPs on chromosome X.

	'X.NP' refers to autosomes on chromosome X.

}
\author{R. Scharpf}
\keyword{manip}
