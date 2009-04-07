\name{computeCopynumber}
\alias{computeCopynumber}
\title{Computes copy number}
\description{
  
  A function that transforms the quantile-normalized fluorescence
  intensities of the polymorphic and nonpolymorphic probes to a copy
  number scale.  
  
}
\usage{
computeCopynumber(chrom, A, B, calls, conf, NP, plate, MIN.OBS=1, 
envir, P, DF.PRIOR = 50, CONF.THR = 0.99, bias.adj=FALSE,
priorProb, gender=NULL, SNR, SNRmin, seed=123, cdfName="genomewidesnp6",
verbose=TRUE, ...)
}
\arguments{
  \item{chrom}{Chromosome (an integer).  Use 23 for X and 24 for Y.}
  \item{A}{The A allele intensities from \code{snprma}}
  \item{B}{The B allele intensities from \code{snprma}}
  \item{calls}{The genotype calls from \code{crlmm}}
  \item{conf}{The genotype confidence scores from \code{crlmm}}
  \item{NP}{The quantile normalized intensities of the nonpolymorphic probes}
  \item{plate}{The bach variable.  Should be the same length as the
    number of columns in A}
  \item{MIN.OBS}{Integer: The minimum number of observations in a genotype
    cluster for which a SNP is deemed complete.}
  \item{envir}{An environment to save intermediate objects}
  \item{P}{Mainly for debugging a particular plate/batch.}
  \item{DF.PRIOR}{The degrees of freedom for the prior.  Higher numbers
    will shrink the variance and correlation more.}
  \item{CONF.THR}{A threshold for the genotype confidence scores.
    Genotypes with scores below the threshold are ignored when computing
  SNP-specific  within-genotype estimates of location and scale.}
  \item{bias.adj}{Logical: whether to adjust the location and scale
    parameters to account for biases due to common copy number
    variants.  This is a SNP-specific adjustment.  Parameters for
    background and slope must have already been estimated and available
    from the environment variable.}
  \item{priorProb}{Numerical vector of length 4.  The prior probability
    of each copy number state (0, 1, 2, 3, and 4). The default is a
    uniform prior.  Ignored if bias.adj=FALSE}
  \item{gender}{Gender of subjects. If not specified, we predict the
    gender from the X chromosome.}
  \item{SNR}{Signal to noise ratio from crlmm.}
  \item{SNRmin}{The minimum value for the SNR -- we suggest 5. Samples
    with SNR below SNRmin are excluded.}
  \item{seed}{Seed used for random samples}
  \item{cdfName}{Annotation package }
  \item{verbose}{Logical: verbose output}
  \item{\dots}{Currently ignored}
}

\details{
  
  This function transforms the intensities to a copy number scale.  We
  assume that the median within-SNP intensity across samples is 2.  We
  make no assumption about the chromosomal copy number. This function is
  useful for detecting rare variants (e.g., variants that would not
  affect the SNP-specific quantile-based estimators of location and
  scale).  A correction for more common variants is coming, as well as
  improved estimates of the uncertainty.
  
}

\value{

All objects created by this function are stored in the environment
  passed to this function.  In addition, each of the elements are
  specific to the chromosome(s) specified by the argument \code{chrom}.
  For instance the element \code{A} is the matrix of quantile-normalized
  intensities for the A-allele on chromosome(s) \code{chrom}.  The
  element of this environment are as follows

  \item{A}{Matrix of quantile-normalized intensities for the A-allele}
  \item{B}{Matrix of quantile-normalized intensities for the A-allele}
  \item{CA}{Copy number estimate for the A-allele (x 100)} 
  \item{CB}{Copy number estimate for the B-allele (x 100)}
  \item{calls}{CRLMM genotype calls (AA=1, AB=2, BB=3)}
  \item{chrom}{Integer(s) indicating the chromosome(s)} 
  \item{cnvs}{Names of the nonpolymorphic probes.  These are the
  rownames of \code{NP} and \code{CT}.}
  \item{conf}{CRLMM confidence scores for the genotypes: 'round(-1000*log2(1-p))'}
  \item{corr}{Correlation of the A and B alleles for genotypes AB}
  \item{corrA.BB}{Correlation of A and B alleles for genotypes BB}
  \item{corrB.AA}{Correlation of A and B alleles for genotypes AA}
  \item{CT}{Copy number estimates for nonpolymorphic probe.  See
  \code{cnvs} for the rownames.}
  \item{CT.sds}{Standard deviation estimates for \code{CT}}
  \item{npflags}{Flags for the nonpolymorphic probes.}
  \item{Ns}{The number of observations for each genotype/plate}
  \item{nuA}{Background/cross-hyb for the A allele (plate- and locus-specific)}
  \item{nuB}{Background/cross-hyb for the B allele (plate- and locus-specific)}
  \item{nuT}{Background for the nonpolymorphic probes (plate- and locus-specific)}
  \item{phiA}{Slope for the A allele (plate- and locus-specific)}
  \item{phiB}{Slope for the B allele (plate- and locus-specific)}
  \item{phiT}{Slope for the nonpolymorphic probes (plate- and locus-specific)}
  \item{plate}{Factor indicating batch (same length as number of cel files)}
  \item{sig2A}{Variance estimate for the A-allele signal (plate- and locus-specific)}
  \item{sig2B}{Variance estimate for the B-allele signal (plate- and locus-specific)}
  \item{sig2T}{Variance estimate for the nonpolymorphic signal (plate- and locus-specific)}
  \item{snpflags}{Flags for polymorphic probes}
  \item{snps}{Rownames for \code{A}, \code{B}, \code{CA}, \code{CB}, ...}
  \item{sns}{ Sample names -- the column names for \code{A}, \code{B}, ...}
  \item{steps}{Steps completed. For internal use.}
  \item{tau2A}{Variance estimate for the B-allele background/cross-hyb (plate- and locus-specific)}
  \item{tau2B}{Variance estimate for the B-allele background/cross-hyb (plate- and locus-specific)}
}
\references{Nothing yet.}
\author{Rob Scharpf}
\keyword{manip}

