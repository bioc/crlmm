\name{plotSNPs}
\alias{plotSNPs}
\alias{plotSamples}

\title{
	Make M vs S plot for SNPs or samples.
}

\description{
	These functions plot the M-values (log-ratios) versus S-values (average intensities)
        for given SNP/(s) or sample/(s) or beanplots for M-values from different samples.
}

\usage{
	plotSNPs(cnSet, row=1, offset=0, xlim=c(9,16), ylim=c(-5,5), verbose=FALSE)
	plotSamples(cnSet, col=1, offset=0, xlim=c(9,16), ylim=c(-5,5), verbose=FALSE, sample=100000, seed=1, type="smoothScatter")
}

\arguments{
  \item{cnSet}{An object of class \code{CNSet}}
  \item{row}{scalar/vector of SNP indexes to plot}
  \item{col}{scalar/vector of sample indexes to plot}
  \item{offset}{numeric, offset to add to intensities in \code{cnSet} 
        before log2-transforming to make log-ratios or average log-intensities}
  \item{xlim}{the x limits of the plot}
  \item{ylim}{the y limits of the plot}
  \item{verbose}{'logical.'  Whether to print descriptive messages during processing}
  \item{sample}{integer indicating the number of SNPs to sample for the plot}
  \item{seed}{integer seed for the random number generator to sample the SNPs}
  \item{type}{character vector specifying the type of sample plot (either 'smoothScatter' or 'beanplot')}
}

\details{
	The \code{plotSNPs} and \code{plotSamples} functions plot the M and S
        values derived from the \code{cnSet} object.
      }

\value{
      One or more M vs S plot for \code{plotSNPs} for a given SNP(/s)
      or either a smoothed scatter plot of M vs S or a beanplot of the M-values 
      for a selected sample(/s) for \code{plotSamples}.
}

\author{
      Matt Ritchie and Cynthia Liu
}

\seealso{
	\code{\link{genotype.Illumina}}
}
\examples{
\dontrun{
   crlmmResult <- genotype.Illumina(sampleSheet=samples[1:10,], path=path,
                                        arrayInfoColNames=list(barcode=NULL,
                                        position="SentrixPosition"),
                                        saveDate=TRUE, cdfName="human370v1c")
   par(mfrow=c(2,2))
   plotSamples(crlmmResult, col=1:4)
   plotSNPs(crlmmResult, row=1:4)
}
}
\keyword{hplot}