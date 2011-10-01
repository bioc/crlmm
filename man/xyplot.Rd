\name{xyplot}
\alias{xyplot}
\alias{xyplot,formula,CNSet-method}
\title{Plot prediction regions and normalized intensities.}
\description{
  Plot prediction regions for integer copy number and normalized intensities.
}
\usage{
xyplot(x, data, ...)
}

\arguments{
  \item{x}{
    A \code{formula}.
  }
  \item{data}{
    A \code{CNSet} object.
  }
  \item{\dots}{
    Additional arguments passed to \code{xyplot} function in lattice.
  }
}

\value{
  A \code{trellis} object.
}

\author{
R. Scharpf
}

\seealso{
  \code{\link[lattice]{xyplot}}, \code{\link{ABpanel}}
}

\examples{
data(sample.CNSet2)
table(batch(sample.CNSet2))
sample.index <- which(batch(sample.CNSet2) == "CUPID")
## A single SNP
pr <- predictionRegion(sample.CNSet2[1, sample.index], copyNumber=0:4)
gt <- as.integer(calls(sample.CNSet2[1, sample.index]))
xyplot(B~A, data=sample.CNSet2[1, sample.index],
       predictRegion=pr,
       panel=ABpanel,
       pch=21,
       fill=c("red", "blue", "green3")[gt],
       xlim=c(6,12), ylim=c(6,12))

## The posterior mean for the above markers
pp <- posteriorProbability(sample.CNSet2[1, sample.index], predictRegion=pr, copyNumber=0:4)
pm <- calculatePosteriorMean(sample.CNSet2[1, sample.index], posteriorProb=pp)

## Muliple SNPs
pr <- predictionRegion(sample.CNSet2[1:10, sample.index], copyNumber=0:4)
gt <- as.integer(calls(sample.CNSet2[1:10, sample.index]))
xyplot(B~A|snpid, data=sample.CNSet2[1:10, sample.index],
       predictRegion=pr,
       panel=ABpanel,
       pch=21,
       fill=c("red", "blue", "green3")[gt],
       xlim=c(6,12), ylim=c(6,12))

## multiple SNPs, prediction regions for 3 batches
\dontrun{
	tab <- table(batch(sample.CNSet2))
	bns <- names(tab)[tab > 50]
	sample.index <- which(batch(sample.CNSet2) %in% bns[1:3])
	pr <- predictionRegion(sample.CNSet2[1:10, sample.index], copyNumber=0:4)
	gt <- as.integer(calls(sample.CNSet2[1:10, sample.index]))
	xyplot(B~A|snpid, data=sample.CNSet2[1:10, sample.index],
	       predictRegion=pr,
	       panel=ABpanel,
	       pch=21,
	       fill=c("red", "blue", "green3")[gt],
	       xlim=c(6,12), ylim=c(6,12))

	## nonpolymorphic markers
	data(sample.CNSet2)
	tab <- table(batch(sample.CNSet2))
	bns <- names(tab)[tab > 50]
	sample.index <- which(batch(sample.CNSet2)%in%bns[1:3])
	np.index <- which(!isSnp(sample.CNSet2))[1:10]
	taus <- tau2(sample.CNSet)[np.index, , , ]
	##trace(predictionRegion, browser, signature=c("CNSet", "integer"))
	pr <- predictionRegion(sample.CNSet2[np.index, sample.index],
			       copyNumber=0:4)
	##trace(posteriorProbability, browser, signature="CNSet")
	pp <- posteriorProbability(sample.CNSet2[np.index, sample.index],
				   predictRegion=pr,
				   copyNumber=0:4)
	pm <- calculatePosteriorMean(sample.CNSet2[np.index, sample.index], posteriorProb=pp)
}
}
\keyword{dplot}
\keyword{hplot}
