\name{readIDAT}
\alias{readIDAT}
\alias{readIDAT}

\title{Low-level function to read idat files}

\description{
  Reads intensity information for each bead type from a single .idat file
  from Infinium II platforms.  This is a low-level function which \code{readIdatFiles}
  is a wrapper to}

\usage{
readIDAT(idatFile)
}

\arguments{
  \item{idatFile}{character string specifying idat file to be read in}
}

\details{

This function returns a list containing summarised intensities and other
information extracted from a single .idat file.

The \code{readIdatFiles} function makes use of \code{readIDAT} to read
the paired Cy3 and Cy5 idats from one or more arrays.

Thanks to Keith Baggerly who providing this code.
}

\value{
  list which includes item \code{Quants} which contains average intensity (\code{Mean}),
  number of beads (\code{NBeads}) and a measure of variability (\code{SD}) for
  each bead type on the array.
}

\references{
  Ritchie ME, Carvalho BS, Hetrick KN, Tavar\'{e} S, Irizarry RA.
  R/Bioconductor software for Illumina's Infinium whole-genome
  genotyping BeadChips. Bioinformatics. 2009 Oct 1;25(19):2621-3.
}

\author{Keith Baggerly, with modifications by Matt Ritchie}

\examples{
#idatdata = readIDAT("4019585367_A_Grn.idat")
#names(idatdata
#idatdata$Quants[1:5,]
}
\keyword{IO}
