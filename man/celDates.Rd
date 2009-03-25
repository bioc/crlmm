\name{celDates}
\alias{celDates}
\title{Extract dates from the cel file header}
\description{
  Extract dates from the cel file header.
}
\usage{
celDates(celfiles)
}
\arguments{
  \item{celfiles}{CEL file names.  Must specify the complete path.}
}
\value{
  date-time class POSIXt
}
\author{R. Scharpf}
\seealso{ \cite{\link[affyio]{read.celfile.header}}, \cite{\link{POSIXt}}}
\keyword{manip}

