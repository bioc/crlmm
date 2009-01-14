\name{list.celfiles}
\alias{list.celfiles}

\title{List CEL files.}
\description{
  Function used to get a list of CEL files.
}
\usage{
list.celfiles(...)
}

\arguments{
  \item{\dots}{Same arguments of \code{\link{list.files}}}
}
\details{
  For the moment, this function returns only uncompressed CEL files (ie,
  no CEL.gz)
}
\value{
  Character vector with filenames.
}

\note{
  Quite often users want to use this function to pass filenames to other
  methods. In this situations, it is safer to use the argument 'full.names=TRUE'.
}
\seealso{\code{\link{list.files}}}
\keyword{IO}
\keyword{utilities}
