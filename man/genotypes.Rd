\name{genotypes}
\alias{genotypes}

\title{
  The possible genotypes for an integer copy number.
}
\description{
  The possible genotypes for an integer copy number (0-4).
}
\usage{
genotypes(copyNumber)
}

\arguments{
  \item{copyNumber}{
    Integer (0-4 allowed).}
}

\value{
  Character vector.
}

\author{
R. Scharpf
}

\examples{

for(i in 0:4) print(genotypes(i))

}

\keyword{manip}