\name{genotypes}
\alias{genotypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  The possible genotypes for an integer copy number.
}
\description{
  The possible genotypes for an integer copy number (0-4).
}
\usage{
genotypes(copyNumber)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{copyNumber}{
    Integer (0-4 allowed).
}
}
\details{

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
x
