\name{ellipse-methods}
\docType{methods}
\alias{ellipse-methods}
\alias{ellipse,CNSet-method}
\title{Method for drawing prediction regions for integer copy number}
\description{

	This method draws the prediction regions for integer copy number. 

}

\section{Methods}{
\describe{

\item{\code{signature(x = "CNSet")}}{Parameters from the linear model
are stored in the featureData slot of the CNSet object.  The ellipse
function extracts parameters for the center and covariance of the
bivariate normal prediction regions.}

}
}

\keyword{methods}
