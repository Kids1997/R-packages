% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GMLHD.R
\name{GMLHD}
\alias{GMLHD}
\title{Generating GMLHDs for Gaussian process models}
\usage{
GMLHD(n, d, X0 = NULL)
}
\arguments{
\item{n}{a positive integer. Sample size.}

\item{d}{a positive integer. Dimension.}

\item{X0}{a real matrix belongs to [0,1]^d. Initial design.}
}
\value{
an n-point d-dimensional GMLHD.
}
\description{
Generating GMLHDs for Gaussian process models
}
\examples{
## e.g.1 An example of 40-point 2-dimensional GMLHD using the default initial design
n = 40
d = 2
X = GMLHD(n,d)
PPlot(X)
## e.g.2 An example of 20-point 2-dimensional GMLHD using initial design X0
n = 20
d = 2
# Initial design is the following 20-point 2-dimensional maximin LHD
X0 = matrix(c(0:19,5,16,11,1,7,17,12,3,9,13,19,0,6,15,10,4,18,2,14,8),ncol=2,byrow=F)/19
X = GMLHD(n,d,X0)
PPlot(X)
}
