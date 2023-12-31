% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mMDesign.R
\name{mMDesign}
\alias{mMDesign}
\title{Generating randomized minimax designs for misspecified linear models}
\usage{
mMDesign(n, d, p, U = NULL)
}
\arguments{
\item{n}{a positive integer. Sample size.}

\item{d}{a positive integer. Dimension.}

\item{p}{a real number belongs to [0,1]. Weight.}

\item{U}{a real matrix belongs to [0,1]^d. Random or quasi-random numbers.}
}
\value{
an n-point d-dimensional randomized minimax design.
}
\description{
Generating randomized minimax designs for misspecified linear models
}
\examples{
## e.g.1 An example of 20-point 2-dimensional random minimax design with p=0.5
n = 20
d = 2
p = 0.5
U = matrix(runif(n*d),ncol=d)
Xm = mMDesign(n,d,p,U)
PPlot(Xm)
## e.g.2 An example of 40-point 3-dimensional quasi-random minimax design with p=0.8
n = 40
d = 3
p = 0.8
Xm = mMDesign(n,d,p)
PPlot(Xm)
}
