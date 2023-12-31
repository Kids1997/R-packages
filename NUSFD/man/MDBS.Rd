% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MDBS.R
\name{MDBS}
\alias{MDBS}
\title{Minimax desing based subsampling algorithm}
\usage{
MDBS(D_full, n, M2, sigma2, flag)
}
\arguments{
\item{D_full}{an N*(d+1) real matrix. Full sample, the first d column
represents covariates and the last column represents response.}

\item{n}{a positive integer. Subsample size.}

\item{M2}{a positive number. L_infinity bound of departures.}

\item{sigma2}{a positive number. Error variance.}

\item{flag}{a Boolean variable. flag=1: minimax design based subsampling method;
flag=0: p-MLHD based subsampling.}
}
\value{
an n-point d-dimensional subsample
}
\description{
Minimax desing based subsampling algorithm
}
\examples{
library(scales)
N = 200
d = 2
n = 20
M2 = 1
sigma2 = 20
const = matrix(rep(1,N),ncol=1)
X0 = matrix(rnorm(N*d),ncol=d)
X <-apply(X0,2,rescale)
e = matrix(rnorm(N,sd=sqrt(sigma2)),ncol=1)
Y = const*2+X\%*\%matrix(c(3,3),ncol=1)+e
D = cbind(X,Y)
## e.g.1 minimax design based subsample
S1 = MDBS(D,n,M2,sigma2,1)
plot(X, main="MDBS")
points(S1[,1:d], col="red", cex=2)
S2 = MDBS(D,n,M2,sigma2,0)
plot(X, main="MDBS-0")
points(S2[,1:d], col="green", cex=2)
}
