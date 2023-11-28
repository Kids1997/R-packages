#' Generating p-MLHDs for Gaussian process models
#'
#' @param n a positive integer. Sample size.
#' @param d a positive integer. Dimension.
#' @param p a real number belongs to [0,1]. Weight.
#' @param X0 a real matrix belongs to [0,1]^d. Initial design.
#'
#' @return an n-point d-dimensional p-MLHD.
#' @export
#'
#' @import lhs
#'
#' @examples
#' ## e.g.1 An example of 40-point 2-dimensional 0.8-MLHD using the default initial design
#' n = 40
#' d = 2
#' p = 0.8
#' X = p_MLHD(n,d,p)
#' PPlot(X)
#' ## e.g.2 An example of 20-point 2-dimensional 0.6-MLHD using initial design X0
#' n = 20
#' d = 2
#' p = 0.6
#' # Initial design is the following 20-point 2-dimensional maximin LHD
#' X0 = matrix(c(0:19,5,16,11,1,7,17,12,3,9,13,19,0,6,15,10,4,18,2,14,8),ncol=2,byrow=F)/19
#' X = p_MLHD(n,d,p,X0)
#' PPlot(X)
p_MLHD <- function(n,d,p,X0=NULL)
{
  if(is.null(X0))
  {
    X0 = maximinLHS(n,d,maxIter = 1000,eps = 0.01)
  }

  Xm = X0
  for(j in 1:d)
  {
    Xm[,j] = Rpinv(p,as.matrix(X0[,j],ncol=1))
  }
  return(Xm)
}
