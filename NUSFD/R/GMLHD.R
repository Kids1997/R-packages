#' Generating GMLHDs for Gaussian process models
#'
#' @param n a positive integer. Sample size.
#' @param d a positive integer. Dimension.
#' @param X0 a real matrix belongs to [0,1]^d. Initial design.
#'
#' @return an n-point d-dimensional GMLHD.
#' @export
#'
#' @examples
#' ## e.g.1 An example of 40-point 2-dimensional GMLHD using the default initial design
#' n = 40
#' d = 2
#' X = GMLHD(n,d)
#' PPlot(X)
#' ## e.g.2 An example of 20-point 2-dimensional GMLHD using initial design X0
#' n = 20
#' d = 2
#' # Initial design is the following 20-point 2-dimensional maximin LHD
#' X0 = matrix(c(0:19,5,16,11,1,7,17,12,3,9,13,19,0,6,15,10,4,18,2,14,8),ncol=2,byrow=F)/19
#' X = GMLHD(n,d,X0)
#' PPlot(X)
GMLHD <- function(n,d,X0=NULL)
{
  if(is.null(X0))
  {
    X0 = maximinLHS(n,d,maxIter = 1000,eps = 0.01)
  }

  Xm = (1-cos(X0*pi))/2
  return(Xm)
}
