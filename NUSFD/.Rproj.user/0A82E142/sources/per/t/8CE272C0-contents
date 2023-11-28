#' Generating randomized minimax designs for misspecified linear models
#'
#' @param n a positive integer. Sample size.
#' @param d a positive integer. Dimension.
#' @param p a real number belongs to [0,1]. Weight.
#' @param U a real matrix belongs to [0,1]^d. Random or quasi-random numbers.
#'
#' @return an n-point d-dimensional randomized minimax design.
#' @export
#'
#' @import UniDOE
#'
#' @examples
#' ## e.g.1 An example of 20-point 2-dimensional random minimax design with p=0.5
#' n = 20
#' d = 2
#' p = 0.5
#' U = matrix(runif(n*d),ncol=d)
#' Xm = mMDesign(n,d,p,U)
#' PPlot(Xm)
#' ## e.g.2 An example of 40-point 3-dimensional quasi-random minimax design with p=0.8
#' n = 40
#' d = 3
#' p = 0.8
#' Xm = mMDesign(n,d,p)
#' PPlot(Xm)
mMDesign <- function(n,d,p,U=NULL)
{
  if(is.null(U)){
    res = GenUD(n=n,s=d,q=n,crit="CD2",maxiter=1000)$final_design
    U = (2*res-1)/(2*n)
  }
  Xm = Rpinv(p,U)
  return(Xm)
}
