#' Calculating the optimal weight p*(r,d) of the minimax design
#'
#' @param r a real number belongs to [0,1]. Relative importance of bias to variance.
#' @param d a positive integer. Dimension.
#'
#' @return the optimal weight p*(r,d) of the minimax design.
#' @export
#'
#' @examples
#' POPT(0.5,1)
#' POPT(0.8,10)
POPT <- function(r,d){
  f<-function(p)
  {
    y <- (1-r)*(1+d/(3-2*p))+r*(1-p)^2*(1+3*d/(3-2*p)^2)
    return(y)
  }
  if (r<= d/(2*d+9))
  {
    return(0)
  }else
  {
    A = optimize(f,c(0,1),tol=1e-5)
    return(A$minimum)
  }
}
