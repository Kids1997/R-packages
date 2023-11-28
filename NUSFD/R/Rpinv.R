#' Minimax transformation from [0,1]^d to [0,1]^d.
#'
#' @param p a real number belongs to [0,1].
#' @param U a real matrix belongs to [0,1]^(n*d).
#'
#' @return a n*d real matrix that represents U after minimax transformation.
#' @export
#'
#' @examples
#' ## An example of 30-point 2-dimensional random minimax design with p=0.5
#' n = 30
#' d = 2
#' p = 0.5
#' U = matrix(runif(n*d),ncol=d)
#' D = Rpinv(p,U)
#' PPlot(D)
Rpinv <- function(p,U){
  Qpinv <- function(p,u){
    n <- length(u)
    a <- matrix(rep(0,n),nrow = n)
    if(p==0)
    {
      a[u>=0.5] = 1
    }else{
      a = (1/p)*(u-(1-p)/2)
    }
    a[a>1] = 1
    a[a<0] = 0
    return(a)
  }
  n <- nrow(U)
  d <- ncol(U)
  X <- matrix(0,nrow = n,ncol = d)
  X[,1] = Qpinv(p,U[,1])
  if(d>1){
    for(j in 2:d)
    {
      cond = matrix(rep(FALSE,n),nrow = n)
      for(i in 1:(j-1))
      {
        cond <- cond | (X[,i]*(X[,i]-1)==0)
      }
      X[cond,j] = Qpinv(0,U[cond,j])
      X[!cond,j] = Qpinv(1,U[!cond,j])
    }
  }
  return(X)
}
