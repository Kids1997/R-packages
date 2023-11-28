#' Minimax desing based subsampling algorithm
#'
#' @param D_full an N*(d+1) real matrix. Full sample, the first d column
#' represents covariates and the last column represents response.
#' @param n a positive integer. Subsample size.
#' @param M2 a positive number. L_infinity bound of departures.
#' @param sigma2 a positive number. Error variance.
#' @param flag a Boolean variable. flag=1: minimax design based subsampling method;
#' flag=0: p-MLHD based subsampling.
#'
#' @return an n-point d-dimensional subsample
#' @export
#'
#' @import nabor
#'
#' @examples
#' library(scales)
#' N = 200
#' d = 2
#' n = 20
#' M2 = 1
#' sigma2 = 20
#' const = matrix(rep(1,N),ncol=1)
#' X0 = matrix(rnorm(N*d),ncol=d)
#' X <-apply(X0,2,rescale)
#' e = matrix(rnorm(N,sd=sqrt(sigma2)),ncol=1)
#' Y = const*2+X%*%matrix(c(3,3),ncol=1)+e
#' D = cbind(X,Y)
#' ## e.g.1 minimax design based subsample
#' S1 = MDBS(D,n,M2,sigma2,1)
#' plot(X, main="MDBS")
#' points(S1[,1:d], col="red", cex=2)
## e.g.2 p-MLHD based subsample
#' S2 = MDBS(D,n,M2,sigma2,0)
#' plot(X, main="MDBS-0")
#' points(S2[,1:d], col="green", cex=2)
MDBS<-function(D_full,n,M2,sigma2,flag)
{
  ### Search the nearest neighbor of the design Xm in the full data D_full
  NNS <- function(D_full,Xm)
  {
    n <- nrow(Xm)
    d <- ncol(Xm)

    NN <- knn(D_full[,1:d],Xm,k=n,searchtype = 2L)$nn.idx
    idx <- c()
    for (i in 1:n)
    {
      if(!(NN[i,1]%in%idx))
      {
        idx<-c(idx,NN[i,1])
      }else{
        idx<-c(idx,NN[i,min(which(!(NN[i,] %in% idx)))])
      }
    }
    S <- matrix(0,nrow = n,ncol = d+1)
    S[1:n,] <- D_full[idx,]
    S<-data.frame(X=S[,1:d],Y=S[,d+1])
    return(S)
  }

  nc <- ncol(D_full)
  d <- nc-1
  r <- M2/(M2+sigma2/n)
  p <- POPT(r,d)

  # MDBS
  if(flag == 1){
    Xm = mMDesign(n,d,p)

  }
  #MDBS-0, p-MLHD based subsampling
  if(flag == 0){
    Xm = p_MLHD(n,d,p)
  }

  S <- NNS(D_full,Xm)
  return(S)
}
