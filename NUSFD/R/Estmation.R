#' Estimating the error variance sigma^2 and the L_infinity bound M^2
#'
#' @param D an N*(d+1) data frame. Full sample, the first d column
#' represents covariates and the last column represents response.
#' @param ntimes a positive integer. Number of repetitions.
#'
#' @return a ntimes*2 matrix. Each row represents an estimate of sigma^2 and M^2.
#' @export
#'
#' @import mlegp
#' @import lhs
#' @import foreach
#' @import iterators
#' @import parallel
#' @import doParallel
#'
#' @examples
#' N = 1000
#' d = 2
#' n = 20
#' M2 = 10
#' sigma2 = 10
#' const = matrix(rep(1,N),ncol=1)
#' X = matrix(rnorm(N*d),ncol=d)
#' e = matrix(rnorm(N,sd=sqrt(sigma2)),ncol=1)
#' Y = const*2+X%*%matrix(c(3,3),ncol=1)+e
#' D = data.frame(X=X,Y=Y)
#' Estmation(D,10)
Estmation <- function(D,ntimes)
{
  remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.025, .975), na.rm = na.rm, ...)
    H <- 1.5* IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    return(na.omit(y))
  }

  n <- nrow(D)
  if(n>1e4)
  {
    ind <- sample(n,1e4)
    D <- D[ind,]
    n <- nrow(D)
  }else{}
  nc <- ncol(D)
  d <- nc-1
  x <- as.matrix(D[,1:d])
  # linear fitting
  linfit <- lm(Y~.,data=D)
  h <- linfit$residuals
  # GP fitting
  cores <- detectCores(logical=F)
  cl <- makeCluster(cores)
  registerDoParallel(cl, cores=cores)

  n0 <- 10*d
  mm <- foreach(i=1:ntimes, .combine='rbind',.packages = c("mlegp","lhs")) %dopar%
    {
      ind <- sample(n,n0)
      xsample <- x[ind,]
      hsample <- h[ind]

      hfit <- mlegp(xsample,hsample,nugget = 10,nugget.known=0)

      E <- h-predict(hfit,x)
      sigma2 <- sum(E^2)/n

      xx <- randomLHS(n,d)
      yhat <- predict(hfit,xx)
      Mhat <- max(abs(remove_outliers(yhat)))
      cc <- c(sigma2,Mhat^2)
      return(cc)
    }
  return(mm)
}
