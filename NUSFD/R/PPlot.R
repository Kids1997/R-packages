#' Pair Plot of a design D
#'
#' @param D a n*d design matrix
#'
#' @return a pair plot.
#' @export
#'
#' @examples
#' D = p_MLHD(30,3,0.8)
#' PPlot(D)
PPlot <- function(D)
{
  d <- ncol(D)
  n <- nrow(D)

  if(n>5000)
  {
    ind <- sample(n,5000)
    Dplot <- D[ind,]
  }else{
    Dplot <- D
  }
  par(mfrow=c(d,d))
  par(mar = c(1,1,1,1))
  for (i in 1:d)
  {
    for (j in 1:d)
    {
      if(i != j)
      {
        plot(Dplot[,i],Dplot[,j],"p")
      }else{
        plot(NULL,xlim=c(0,1),ylim=c(0,1))
        text(x=0.5,y=0.5,paste("x",i,sep=""),cex=2)
      }
    }
  }
}
