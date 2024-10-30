## Auxiliary functions  ====
CRE<-function(n,n1){
  return(sample(n,n1))
}
ReM<-function(X,n1,alpha){
  n = nrow(X)
  p = ncol(X)
  cci = solve(cov(X))
  pw = n1/n
  imb = Inf
  while (imb>qchisq(alpha,p)) {
    idx = sample(n,n1)
    dmX = apply(X[idx,,drop=F],2,mean)-apply(X[-idx,,drop=F],2,mean)
    imb = n*pw*(1-pw)*t(dmX)%*%cci%*%dmX
  }
  return(idx)
}
BODE<-function(X,n1,C){
  M = 500
  n = nrow(X)
  s = 0
  valstar = Inf
  while (s<=M) {
    idx = sample(n,n1)
    w1 = matrix(rep(1,n),ncol=1)
    w1[idx,]  = (1/n1-1/n)*matrix(1,nrow = n1)
    w1[-idx,] = (-1/n)*matrix(1,nrow = n-n1)
    w0 = matrix(rep(1,n),ncol=1)
    w0[idx]  = (1/n)*matrix(1,nrow = n1)
    w0[-idx] = (-1/n1+1/n)*matrix(1,nrow = n-n1)
    w = rbind(w1,w0)
    valnew= t(w)%*%C%*%w
    if(valnew<valstar){
      valstar = valnew
      idxstar = idx
    }
    s = s+1
  }
  return(idxstar)
}
MODE<-function(X,n1,flag,DK=NULL){
  n = nrow(X)
  p = ncol(X)
  if(flag==1){
    r = ceiling(n/n1)
    idxstar = twin(X,r)
  }else if(flag==2){
    idxstar = PDDS(X,n1,partperd= c(5,2))
  }else{
    M = 500
    s = 0
    valstar = Inf
    while (s<=M) {
      print(s)
      idx = sample(n,n1)
      X1 = X[idx,]
      X0 = X[-idx,]
      valnew = DK(X1,X)
      if(valnew<valstar){
        valstar = valnew
        idxstar = idx
      }
      s = s+1
    }
  }
  return(idxstar)
}

Calculating<-function(U,V){
  n = nrow(U)
  C = matrix(0,nrow = n,ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      C[i,j] = 10*exp(-sum((U[i,]-V[j,])^2)/10)
    }
  }
  C = C + t(C)+ 10 *diag(n)
  return(C)
}
CalC<-function(X){
  KED<-function(x,y){
    res = -(sum((x-y)^2))^{1/2}
    return(res)
  }
  Xs = scale(X)
  n = nrow(Xs)
  C = 0
  for (i in 1:n) {
    for (j in 1:n) {
      C = C+KED(Xs[i,],Xs[j,])
    }
  }
  return(C/(n^2))
}
DK<-function(X1,X){
  GKer<-function(x,y){
    return(exp(-sum((x-y)^2)))
  }
  n1 = nrow(X1)
  n = nrow(X)
  val = 0
  for (i in 1:n1) {
    for (j in 1:n) {
      val = val - 2*GKer(X1[i,],X[j,])/(n*n1)
    }
    for (j in 1:n1) {
      val = val + GKer(X1[i,],X1[j,])/(n1*n1)
    }
  }
  return(val)
}
POC<-function(X,flag1,flag2){
  n = nrow(X)
  p = ncol(X)
  v = matrix(rnorm(n),ncol = 1)
  u = matrix(rnorm(n),ncol = 1)
  mx = sqrt(p)*apply(X,1,mean)
  if(flag1==0){
    f1X = exp(mx)
    f0X = 1+mx
  }else{
    f1X = exp(mx)+exp(mx/2)
    f0X = exp(mx)-exp(mx/2)
  }
  if(flag2==0){
    c = 0
  }else{
    c = 1
  }
  Y1 = f1X+c*u
  Y0 = f0X+c*v
  return(cbind(Y1,Y0))
}
Mdist<-function(X1,X0,Cinv){
  n1 = nrow(X1)
  n0 = nrow(X0)
  n = n1+n0
  pw = n1/n
  data = apply(X1,2,mean)-apply(X0,2,mean)
  dmX = matrix(data,ncol=1)
  Md = n*pw*(1-pw)*t(dmX)%*%Cinv%*%dmX
  return(Md)
}
Edist<-function(X1,X0,Const){
  X = rbind(X1,X0)
  n1 = nrow(X1)
  n0 = nrow(X0)
  Ed = energy(X,X1)
  res = ((n1+n0)/n0)*(Ed+Const)
  return(res)
}

## R packages====

library(foreach)
library(iterators)
library(doParallel)
library(doSNOW)
library(ggplot2)
library(gridExtra)
library(twinning)
library(ggpubr)

cores <- detectCores(logical=F)
cl <- makeCluster(cores)
registerDoSNOW(cl)

# Before running the following codes, 
# you need to specify a path to store the calculation results

path = "*******"

## Covariate imbalance measures ====
N = c(1000,2000)
P = c(5,10)
nloop = 1000
for (kk in 1:2) {
  n = N[kk]
  p = P[kk]
  n1 = n/2
  set.seed(88888888)
  
  for(f in 1:2){
    if(f==1){
      X = matrix(-2+4*runif(n*p),ncol=p)
    }else{
      X = matrix(rnorm(n*p,0,sqrt(3)),ncol=p)
    }
    X = scale(X)
    
    Z1 = cbind(X,p*matrix(rep(1,n),ncol = 1))
    Z0 = cbind(X,p*matrix(rep(0,n),ncol = 1))
    C11 = Calculating(Z1,Z1)
    C10 = Calculating(Z1,Z0)
    C00 = Calculating(Z0,Z0)
    C = rbind(cbind(C11,C10),cbind(t(C10),C00))
    
    Cinv = solve(cov(X))
    Const = CalC(X)
    
    print(paste("***n =",n,"p =", p,"case",f,"***"))
    pb <- txtProgressBar(max = nloop, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    allrep <- foreach(i=1:nloop, .combine='rbind', 
                      .packages = c("PDDS","twinning"),
                      .options.snow = opts) %dopar%
      {
        idx1 = CRE(n,n1)
        idx2 = ReM(X,n1,0.2)
        idx3 = BODE(X,n1,C)
        idx4 = MODE(X,n1,1)
        
        md1 = Mdist(X[idx1,],X[-idx1,],Cinv)
        md2 = Mdist(X[idx2,],X[-idx2,],Cinv)
        md3 = Mdist(X[idx3,],X[-idx3,],Cinv)
        md4 = Mdist(X[idx4,],X[-idx4,],Cinv)
        
        
        ed1 = Edist(X[idx1,],X[-idx1,],Const)
        ed2 = Edist(X[idx2,],X[-idx2,],Const)
        ed3 = Edist(X[idx3,],X[-idx3,],Const)
        ed4 = Edist(X[idx4,],X[-idx4,],Const)
        
        mdall = c(md1,md2,md3,md4)
        edall = c(ed1,ed2,ed3,ed4)
        return(c(edall,mdall))
      }
    close(pb)
    FILE = paste(path,"allimb-",f,"-",kk,".txt")
    write.table (allrep, file =FILE, sep =" ", row.names =F, col.names =F)
  }
}
# Visualization of imbalance measures
LABs = matrix(c(1,1,2,1,1,2,2,2),ncol=2)
ED = c()
MD = c()
for (lab in 1:4) {
  FILE = paste(path,"allimb-",LABs[lab,1],"-", LABs[lab,2],".txt")
  res = as.matrix(read.table(FILE))
  ED = c(ED,as.vector(res[,1:4]))
  MD = c(MD,as.vector(res[,5:8]))
}

M4 <- c("CRE","ReM","BODE","MODE")
SD <- c("n=1000,p=5","n=2000,p=10")
Dist <-c("unifrom-X","nrom-X")
methods = rep(rep(M4,each=nloop),length(SD)*length(Dist))
cases = rep(SD,each =nloop*length(M4)*length(Dist))
dists = rep(rep(Dist, each = nloop*length(M4)),length(Dist))
DF = data.frame(
  method = factor(methods,levels = M4),
  case = factor(cases,levels= SD),
  dist = factor(dists,levels = Dist),
  ed = ED,
  md = MD)

# Mahalanobis distance
picM = list(NULL)
length(picM) = 2
for (i in 1:2) {
  picM[[i]]<-ggplot(subset(DF,case==SD[i]),aes(x=dist,y=md,fill=method))+
    geom_boxplot(aes(color=method))+
    labs(x="",y="Mahalanobis distance",title = SD[i])+
    theme_bw()+ 
    theme(legend.title=element_blank(),
          legend.position = "top",#c(0.64,0.95), #
          legend.direction = "horizontal",
          plot.title = element_text(hjust = 0.5))+
    guides(fill = guide_legend( nrow = 1, byrow = TRUE))
}
# Energy distance
picE = list(NULL)
length(picE) = 2
for (i in 1:2) {
  picE[[i]]<-ggplot(subset(DF,case==SD[i]),aes(x=dist,y=ed,fill=method))+
    geom_boxplot(aes(color=method))+
    labs(x="",y="Energy distance",title = SD[i])+
    theme_bw()+ 
    theme(legend.title=element_blank(),
          legend.position = "top",#c(0.64,0.95), #
          legend.direction = "horizontal",
          plot.title = element_text(hjust = 0.5))+
    guides(fill = guide_legend( nrow = 1, byrow = TRUE))
}

ggarrange(picM[[1]],picM[[2]],picE[[1]],picE[[2]], common.legend = TRUE, legend="top")


## Estimation precision=====  
N = c(250,1000,500,2000)
P = c(5,5,10,10)
nloop = 1000
FLAG = matrix(c(0,0,0,1,1,0,1,1),ncol = 2,byrow = T)

# Calculating mse
for (kk in 1:4) {
  n = N[kk]
  p = P[kk]
  n1 = n/2
  set.seed(88888888)
  X = matrix(-2+4*runif(n*p),ncol=p)
  
  Z1 = cbind(X,p*matrix(rep(1,n),ncol = 1))
  Z0 = cbind(X,p*matrix(rep(0,n),ncol = 1))
  C11 = Calculating(Z1,Z1)
  C10 = Calculating(Z1,Z0)
  C00 = Calculating(Z0,Z0)
  C = rbind(cbind(C11,C10),cbind(t(C10),C00))
  
  for(f in 1:4){
    print(paste("***n =",n," p =", p, " case",f,"***"))
    flag1 = FLAG[f,1]
    flag2 = FLAG[f,2]
    pb <- txtProgressBar(max = nloop, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    allrep <- foreach(i=1:nloop, .combine='rbind', 
                      .packages = c("PDDS","twinning"),
                      .options.snow = opts) %dopar%
      {
        L = POC(X,flag1,flag2)
        Y1 = L[,1]
        Y0 = L[,2]
        taui = mean(Y1-Y0)
        
        idx1 = CRE(n,n1)
        idx2 = ReM(X,n1,0.2)
        idx3 = BODE(X,n1,C)
        idx4 = MODE(X,n1,1)
        
        IDX = rbind(idx1,idx2,idx3,idx4)
        tauall = c()
        ks1all = c()
        ks0all = c()
        for (ii in 1:4) {
          idxii = IDX[ii,]
          res = mean(Y1[idxii])-mean(Y0[-idxii])
          tauall = c(tauall,res)
          res1 = ks.test(Y1,Y1[idxii])$statistic
          names(res1) = NULL
          ks1all = c(ks1all,res1)
          res0 = ks.test(Y0,Y0[idxii])$statistic
          names(res0) = NULL
          ks0all = c(ks0all,res0)
        }
        sumall = c(taui,tauall,ks1all,ks0all)
        return(sumall)
      }
    close(pb)
    FILE = paste(path,"allsummary-",f,"-",kk,".txt")
    write.table (allrep, file =FILE, sep =" ", row.names =F, col.names =F)
  }
}

# MSE tables for various experiments
for(kk in 1:4){
  MSE = matrix(0,ncol=4,nrow = 4)
  for (f in 1:4) {
    FILE = paste(path,"allsummary-",f,"-",kk,".txt")
    allrep = read.table(FILE)
    mse = (allrep[,2:5]-allrep[,1])^2
    MSE[f,] = apply(mse,2,mean)
  }
  PR = (MSE[,1]-MSE[,4])/MSE[,1]
  perf = cbind(round(MSE*10^4,0),100*round(PR,2))
  print(paste("MSEs*10^4 of of different experiments under n =",N[kk],"p =",P[kk]))
  print(perf)
}

# KS tables for various experiments
for(kk in 1:4){
  MSE = matrix(0,ncol=4,nrow = 4)
  KS1 = matrix(0,ncol=4,nrow = 4)
  KS0 = matrix(0,ncol=4,nrow = 4)
  for (f in 1:4) {
    FILE = paste(path,"allsummary-",f,"-",kk,".txt")
    allrep = read.table(FILE)
    KS1[f,] = apply(allrep[,5:8],2,mean)
    KS0[f,] = apply(allrep[,9:12],2,mean)
  }
  print(paste("KS1*10^4 of of different experiments under n =",N[kk],"p =",P[kk]))
  print(round(KS1*1e2,2))
  # print(paste("KS0*10^4 of of different experiments under n =",N[kk],"p =",P[kk]))
  # print(KS0)
}