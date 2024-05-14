setwd("/panfs/jay/groups/20/panwei/wan01299/EL_cusum/reviewer_comments/")
library(doParallel)
myCluster <- makeCluster(80, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

repetition=1000
d=1

 
threshold_el_mean=6.74
m=1000
nmax=500
winlen=100
mu0=0
mu1=0.5
Sigma=1 
stopfun <- function(x, th) {
  idx <- which(x >= th)
  if (length(idx) == 0) return(length(x))
  return(idx[1])
}  

results=foreach (i = 1:repetition, .combine='cbind' ) %dopar% {

 
library(el.convex)
 

 
train0=matrix(rnorm(n=m*d,mean=mu0,sd=Sigma),nrow=m,ncol=d)
train1=matrix(rnorm(n=m*d,mean=mu1,sd=Sigma),nrow=m,ncol=d)
mu0hat=colMeans(train0)
mu1hat=colMeans(train1)

  x=rnorm(nmax, mu0, Sigma)	
  g <- numeric(nmax)
  for (n in 2:nmax) {
    # if (n %% 50 == 0) print(n)
    start <- seq(max(1, n-winlen), n-1)
    lr <- numeric(n)
    lam0 <- lam1 <- 0
    for (j in start) {
      l1 <- el.test.newton(x[j:n], mu1hat, lam1)
      if(length(na.omit(na.omit(l1$wts)))==0) next # NAs may appear in some repetitions, we need to remove NAs otherwise system returns an error  
      if(length(na.omit(na.omit(l1$mu)))==0) next
      if (abs(sum(na.omit(l1$wts))-1) >= .01) next
      if (abs(na.omit(l1$mu)) >= .01) next
      l0 <- el.test.newton(x[j:n], mu0hat, lam0)
      if(length(na.omit(na.omit(l0$wts)))==0) next
      if(length(na.omit(na.omit(l0$mu)))==0) next
      if (abs(sum(na.omit(l0$wts))-1) >= .01) next
      if (abs(na.omit(l0$mu)) >= .01) next
      lr[j] <- (l0[["-2LLR"]] - l1[["-2LLR"]]) / 2
      lam0 <- l0$lambda
      lam1 <- l1$lambda
    }
    g[n] <- max(lr[start])
  }
  g
 
 
}  

rl0=apply(results, 2, stopfun, th = threshold_el_mean)
ARL0_el_mean=data.frame(ARL0=mean(rl0),SD=sd(rl0) / sqrt(repetition))

ARL0=list(ARL0_el_mean=ARL0_el_mean)
save(ARL0,file = "ARL0_EL_Table1_1000.RData")