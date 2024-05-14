setwd("/panfs/jay/groups/20/panwei/wan01299/EL_cusum/reviewer_comments/")
library(doParallel)
myCluster <- makeCluster(80, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

repetition=1000
d=1

threshold_opt=2.20
threshold_app=2.15
threshold_ks=2.10
threshold_el_mean=6.74
m=1000
maxobs=500
winlen=100
mu0=0
mu1=0.5 

results=foreach (i = 1:repetition, .combine='c', .multicombine=FALSE) %dopar% {

library(ks)
library(el.convex)
library(mvtnorm)
if(d==1){

Sigma=1
library(el.convex)
train0=matrix(rnorm(n=m*d,mean=mu0,sd=Sigma),nrow=m,ncol=d)
train1=matrix(rnorm(n=m*d,mean=mu1,sd=Sigma),nrow=m,ncol=d)
mu0hat=colMeans(train0)
mu1hat=colMeans(train1)
cov0hat=cov(train0)
cov1hat=cov(train1)
var0hat=diag(cov0hat)
var1hat=diag(cov1hat)
lap0hat=colMeans(cbind(exp(-train0%*%rep(-0.5,times=d)),exp(-train0%*%rep(0.5,times=d))))
lap1hat=colMeans(cbind(exp(-train1%*%rep(-0.5,times=d)),exp(-train1%*%rep(0.5,times=d))))
obs=matrix(rnorm(n=maxobs*d,mean=mu0,sd=Sigma),nrow=maxobs,ncol=d)

lapobs=cbind(exp(-obs%*%rep(-0.5,times=d)),exp(-obs%*%rep(0.5,times=d)))
used_mean_obs0=matrix(0,nrow=maxobs,ncol=d)
used_mean_obs1=matrix(0,nrow=maxobs,ncol=d)
used_var_obs0=matrix(0,nrow=maxobs,ncol=d)
used_var_obs1=matrix(0,nrow=maxobs,ncol=d)
used_lap_obs0=matrix(0,nrow=maxobs,ncol=2)
used_lap_obs1=matrix(0,nrow=maxobs,ncol=2)
used_mean_and_var_obs0=matrix(0,nrow=maxobs,ncol=2*d)
used_mean_and_var_obs1=matrix(0,nrow=maxobs,ncol=2*d)

}
if(d>1){
  mu0=rep(0,times=d)
  mu1=rep(0.5,times=d)
  Sigma=diag(rep(1,times=d),nrow = d,ncol = d)
   
  train0=rmvnorm(n=m,mean = mu0,sigma = Sigma)
  train1=rmvnorm(n=m,mean = mu0,sigma = Sigma)
  mu0hat=colMeans(train0)
  mu1hat=colMeans(train1)
  cov0hat=cov(train0)
  cov1hat=cov(train1)
  var0hat=diag(cov0hat)
  var1hat=diag(cov1hat)
  lap0hat=colMeans(cbind(exp(-train0%*%rep(-0.5,times=d)),exp(-train0%*%rep(0.5,times=d))))
  lap1hat=colMeans(cbind(exp(-train1%*%rep(-0.5,times=d)),exp(-train1%*%rep(0.5,times=d))))
  obs=rmvnorm(n=maxobs,mean = mu0,sigma = Sigma)
  
  lapobs=cbind(exp(-obs%*%rep(-0.5,times=d)),exp(-obs%*%rep(0.5,times=d)))
  used_mean_obs0=matrix(0,nrow=maxobs,ncol=d)
  used_mean_obs1=matrix(0,nrow=maxobs,ncol=d)
  used_var_obs0=matrix(0,nrow=maxobs,ncol=d)
  used_var_obs1=matrix(0,nrow=maxobs,ncol=d)
  used_lap_obs0=matrix(0,nrow=maxobs,ncol=2)
  used_lap_obs1=matrix(0,nrow=maxobs,ncol=2)
used_mean_and_var_obs0=matrix(0,nrow=maxobs,ncol=2*d)
used_mean_and_var_obs1=matrix(0,nrow=maxobs,ncol=2*d)

}

  for(i in 1:maxobs){
    used_mean_obs0[i,]=obs[i,]-mu0hat
    used_mean_obs1[i,]=obs[i,]-mu1hat
    used_var_obs0[i,]=obs[i,]^2-mu0hat^2-var0hat
    used_var_obs1[i,]=obs[i,]^2-mu1hat^2-var1hat
    used_lap_obs0[i,]=lapobs[i,]-lap0hat
    used_lap_obs1[i,]=lapobs[i,]-lap1hat
    used_mean_and_var_obs0[i,]=c(used_mean_obs0[i,],used_var_obs0[i,])
    used_mean_and_var_obs1[i,]=c(used_mean_obs1[i,],used_var_obs1[i,])
  }
  
  
el_mean_stop=rep(0,times=length(threshold_el_mean))
for(i1 in 1:length(threshold_el_mean)){
  g_el_mean=numeric(maxobs)
  for (n in (d+1): maxobs) {
  start1=seq(max(1, n-winlen), n-d)
  lr=numeric(length(start1))
  for (j in  start1 ) {
    l1=el.test.newton(x=used_mean_obs1[j:n,],mu=0)
    if (abs(sum(na.omit(c(l1$wts,0)))-1) >= .01) next
    if (max(abs(l1$mu)) >= .01) next
    l0=el.test.newton(x=used_mean_obs0[j:n,],mu=0)
    if (abs(sum(na.omit(c(l0$wts,0)))-1) >= .01) next
    if (max(abs(na.omit(c(l0$mu,0)))) >= .01) next
    lr[which(start1==j)]=(l0[["-2LLR"]] - l1[["-2LLR"]]) / 2
  }
  g_el_mean[n] <- max(c(na.omit(lr[start1]),0))
  print(g_el_mean[n])
  print(n)
  if (g_el_mean[n] >= threshold_el_mean[i1]) break
  if (n >= maxobs) break
  } 
  el_mean_stop[i1]=n
}



opt_stop=rep(0,times=length(threshold_opt))
for(i1 in 1:length(threshold_opt)){
lr=numeric(maxobs)
for(i in 1:maxobs){
  if(d==1){
    lr[i]=sum(dmvnorm(x=obs[i,],mean=mu1,sigma = diag(Sigma,nrow=d,ncol=d),log = TRUE)-dmvnorm(x=obs[i,],mean=mu0,sigma = diag(Sigma,nrow=d,ncol=d),log = TRUE))
  }
  if(d>1){
    lr[i]=sum(dmvnorm(x=obs[i,],mean=mu1,sigma = Sigma,log = TRUE)-dmvnorm(x=obs[i,],mean=mu0,sigma =Sigma,log = TRUE))
  }
}
lr=na.omit(lr)
g_opt=numeric(length(lr))
for(n in 1:length(lr)){
  if(n==1){
    g_opt[n]=max(0,lr[n])
  }
  if(n>1){
    g_opt[n]=max(0,lr[n]+g_opt[n-1])
  }
  if (n >= maxobs) break
  if(g_opt[n]>threshold_opt[i1]) break
}
opt_stop[i1]=n
}

app_stop=rep(0,times=length(threshold_app))
for(i1 in 1:length(threshold_app)){
lr=numeric(maxobs)
for(i in 1:maxobs){
   
    lr[i]=sum(dmvnorm(x=obs[i,],mean=mu1hat,sigma = cov1hat,log = TRUE)-dmvnorm(x=obs[i,],mean=mu0hat,sigma =cov0hat,log = TRUE))
   
}
lr=na.omit(lr)
g_app=numeric(length(lr))
for(n in 1:length(lr)){
  if(n==1){
    g_app[n]=max(0,lr[n])
  }
  if(n>1){
    g_app[n]=max(0,lr[n]+g_opt[n-1])
  }
  if (n >= maxobs) break
  if(g_app[n]>threshold_app[i1]) break
}
app_stop[i1]=n
}

ks_stop=rep(0,times=length(threshold_ks))
for(i1 in 1:length(threshold_ks)){
lr=numeric(maxobs)
if(d==1){
  bd0=hpi(train0)
  bd1=hpi(train1)
  f0=kde(train0, h = bd0, eval.points = obs)$estimate
  f1=kde(train1, h = bd1, eval.points = obs)$estimate
  lr=na.omit(log(f1/f0))
}
if(d>1){
  bd0=numeric(d)
  bd1=numeric(d)
  for(i in 1:d){
    bd0[i]=hpi(train0[,i])
    bd1[i]=hpi(train1[,i])
    h0=diag(bd0,nrow=d,ncol=d)
    h1=diag(bd1,nrow=d,ncol=d)
  }
  f0=kde(train0, h = h0, eval.points = obs)$estimate
  f1=kde(train1, h = h1, eval.points = obs)$estimate
  lr=na.omit(log(f1/f0))
}
g_ks=numeric(length(lr))
for(n in 1:length(lr)){
  if(n==1){
    g_ks[n]=max(0,lr[n])
  }
  if(n>1){
    g_ks[n]=max(0,lr[n]+g_opt[n-1])
  }
  if (n >= maxobs) break
  if(g_ks[n]>threshold_ks[i1]) break
}
ks_stop[i1]=n
}
stop_all=list(el_mean_stop=el_mean_stop,
              
              opt_stop=opt_stop,
              app_stop=app_stop,
              ks_stop=ks_stop)
              
list(stop_all)              
}              
              
ARL0_el_mean=data.frame(a=threshold_el_mean,ARL0=rep(0,times=length(threshold_el_mean)),SD=rep(0,times=length(threshold_el_mean)) )
ARL0_opt=data.frame(a=threshold_opt,ARL0=rep(0,times=length(threshold_opt)),SD=rep(0,times=length(threshold_opt)) )
ARL0_app=data.frame(a=threshold_app,ARL0=rep(0,times=length(threshold_app)),SD=rep(0,times=length(threshold_app)) )
ARL0_ks=data.frame(a=threshold_ks,ARL0=rep(0,times=length(threshold_ks)),SD=rep(0,times=length(threshold_ks)) )

for(i in 1:length(threshold_el_mean)){
   stop_rep=rep(0,times=repetition)
   for(j in 1:repetition){
     stop_rep[j]=results[[j]]$el_mean_stop[i]
   }
   ARL0_el_mean$ARL0[i]=mean(stop_rep)
   ARL0_el_mean$SD[i]=sd(stop_rep)
}
 
for(i in 1:length(threshold_opt)){
  stop_rep=rep(0,times=repetition)
  for(j in 1:repetition){
    stop_rep[j]=results[[j]]$opt_stop[i]
  }
  ARL0_opt$ARL0[i]=mean(stop_rep)
  ARL0_opt$SD[i]=sd(stop_rep)
}
for(i in 1:length(threshold_app)){
  stop_rep=rep(0,times=repetition)
  for(j in 1:repetition){
    stop_rep[j]=results[[j]]$app_stop[i]
  }
  ARL0_app$ARL0[i]=mean(stop_rep)
  ARL0_app$SD[i]=sd(stop_rep)
}
for(i in 1:length(threshold_ks)){
  stop_rep=rep(0,times=repetition)
  for(j in 1:repetition){
    stop_rep[j]=results[[j]]$ks_stop[i]
  }
  ARL0_ks$ARL0[i]=mean(stop_rep)
  ARL0_ks$SD[i]=sd(stop_rep)
}  

ARL0=list(ARL0_el_mean=ARL0_el_mean,
         
         ARL0_opt=ARL0_opt,
         ARL0_app=ARL0_app,
         ARL0_ks=ARL0_ks
         )

save(ARL0,file = "ARL0_Table1_1000.RData")