setwd("/panfs/jay/groups/20/panwei/wan01299/EL_cusum/reviewer_comments/")
library(doParallel)
myCluster <- makeCluster(80, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

repetition=1000
d=1
threshold_opt=2.2
threshold_app=2.15
threshold_ks=2.1
threshold_el_mean=6.74
m=1000
maxobs=500 
winlen=100

 
results=foreach (i = 1:repetition, .combine='c', .multicombine=FALSE) %dopar% {

library(ks)
library(el.convex)
library(mvtnorm)
if(d==1){
mu0=0
mu1=0.5
Sigma=1

 
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
obs1=matrix(rnorm(n=49*d,mean=mu0,sd=Sigma),nrow=49,ncol=d)
obs2=matrix(rnorm(n=(maxobs-49)*d,mean=mu1,sd=Sigma),nrow=maxobs-49,ncol=d)
obs=rbind(obs1,obs2)
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
  obs1=rmvnorm(n=49,mean = mu0,sigma = Sigma)
  obs2=rmvnorm(n=maxobs-49,mean = mu0,sigma = Sigma)
  obs=rbind(obs1,obs2)
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

 
  g_el_mean=numeric(maxobs)
  for (n in (d+1): maxobs) {
  start1=seq(max(1, n-winlen), n-d)
  lr=numeric(length(start1))
  for (j in  start1 ) {
    l1=el.test.newton(x=used_mean_obs1[j:n,],mu=rep(0,times=d))
    if (abs(sum(na.omit(c(l1$wts,0)))-1) >= .01) next
    if (max(abs(l1$mu)) >= .01) next
    l0=el.test.newton(x=used_mean_obs0[j:n,],mu=rep(0,times=d))
    if (abs(sum(na.omit(c(l0$wts,0)))-1) >= .01) next
    if (max(abs(na.omit(c(l0$mu,0)))) >= .01) next
    lr[which(start1==j)]=(l0[["-2LLR"]] - l1[["-2LLR"]]) / 2
  }
  g_el_mean[n] <- max(c(na.omit(lr),0))
  if (g_el_mean[n] >= threshold_el_mean ) break
  if (n >= maxobs) break
  } 
  el_mean_stop=n-50
 
 
 
 
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
  if(g_opt[n]>threshold_opt ) break
}
opt_stop=n-50
 

 
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
    g_app[n]=max(0,lr[n]+g_app[n-1])
  }
  if (n >= maxobs) break
  if(g_app[n]>threshold_app ) break
}
app_stop=n-50
 

 
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
    g_ks[n]=max(0,lr[n]+g_ks[n-1])
  }
  if (n >= maxobs) break
  if(g_ks[n]>threshold_ks) break
}
ks_stop=n-50
 
stop_all=list(el_mean_stop=el_mean_stop,
              opt_stop=opt_stop,
              app_stop=app_stop,
              ks_stop=ks_stop)
              
list(stop_all)              
}              

 
ARL1_el_mean=data.frame(a=threshold_el_mean,ARL0=rep(0,times=length(threshold_el_mean)),
                        SD=rep(0,times=length(threshold_el_mean)),FD=rep(0,times=length(threshold_el_mean)),
                        ND=rep(0,times=length(threshold_el_mean)))
ARL1_opt=data.frame(a=threshold_opt,ARL0=rep(0,times=length(threshold_opt)),
                    SD=rep(0,times=length(threshold_opt)),FD=rep(0,times=length(threshold_opt)),
                    ND=rep(0,times=length(threshold_opt)))
ARL1_app=data.frame(a=threshold_app,ARL0=rep(0,times=length(threshold_app)),
                    SD=rep(0,times=length(threshold_app)),FD=rep(0,times=length(threshold_app)),
                    ND=rep(0,times=length(threshold_app)))
ARL1_ks=data.frame(a=threshold_ks,ARL0=rep(0,times=length(threshold_ks)),
                   SD=rep(0,times=length(threshold_ks)),FD=rep(0,times=length(threshold_ks)),
                   ND=rep(0,times=length(threshold_ks)))

for(i in 1:length(threshold_el_mean)){
  stop_rep=rep(0,times=repetition)
  for(j in 1:repetition){
    stop_rep[j]=results[[j]]$el_mean_stop[i]
  }
  ARL1_el_mean$ARL0[i]=mean(stop_rep[stop_rep>=0])
  ARL1_el_mean$SD[i]=sd(stop_rep[stop_rep>=0])
  ARL1_el_mean$ND[i]=sum(stop_rep=maxobs-50)
  ARL1_el_mean$ND[i]=sum(stop_rep<0)
}

for(i in 1:length(threshold_opt)){
  stop_rep=rep(0,times=repetition)
  for(j in 1:repetition){
    stop_rep[j]=results[[j]]$opt_stop[i]
  }
  ARL1_opt$ARL0[i]=mean(stop_rep[stop_rep>=0])
  ARL1_opt$SD[i]=sd(stop_rep[stop_rep>=0])
  ARL1_opt$ND[i]=sum(stop_rep=maxobs-50)
  ARL1_opt$FD[i]=sum(stop_rep<0)
}

for(i in 1:length(threshold_app)){
  stop_rep=rep(0,times=repetition)
  for(j in 1:repetition){
    stop_rep[j]=results[[j]]$app_stop[i]
  }
  ARL1_app$ARL0[i]=mean(stop_rep[stop_rep>=0])
  ARL1_app$SD[i]=sd(stop_rep[stop_rep>=0])
  ARL1_app$ND[i]=sum(stop_rep=maxobs-50)
  ARL1_app$FD[i]=sum(stop_rep<0)
}

for(i in 1:length(threshold_ks)){
  stop_rep=rep(0,times=repetition)
  for(j in 1:repetition){
    stop_rep[j]=results[[j]]$ks_stop[i]
  }
  ARL1_ks$ARL0[i]=mean(stop_rep[stop_rep>=0])
  ARL1_ks$SD[i]=sd(stop_rep[stop_rep>=0])
  ARL1_ks$ND[i]=sum(stop_rep=maxobs-50)
  ARL1_ks$FD[i]=sum(stop_rep<0)
}

 
 
ARL1=list(ARL1_el_mean=ARL1_el_mean,
         ARL1_opt=ARL1_opt,
         ARL1_app=ARL1_app,
         ARL1_ks=ARL1_ks
         )

save(ARL1,file = "ARL1_Table1_1000.RData")