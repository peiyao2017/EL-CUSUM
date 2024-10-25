setwd("/panfs/jay/groups/20/panwei/wan01299/EL_cusum/table7/")
library(doParallel)
myCluster <- makeCluster(80, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

repetition=1000
d=2


threshold_opt=3.3
threshold_app=4.5
threshold_ks=5
threshold_el_var=29
threshold_el_laplace=57
threshold_ael_var=3.4
threshold_ael_laplace=5.1
threshold_tel_var=14
threshold_tel_laplace=25
threshold_tael_var=2
threshold_tael_laplace=2.4
 
 

lambda0=4
lambda1=4 
m=1000
maxobs=500
winlen=100  

results=foreach (i = 1:repetition, .combine='c', .multicombine=FALSE) %dopar% {
  
  
  
  library(ks)
  library(mixsmsn)  
  library(el.convex)
  library(mvtnorm)
  library(sn) 
  
 
  
  train0=matrix(rpois(n=d*m,lambda =lambda0 ),nrow=m,ncol=d)
  train1=matrix(rpois(n=d*m,lambda =lambda1 ),nrow=m,ncol=d)
  mu0hat=colMeans(train0)
  mu1hat=colMeans(train1)
  cov0hat=cov(train0)
  cov1hat=cov(train1)
  var0hat=diag(cov0hat)
  var1hat=diag(cov1hat)
  lap0hat=colMeans(cbind(exp(-train0%*%rep(0.2,times=d)),exp(-train0%*%rep(0.6,times=d))))
  lap1hat=colMeans(cbind(exp(-train1%*%rep(0.2,times=d)),exp(-train1%*%rep(0.6,times=d))))
  obs1=matrix(rpois(n=d*49,lambda = lambda0),nrow=49,ncol=d)
  obs2=matrix(rpois(n=d*(maxobs-49),lambda = lambda1),nrow=maxobs-49,ncol=d)
  obs=rbind(obs1,obs2)
  lapobs=cbind(exp(-obs%*%rep(0.2,times=d)),exp(-obs%*%rep(0.6,times=d)))
  used_mean_obs0=matrix(0,nrow=maxobs,ncol=d)
  used_mean_obs1=matrix(0,nrow=maxobs,ncol=d)
  used_var_obs0=matrix(0,nrow=maxobs,ncol=d)
  used_var_obs1=matrix(0,nrow=maxobs,ncol=d)
  used_lap_obs0=matrix(0,nrow=maxobs,ncol=2)
  used_lap_obs1=matrix(0,nrow=maxobs,ncol=2)
  used_mean_and_var_obs0=matrix(0,nrow=maxobs,ncol=d*2)
  used_mean_and_var_obs1=matrix(0,nrow=maxobs,ncol=d*2) 
  
  
  for(i in 1:maxobs){
    used_mean_obs0[i,]=obs[i,]-mu0hat
    used_mean_obs1[i,]=obs[i,]-mu1hat
    used_var_obs0[i,]=(obs[i,]-mu0hat)^2-var0hat
    used_var_obs1[i,]=(obs[i,]-mu1hat)^2-var1hat
    used_lap_obs0[i,]=lapobs[i,]-lap0hat
    used_lap_obs1[i,]=lapobs[i,]-lap1hat
    used_mean_and_var_obs0[i,]=c(used_mean_obs0[i,],used_var_obs0[i,])
    used_mean_and_var_obs1[i,]=c(used_mean_obs1[i,],used_var_obs1[i,])
  }
  
  el_var_stop=rep(0,times=length(threshold_el_var))
  for(i1 in 1:length(threshold_el_var)){
    g_el_var=numeric(maxobs)
    for (n in (2*d+1): maxobs) {
      start1=seq(max(1, n-winlen), n-2*d)
      lr=numeric(length(start1))
      for (j in start1) {
        l1=el.test.newton(x=used_mean_and_var_obs1[j:n,],mu=rep(0,times=2*d))
        if (abs(sum(na.omit(c(l1$wts,1e-10)))-1) >= .01) next
        if (max(abs(na.omit(c(l1$mu,0)))) >= .01) next
        l0=el.test.newton(x=used_mean_and_var_obs0[j:n,],mu=rep(0,times=2*d))
        if (abs(sum(na.omit(c(l0$wts,1e-10)))-1) >= .01) next
        if (max(abs(na.omit(c(l0$mu,0)))) >= .01) next
        lr[which(start1==j)]=(l0[["-2LLR"]] - l1[["-2LLR"]]) / 2
      }
      g_el_var[n] <- max(c(na.omit(lr),0))
      if (g_el_var[n] >= threshold_el_var[i1]) break
      if (n>=maxobs) break
    } 
    el_var_stop[i1]=n
  }
  
  el_laplace_stop=rep(0,times=length(threshold_el_laplace))
  for(i1 in 1:length(threshold_el_laplace)){
    g_el_laplace=numeric(maxobs)
    for (n in (2+1): maxobs) {
      start1=seq(max(1, n-winlen), n-2)
      lr=numeric(length(start1))
      for (j in start1) {
        l1=el.test.newton(x=used_lap_obs1[j:n,],mu=rep(0,times=2))
        if (abs(sum(na.omit(c(l1$wts,1e-10)))-1) >= .01) next
        if (abs(max(na.omit(c(l1$mu,0)))) >= .01) next
        l0=el.test.newton(x=used_lap_obs0[j:n,],mu=rep(0,times=2))
        if (abs(sum(na.omit(c(l0$wts,1e-10)))-1) >= .01) next
        if (abs(max(na.omit(c(l0$mu,0)))) >= .01) next
        lr[which(start1==j)] <- (l0[["-2LLR"]] - l1[["-2LLR"]]) / 2
      }
      g_el_laplace[n] <- max(c(na.omit(lr),0))
      if (g_el_laplace[n] >= threshold_el_laplace[i1]) break
      if (n>=maxobs) break
    }
    el_laplace_stop[i1]=n
  }
  
  ael_var_stop=rep(0,times=length(threshold_ael_var))
  for(i1 in 1:length(threshold_ael_var)){
    g_ael_var=numeric(maxobs)
    for (n in (2*d+1): maxobs) {
      start1=seq(max(1, n-winlen), n-2*d)
      lr=numeric(length(start1))
      for (j in start1) {
        l1=el.test.newton(x=rbind(as.matrix(used_mean_and_var_obs1[j:n,]),-max(1,log(n)/2)*colMeans(as.matrix(used_mean_and_var_obs1[j:n,]))),mu=rep(0,times=2*d))
        
        l0=el.test.newton(x=rbind(as.matrix(used_mean_and_var_obs0[j:n,]),-max(1,log(n)/2)*colMeans(as.matrix(used_mean_and_var_obs0[j:n,]))),mu=rep(0,times=2*d))
        
        lr[which(start1==j)]=(l0[["-2LLR"]] - l1[["-2LLR"]]) / 2
      }
      g_ael_var[n] <- max(c(na.omit(lr),0))
      print(n)
      print(length(lr))
      if (g_ael_var[n] >= threshold_ael_var[i1]) break
      if (n>=maxobs) break
    }
    ael_var_stop[i1]=n
  }
  
  ael_laplace_stop=rep(0,times=length(threshold_ael_laplace))
  for(i1 in 1:length(threshold_ael_laplace)){
    g_ael_laplace=numeric(maxobs)
    for (n in (2+1): maxobs) {
      start1=seq(max(1, n-winlen), n-2)
      lr=numeric(length(start1))
      for (j in start1) {
        l1=el.test.newton(x=rbind(used_lap_obs1[j:n,],-max(1,log(n)/2)*colMeans(used_lap_obs1[j:n,])),mu=rep(0,times=2),tol = 1e-04)
        
        l0=el.test.newton(x=rbind(used_lap_obs0[j:n,],-max(1,log(n)/2)*colMeans(used_lap_obs0[j:n,])),mu=rep(0,times=2),tol = 1e-04)
        
        lr[which(start1==j)] <- (l0[["-2LLR"]] - l1[["-2LLR"]]) / 2
      }
      g_ael_laplace[n] <- max(c(na.omit(lr),0))
      if (g_ael_laplace[n] >= threshold_ael_laplace[i1]) break
      if (n>=maxobs) break
      print(n)
      print(length(lr))
    }
    ael_laplace_stop[i1]=n
  }
  
  tel_var_stop=rep(0,times=length(threshold_tel_var))
  for(i1 in 1:length(threshold_tel_var)){
    g_tel_var=numeric(maxobs)
    for (n in (2*d+1): maxobs) {
      start1=seq(max(1, n-winlen), n-2*d)
      lr=numeric(length(start1))
      for (j in start1) {
        l1=el.test.newton(x=used_mean_and_var_obs1[j:n,],mu=rep(0,times=2*d))
        if (abs(sum(na.omit(c(l1$wts,1e-10)))-1) >= .01) next
        if (max(abs(na.omit(c(l1$mu,0)))) >= .01) next
        l0=el.test.newton(x=used_mean_and_var_obs0[j:n,],mu=rep(0,times=2*d))
        if (abs(sum(na.omit(c(l0$wts,1e-10)))-1) >= .01) next
        if (max(abs(na.omit(c(l0$mu,0)))) >= .01) next
        l0t=l0[["-2LLR"]]*max(1-l0[["-2LLR"]]/length(j:n),0.5)
        l1t=l1[["-2LLR"]]*max(1-l1[["-2LLR"]]/length(j:n),0.5)
        lr[which(start1==j)] <- (l0t - l1t) / 2
      }
      g_tel_var[n] <- max(c(na.omit(lr),0))
      if (g_tel_var[n] >= threshold_tel_var[i1]) break
      if (n>=maxobs) break
    }
    tel_var_stop[i1]=n
  }
  
  tel_laplace_stop=rep(0,times=length(threshold_tel_laplace))
  for(i1 in 1:length(threshold_tel_laplace)){
    g_tel_laplace=numeric(maxobs)
    for (n in (2+1): maxobs) {
      start1=seq(max(1, n-winlen), n-2)
      lr=numeric(length(start1))
      for (j in start1) {
        l1=el.test.newton(x=used_lap_obs1[j:n,],mu=rep(0,times=2))
        if (abs(sum(na.omit(c(l1$wts,1e-10)))-1) >= .01) next
        if (abs(max(na.omit(c(l1$mu,0)))) >= .01) next
        l0=el.test.newton(x=used_lap_obs0[j:n,],mu=rep(0,times=2))
        if (abs(sum(na.omit(c(l0$wts,1e-10)))-1) >= .01) next
        if (abs(max(na.omit(c(l0$mu,0)))) >= .01) next
        l0t=l0[["-2LLR"]]*max(1-l0[["-2LLR"]]/length(j:n),0.5)
        l1t=l1[["-2LLR"]]*max(1-l1[["-2LLR"]]/length(j:n),0.5)
        lr[which(start1==j)] <- (l0t - l1t) / 2
      }
      g_tel_laplace[n] <- max(c(na.omit(lr),0))
      if (g_tel_laplace[n] >= threshold_tel_laplace[i1]) break
      if (n>=maxobs) break
    }
    tel_laplace_stop[i1]=n
  }
  
  tael_var_stop=rep(0,times=length(threshold_tael_var))
  for(i1 in 1:length(threshold_tael_var)){
    g_tael_var=numeric(maxobs)
    for (n in (2*d+1): maxobs) {
      start1=seq(max(1, n-winlen), n-2*d)
      lr=numeric(length(start1))
      for (j in start1) {
        l1=el.test.newton(x=rbind(as.matrix(used_mean_and_var_obs1[j:n,]),-max(1,log(n)/2)*colMeans(as.matrix(used_mean_and_var_obs1[j:n,]))),mu=rep(0,times=2*d))
        
        l0=el.test.newton(x=rbind(as.matrix(used_mean_and_var_obs0[j:n,]),-max(1,log(n)/2)*colMeans(as.matrix(used_mean_and_var_obs0[j:n,]))),mu=rep(0,times=2*d))
        l0t=l0[["-2LLR"]]*max(1-l0[["-2LLR"]]/length(j:n),0.5)
        l1t=l1[["-2LLR"]]*max(1-l1[["-2LLR"]]/length(j:n),0.5)
        lr[which(start1==j)] <- (l0t - l1t) / 2
      }
      g_tael_var[n] <- max(c(na.omit(lr),0))
      print(n)
      print(length(lr))
      if (g_tael_var[n] >= threshold_tael_var[i1]) break
      if (n>=maxobs) break
    }
    tael_var_stop[i1]=n
  }
  
  tael_laplace_stop=rep(0,times=length(threshold_tael_laplace))
  for(i1 in 1:length(threshold_tael_laplace)){
    g_tael_laplace=numeric(maxobs)
    for (n in (2+1): maxobs) {
      start1=seq(max(1, n-winlen), n-2)
      lr=numeric(length(start1))
      for (j in start1) {
        l1=el.test.newton(x=rbind(used_lap_obs1[j:n,],-max(1,log(n)/2)*colMeans(used_lap_obs1[j:n,])),mu=rep(0,times=2),tol = 1e-04)
        
        l0=el.test.newton(x=rbind(used_lap_obs0[j:n,],-max(1,log(n)/2)*colMeans(used_lap_obs0[j:n,])),mu=rep(0,times=2),tol = 1e-04)
        
        l0t=l0[["-2LLR"]]*max(1-l0[["-2LLR"]]/length(j:n),0.5)
        l1t=l1[["-2LLR"]]*max(1-l1[["-2LLR"]]/length(j:n),0.5)
        lr[which(start1==j)] <- (l0t - l1t) / 2
      }
      g_tael_laplace[n] <- max(c(na.omit(lr),0))
      if (g_tael_laplace[n] >= threshold_tael_laplace[i1]) break
      if (n>=maxobs) break
      print(n)
      print(length(lr))
    }
    tael_laplace_stop[i1]=n
  }
  
  opt_stop=rep(0,times=length(threshold_opt))
  for(i1 in 1:length(threshold_opt)){
    lr =  rowSums(dpois(x=obs ,lambda = lambda1,log = TRUE))-rowSums(dpois(x=obs ,lambda = lambda0,log = TRUE))
    
    lr=na.omit(lr)
    g_opt=numeric(length(lr))
    for(n in 1:length(lr)){
      if(n==1){
        g_opt[n]=lr[n]
      }
      if(n>1){
        g_opt[n]=max(0,g_opt[n-1])+lr[n]
      }
      if(n>=maxobs) break
      if(g_opt[n]>threshold_opt[i1]) break
    }
    opt_stop[i1]=n
  }
  
  app_stop=rep(0,times=length(threshold_app))
  for(i1 in 1:length(threshold_app)){
      lr = dmvnorm(x=obs ,mean=mu1hat,sigma = cov1hat,log = TRUE)-dmvnorm(x=obs ,mean=mu0hat,sigma =cov0hat,log = TRUE) 
      
   
    lr=na.omit(lr)
    g_app=numeric(length(lr))
    for(n in 1:length(lr)){
      if(n==1){
        g_app[n]=lr[n]
      }
      if(n>1){
        g_app[n]=max(0,g_app[n-1])+lr[n]
      }
      if(n>=maxobs) break
      if(g_app[n]>threshold_app[i1]) break
    }
    app_stop[i1]=n
  }
  
  ks_stop=rep(0,times=length(threshold_ks))
  for(i1 in 1:length(threshold_ks)){
     
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
        g_ks[n]=lr[n]
      }
      if(n>1){
        g_ks[n]=max(0,g_ks[n-1])+lr[n]
      }
      if(n>=maxobs) break
      if(g_ks[n]>threshold_ks[i1]) break
    }
    ks_stop[i1]=n
  }
  stop_all=list(el_var_stop=el_var_stop,
                el_laplace_stop=el_laplace_stop,
                ael_var_stop=ael_var_stop,
                ael_laplace_stop=ael_laplace_stop,
                tel_var_stop=tel_var_stop,
                tel_laplace_stop=tel_laplace_stop,
                tael_var_stop=tael_var_stop,
                tael_laplace_stop=tael_laplace_stop,
                opt_stop=opt_stop,
                app_stop=app_stop,
                ks_stop=ks_stop)
  
  list(stop_all)              
}              

ARL0_el_var=data.frame(a=threshold_el_var,ARL0=rep(0,times=length(threshold_el_var)),SD=rep(0,times=length(threshold_el_var)))
ARL0_el_laplace=data.frame(a=threshold_el_laplace,ARL0=rep(0,times=length(threshold_el_laplace)),SD=rep(0,times=length(threshold_el_laplace)))
ARL0_ael_var=data.frame(a=threshold_ael_var,ARL0=rep(0,times=length(threshold_ael_var)),SD=rep(0,times=length(threshold_ael_var)))
ARL0_ael_laplace=data.frame(a=threshold_ael_laplace,ARL0=rep(0,times=length(threshold_ael_laplace)),SD=rep(0,times=length(threshold_ael_laplace)))
ARL0_tel_var=data.frame(a=threshold_tel_var,ARL0=rep(0,times=length(threshold_tel_var)),SD=rep(0,times=length(threshold_tel_var)))
ARL0_tel_laplace=data.frame(a=threshold_tel_laplace,ARL0=rep(0,times=length(threshold_tel_laplace)),SD=rep(0,times=length(threshold_tel_laplace)))
ARL0_tael_var=data.frame(a=threshold_tael_var,ARL0=rep(0,times=length(threshold_tael_var)),SD=rep(0,times=length(threshold_tael_var)))
ARL0_tael_laplace=data.frame(a=threshold_tael_laplace,ARL0=rep(0,times=length(threshold_tael_laplace)),SD=rep(0,times=length(threshold_tael_laplace)))
ARL0_opt=data.frame(a=threshold_opt,ARL0=rep(0,times=length(threshold_opt)),SD=rep(0,times=length(threshold_opt)))
ARL0_app=data.frame(a=threshold_app,ARL0=rep(0,times=length(threshold_app)),SD=rep(0,times=length(threshold_app)))
ARL0_ks=data.frame(a=threshold_ks,ARL0=rep(0,times=length(threshold_ks)),SD=rep(0,times=length(threshold_ks)))

for(i in 1:length(threshold_el_var)){
   stop_rep=rep(0,times=repetition)
   for(j in 1:repetition){
     stop_rep[j]=results[[j]]$el_var_stop[i]
   }
   ARL0_el_var$ARL0[i]=mean(stop_rep)
   ARL0_el_var$SD[i]=sd(stop_rep)
}
for(i in 1:length(threshold_el_laplace)){
  stop_rep=rep(0,times=repetition)
  for(j in 1:repetition){
    stop_rep[j]=results[[j]]$el_laplace_stop[i]
  }
  ARL0_el_laplace$ARL0[i]=mean(stop_rep)
  ARL0_el_laplace$SD[i]=sd(stop_rep)
}
for(i in 1:length(threshold_ael_var)){
  stop_rep=rep(0,times=repetition)
  for(j in 1:repetition){
    stop_rep[j]=results[[j]]$ael_var_stop[i]
  }
  ARL0_ael_var$ARL0[i]=mean(stop_rep)
  ARL0_ael_var$SD[i]=sd(stop_rep)
}
for(i in 1:length(threshold_ael_laplace)){
  stop_rep=rep(0,times=repetition)
  for(j in 1:repetition){
    stop_rep[j]=results[[j]]$ael_laplace_stop[i]
  }
  ARL0_ael_laplace$ARL0[i]=mean(stop_rep)
  ARL0_ael_laplace$SD[i]=sd(stop_rep)
}  
for(i in 1:length(threshold_tel_var)){
  stop_rep=rep(0,times=repetition)
  for(j in 1:repetition){
    stop_rep[j]=results[[j]]$tel_var_stop[i]
  }
  ARL0_tel_var$ARL0[i]=mean(stop_rep)
  ARL0_tel_var$SD[i]=sd(stop_rep)
}
for(i in 1:length(threshold_tel_laplace)){
  stop_rep=rep(0,times=repetition)
  for(j in 1:repetition){
    stop_rep[j]=results[[j]]$tel_laplace_stop[i]
  }
  ARL0_tel_laplace$ARL0[i]=mean(stop_rep)
  ARL0_tel_laplace$SD[i]=sd(stop_rep)
}
for(i in 1:length(threshold_tael_var)){
  stop_rep=rep(0,times=repetition)
  for(j in 1:repetition){
    stop_rep[j]=results[[j]]$tael_var_stop[i]
  }
  ARL0_tael_var$ARL0[i]=mean(stop_rep)
  ARL0_tael_var$SD[i]=sd(stop_rep)
}
for(i in 1:length(threshold_tael_laplace)){
  stop_rep=rep(0,times=repetition)
  for(j in 1:repetition){
    stop_rep[j]=results[[j]]$tael_laplace_stop[i]
  }
  ARL0_tael_laplace$ARL0[i]=mean(stop_rep)
  ARL0_tael_laplace$SD[i]=sd(stop_rep)
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

ARL0=list(ARL0_el_var=ARL0_el_var,
         ARL0_el_laplae=ARL0_el_laplace,
         ARL0_ael_var=ARL0_ael_var,
         ARL0_ael_laplae=ARL0_ael_laplace,
         ARL0_tel_var=ARL0_tel_var,
         ARL0_tel_laplae=ARL0_tel_laplace,
         ARL0_tael_var=ARL0_tael_var,
         ARL0_tael_laplae=ARL0_tael_laplace,
         ARL0_opt=ARL0_opt,
         ARL0_app=ARL0_app,
         ARL0_ks=ARL0_ks
         )

save(ARL0,file = "ARL0_Table7_1000.RData")