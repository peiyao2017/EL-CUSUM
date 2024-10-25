setwd("/panfs/jay/groups/20/panwei/wan01299/EL_cusum/table6/")
library(doParallel)
myCluster <- makeCluster(80, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

repetition=1000
d=4

threshold_opt=5 
threshold_app=38
threshold_ks=217
threshold_el_var=38
threshold_el_laplace=28
threshold_ael_var=4
threshold_ael_laplace=9
threshold_tel_var=19
threshold_tel_laplace=13
threshold_tael_var=1.3
threshold_tael_laplace=4.4
v0=list(c(-1.5,1,0.9),c(1.5,1,-0.9))
v1=list(c(-2,1,0.9),c(2,1,-0.9)) 

m=20
maxobs=200 
winlen=100  

c=0

repeat{

results=foreach (i = 1:repetition, .combine='c', .multicombine=FALSE) %dopar% {
  
  
  
  library(ks)
  library(mixsmsn)  
  library(el.convex)
  library(mvtnorm)
  library(sn) 
  
  
  
  train0=matrix(rmix(n=d*m,pii=c(0.5,0.5),family = "Skew.normal",arg = v0),nrow=m,ncol=d)
  train1=matrix(rmix(n=d*m,pii=c(0.5,0.5),family = "Skew.normal",arg = v1),nrow=m,ncol=d)
  mu0hat=colMeans(train0)
  mu1hat=colMeans(train1)
  cov0hat=cov(train0)
  cov1hat=cov(train1)
  var0hat=diag(cov0hat)
  var1hat=diag(cov1hat)
  lap0hat=colMeans(cbind(exp(-train0%*%rep(0.4,times=d)),exp(-train0%*%rep(0.8,times=d))))
  lap1hat=colMeans(cbind(exp(-train1%*%rep(0.4,times=d)),exp(-train1%*%rep(0.8,times=d))))
  obs1=matrix(rmix(n=d*49,pii=c(0.5,0.5),family = "Skew.normal",arg = v0),nrow=49,ncol=d)
  obs2=matrix(rmix(n=d*(maxobs-49),pii=c(0.5,0.5),family = "Skew.normal",arg = v1),nrow=maxobs-49,ncol=d)
  obs=rbind(obs1,obs2)
  lapobs=cbind(exp(-obs%*%rep(0.4,times=d)),exp(-obs%*%rep(0.8,times=d)))
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
  
  g_el_var=numeric(maxobs)
  for (n in (d+1): maxobs) {
    start1=seq(max(1, n-winlen), n-d)
    lr=numeric(length(start1))
    for (j in start1) {
      l1=el.test.newton(x=used_var_obs1[j:n,],mu=rep(0,times=d))
      if (abs(sum(na.omit(c(l1$wts,0)))-1) >= .01) next
      if (max(abs(l1$mu)) >= .01) next
      l0=el.test.newton(x=used_var_obs0[j:n,],mu=rep(0,times=d))
      if (abs(sum(na.omit(c(l0$wts,0)))-1) >= .01) next
      if (max(abs(na.omit(c(l0$mu,0)))) >= .01) next
      lr[which(start1==j)]=(l0[["-2LLR"]] - l1[["-2LLR"]]) / 2
    }
    print(n)
    
    g_el_var[n] <- max(c(na.omit(lr),0))
    print(g_el_var[n])
    if (g_el_var[n] >= threshold_el_var ) break
    if (n >= maxobs) break
  } 
  el_var_stop=n-50
  
  
  g_el_laplace=numeric(maxobs)
  for (n in (2+1): maxobs) {
    start1=seq(max(1, n-winlen), n-2)
    lr=numeric(length(start1))
    for (j in start1) {
      l1=el.test.newton(x=used_lap_obs1[j:n,],mu=rep(0,times=2))
      if (abs(sum(na.omit(c(l1$wts,0)))-1) >= .01) next
      if (abs(max(l1$mu)) >= .01) next
      l0=el.test.newton(x=used_lap_obs0[j:n,],mu=rep(0,times=2))
      if (abs(sum(na.omit(c(l0$wts,0)))-1) >= .01) next
      if (abs(max(na.omit(c(l0$mu,0)))) >= .01) next
      lr[which(start1==j)] <- (l0[["-2LLR"]] - l1[["-2LLR"]]) / 2
    }
    g_el_laplace[n] <- max(c(na.omit(lr),0))
    if (g_el_laplace[n] >= threshold_el_laplace) break
    if (n >= maxobs) break
  }
  el_laplace_stop=n-50
  
  
  g_ael_var=numeric(maxobs)
  for (n in (d+1): maxobs) {
    start1=seq(max(1, n-winlen), n-d)
    lr=numeric(length(start1))
    for (j in start1) {
      l1=el.test.newton(x=rbind(as.matrix(used_var_obs1[j:n,]),-max(1,log(n)/2)*colMeans(as.matrix(used_var_obs1[j:n,]))),mu=rep(0,times=d))
      
      l0=el.test.newton(x=rbind(as.matrix(used_var_obs0[j:n,]),-max(1,log(n)/2)*colMeans(as.matrix(used_var_obs0[j:n,]))),mu=rep(0,times=d))
      
      lr[which(start1==j)]=(l0[["-2LLR"]] - l1[["-2LLR"]]) / 2
    }
    g_ael_var[n] <- max(c(na.omit(lr),0))
    print(n)
    print(length(lr))
    if (g_ael_var[n] >= threshold_ael_var ) break
    if (n >= maxobs) break
  }
  ael_var_stop =n-50
  
  
  
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
    if (g_ael_laplace[n] >= threshold_ael_laplace ) break
    if (n >= maxobs) break
    print(n)
    print(length(lr))
  }
  ael_laplace_stop=n-50
  
  
  
  g_tel_var=numeric(maxobs)
  for (n in (d+1): maxobs) {
    start1=seq(max(1, n-winlen), n-d)
    lr=numeric(length(start1))
    for (j in start1) {
      l1=el.test.newton(x=used_var_obs1[j:n,],mu=rep(0,times=d))
      if (abs(sum(na.omit(c(l1$wts,0)))-1) >= .01) next
      if (max(abs(l1$mu)) >= .01) next
      l0=el.test.newton(x=used_var_obs0[j:n,],mu=rep(0,times=d))
      if (abs(sum(na.omit(c(l0$wts,0)))-1) >= .01) next
      if (max(abs(na.omit(c(l0$mu,0)))) >= .01) next
      l0t=l0[["-2LLR"]]*max(1-l0[["-2LLR"]]/length(j:n),0.5)
      l1t=l1[["-2LLR"]]*max(1-l1[["-2LLR"]]/length(j:n),0.5)
      lr[which(start1==j)] <- (l0t - l1t) / 2
    }
    g_tel_var[n] <- max(c(na.omit(lr),0))
    if (g_tel_var[n] >= threshold_tel_var) break
    if (n >= maxobs) break
  }
  tel_var_stop=n-50
  
  
  
  g_tel_laplace=numeric(maxobs)
  for (n in (2+1): maxobs) {
    start1=seq(max(1, n-winlen), n-2)
    lr=numeric(length(start1))
    for (j in start1) {
      l1=el.test.newton(x=used_lap_obs1[j:n,],mu=rep(0,times=2))
      if (abs(sum(na.omit(c(l1$wts,0)))-1) >= .01) next
      if (abs(max(l1$mu)) >= .01) next
      l0=el.test.newton(x=used_lap_obs0[j:n,],mu=rep(0,times=2))
      if (abs(sum(na.omit(c(l0$wts,0)))-1) >= .01) next
      if (abs(max(na.omit(c(l0$mu,0)))) >= .01) next
      l0t=l0[["-2LLR"]]*max(1-l0[["-2LLR"]]/length(j:n),0.5)
      l1t=l1[["-2LLR"]]*max(1-l1[["-2LLR"]]/length(j:n),0.5)
      lr[which(start1==j)] <- (l0t - l1t) / 2
    }
    g_tel_laplace[n] <- max(c(na.omit(lr),0))
    if (g_tel_laplace[n] >= threshold_tel_laplace ) break
    if (n >= maxobs) break
  }
  tel_laplace_stop=n-50
  
  
  g_tael_var=numeric(maxobs)
  for (n in (d+1): maxobs) {
    start1=seq(max(1, n-winlen), n-d)
    lr=numeric(length(start1))
    for (j in start1) {
      l1=el.test.newton(x=rbind(as.matrix(used_var_obs1[j:n,]),-max(1,log(n)/2)*colMeans(as.matrix(used_var_obs1[j:n,]))),mu=rep(0,times=d))
      
      l0=el.test.newton(x=rbind(as.matrix(used_var_obs0[j:n,]),-max(1,log(n)/2)*colMeans(as.matrix(used_var_obs0[j:n,]))),mu=rep(0,times=d))
      l0t=l0[["-2LLR"]]*max(1-l0[["-2LLR"]]/length(j:n),0.5)
      l1t=l1[["-2LLR"]]*max(1-l1[["-2LLR"]]/length(j:n),0.5)
      lr[which(start1==j)] <- (l0t - l1t) / 2
    }
    g_tael_var[n] <- max(c(na.omit(lr),0))
    print(n)
    print(length(lr))
    if (g_tael_var[n] >= threshold_tael_var ) break
    if (n >= maxobs) break
  }
  tael_var_stop=n-50
  
  
  
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
    if (g_tael_laplace[n] >= threshold_tael_laplace) break
    if (n >= maxobs) break
    print(n)
    print(length(lr))
  }
  tael_laplace_stop=n-50
  
  
  
  lr=numeric(maxobs)
  for(i in 1:maxobs){
    lr[i]=sum(log((0.5*dsn(obs[i,],xi=v1[[1]][1],omega=v1[[1]][2], alpha=v1[[1]][3])+0.5*dsn(obs[i,],xi=v1[[2]][1],omega=v1[[2]][2], alpha=v1[[2]][3]))/(0.5*dsn(obs[i,],xi=v0[[1]][1],omega=v0[[1]][2], alpha=v0[[1]][3])+0.5*dsn(obs[i,],xi=v0[[2]][1],omega=v0[[2]][2], alpha=v0[[2]][3]))))
  }
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
      if(g_opt[n]>threshold_opt ) break
    }
    opt_stop=n-50
 
  
     lr = dmvnorm(x=obs ,mean=mu1hat,sigma = cov1hat,log = TRUE)-dmvnorm(x=obs ,mean=mu0hat,sigma =cov0hat,log = TRUE) 
    
 
  lr=na.omit(lr)
  g_app=numeric(length(lr))
  for(n in 1:length(lr)){
    if(n==1){
      g_app[n]= lr[n] 
    }
    if(n>1){
      g_app[n]=max(0,g_app[n-1])+lr[n]
    }
    if (n >= maxobs) break
    if(g_app[n]>threshold_app ) break
  }
  app_stop=n-50
  
  
  
   
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
      g_ks[n]= lr[n] 
    }
    if(n>1){
      g_ks[n]=max(0,g_ks[n-1])+lr[n]
    }
    if (n >= maxobs) break
    if(g_ks[n]>threshold_ks ) break
  }
  ks_stop=n-50
  
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

stop_rep_el_var=rep(0,times=repetition)
stop_rep_el_laplace=rep(0,times=repetition)
stop_rep_ael_var=rep(0,times=repetition)
stop_rep_ael_laplace=rep(0,times=repetition)
stop_rep_tel_var=rep(0,times=repetition)
stop_rep_tel_laplace=rep(0,times=repetition)
stop_rep_tael_var=rep(0,times=repetition)
stop_rep_tael_laplace=rep(0,times=repetition)
stop_rep_opt=rep(0,times=repetition)
stop_rep_app=rep(0,times=repetition)
stop_rep_ks=rep(0,times=repetition)

for(i in 1:repetition){
  stop_rep_opt[i]=results[[i]]$opt_stop
  stop_rep_app[i]=results[[i]]$app_stop 
  stop_rep_ks[i]=results[[i]]$ks_stop 
  stop_rep_el_var[i]=results[[i]]$el_var_stop 
  stop_rep_el_laplace[i]=results[[i]]$el_laplace_stop
  stop_rep_ael_laplace[i]=results[[i]]$ael_laplace_stop
  stop_rep_ael_var[i]=results[[i]]$ael_var_stop
  stop_rep_tael_var[i]=results[[i]]$tael_var_stop 
  stop_rep_tael_laplace[i]=results[[i]]$tael_laplace_stop 
  stop_rep_tel_var[i]=results[[i]]$tel_var_stop 
  stop_rep_tel_laplace[i]=results[[i]]$tel_laplace_stop 
}


if(c==0){
  stop_rep_el_var_total=stop_rep_el_var 
  stop_rep_el_laplace_total=stop_rep_el_laplace 
  stop_rep_ael_var_total=stop_rep_ael_var 
  stop_rep_ael_laplace_total=stop_rep_ael_laplace 
  stop_rep_tel_var_total=stop_rep_tel_var 
  stop_rep_tel_laplace_total=stop_rep_tel_laplace 
  stop_rep_tael_var_total=stop_rep_tael_var 
  stop_rep_tael_laplace_total=stop_rep_tael_laplace 
  stop_rep_opt_total=stop_rep_opt 
  stop_rep_app_total=stop_rep_app
  stop_rep_ks_total=stop_rep_ks
}

if(c>0){
  stop_rep_el_var_total=c(stop_rep_el_var_total,stop_rep_el_var) 
  stop_rep_el_laplace_total=c(stop_rep_el_laplace_total,stop_rep_el_laplace) 
  stop_rep_ael_var_total=c(stop_rep_ael_var_total,stop_rep_ael_var) 
  stop_rep_ael_laplace_total=c(stop_rep_ael_laplace_total,stop_rep_ael_laplace) 
  stop_rep_tel_var_total=c(stop_rep_tel_var_total,stop_rep_tel_var) 
  stop_rep_tel_laplace_total=c(stop_rep_tel_laplace_total,stop_rep_tel_laplace) 
  stop_rep_tael_var_total=c(stop_rep_tael_var_total,stop_rep_tael_var) 
  stop_rep_tael_laplace_total=c(stop_rep_tael_laplace_total,stop_rep_tael_laplace) 
  stop_rep_opt_total=c(stop_rep_opt_total,stop_rep_opt) 
  stop_rep_app_total=c(stop_rep_app_total,stop_rep_app)
  stop_rep_ks_total=c(stop_rep_ks_total,stop_rep_ks)
}

stop_rep_el_var_total=stop_rep_el_var_total[stop_rep_el_var_total>=0]
stop_rep_el_laplace_total=stop_rep_el_laplace_total[stop_rep_el_laplace_total>=0]
stop_rep_ael_var_total=stop_rep_ael_var_total[stop_rep_ael_var_total>=0]
stop_rep_ael_laplace_total=stop_rep_ael_laplace_total[stop_rep_ael_laplace_total>=0]
stop_rep_tel_var_total=stop_rep_tel_var_total[stop_rep_tel_var_total>=0]
stop_rep_tel_laplace_total=stop_rep_tel_laplace_total[stop_rep_tel_laplace_total>=0]
stop_rep_tael_var_total=stop_rep_tael_var_total[stop_rep_tael_var_total>=0]
stop_rep_tel_laplace_total=stop_rep_tael_laplace_total[stop_rep_tael_laplace_total>=0]
stop_rep_opt_total=stop_rep_opt_total[stop_rep_opt_total>=0]
stop_rep_app_total=stop_rep_app_total[stop_rep_app_total>=0]
stop_rep_ks_total=stop_rep_ks_total[stop_rep_ks_total>=0]

c=c+1

if(length(stop_rep_el_var_total)>repetition&
   length(stop_rep_el_laplace_total)>repetition&
   length(stop_rep_ael_var_total)>repetition&
   length(stop_rep_ael_laplace_total)>repetition&
   length(stop_rep_tel_var_total)>repetition&
   length(stop_rep_tel_laplace_total)>repetition&
   length(stop_rep_tael_var_total)>repetition&
   length(stop_rep_tael_laplace_total)>repetition&
   length(stop_rep_ks_total)>repetition&
   length(stop_rep_opt_total)>repetition&
   length(stop_rep_app_total)>repetition
   ){
  break
}

}


ARL1_el_var=data.frame(a=threshold_el_var,ARL1=mean(stop_rep_el_var_total[1:repetition]),SD=sd(stop_rep_el_var_total[1:repetition]),ND=sum(stop_rep_el_var_total[1:repetition]>=maxobs-50))
ARL1_el_laplace=data.frame(a=threshold_el_laplace,ARL1=mean(stop_rep_el_laplace_total[1:repetition]),SD=sd(stop_rep_el_laplace_total[1:repetition]),ND=sum(stop_rep_el_laplace_total[1:repetition]>=maxobs-50))
ARL1_ael_var=data.frame(a=threshold_ael_var,ARL1=mean(stop_rep_ael_var_total[1:repetition]),SD=sd(stop_rep_ael_var_total[1:repetition]),ND=sum(stop_rep_ael_var_total[1:repetition]>=maxobs-50))
ARL1_ael_laplace=data.frame(a=threshold_ael_laplace,ARL1=mean(stop_rep_ael_laplace_total[1:repetition]),SD=sd(stop_rep_ael_laplace_total[1:repetition]),ND=sum(stop_rep_ael_laplace_total[1:repetition]>=maxobs-50))
ARL1_tel_var=data.frame(a=threshold_tel_var,ARL1=mean(stop_rep_tel_var_total[1:repetition]),SD=sd(stop_rep_tel_var_total[1:repetition]),ND=sum(stop_rep_tel_var_total[1:repetition]>=maxobs-50))
ARL1_tel_laplace=data.frame(a=threshold_tel_laplace,ARL1=mean(stop_rep_tel_laplace_total[1:repetition]),SD=sd(stop_rep_tel_laplace_total[1:repetition]),ND=sum(stop_rep_tel_laplace_total[1:repetition]>=maxobs-50))
ARL1_tael_var=data.frame(a=threshold_tael_var,ARL1=mean(stop_rep_tael_var_total[1:repetition]),SD=sd(stop_rep_tael_var_total[1:repetition]),ND=sum(stop_rep_tael_var_total[1:repetition]>=maxobs-50))
ARL1_tael_laplace=data.frame(a=threshold_tael_laplace,ARL1=mean(stop_rep_tael_laplace_total[1:repetition]),SD=sd(stop_rep_tael_laplace_total[1:repetition]),ND=sum(stop_rep_tael_laplace_total[1:repetition]>=maxobs-50))
ARL1_opt=data.frame(a=threshold_opt,ARL1=mean(stop_rep_opt_total[1:repetition]),SD=sd(stop_rep_opt_total[1:repetition]),ND=sum(stop_rep_opt_total[1:repetition]>=maxobs-50))
ARL1_app=data.frame(a=threshold_app,ARL1=mean(stop_rep_app_total[1:repetition]),SD=sd(stop_rep_app_total[1:repetition]),ND=sum(stop_rep_app_total[1:repetition]>=maxobs-50))
ARL1_ks=data.frame(a=threshold_ks,ARL1=mean(stop_rep_ks_total[1:repetition]),SD=sd(stop_rep_ks_total[1:repetition]),ND=sum(stop_rep_ks_total[1:repetition]>=maxobs-50))


ARL1=list(ARL1_el_var=ARL1_el_var,
          ARL1_el_laplae=ARL1_el_laplace,
          ARL1_ael_var=ARL1_ael_var,
          ARL1_ael_laplae=ARL1_ael_laplace,
          ARL1_tel_var=ARL1_tel_var,
          ARL1_tel_laplae=ARL1_tel_laplace,
          ARL1_tael_var=ARL1_tael_var,
          ARL1_tael_laplae=ARL1_tael_laplace,
          ARL1_opt=ARL1_opt,
          ARL1_app=ARL1_app,
          ARL1_ks=ARL1_ks
)

save(ARL1,file = "ARL1_Table6_20.RData")