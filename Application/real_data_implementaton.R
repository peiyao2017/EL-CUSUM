setwd("D:/")
threshold_ael_mean=4
d=1
winlen=100
load("D:/data_clean.RData")
location_fixed=data_clean[[2]]
data=data_clean[[1]]
data_use=data[1:(location_fixed[2]-1)]

data_use1=data_use[1:(location_fixed[1]-1)]
data_use2=data_use[location_fixed[1]:length(data_use)]
 
load("ARL1_real_data.RData")

data_use1=data_use1[!data_use1%in%c(ARL1$train0,ARL1$train1)]
data_use2=data_use2[!data_use2%in%c(ARL1$train0,ARL1$train1)]
data_use=c(data_use1,data_use2)



mu0hat=ARL1$mean0
mu1hat=ARL1$mean1
obs=data_use
maxobs=length(obs)

used_mean_obs0=matrix(obs-mu0hat,nrow=length(obs),ncol=d)
used_mean_obs1=matrix(obs-mu1hat,nrow=length(obs),ncol=d)
library(el.convex)
g_ael_mean=numeric(maxobs)
for (n in (d+1): maxobs) {
  start1=seq(max(1, n-winlen), n-d)
  lr=numeric(length(start1))
  for (j in  start1 ) {
    l1=el.test.newton(x=rbind(as.matrix(used_mean_obs1[j:n,]),-max(1,log(n)/2)*colMeans(as.matrix(used_mean_obs1[j:n,]))),mu=rep(0,times=d))
    
    l0=el.test.newton(x=rbind(as.matrix(used_mean_obs0[j:n,]),-max(1,log(n)/2)*colMeans(as.matrix(used_mean_obs0[j:n,]))),mu=rep(0,times=d))
    
    lr[which(start1==j)]=(l0[["-2LLR"]] - l1[["-2LLR"]]) / 2
  }
  g_ael_mean[n] <- max(c(na.omit(lr),0))
  print(n)
  print(g_ael_mean[n])
  if (g_ael_mean[n] >= threshold_ael_mean ) break
  if (n >= maxobs) break
}
stop=n 
change=n-9

load("D:/data_clean.RData")
location_fixed=data_clean[[2]]
data=data_clean[[1]]
data_use=data[1:(location_fixed[2]-1)]

data_use1=data_use[1:(location_fixed[1]-1)]
data_use2=data_use[location_fixed[1]:length(data_use)] 
 
plot(x=c(1:length(data_use)),y=data_use,main="Observations before and after change",xlab="time",ylab="nuclear-magnetic response",type="n",ylim=c(min(data_use)-10000,max(data_use)+10000))
points(x=c(1:length(data_use1)),y=data_use1,col="#FFCC00", cex=0.8,pch=20)
points(x=c((1+length(data_use1)):(length(data_use1)+length(data_use2))),y=data_use2,col="#3399FF", cex=0.8,pch=20)
legend("bottomleft",legend=c("before change","after change" ),cex=c(0.8,0.8 ),col=c("#FFCC00","#3399FF" ),pch=c(20,20 ) )

data_use1=data_use1[!data_use1%in%c(ARL1$train0,ARL1$train1)]
data_use2=data_use2[!data_use2%in%c(ARL1$train0,ARL1$train1)]
data_use=c(data_use1,data_use2)
plot(x=c(1:length(data_use)),y=data_use,main="Chang-point detection for test data",xlab="time",ylab="nuclear-magnetic response",type="n")
abline(v=stop-10,col="#3399FF",lwd=2)
points(x=c(1:length(data_use1)),y=data_use1,col="#666666", cex=0.8,pch=20)
points(x=c((length(data_use1)+1):(length(data_use1)+length(data_use2))),y=data_use2,col="#666666", cex=0.8,pch=20)