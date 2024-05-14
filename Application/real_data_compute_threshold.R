setwd("/panfs/jay/groups/20/panwei/wan01299/EL_cusum/real_data/")
a1=read.table(file="well_log.dat",skip = 7)
a1=a1$V1
data_use=a1[1:1500]
 
data_use1=data_use[1:1070]
data_use2=data_use[1071:length(data_use)]
data_use1=data_use1[data_use1>mean(data_use1)-3*sd(data_use1)&data_use1<mean(data_use1)+3*sd(data_use1)]
data_use2=data_use2[data_use2>mean(data_use2)-3*sd(data_use2)&data_use2<mean(data_use2)+3*sd(data_use2)]
data_use=c(data_use1,data_use2)
 
 
library(doParallel)
myCluster <- makeCluster(80, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

repetition=1000
d=1
maxobs=40 
threshold_ael_mean=seq(from=0,to=2,by=0.1)
 
winlen=100 

train0=sample(data_use1,size = 100,replace = TRUE)
train1=sample(data_use2,size = 100,replace = TRUE)
mean0=mean(train0)
mean1=mean(train1)
  
results=foreach (i = 1:repetition, .combine='c', .multicombine=FALSE) %dopar% {
  
obs0=sample(train0,size = 20,replace = TRUE)
obs1=sample(train1,size =maxobs-20,replace = TRUE)
  library(el.convex)
  
  
  
  obs=c(obs0,obs1)
  used_mean_obs1=as.matrix(obs-mean1)
  used_mean_obs0=as.matrix(obs-mean0)


  ael_mean_stop=rep(0,times=length(threshold_ael_mean))
  for(i1 in 1:length(threshold_ael_mean)){
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
      if (g_ael_mean[n] >= threshold_ael_mean[i1]) break
      if (n >= maxobs) break
    }
    ael_mean_stop[i1]=n-20
  }
  
  stop_all=list( ael_mean_stop=ael_mean_stop )
  
  
  list(stop_all)              
}              


ARL1=data.frame(a=threshold_ael_mean,ARL1=rep(0,times=length(threshold_ael_mean)),SD=rep(0,times=length(threshold_ael_mean)))


for(i in 1:length(threshold_ael_mean)){
  stop_rep=rep(0,times=repetition)
  for(j in 1:repetition){
    stop_rep[j]=results[[j]]$ael_mean_stop[i]
  }
  ARL1$ARL1[i]=mean(stop_rep[stop_rep>=0])
  ARL1$SD[i]=sd(stop_rep[stop_rep>=0])
}

ARL1=list( ARL1=ARL1,mean0=mean0,mean1=mean1,train0=train0,train1=train1)

save(ARL1,file = "ARL1_real_data.RData")