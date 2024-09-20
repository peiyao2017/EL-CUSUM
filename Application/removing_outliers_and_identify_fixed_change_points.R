setwd("D:/")
library(changepoint.np)
a=read.table("D:/well_log.dat",header = FALSE)
a=a[,1]
n=length(a)
x=1:n
plot(x,a,main="Well-log data",xlab="time",ylab="nuclear-magnetic response",cex=0.5,pch=20,col="#666666")


h=75

move_median<-rep(0,n)
move_mad<-rep(0,n)
for(t in 1:n){
  move_median[t]=median(a[max(1,t-h):min(n,t+h)])
  move_mad[t]=mad(a[max(1,t-h):min(n,t+h)])
}

residuals=a-move_median

##remove outliers

data_clean=a[abs(residuals)<2*move_mad]
plot(x=c(1:length(data_clean)),data_clean,main="Well-log data after removing outliers",xlab="time",ylab="nuclear-magnetic response",cex=0.5,pch=20,col="#666666")


change_fixed=cpt.np(data_clean,class = FALSE,minseglen=1,penalty="CROPS",pen.value = c(20,1000))

number_of_cp=change_fixed$cpt.out[2,]
p_cost=change_fixed$cpt.out[3,]

plot( number_of_cp , p_cost,xlab="number of change points",ylab="unpenalized cost",main="Number of change points vs penalized cost" )
abline(v=13)
location_fixed=change_fixed$changepoints[[13]]
length(location_fixed)

plot(x=c(1:length(data_clean)),y=data_clean,pch=20,cex=0.5,xlab="time",ylab="nuclear magnetic response", col="#666666",main="Well-log data with identified change points")
for(i in 1:length(location_fixed)){
  abline(v=location_fixed[i],col="#3399FF")
}

data_clean=list(data_clean,location_fixed)
save(data_clean,file = "data_clean.RData")