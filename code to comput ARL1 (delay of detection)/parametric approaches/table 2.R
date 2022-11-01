arl11=numeric()
arl12=numeric()
arl13=numeric()

repeat{
  a=foreach(i=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    opt1_cusum(threshold =7)
  }
  arl11=c(arl11,a[a!="false"])
  if(length(arl11)>2000){
    break
  }
}

repeat{
  a=foreach(i=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    app1_cusum(threshold =3005.5,m=20)
  }
  arl12=c(arl12,a[a!="false"])
  if(length(arl12)>2000){
    break
  }
}

repeat{
  a=foreach(i=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    library(ks)
    ks1_cusum(threshold =192246.5,m=20)
  }
  arl13=c(arl13,a[a!="false"])
  if(length(arl13)>2000){
    break
  }
}






x=data.frame(arl11=arl11[1:2000],arl12=arl12[1:2000],arl13=arl13[1:2000])

for(i in 1:3){
  x[,i]=as.numeric(x[,i])
}

delay=matrix(0,nrow=3,ncol=1000)
nodetect=matrix(0,nrow=3,ncol=1000)
for(i in 1:3){
  for(j in 1:1000){
    delay[i,j]=x[2*j-1,i]
    nodetect[i,j]=x[2*j,i]
  }
}


mean(delay[1,])
sd(delay[1,])/sqrt(1000)
sum(nodetect[1,])


mean(delay[2,])
sd(delay[2,])/sqrt(1000)
sum(nodetect[2,])


mean(delay[3,])
sd(delay[3,])/sqrt(1000)
sum(nodetect[3,])


app1_cusum=function(m=100,d=4,threshold=200,mu0=rep(0,times=d),mu1=rep(0.5,times=d),covar=diag(1,d,d)){
  false=0
  c=0
  y0=mvrnorm(n=m,mu=mu0,Sigma=covar)
  y1=mvrnorm(n=m,mu=mu1,Sigma=covar)
  estmu0=colMeans(y0)
  estmu1=colMeans(y1)
  estcovar0=cov(y0)
  estcovar1=cov(y1)
  repeat{
    repeat{
      x=mvrnorm(n=1,mu=mu0,Sigma = covar)
      c=c+1
      print(c)
      x1=numeric()
      x2=numeric()
      c1=0
      repeat{
        if(c1==0){
          x2[1]=log(dmvnorm(x ,mean = estmu1,sigma = estcovar1)/dmvnorm(x ,mean = estmu0,sigma = estcovar0))
          x1[1]=x2[1]
        }
        print(max(x1))
        if(max(x1)>threshold){
          false=false+1
          break
        }
        if(max(x1)<0){
          break
        }
        if(c1>=1){
          x2[length(x2)+1]=log(dmvnorm(x[nrow(x),],mean = estmu1,sigma = estcovar1)/dmvnorm(x[nrow(x),],mean = estmu0,sigma = estcovar0))
          for(i in 1:length(x2)){
            x1[i]=sum(x2[i:length(x2)])
          }
        }
        print(max(x1))
        if(max(x1)>threshold){
          false=false+1
          break
        }
        if(max(x1)<0){
          break
        }
        x=rbind(x,mvrnorm(n=1,mu=mu0,Sigma = covar))
        c=c+1
        c1=c1+1
        print(c) 
        if(c>=49){
          break
        }
        if(false>=10){
          break
        }
      }
      if(max(x1)>threshold){
        break
      }
      if(c>=49){
        break
      }
      if(false>=10){
        break
      }
    }
    if(c>=49){
      break
    }
    if(false>=10){
      break
    }
  }
  if(false<10){
    nodetect=0
    repeat{
      x=rbind(x,mvrnorm(n=1,mu=mu1,Sigma = covar))
      c=c+1
      x2[length(x2)+1]=log(dmvnorm(x[nrow(x),],mean = estmu1,sigma = estcovar1)/dmvnorm(x[nrow(x),],mean = estmu0,sigma = estcovar0))
      for(i in 1:length(x2)){
        x1[i]=sum(x2[i:length(x2)])
      }
      if(max(x1)>threshold){
        break
      }
      if(max(x1)<0){
        break
      }
      if(c>1000){
        nodetect=1
        break
      }
    }
    
    if(max(x1)<=threshold){
      repeat{
        x=mvrnorm(n=1,mu=mu1,Sigma = covar)
        c=c+1
        print(c)
        x1=numeric()
        x2=numeric()
        c1=0
        repeat{
          if(c1==0){
            x2[1]=log(dmvnorm(x ,mean = estmu1,sigma = estcovar1)/dmvnorm(x ,mean = estmu0,sigma = estcovar0))
            x1[1]=x2[1]
          }
          print(max(x1))
          if(max(x1)>threshold){
            break
          }
          if(max(x1)<0){
            break
          }
          if(c1>=1){
            x2[length(x2)+1]=log(dmvnorm(x[nrow(x),],mean = estmu1,sigma = estcovar1)/dmvnorm(x[nrow(x),],mean = estmu0,sigma = estcovar0))
            for(i in 1:length(x2)){
              x1[i]=sum(x2[i:length(x2)])
            }
          }
          print(max(x1))
          if(max(x1)>threshold){
            break
          }
          if(max(x1)<0){
            break
          }
          x=rbind(x,mvrnorm(n=1,mu=mu1,Sigma = covar))
          c=c+1
          c1=c1+1
          print(c) 
          if(c>1000){
            nodetect=1
            break
          }
        }
        if(max(x1)>threshold){
          break
        }
        if(c>1000){
          nodetect=1
          break
        }
      }
    }
    return(c(c-50,nodetect))
  }
  if(false>=10){
    return("false")
  }
}

ks1_cusum=function(m=100,d=4,threshold=100,mu0=rep(0,times=d),mu1=rep(0.5,times=d),covar=diag(1,d,d)){
  false=0
  c=0
  y0=mvrnorm(n=m,mu=mu0,Sigma=covar)
  y1=mvrnorm(n=m,mu=mu1,Sigma=covar)
  wh0=numeric()
  wh1=numeric()
  for(i in 1:d){
    wh0[i]=hpi(y0[,i],deriv.order =0)
    wh1[i]=hpi(y1[,i],deriv.order =0)
  }
  f0=function(x=rep(0,times=d)){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(sqrt(2*pi)*exp(-((x-y0[i,])/wh0)^2/2)/wh0)
    }
    return(mean(p))
  }
  f1=function(x=rep(0,times=d)){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(sqrt(2*pi)*exp(-((x-y1[i,])/wh1)^2/2)/wh1)
    }
    return(mean(p))
  }
  repeat{
    repeat{
      x=mvrnorm(n=1,mu=mu0,Sigma = covar)
      c=c+1
      print(c)
      x1=numeric()
      x2=numeric()
      c1=0
      repeat{
        if(c1==0){
          x2[1]=log(f1(x)/f0(x))
          x1[1]=x2[1]
        }
        print(max(x1))
        if(max(x1)>threshold){
          false=false+1
          break
        }
        if(max(x1)<0){
          break
        }
        if(c1>=1){
          x2[length(x2)+1]=log(f1(x[nrow(x),])/f0(x[nrow(x),]))
          for(i in 1:length(x2)){
            x1[i]=sum(x2[i:length(x2)])
          }
        }
        print(max(x1))
        if(max(x1)>threshold){
          false=false+1
          break
        }
        if(max(x1)<0){
          break
        }
        x=rbind(x,mvrnorm(n=1,mu=mu0,Sigma = covar))
        c=c+1
        c1=c1+1
        print(c) 
        if(c>=49){
          break
        }
        if(false>=10){
          break
        }
      }
      if(max(x1)>threshold){
        break
      }
      if(c>=49){
        break
      }
      if(false>=10){
        break
      }
    }
    if(c>=49){
      break
    }
    if(false>=10){
      break
    }
  }
  if(false<10){
    nodetect=0
    repeat{
      x=rbind(x,mvrnorm(n=1,mu=mu1,Sigma = covar))
      c=c+1
      x2[length(x2)+1]=log(f1(x[nrow(x),])/f0(x[nrow(x),]))
      for(i in 1:length(x2)){
        x1[i]=sum(x2[i:length(x2)])
      }
      if(max(x1)>threshold){
        break
      }
      if(max(x1)<0){
        break
      }
      if(c>1000){
        nodetect=1
        break
      }
    }
    
    if(max(x1)<=threshold){
      repeat{
        x=mvrnorm(n=1,mu=mu1,Sigma = covar)
        c=c+1
        print(c)
        x1=numeric()
        x2=numeric()
        c1=0
        repeat{
          if(c1==0){
            x2[1]=log(f1(x)/f0(x))
            x1[1]=x2[1]
          }
          print(max(x1))
          if(max(x1)>threshold){
            break
          }
          if(max(x1)<0){
            break
          }
          if(c1>=1){
            x2[length(x2)+1]=log(f1(x[nrow(x),])/f0(x[nrow(x),]))
            for(i in 1:length(x2)){
              x1[i]=sum(x2[i:length(x2)])
            }
          }
          print(max(x1))
          if(max(x1)>threshold){
            break
          }
          if(max(x1)<0){
            break
          }
          x=rbind(x,mvrnorm(n=1,mu=mu1,Sigma = covar))
          c=c+1
          c1=c1+1
          print(c) 
          if(c>1000){
            nodetect=1
            break
          }
        }
        if(max(x1)>threshold){
          break
        }
        if(c>1000){
          nodetect=1
          break
        }
      }
    }
    return(c(c-50,nodetect))
  }
  if(false>=10){
    return("false")
  }
}

opt1_cusum=function(d=4,threshold=1,mu0=rep(0,times=d),mu1=rep(0.5,times=d),covar=diag(1,d,d)){
  false=0
  c=0
  repeat{
    repeat{
      x=mvrnorm(n=1,mu=mu0,Sigma = covar)
      c=c+1
      print(c)
      x1=numeric()
      x2=numeric()
      c1=0
      repeat{
        if(c1==0){
          x2[1]=(t(x-mu0)%*%solve(covar)%*%(x-mu0)-t(x-mu1)%*%solve(covar)%*%(x-mu1))*0.5
          x1[1]=x2[1]
        }
        print(max(x1))
        if(max(x1)>threshold){
          false=false+1
          break
        }
        if(max(x1)<0){
          break
        }
        if(c1>=1){
          x2[length(x2)+1]=(t(x[nrow(x),]-mu0)%*%solve(covar)%*%(x[nrow(x),]-mu0)-t(x[nrow(x),]-mu1)%*%solve(covar)%*%(x[nrow(x),]-mu1))*0.5
          for(i in 1:length(x2)){
            x1[i]=sum(x2[i:length(x2)])
          }
        }
        print(max(x1))
        if(max(x1)>threshold){
          false=false+1
          break
        }
        if(max(x1)<0){
          break
        }
        x=rbind(x,mvrnorm(n=1,mu=mu0,Sigma = covar))
        c=c+1
        c1=c1+1
        print(c) 
        if(c>=49){
          break
        }
        if(false>=10){
          break
        }
      }
      if(max(x1)>threshold){
        break
      }
      if(c>=49){
        break
      }
      if(false>=10){
        break
      }
    }
    if(c>=49){
      break
    }
    if(false>=10){
      break
    }
  }
  if(false<10){
    nodetect=0
    repeat{
      x=rbind(x,mvrnorm(n=1,mu=mu1,Sigma = covar))
      c=c+1
      x2[length(x2)+1]=(t(x[nrow(x),]-mu0)%*%solve(covar)%*%(x[nrow(x),]-mu0)-t(x[nrow(x),]-mu1)%*%solve(covar)%*%(x[nrow(x),]-mu1))*0.5
      for(i in 1:length(x2)){
        x1[i]=sum(x2[i:length(x2)])
      }
      if(max(x1)>threshold){
        break
      }
      if(max(x1)<0){
        break
      }
      if(c>1000){
        nodetect=1
        break
      }
    }
    
    if(max(x1)<=threshold){
      repeat{
        x=mvrnorm(n=1,mu=mu1,Sigma = covar)
        c=c+1
        print(c)
        x1=numeric()
        x2=numeric()
        c1=0
        repeat{
          if(c1==0){
            x2[1]=(t(x-mu0)%*%solve(covar)%*%(x-mu0)-t(x-mu1)%*%solve(covar)%*%(x-mu1))*0.5
            x1[1]=x2[1]
          }
          print(max(x1))
          if(max(x1)>threshold){
            break
          }
          if(max(x1)<0){
            break
          }
          if(c1>=1){
            x2[length(x2)+1]=(t(x[nrow(x),]-mu0)%*%solve(covar)%*%(x[nrow(x),]-mu0)-t(x[nrow(x),]-mu1)%*%solve(covar)%*%(x[nrow(x),]-mu1))*0.5
            for(i in 1:length(x2)){
              x1[i]=sum(x2[i:length(x2)])
            }
          }
          print(max(x1))
          if(max(x1)>threshold){
            break
          }
          if(max(x1)<0){
            break
          }
          x=rbind(x,mvrnorm(n=1,mu=mu1,Sigma = covar))
          c=c+1
          c1=c1+1
          print(c)
          if(c>1000){
            nodetect=1
            break
          }
        }
        if(max(x1)>threshold){
          break
        }
        if(c>1000){
          nodetect=1
          break
        }
      }
    }
    return(c(c-50,nodetect))
  }
  if(false>=10){
    return("false")
  }
}


