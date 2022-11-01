install.packages('doParallel')
library(doParallel)

myCluster <- makeCluster(2, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)

 

library(el.convex)

 
tel_laplace=function(m=20,d=2,mu0=rep(0,times=d),mu1=rep(0.5,times=d),threshold=d){
  x=mvrnorm(n = m, mu=mu0, Sigma=diag(1,d,d))
  x1=matrix(0,ncol=1,nrow=m)
  x2=matrix(0,ncol=1,nrow=m)
  y1=rep(0.5,d)
  y2=rep(-0.5,d)
  for(i in 1:m){
    x1[i,]=exp(-y1%*%x[i,])
    x2[i,]=exp(-y2%*%x[i,])
  }
  x3=colMeans(x1)
  x4=colMeans(x2)
  x5=cbind(x3,x4)
  x6=mvrnorm(n = m, mu=mu1, Sigma=diag(1,d,d))
  x7=matrix(0,ncol=2,nrow=m)
  for(i in 1:m){
    x7[i,]=c(exp(-y1%*%x6[i,]),exp(-y2%*%x6[i,]))
  }
  x8=colMeans(x7)
  c=d+1
  x9=matrix(0,c,d)
  h0=matrix(0,c,2)
  h1=matrix(0,c,2)
  for(i in 1:c){
    x9[i,]=mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d))
  }
  for(i in 1:c){
    h0[i,]=c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x5
    h1[i,]=c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x8
  }
  z=numeric()
  repeat{
    for(i in 1:(c-d)){
      x10=h0[i:c,]
      x11=h1[i:c,]
      l0=el.test.newton(x=x10,mu=rep(0,2)) 
      l1=el.test.newton(x=x11,mu=rep(0,2)) 
      if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
        tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
        tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
        z[i]=0.5*(tl0-tl1)
      }
    }
    z=na.omit(z)
    if(length(z)==0){
      x9=rbind(x9,mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d)))
      c=c+1
      print(c)
      h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
      h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
    }
    if(length(z)>0){
      if(max(z)<=threshold){
        x9=rbind(x9,mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d)))
        c=c+1
        print(c)
        h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
        h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
      }
    }
    if(length(z)>0){
      if(max(z)>threshold){
        break
      }
    }
    if(c>100){
      break
    }
  }
  if(length(z)==0){
    z=c(0)
  }
  if(max(z)<=threshold){
    repeat{
      current0=h0[(nrow(h0)-99):nrow(h0),]
      current1=h1[(nrow(h1)-99):nrow(h1),]
      z=numeric()
      for(i in 1:(nrow(current0)-d)){
        x10=current0[i:nrow(current0),]
        x11=current1[i:nrow(current0),]
        l0=el.test.newton(x=x10,mu=rep(0,2)) 
        l1=el.test.newton(x=x11,mu=rep(0,2))
        if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
          tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
          tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
          z[i]=0.5*(tl0-tl1)
        }
      }
      z=na.omit(z)
      if(length(z)==0){
        x9=rbind(x9,mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d)))
        c=c+1
        print(c)
        h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
        h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
      }
      if(length(z)>0){
        if(max(z)>threshold){
          print(max(z))
          break
        }
      }
      if(length(z)>0){
        if(max(z)<=threshold){
          x9=rbind(x9,mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d)))
          c=c+1
          print(c)
          h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
          h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
        }
      }
      if(c>=1000){
        break
      }
    }
  }
  return(c)
}




tel1_laplace=function(m=20,d=2,mu0=rep(0,times=d),mu1=rep(0.5,times=d),threshold=2){
  false=0
  x=mvrnorm(n = m, mu=mu0, Sigma=diag(1,d,d))
  x1=matrix(0,ncol=1,nrow=m)
  x2=matrix(0,ncol=1,nrow=m)
  y1=rep(0.4,times=d)
  y2=rep(0.8,times=d)
  for(i in 1:m){
    x1[i,]=exp(-y1%*%x[i,])
    x2[i,]=exp(-y2%*%x[i,])
  }
  x3=colMeans(x1)
  x4=colMeans(x2)
  x5=cbind(x3,x4)
  x6=mvrnorm(n = m, mu=mu1, Sigma=diag(1,d,d))
  x7=matrix(0,ncol=2,nrow=m)
  for(i in 1:m){
    x7[i,]=c(exp(-y1%*%x6[i,]),exp(-y2%*%x6[i,]))
  }
  x8=colMeans(x7)
  repeat{
    c=d+1
    repeat{
      x9=matrix(0,c,d)
      h0=matrix(0,c,2)
      h1=matrix(0,c,2)
      for(i in 1:c){
        x9[i,]=mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d))
      }
      for(i in 1:c){
        h0[i,]=c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x5
        h1[i,]=c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x8
      }
      z=numeric()
      repeat{
        for(i in 1:(c-d)){
          x10=h0[i:c,]
          x11=h1[i:c,]
          l0=el.test.newton(x=x10,mu=rep(0,2)) 
          l1=el.test.newton(x=x11,mu=rep(0,2)) 
          if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
            tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
            tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
            z[i]=0.5*(tl0-tl1)
          }
        }
        z=na.omit(z)
        if(length(z)==0){
          x9=rbind(x9,mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d)))
          c=c+1
          print(c)
          h0=rbind(h0,c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x5)
          h1=rbind(h1,c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x8)
        }
        if(length(z)>0){
          if(max(z)>threshold){
            false=false+1
            break
          }
        }
        if(length(z)>0){
          if(max(z)<=threshold){
            x9=rbind(x9,mvrnorm(n=1,mu=mu0,Sigma = diag(1,d,d)))
            c=c+1
            print(c)
            h0=rbind(h0,c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x5)
            h1=rbind(h1,c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x8)
          }
        }
        if(c>=49){
          break
        }
        if(false>=10){
          break
        }
      }
      if(max(z)>threshold){
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
    x9=rbind(x9,mvrnorm(n=1,mu=mu1,Sigma = diag(1,d,d)))
    c=c+1
    h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
    h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
    z=numeric()
    for(i in 1:(c-d)){
      x10=h0[i:c,]
      x11=h1[i:c,]
      l0=el.test.newton(x=x10,mu=rep(0,2)) 
      l1=el.test.newton(x=x11,mu=rep(0,2)) 
      if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
        tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
        tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
        z[i]=0.5*(tl0-tl1)
      }
    }
    z=na.omit(z)
    if(length(z)==0){
      x9=rbind(x9,mvrnorm(n=1,mu=mu1,Sigma = diag(1,d,d)))
      c=c+1
      print(c)
      h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
      h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
    }
    if(length(z)>0){
      if(max(z)>threshold){
        break
      }
    }
    if(length(z)>0){
      if(max(z)<=threshold){
        x9=rbind(x9,mvrnorm(n=1,mu=mu1,Sigma = diag(1,d,d)))
        c=c+1
        print(c)
        h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
        h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
      }
    }
    if(c>=101){
      break
    }
  }
  if(length(z)==0){
    z=c(0)
  }
  if(max(z)<=threshold){
    repeat{
      z=numeric()
      x9=rbind(x9,mvrnorm(n=1,mu=mu1,Sigma = diag(1,d,d)))
      h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
      h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
      current0=h0[(nrow(h0)-99):nrow(h0),]
      current1=h1[(nrow(h1)-99):nrow(h1),]
      z=numeric()
      for(i in 1:(nrow(current0)-d)){
        x10=current0[i:nrow(current0),]
        x11=current1[i:nrow(current0),]
        l0=el.test.newton(x=x10,mu=rep(0,2)) 
        l1=el.test.newton(x=x11,mu=rep(0,2)) 
        if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
          tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
          tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
          z[i]=0.5*(tl0-tl1)
        }
      }
      z=na.omit(z)
      if(length(z)==0){
        x9=rbind(x9,mvrnorm(n=1,mu=mu1,Sigma = diag(1,d,d)))
        c=c+1
        print(c)
        h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
        h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
      }
      if(length(z)>0){
        if(max(z)>threshold){
          break
        }
      }
      if(length(z)>0){
        if(max(z)<=threshold){
          x9=rbind(x9,mvrnorm(n=1,mu=mu1,Sigma = diag(1,d,d)))
          c=c+1
          print(c)
          h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
          h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
        }
      }
      if(c>1000){
        nodetect=nodetect+1
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



