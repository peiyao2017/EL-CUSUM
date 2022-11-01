
library(el.convex)


el_laplace=function(m=20,d=2,threshod=200){
  x=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x[i,]=rgamma(n=d,shape=10,rate=2)
  }
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
  x7=matrix(0,nrow=m,ncol=2)
  x71=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x71[i,]=rgamma(n=d,shape=5,rate=1)
  }
  for(i in 1:m){
    x7[i,]=c(exp(-y1%*%x71[i,]),exp(-y2%*%x71[i,]))
  }
  x8=colMeans(x7)
  c=2+1
  x9=matrix(0,c,d)
  h0=matrix(0,c,2)
  h1=matrix(0,c,2)
  for(i in 1:c){
    x9[i,]=rgamma(n=d,shape=10,rate = 2)
  }
  for(i in 1:c){
    h0[i,]=c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x5
    h1[i,]=c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x8
  }
  z=numeric()
  repeat{
    for(i in 1:(c-2)){
      x10=h0[i:c,]
      x11=h1[i:c,]
      l0=el.test.newton(x=x10,mu=rep(0,2))
      l1=el.test.newton(x=x11,mu=rep(0,2))
      if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
        z[i]=0.5*(l0$`-2LLR`-l1$`-2LLR`)
      }
    }
    z=na.omit(z)
    if(length(z)==0){
      x9=rbind(x9,rgamma(n=d,shape=10,rate = 2))
      c=c+1
      print(c)
      h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
      h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
    }
    if(length(z)>0){
      if(max(z)<=threshod){
        x9=rbind(x9,rgamma(n=d,shape=10,rate = 2))
        c=c+1
        print(c)
        h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
        h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
      }
    }
    if(length(z)>0){
      if(max(z)>threshod){
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
  if(max(z)<=threshod){
    repeat{
      current0=h0[(nrow(h0)-99):nrow(h0),]
      current1=h1[(nrow(h1)-99):nrow(h1),]
      z=numeric()
      for(i in 1:(nrow(current0)-2)){
        x10=current0[i:nrow(current0),]
        x11=current1[i:nrow(current0),]
        l0=el.test.newton(x=x10,mu=rep(0,2))
        l1=el.test.newton(x=x11,mu=rep(0,2))
        if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
          z[i]=0.5*(l0$`-2LLR`-l1$`-2LLR`)
        }
      }
      z=na.omit(z)
      if(length(z)==0){
        x9=rbind(x9,rgamma(n=d,shape=10,rate = 2))
        c=c+1
        print(c)
        h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
        h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
      }
      if(length(z)>0){
         print(max(z))
        if(max(z)>threshod){
          break
        }
      }
      if(length(z)>0){
        print(max(z))
        if(max(z)<=threshod){
          x9=rbind(x9,rgamma(n=d,shape=10,rate = 2))
          c=c+1
          print(c)
          h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
          h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
        }
      }
      if(c>1000){
        break
      }
    }
  }
  return(c)
}



 

el1_laplace=function(m=20,d=2,threshold=2.1){
  false=0
  x=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x[i,]=rgamma(n=d,shape=10,rate=2)
  }
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
  x7=matrix(0,nrow=m,ncol=2)
  x71=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x71[i,]=rgamma(n=d,shape=5,rate=1)
  }
  for(i in 1:m){
    x7[i,]=c(exp(-y1%*%x71[i,]),exp(-y2%*%x71[i,]))
  }
  x8=colMeans(x7)
  repeat{
    c=2+1
    repeat{
      x9=matrix(0,c,d)
      h0=matrix(0,c,2)
      h1=matrix(0,c,2)
      for(i in 1:c){
        x9[i,]=rgamma(n=d,shape=10,rate = 2)
      }
      for(i in 1:c){
        h0[i,]=c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x5
        h1[i,]=c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x8
      }
      z=numeric()
      repeat{
        for(i in 1:(c-2)){
          x10=h0[i:c,]
          x11=h1[i:c,]
          l0=el.test.newton(x=x10,mu=rep(0,2))
          l1=el.test.newton(x=x11,mu=rep(0,2))
          if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
            z[i]=0.5*(l0$`-2LLR`-l1$`-2LLR`)
          }
        }
        z=na.omit(z)
        if(length(z)==0){
          x9=rbind(x9,rgamma(n=d,shape=10,rate = 2))
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
            x9=rbind(x9,rgamma(n=d,shape=10,rate = 2))
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
    x9=rbind(x9,rgamma(n=d,shape=5,rate = 1))
    c=c+1
    h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
    h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
    z=numeric()
    for(i in 1:(c-2)){
      x10=h0[i:c,]
      x11=h1[i:c,]
      l0=el.test.newton(x=x10,mu=rep(0,2))
      l1=el.test.newton(x=x11,mu=rep(0,2))
      if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
        z[i]=0.5*(l0$`-2LLR`-l1$`-2LLR`)
      }
    }
    z=na.omit(z)
    if(length(z)==0){
      x9=rbind(x9,rgamma(n=d,shape=5,rate = 1))
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
        x9=rbind(x9,rgamma(n=d,shape=5,rate = 1))
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
      x9=rbind(x9,rgamma(n=d,shape=5,rate = 1))
      h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
      h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
      current0=h0[(nrow(h0)-99):nrow(h0),]
      current1=h1[(nrow(h1)-99):nrow(h1),]
      z=numeric()
      for(i in 1:(nrow(current0)-2)){
        x10=current0[i:nrow(current0),]
        x11=current1[i:nrow(current0),]
        l0=el.test.newton(x=x10,mu=rep(0,2))
        l1=el.test.newton(x=x11,mu=rep(0,2))
        if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
          z[i]=0.5*(l0$`-2LLR`-l1$`-2LLR`)
        }
      }
      z=na.omit(z)
      if(length(z)==0){
        x9=rbind(x9,rgamma(n=d,shape=5,rate = 1))
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
        print(max(z))
        if(max(z)<=threshold){
          x9=rbind(x9,rgamma(n=d,shape=5,rate = 1))
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
