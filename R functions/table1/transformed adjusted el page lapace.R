


library(el.convex)

 

tael_laplace=function(m=100,mu0=0,mu1=0.5,threshold=200){
  x=rnorm(n = m, mean=mu0,sd=1)
  x1=rep(0,times=m)
  x2=rep(0,times=m)
  y1=0.5
  y2=-0.5
  for(i in 1:m){
    x1[i]=exp(-y1*x[i])
    x2[i]=exp(-y2*x[i])
  }
  x3=mean(x1)
  x4=mean(x2)
  x5=c(x3,x4)
  x6=rnorm(n = m, mean=mu1,sd=1)
  x7=matrix(0,ncol=2,nrow=m)
  for(i in 1:m){
    x7[i,]=c(exp(-y1*x6[i]),exp(-y2*x6[i]))
  }
  x8=colMeans(x7)
  c=3
  if(c<=100){
    x9=numeric()
    h0=matrix(0,c,2)
    h1=matrix(0,c,2)
    for(i in 1:c){
      x9[i]=rnorm(n=1,mean=mu0,sd=1)
    }
    for(i in 1:c){
      h0[i,]=c(exp(-y1*x9[i]),exp(-y2*x9[i]))-x5
      h1[i,]=c(exp(-y1*x9[i]),exp(-y2*x9[i]))-x8
    }
    z=numeric()
    repeat{
      for(i in 1:(c-1)){
        x10=h0[i:c,]
        x11=h1[i:c,]
        if(i==c){
          x10s=x10
          x11s=x11
        }else{
          x10s=colMeans(x10)
          x11s=colMeans(x11)
        }
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,2))$`-2LLR`
        l1=el.test.newton(x=x11,mu=rep(0,2))$`-2LLR`
        tl0=l0*max(1-l0/nrow(x10),0.5)
        tl1=l1*max(1-l1/nrow(x11),0.5)
        z[i]=0.5*(tl0-tl1)
      }
      if(max(z)>threshold){
        break
      }
      x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
      c=c+1
      print(c)
      h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
      h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
      if(c>100){
        break
      }
    }
  }
  if(max(z)<=threshold){
    repeat{
      current0=h0[(nrow(h0)-99):nrow(h0),]
      current1=h1[(nrow(h1)-99):nrow(h1),]
      z=numeric()
      for(i in 1:(nrow(current0)-1)){
        x10=current0[i:nrow(current0),]
        x11=current1[i:nrow(current0),]
        if(i==nrow(current0)){
          x10s=x10
          x11s=x11
        }else{
          x10s=colMeans(x10)
          x11s=colMeans(x11)
        }
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,2))$`-2LLR`
        l1=el.test.newton(x=x11,mu=rep(0,2))$`-2LLR`
        tl0=l0*max(1-l0/nrow(x10),0.5)
        tl1=l1*max(1-l1/nrow(x11),0.5)
        z[i]=0.5*(tl0-tl1)
      }
      print(max(z))
      if(max(z)>threshold){
        break
      }
      x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
      c=c+1
      print(c)
      h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
      h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
      if(c>1000){
        break
      }
    }
    
  }
  return(c)
}




tael1_laplace=function(m=20,mu0=0,mu1=0.5,threshold=200){
  false=0
  x=rnorm(n = m, mean=mu0,sd=1)
  x1=rep(0,times=m)
  x2=rep(0,times=m)
  y1=0.5
  y2=-0.5
  for(i in 1:m){
    x1[i]=exp(-y1*x[i])
    x2[i]=exp(-y2*x[i])
  }
  x3=mean(x1)
  x4=mean(x2)
  x5=c(x3,x4)
  x6=rnorm(n = m, mean=mu1,sd=1)
  x7=matrix(0,ncol=2,nrow=m)
  for(i in 1:m){
    x7[i,]=c(exp(-y1*x6[i]),exp(-y2*x6[i]))
  }
  x8=colMeans(x7)
  repeat{
    c=3
    repeat{
      x9=numeric()
      h0=matrix(0,c,2)
      h1=matrix(0,c,2)
      for(i in 1:c){
        x9[i]=rnorm(n=1,mean=mu0,sd=1)
      }
      for(i in 1:c){
        h0[i,]=c(exp(-y1*x9[i]),exp(-y2*x9[i]))-x5
        h1[i,]=c(exp(-y1*x9[i]),exp(-y2*x9[i]))-x8
      }
      z=numeric()
      repeat{
        for(i in 1:(c-1)){
          x10=h0[i:c,]
          x11=h1[i:c,]
          x10s=colMeans(x10)
          x11s=colMeans(x11)
          x10=rbind(x10,-max(1,log(c)/2)*x10s)
          x11=rbind(x11,-max(1,log(c)/2)*x11s)
          l0=el.test.newton(x=x10,mu=rep(0,2))$`-2LLR`
          l1=el.test.newton(x=x11,mu=rep(0,2))$`-2LLR`
          tl0=l0*max(1-l0/nrow(x10),0.5)
          tl1=l1*max(1-l1/nrow(x11),0.5)
          z[i]=0.5*(tl0-tl1)
        }
        if(max(z)>threshold){
          false=false+1
          break
        }
        c=c+1
        if(c>=50){
          break
        }
        if(false>=10){
          break
        }
        x9=c(x9,rnorm(n=1,mean=mu0,sd=1))
        h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
        h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
      }
      if(max(z)>threshold){
        break
      }
      if(c>=50){
        break
      }
      if(false>=10){
        break
      }
    }
    if(c>=50){
      break
    }
    if(false>=10){
      break
    }
  }
  if(false<10){
  nodetect=0
  repeat{
    x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
    h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
    h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
    z=numeric()
    for(i in 1:(c-1)){
      x10=h0[i:c,]
      x11=h1[i:c,]
      x10s=colMeans(x10)
      x11s=colMeans(x11)
      x10=rbind(x10,-max(1,log(c)/2)*x10s)
      x11=rbind(x11,-max(1,log(c)/2)*x11s)
      l0=el.test.newton(x=x10,mu=rep(0,2))$`-2LLR`
      l1=el.test.newton(x=x11,mu=rep(0,2))$`-2LLR`
      tl0=l0*max(1-l0/nrow(x10),0.5)
      tl1=l1*max(1-l1/nrow(x11),0.5)
      z[i]=0.5*(tl0-tl1)
    }
    c=c+1
    if(max(z)>threshold){
      break
    }
    if(c>=101){
      break
    }
  }
  if(max(z)<=threshold){
    repeat{
      x9=c(x9,rnorm(n=1,mean=mu1,sd=1))
      h0=rbind(h0,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x5)
      h1=rbind(h1,c(exp(-y1*x9[c]),exp(-y2*x9[c]))-x8)
      current0=h0[(nrow(h0)-99):nrow(h0),]
      current1=h1[(nrow(h1)-99):nrow(h1),]
      z=numeric()
      for(i in 1:(nrow(current0)-1)){
        x10=current0[i:nrow(current0),]
        x11=current1[i:nrow(current0),]
        x10s=colMeans(x10)
        x11s=colMeans(x11)
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,2))$`-2LLR`
        l1=el.test.newton(x=x11,mu=rep(0,2))$`-2LLR`
        tl0=l0*max(1-l0/nrow(x10),0.5)
        tl1=l1*max(1-l1/nrow(x11),0.5)
        z[i]=0.5*(tl0-tl1)
      }
      print(max(z))
      if(max(z)>threshold){
        break
      }
      c=c+1
      print(c)
      if(max(z)>threshold){
        break
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

