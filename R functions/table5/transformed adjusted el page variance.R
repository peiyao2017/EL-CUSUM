


library(el.convex)
tael_var=function(m=100,d=4,threshold=200,min0=-3,max0=3,min1=-3.5,max1=3.5){
  x=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x[i,]=runif(n=d,min=min0,max=max0)
  }
  
  x5=apply(x,MARGIN =2,var)
  x51=colMeans(x)
  x7=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x7[i,]=runif(n=d,min=min1,max=max1)
  }
  x81=colMeans(x7)
  x8=apply(x7,MARGIN =2,var)
  c=d+1
  if(c<=100){
    x9=matrix(0,c,d)
    h0=matrix(0,c,d)
    h1=matrix(0,c,d)
    for(i in 1:c){
       x9[i,]=runif(n=d,min=min0,max=max0)
    }
    for(i in 1:c){
      h0[i,]=(x9[i,]-x51)^2-x5
      h1[i,]=(x9[i,]-x81)^2-x8
    }
    z=numeric()
    repeat{
      for(i in 1:(c-(d-1))){
        x10=h0[i:c,]
        x11=h1[i:c,]
        x10s=colMeans(x10)
        x11s=colMeans(x11)
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,d))$`-2LLR`
        l1=el.test.newton(x=x11,mu=rep(0,d))$`-2LLR`
        tl0=l0*max(1-l0/nrow(x10),0.5)
        tl1=l1*max(1-l1/nrow(x11),0.5)
        z[i]=0.5*(tl0-tl1)
      }
      if(max(z)>threshold){
        break
      }
      x9=rbind(x9,runif(n=d,min=min0,max=max0))
      c=c+1
      print(c)
      h0=rbind(h0,(x9[c,]-x51)^2-x5)
      h1=rbind(h1,(x9[c,]-x81)^2-x8)
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
      for(i in 1:(nrow(current0)-(d-1))){
        x10=current0[i:nrow(current0),]
        x11=current1[i:nrow(current0),]
        x10s=colMeans(x10)
        x11s=colMeans(x11)
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,d))$`-2LLR`
        l1=el.test.newton(x=x11,mu=rep(0,d))$`-2LLR`
        tl0=l0*max(1-l0/nrow(x10),0.5)
        tl1=l1*max(1-l1/nrow(x11),0.5)
        z[i]=0.5*(tl0-tl1)
      }
      print(max(z))
      if(max(z)>threshold){
        break
      }
      x9=rbind(x9,runif(n=d,min=min0,max=max0))
      c=c+1
      print(c)
      h0=rbind(h0,(x9[c,]-x51)^2-x5)
      h1=rbind(h1,(x9[c,]-x81)^2-x8)
      if(c>=1000){
        break
      }
    }
  }
  return(c)
}



tael1_var=function(m=20,d=4,threshold=200,min0=-3,max0=3,min1=-3.5,max1=3.5){
  x=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x[i,]=runif(n=d,min=min0,max=max0)
  }
  
  x5=apply(x,MARGIN =2,var)
  x51=colMeans(x)
  x7=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x7[i,]=runif(n=d,min=min1,max=max1)
  }
  x81=colMeans(x7)
  x8=apply(x7,MARGIN =2,var)
  repeat{
    c=d+1
    repeat{
      x9=matrix(0,c,d)
      h0=matrix(0,c,d)
      h1=matrix(0,c,d)
      for(i in 1:c){
        x9[i,]=runif(n=d,min=min0,max=max0)
      }
      for(i in 1:c){
        h0[i,]=(x9[i,]-x51)^2-x5
        h1[i,]=(x9[i,]-x81)^2-x8
      }
      z=numeric()
      repeat{
        for(i in 1:(c-(d-1))){
          x10=h0[i:c,]
          x11=h1[i:c,]
          x10s=colMeans(x10)
          x11s=colMeans(x11)
          x10=rbind(x10,-max(1,log(c)/2)*x10s)
          x11=rbind(x11,-max(1,log(c)/2)*x11s)
          l0=el.test.newton(x=x10,mu=rep(0,d))$`-2LLR`
          l1=el.test.newton(x=x11,mu=rep(0,d))$`-2LLR`
          tl0=l0*max(1-l0/nrow(x10),0.5)
          tl1=l1*max(1-l1/nrow(x11),0.5)
          z[i]=0.5*(tl0-tl1)
        }
        if(max(z)>threshold){
          break
        }
        c=c+1
        print(c)
        if(c>=50){
          break
        }
        x9=rbind(x9,runif(n=d,min=min0,max=max0))
        h0=rbind(h0,(x9[c,]-x51)^2-x5)
        h1=rbind(h1,(x9[c,]-x81)^2-x8)
      }
      if(max(z)>threshold){
        break
      }
      if(c>=50){
        break
      }
    }
    if(c>=50){
      break
    }
  }
  nodetect=0
  repeat{
    x9=rbind(x9,runif(n=d,min=min1,max=max1))
    h0=rbind(h0,(x9[c,]-x51)^2-x5)
    h1=rbind(h1,(x9[c,]-x81)^2-x8)
    z=numeric()
    for(i in 1:(c-(d-1))){
      x10=h0[i:c,]
      x11=h1[i:c,]
      x10s=colMeans(x10)
      x11s=colMeans(x11)
      x10=rbind(x10,-max(1,log(c)/2)*x10s)
      x11=rbind(x11,-max(1,log(c)/2)*x11s)
      l0=el.test.newton(x=x10,mu=rep(0,d))$`-2LLR`
      l1=el.test.newton(x=x11,mu=rep(0,d))$`-2LLR`
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
      x9=rbind(x9,runif(n=d,min=min1,max=max1))
      h0=rbind(h0,(x9[c,]-x51)^2-x5)
      h1=rbind(h1,(x9[c,]-x81)^2-x8)
      current0=h0[(nrow(h0)-99):nrow(h0),]
      current1=h1[(nrow(h1)-99):nrow(h1),]
      z=numeric()
      for(i in 1:(nrow(current0)-(d-1))){
        x10=current0[i:nrow(current0),]
        x11=current1[i:nrow(current0),]
        x10s=colMeans(x10)
        x11s=colMeans(x11)
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,d))$`-2LLR`
        l1=el.test.newton(x=x11,mu=rep(0,d))$`-2LLR`
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
