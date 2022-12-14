arl1=numeric()
arl2=numeric()
arl3=numeric()
arl4=numeric()
arl5=numeric()
arl6=numeric()
arl7=numeric()
arl8=numeric()

repeat{
  a=foreach(i=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    el1_var(threshold =20,m=20)
  }
  arl1=c(arl1,a[a!="false"])
  if(length(arl1)>2000){
    break
  }
}

repeat{
  a=foreach(i=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    el1_laplace(threshold =20,m=20)
  }
  arl2=c(arl2,a[a!="false"])
  if(length(arl2)>2000){
    break
  }
}

repeat{
  a=foreach(i=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    ael1_var(threshold =20,m=20)
  }
  arl3=c(arl3,a[a!="false"])
  if(length(arl3)>2000){
    break
  }
}


repeat{
  a=foreach(i=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    ael1_laplace(threshold =20,m=20)
  }
  arl4=c(arl4,a[a!="false"])
  if(length(arl4)>2000){
    break
  }
}


repeat{
  a=foreach(i=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    tel1_var(threshold =20,m=20)
  }
  arl5=c(arl5,a[a!="false"])
  if(length(arl5)>2000){
    break
  }
}

repeat{
  a=foreach(i=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    tel1_laplace(threshold =20,m=20)
  }
  arl6=c(arl6,a[a!="false"])
  if(length(arl6)>2000){
    break
  }
}

repeat{
  a=foreach(i=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    tael1_var(threshold =20,m=20)
  }
  arl7=c(arl7,a[a!="false"])
  if(length(arl7)>2000){
    break
  }
}

repeat{
  a=foreach(i=1:1000, .combine="c")%dopar%{
    library(el.convex)
    library(MASS)
    tael1_laplace(threshold =20,m=20)
  }
  arl8=c(arl8,a[a!="false"])
  if(length(arl8)>2000){
    break
  }
}

x=data.frame(arl1=arl1[1:2000],arl2=arl2[1:2000],arl3=arl3[1:2000],arl4=arl4[1:2000],arl5=arl5[1:2000],arl6=arl6[1:2000],arl7=arl7[1:2000],arl8=arl8[1:2000])

for(i in 1:8){
  x[,i]=as.numeric(x[,i])
}




delay=matrix(0,nrow=8,ncol=1000)
nodetect=matrix(0,nrow=8,ncol=1000)
for(i in 1:8){
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

mean(delay[4,])
sd(delay[4,])/sqrt(1000)
sum(nodetect[4,])


mean(delay[5,])
sd(delay[5,])/sqrt(1000)
sum(nodetect[5,])

mean(delay[6,])
sd(delay[6,])/sqrt(1000)
sum(nodetect[6,])


mean(delay[7,])
sd(delay[7,])/sqrt(1000)
sum(nodetect[7,])

mean(delay[8,])
sd(delay[8,])/sqrt(1000)
sum(nodetect[8,])



ael1_laplace=function(m=20,d=4,threshold=200,min0=-3,max0=3,min1=-3.5,max1=3.5){
  false=0
  x=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x[i,]=runif(n=d,min=min0,max=max0)
  }
  x1=matrix(0,ncol=1,nrow=m)
  x2=matrix(0,ncol=1,nrow=m)
  y1=rep(-0.5,times=d)
  y2=rep(0.5,times=d)
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
    x71[i,]=runif(n=d,min=min1,max=max1)
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
        x9[i,]=runif(n=d,min=min0,max=max0)
      }
      for(i in 1:c){
        h0[i,]=c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x5
        h1[i,]=c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x8
      }
      z=numeric()
      repeat{
        for(i in 1:(c-(2-1))){
          x10=h0[i:c,]
          x11=h1[i:c,]
          x10s=colMeans(x10)
          x11s=colMeans(x11)
          x10=rbind(x10,-max(1,log(c)/2)*x10s)
          x11=rbind(x11,-max(1,log(c)/2)*x11s)
          l0=el.test.newton(x=x10,mu=rep(0,2))
          l1=el.test.newton(x=x11,mu=rep(0,2))
          l=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
          z[i]=log(l)
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
        x9=rbind(x9,runif(n=d,min=min0,max=max0))
        h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
        h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
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
      x9=rbind(x9,runif(n=d,min=min1,max=max1))
      h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
      h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
      z=numeric()
      for(i in 1:(c-(2-1))){
        x10=h0[i:c,]
        x11=h1[i:c,]
        x10s=colMeans(x10)
        x11s=colMeans(x11)
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,2))
        l1=el.test.newton(x=x11,mu=rep(0,2))
        l=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
        z[i]=log(l)
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
        h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
        h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
        current0=h0[(nrow(h0)-99):nrow(h0),]
        current1=h1[(nrow(h1)-99):nrow(h1),]
        z=numeric()
        for(i in 1:(nrow(current0)-(2-1))){
          x10=current0[i:nrow(current0),]
          x11=current1[i:nrow(current0),]
          x10s=colMeans(x10)
          x11s=colMeans(x11)
          x10=rbind(x10,-max(1,log(c)/2)*x10s)
          x11=rbind(x11,-max(1,log(c)/2)*x11s)
          l0=el.test.newton(x=x10,mu=rep(0,2))
          l1=el.test.newton(x=x11,mu=rep(0,2))
          l=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
          z[i]=log(l)
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

ael1_var=function(m=20,d=4,min0=-3,max0=3,min1=-3.5,max1=3.5,threshold=200){
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
          l0=el.test.newton(x=x10,mu=rep(0,d))
          l1=el.test.newton(x=x11,mu=rep(0,d))
          l=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
          z[i]=log(l)
        }
        if(max(z)>threshold){
          break
        }
        c=c+1
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
      l0=el.test.newton(x=x10,mu=rep(0,d))
      l1=el.test.newton(x=x11,mu=rep(0,d))
      l=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
      z[i]=log(l)
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
        l0=el.test.newton(x=x10,mu=rep(0,d))
        l1=el.test.newton(x=x11,mu=rep(0,d))
        l=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
        z[i]=log(l)
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

el1_laplace=function(m=20,d=4,threshold=2.1,min0=-3,max0=3,min1=-3.5,max1=3.5){
  false=0
  x=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x[i,]=runif(n=d,min=min0,max=max0)
  }
  x1=matrix(0,ncol=1,nrow=m)
  x2=matrix(0,ncol=1,nrow=m)
  y1=rep(-0.5,times=d)
  y2=rep(0.5,times=d)
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
    x71[i,]=runif(n=d,min=min1,max=max1)
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
        x9[i,]=runif(n=d,min=min0,max=max0)
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
          x9=rbind(x9,runif(n=d,min=min0,max=max0))
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
            x9=rbind(x9,runif(n=d,min=min0,max=max0))
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
      x9=rbind(x9,runif(n=d,min=min1,max=max1))
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
        x9=rbind(x9,runif(n=d,min=min1,max=max1))
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
          x9=rbind(x9,runif(n=d,min=min1,max=max1))
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
        x9=rbind(x9,runif(n=d,min=min1,max=max1))
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
          x9=rbind(x9,runif(n=d,min=min1,max=max1))
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
            x9=rbind(x9,runif(n=d,min=min1,max=max1))
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

el1_var=function(m=20,d=4,threshold=200,min0=-3,max0=3,min1=-3.5,max1=3.5){
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
      repeat{
        z=numeric()
        for(i in 1:(c-d)){
          x10=h0[i:c,]
          x11=h1[i:c,]
          l0=el.test.newton(x=x10,mu=rep(0,d))
          l1=el.test.newton(x=x11,mu=rep(0,d))
          l=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
          if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
            z[i]=log(l)
          }
        }
        z=na.omit(z)
        if(length(z)==0){
          x9=rbind(x9,runif(n=d,min=min0,max=max0))
          c=c+1
          print(c)
          h0=rbind(h0,(x9[c,]-x51)^2-x5)
          h1=rbind(h1,(x9[c,]-x81)^2-x8)
        }
        if(length(z)>0){
          if(max(z)>threshold){
            break
          }
        }
        if(length(z)>0){
          if(max(z)<=threshold){
            x9=rbind(x9,runif(n=d,min=min0,max=max0))
            c=c+1
            print(c)
            h0=rbind(h0,(x9[c,]-x51)^2-x5)
            h1=rbind(h1,(x9[c,]-x81)^2-x8)
          }
        }
        if(c>=49){
          break
        }
      }
      if(length(z)>0){
        if(max(z)>threshold){
          break
        }
      }
      if(c>=49){
        break
      }
    }
    if(c>=49){
      break
    }
  }
  nodetect=0
  repeat{
    z=numeric()
    x9=rbind(x9,runif(n=d,min=min1,max=max1))
    c=c+1
    print(c)
    h0=rbind(h0,(x9[c,]-x51)^2-x5)
    h1=rbind(h1,(x9[c,]-x81)^2-x8)
    z=numeric()
    for(i in 1:(c-d)){
      x10=h0[i:c,]
      x11=h1[i:c,]
      l0=el.test.newton(x=x10,mu=rep(0,d))
      l1=el.test.newton(x=x11,mu=rep(0,d))
      l=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
      if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
        z[i]=log(l)
      }
    }
    z=na.omit(z)
    if(length(z)==0){
      x9=rbind(x9,runif(n=d,min=min1,max=max1))
      c=c+1
      print(c)
      h0=rbind(h0,(x9[c,]-x51)^2-x5)
      h1=rbind(h1,(x9[c,]-x81)^2-x8)
    }
    if(length(z)>0){
      if(max(z)>threshold){
        break
      }
    }
    if(length(z)>0){
      if(max(z)<=threshold){
        x9=rbind(x9,runif(n=d,min=min1,max=max1))
        c=c+1
        print(c)
        h0=rbind(h0,(x9[c,]-x51)^2-x5)
        h1=rbind(h1,(x9[c,]-x81)^2-x8)
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
      x9=rbind(x9,runif(n=d,min=min1,max=max1))
      c=c+1
      print(c)
      h0=rbind(h0,(x9[c,]-x51)^2-x5)
      h1=rbind(h1,(x9[c,]-x81)^2-x8)
      current0=h0[(nrow(h0)-99):nrow(h0),]
      current1=h1[(nrow(h1)-99):nrow(h1),]
      z=numeric()
      for(i in 1:(nrow(current0)-d)){
        x10=current0[i:nrow(current0),]
        x11=current1[i:nrow(current0),]
        l0=el.test.newton(x=x10,mu=rep(0,d))
        l1=el.test.newton(x=x11,mu=rep(0,d))
        l=exp(-0.5*l1$`-2LLR`+0.5*l0$`-2LLR`)
        if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
          z[i]=log(l)
        }
      }
      z=na.omit(z)
      if(length(z)==0){
        x9=rbind(x9,runif(n=d,min=min1,max=max1))
        c=c+1
        print(c)
        h0=rbind(h0,(x9[c,]-x51)^2-x5)
        h1=rbind(h1,(x9[c,]-x81)^2-x8)
      }
      if(length(z)>0){
        if(max(z)>threshold){
          break
        }
      }
      if(length(z)>0){
        print(max(z))
        if(max(z)<=threshold){
          x9=rbind(x9,runif(n=d,min=min1,max=max1))
          c=c+1
          print(c)
          h0=rbind(h0,(x9[c,]-x51)^2-x5)
          h1=rbind(h1,(x9[c,]-x81)^2-x8)
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

tael1_laplace=function(m=20,d=4,threshold=200,min0=-3,max0=3,min1=-3.5,max1=3.5){
  false=0
  x=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x[i,]=runif(n=d,min=min0,max=max0)
  }
  x1=matrix(0,ncol=1,nrow=m)
  x2=matrix(0,ncol=1,nrow=m)
  y1=rep(-0.5,times=d)
  y2=rep(0.5,times=d)
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
    x71[i,]=runif(n=d,min=min1,max=max1)
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
        x9[i,]=runif(n=d,min=min0,max=max0)
      }
      for(i in 1:c){
        h0[i,]=c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x5
        h1[i,]=c(exp(-y1%*%x9[i,]),exp(-y2%*%x9[i,]))-x8
      }
      z=numeric()
      repeat{
        for(i in 1:(c-(2-1))){
          x10=h0[i:c,]
          x11=h1[i:c,]
          x10s=colMeans(x10)
          x11s=colMeans(x11)
          x10=rbind(x10,-max(1,log(c)/2)*x10s)
          x11=rbind(x11,-max(1,log(c)/2)*x11s)
          l0=el.test.newton(x=x10,mu=rep(0,2))
          l1=el.test.newton(x=x11,mu=rep(0,2))
          l0=l0$`-2LLR`
          l1=l1$`-2LLR`
          tl0=l0*max(1-l0/nrow(x10),0.5)
          tl1=l1*max(1-l1/nrow(x11),0.5)
          z[i]=0.5*(tl0-tl1)
        }
        if(max(z)>threshold){
          false=false+1
          break
        }
        c=c+1
        print(c)
        if(c>=50){
          break
        }
        x9=rbind(x9,runif(n=d,min=min0,max=max0))
        h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
        h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
        if(false>=10){
          break
        }
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
      x9=rbind(x9,runif(n=d,min=min1,max=max1) )
      h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
      h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
      z=numeric()
      for(i in 1:(c-(2-1))){
        x10=h0[i:c,]
        x11=h1[i:c,]
        x10s=colMeans(x10)
        x11s=colMeans(x11)
        x10=rbind(x10,-max(1,log(c)/2)*x10s)
        x11=rbind(x11,-max(1,log(c)/2)*x11s)
        l0=el.test.newton(x=x10,mu=rep(0,2))
        l1=el.test.newton(x=x11,mu=rep(0,2))
        l0=l0$`-2LLR`
        l1=l1$`-2LLR`
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
        h0=rbind(h0,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x5)
        h1=rbind(h1,c(exp(-y1%*%x9[c,]),exp(-y2%*%x9[c,]))-x8)
        current0=h0[(nrow(h0)-99):nrow(h0),]
        current1=h1[(nrow(h1)-99):nrow(h1),]
        z=numeric()
        for(i in 1:(nrow(current0)-(2-1))){
          x10=current0[i:nrow(current0),]
          x11=current1[i:nrow(current0),]
          x10s=colMeans(x10)
          x11s=colMeans(x11)
          x10=rbind(x10,-max(1,log(c)/2)*x10s)
          x11=rbind(x11,-max(1,log(c)/2)*x11s)
          l0=el.test.newton(x=x10,mu=rep(0,2))
          l1=el.test.newton(x=x11,mu=rep(0,2))
          l0=l0$`-2LLR`
          l1=l1$`-2LLR`
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


tel1_laplace=function(m=20,d=4,threshold=200,min0=-3,max0=3,min1=-3.5,max1=3.5){
  false=0
  x=matrix(0,nrow=m,ncol=d)
  for(i in 1:m){
    x[i,]=runif(n=d,min=min0,max=max0)
  }
  x1=matrix(0,ncol=1,nrow=m)
  x2=matrix(0,ncol=1,nrow=m)
  y1=rep(0.5,times=d)
  y2=rep(-0.5,times=d)
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
    x71[i,]=runif(n=d,min=min1,max=max1)
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
        x9[i,]=runif(n=d,min=min0,max=max0)
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
            tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
            tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
            z[i]=0.5*(tl0-tl1)
          }
        }
        z=na.omit(z)
        if(length(z)==0){
          x9=rbind(x9,runif(n=d,min=min0,max=max0) )
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
            x9=rbind(x9, runif(n=d,min=min0,max=max0))
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
      x9=rbind(x9, runif(n=d,min=min1,max=max1) )
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
          tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
          tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
          z[i]=0.5*(tl0-tl1)
        }
      }
      z=na.omit(z)
      if(length(z)==0){
        x9=rbind(x9,runif(n=d,min=min1,max=max1) )
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
          x9=rbind(x9, runif(n=d,min=min1,max=max1))
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
        x9=rbind(x9, runif(n=d,min=min1,max=max1))
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
            tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
            tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
            z[i]=0.5*(tl0-tl1)
          }
        }
        z=na.omit(z)
        if(length(z)==0){
          x9=rbind(x9, runif(n=d,min=min1,max=max1))
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
            x9=rbind(x9,runif(n=d,min=min1,max=max1))
            c=c+1
            print(c)
            print(max(z))
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


tel1_var=function(m=20,d=4,threshold=200,min0=-3,max0=3,min1=-3.5,max1=3.5){
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
      repeat{
        z=numeric()
        for(i in 1:(c-d)){
          x10=h0[i:c,]
          x11=h1[i:c,]
          l0=el.test.newton(x=x10,mu=rep(0,d)) 
          l1=el.test.newton(x=x11,mu=rep(0,d)) 
          if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
            tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
            tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
            z[i]=0.5*(tl0-tl1)
          }
        }
        z=na.omit(z)
        if(length(z)==0){
          x9=rbind(x9,runif(n=d,min=min0,max=max0))
          c=c+1
          print(c)
          h0=rbind(h0,(x9[c,]-x51)^2-x5)
          h1=rbind(h1,(x9[c,]-x81)^2-x8)
        }
        if(length(z)>0){
          if(max(z)>threshold){
            break
          }
        }
        if(length(z)>0){
          if(max(z)<=threshold){
            x9=rbind(x9,runif(n=d,min=min0,max=max0))
            c=c+1
            print(c)
            h0=rbind(h0,(x9[c,]-x51)^2-x5)
            h1=rbind(h1,(x9[c,]-x81)^2-x8)
          }
        }
        if(c>=49){
          break
        }
      }
      if(length(z)>0){
        if(max(z)>threshold){
          break
        }
      }
      if(c>=49){
        break
      }
    }
    if(c>=49){
      break
    }
  }
  nodetect=0
  repeat{
    z=numeric()
    x9=rbind(x9,runif(n=d,min=min1,max=max1))
    c=c+1
    print(c)
    h0=rbind(h0,(x9[c,]-x51)^2-x5)
    h1=rbind(h1,(x9[c,]-x81)^2-x8)
    z=numeric()
    for(i in 1:(c-d)){
      x10=h0[i:c,]
      x11=h1[i:c,]
      l0=el.test.newton(x=x10,mu=rep(0,d)) 
      l1=el.test.newton(x=x11,mu=rep(0,d)) 
      if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
        tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
        tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
        z[i]=0.5*(tl0-tl1)
      }
    }
    z=na.omit(z)
    if(length(z)==0){
      x9=rbind(x9,runif(n=d,min=min1,max=max1))
      c=c+1
      print(c)
      h0=rbind(h0,(x9[c,]-x51)^2-x5)
      h1=rbind(h1,(x9[c,]-x81)^2-x8)
    }
    if(length(z)>0){
      if(max(z)>threshold){
        break
      }
    }
    if(length(z)>0){
      if(max(z)<=threshold){
        x9=rbind(x9,runif(n=d,min=min1,max=max1))
        c=c+1
        print(c)
        h0=rbind(h0,(x9[c,]-x51)^2-x5)
        h1=rbind(h1,(x9[c,]-x81)^2-x8)
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
      x9=rbind(x9,runif(n=d,min=min1,max=max1))
      c=c+1
      print(c)
      h0=rbind(h0,(x9[c,]-x51)^2-x5)
      h1=rbind(h1,(x9[c,]-x81)^2-x8)
      current0=h0[(nrow(h0)-99):nrow(h0),]
      current1=h1[(nrow(h1)-99):nrow(h1),]
      z=numeric()
      for(i in 1:(nrow(current0)-d)){
        x10=current0[i:nrow(current0),]
        x11=current1[i:nrow(current0),]
        l0=el.test.newton(x=x10,mu=rep(0,d)) 
        l1=el.test.newton(x=x11,mu=rep(0,d)) 
        if(abs(sum(l0$wts)-1)<0.01&abs(sum(l1$wts)-1)<0.01){
          tl0=l0$`-2LLR`*max(1-l0$`-2LLR`/nrow(x10),0.5)
          tl1=l1$`-2LLR`*max(1-l1$`-2LLR`/nrow(x11),0.5)
          z[i]=0.5*(tl0-tl1)
        }
      }
      z=na.omit(z)
      if(length(z)==0){
        x9=rbind(x9,runif(n=d,min=min1,max=max1))
        c=c+1
        print(c)
        h0=rbind(h0,(x9[c,]-x51)^2-x5)
        h1=rbind(h1,(x9[c,]-x81)^2-x8)
      }
      if(length(z)>0){
        print(max(z))
        if(max(z)>threshold){
          break
        }
      }
      if(length(z)>0){
        print(max(z))
        if(max(z)<=threshold){
          x9=rbind(x9,runif(n=d,min=min1,max=max1))
          c=c+1
          print(c)
          h0=rbind(h0,(x9[c,]-x51)^2-x5)
          h1=rbind(h1,(x9[c,]-x81)^2-x8)
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
