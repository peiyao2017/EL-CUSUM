ks_cusum=function(threshold=1,m=100,mu0=0,mu1=0.5, covar=1){
  c=0
  y0=rnorm(n=m,mean=mu0,sd=covar)
  y1=rnorm(n=m,mean=mu1,sd=covar)
  wh0=hpi(y0,deriv.order =0)
  wh1=hpi(y1,deriv.order =0)
  f0=function(x=0){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y0[i])/wh0)^2/2)/wh0)
    }
    return(mean(p))
  }
  f1=function(x=0){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y1[i])/wh1)^2/2)/wh1)
    }
    return(mean(p))
  }
  repeat{
    x=rnorm(n=1,mean=mu0,sd=covar)
    c=c+1
    print(c)
    x1=numeric()
    x2=numeric()
    repeat{
      x2[length(x2)+1]=log(f1(x)/f0(x))
      for(i in 1:length(x2)){
        x1[i]=sum(x2[i:length(x2)])
      }
      print(max(x1))
      if(max(x1)>threshold){
        break
      }
      if(max(x1)<0){
        break
      }
      x=c(x,rnorm(n=1,mean=mu0,sd = covar))
      c=c+1
      print(c) 
    }
    if(max(x1)>threshold){
      break
    }
    if(c>1000){
      break
    }
  }
  return(c)
}



ks1_cusum=function(m=100,threshold=1,mu0=0,mu1=0.5,covar=1){
  c=0
  false=0
  y0=rnorm(n=m,mean=mu0,sd=covar)
  y1=rnorm(n=m,mean=mu1,sd=covar)
  wh0=hpi(y0,deriv.order =0)
  wh1=hpi(y1,deriv.order =0)
  f0=function(x=0){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y0[i])/wh0)^2/2)/wh0)
    }
    return(mean(p))
  }
  f1=function(x=0){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y1[i])/wh1)^2/2)/wh1)
    }
    return(mean(p))
  }
  repeat{
    repeat{
      x=rnorm(n=1,mean=mu0,sd = covar)
      c=c+1
      print(c)
      x1=numeric()
      x2=numeric()
      repeat{
        x2[length(x2)+1]=log(f1(x)/f0(x))
        for(i in 1:length(x2)){
          x1[i]=sum(x2[i:length(x2)])
        }
        print(max(x1))
        if(max(x1)>threshold){
          false=false+1
          break
        }
        if(max(x1)<0){
          break
        }
        x=c(x,rnorm(n=1,mean=mu0,sd = covar))
        c=c+1
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
    x=rbind(x,rnorm(n=1,mean=mu1,sd = covar))
    c=c+1
    x2[length(x2)+1]=log(f1(x)/f0(x))
    for(i in 1:length(x2)){
      x1[i]=sum(x2[i:length(x2)])
    }
    if(max(x1)>threshold){
      break
    }
    if(max(x1)<0){
      break
    }
    if(c>=1000){
      nondetect=1
      break
    }
  }
  if(max(x1)<=threshold){
    repeat{
      x=rnorm(n=1,mean=mu1,sd = covar)
      c=c+1
      print(c)
      x1=numeric()
      x2=numeric()
      repeat{
        x2[length(x2)+1]=log(f1(x)/f0(x))
        for(i in 1:length(x2)){
          x1[i]=sum(x2[i:length(x2)])
        }
        print(max(x1))
        if(max(x1)>threshold){
          break
        }
        if(max(x1)<0){
          break
        }
        x=c(x,rnorm(n=1,mean=mu1,sd = covar))
        c=c+1
        print(c) 
        if(c>=1000){
          nodetect=1
          break
        }
      }
      if(max(x1)>threshold){
        break
      }
      if(c>=1000){
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

