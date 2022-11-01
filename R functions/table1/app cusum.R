
app_cusum=function(m=100,threshold=200,mu0=0,mu1=0.5,covar=1){
  c=0
  y0=rnorm(n=m,mean=mu0,sd=covar)
  y1=rnorm(n=m,mean=mu1,sd=covar)
  estmu0=mean(y0)
  estmu1=mean(y1)
  estcovar0=sd(y0)
  estcovar1=sd(y1)
  repeat{
    x=rnorm(n=1,mean=mu0,sd=covar)
    c=c+1
    print(c)
    x1=numeric()
    x2=numeric()
    repeat{
      x2[length(x2)+1]=log(dnorm(x[length(x)],mean = estmu1,sd = estcovar1)/dnorm(x[length(x)],mean = estmu0,sd= estcovar0))
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



app1_cusum=function(m=100,threshold=200,mu0=0,mu1=0.5,covar=1){
  false=0
  c=0
  y0=rnorm(n=m,mean=mu0,sd=covar)
  y1=rnorm(n=m,mean=mu1,sd=covar)
  estmu0=mean(y0)
  estmu1=mean(y1)
  estcovar0=sd(y0)
  estcovar1=sd(y1)
  repeat{
    repeat{
      x=rnorm(n=1,mean=mu0,sd = covar)
      c=c+1
      print(c)
      x1=numeric()
      x2=numeric()
      repeat{
        x2[length(x2)+1]=log(dnorm(x[length(x)],mean = estmu1,sd= estcovar1)/dnorm(x[length(x)],mean = estmu0,sd= estcovar0))
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
    x2[length(x2)+1]=log(dnorm(x[length(x)],mean = estmu1,sd= estcovar1)/dnorm(x[length(x)],mean = estmu0,sd= estcovar0))
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
      nodetect=1
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
        x2[length(x2)+1]=log(dnorm(x[length(x)],mean = estmu1,sd= estcovar1)/dnorm(x[length(x)],mean = estmu0,sd= estcovar0))
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
