





ks_cusum=function(m=100,d=2,threshold=1, mu1=rep(0.5,d)){
  c=0
  y0=matrix(0,m,d)
  y1=matrix(0,m,d)
  for(i in 1:m){
    y0[i,]=rgamma(n=d,shape = 9,rate = 2)
    y1[i,]=y0[i,]+mu1
    
  }
  wh0=numeric()
  wh1=numeric()
  for(i in 1:d){
    wh0[i]=hpi(y0[,i],deriv.order =0)
    wh1[i]=hpi(y1[,i],deriv.order =0)
  }
  f0=function(x=rep(0,times=d)){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y0[i,])/wh0)^2/2)/wh0)
    }
    return(mean(p))
  }
  f1=function(x=rep(0,times=d)){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y1[i,])/wh1)^2/2)/wh1)
    }
    return(mean(p))
  }
  repeat{
    x=rgamma(n=d,shape = 9,rate = 2)
    c=c+1
    
    x1=numeric()
    x2=numeric()
    c1=0
    repeat{
      if(c1==0){
        x2[1]=log(f1(x)/f0(x))
        x1[1]=x2[1]
      }
      
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
      if(max(x1)>threshold){
        break
      }
      if(max(x1)<0){
        break
      }
      x=rbind(x,rgamma(n=d,shape = 9,rate = 2))
      c=c+1
      c1=c1+1
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



ks1_cusum=function(m=100,d=2,threshold=1,mu1=rep(0.5,times=d)){
  false=0
  c=0
  y0=matrix(0,m,d)
  y1=matrix(0,m,d)
  for(i in 1:m){
    y0[i,]=rgamma(n=d,shape = 9,rate = 2)
    y1[i,]=y0[i,]+mu1
    
  }
  wh0=numeric()
  wh1=numeric()
  for(i in 1:d){
    wh0[i]=hpi(y0[,i],deriv.order =0)
    wh1[i]=hpi(y1[,i],deriv.order =0)
  }
  f0=function(x=rep(0,times=d)){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y0[i,])/wh0)^2/2)/wh0)
    }
    return(mean(p))
  }
  f1=function(x=rep(0,times=d)){
    p=numeric()
    for(i in 1:m){
      p[i]=prod(1/sqrt(2*pi)*exp(-((x-y1[i,])/wh1)^2/2)/wh1)
    }
    return(mean(p))
  }
  repeat{
    repeat{
      x=rgamma(n=d,shape = 9,rate = 2)
      c=c+1
      
      x1=numeric()
      x2=numeric()
      c1=0
      repeat{
        if(c1==0){
          x2[1]=log(f1(x)/f0(x))
          x1[1]=x2[1]
        }
        
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
        
        if(max(x1)>threshold){
          false=false+1
          break
        }
        if(max(x1)<0){
          break
        }
        x=rbind(x,rgamma(n=d,shape = 9,rate = 2))
        c=c+1
        c1=c1+1
        
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
  nondetect=0
  repeat{
    x=rbind(x,rgamma(n=d,shape = 9,rate = 2)+mu1)
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
      nondetect=1
      break
    }
  }
  
  if(max(x1)<=threshold){
    repeat{
      x=rgamma(n=d,shape = 9,rate = 2)+mu1
      c=c+1
      
      x1=numeric()
      x2=numeric()
      c1=0
      repeat{
        if(c1==0){
          x2[1]=log(f1(x)/f0(x))
          x1[1]=x2[1]
        }
        
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
        
        if(max(x1)>threshold){
          break
        }
        if(max(x1)<0){
          break
        }
        x=rbind(x,rgamma(n=d,shape = 9,rate = 2)+mu1)
        c=c+1
        c1=c1+1
        if(c>1000){
          nondetect=1
          break
        }
      }
      if(max(x1)>threshold){
        break
      }
      if(c>1000){
        nondetect=1
        break
      }
    }
  }
  return(c(c-50,nondetect))
  }
  if(false>=10){
    return("false")
  }
}
