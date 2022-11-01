opt_cusum=function(threshold=1, d=2){
  
  c=0
  repeat{
    c1=0
    l=numeric()
    l1=numeric()
    repeat{
      x0=rgamma(d,shape = 9,rate=2)
      c=c+1
      c1=c1+1
      print(c)
      l[c1]=prod(dgamma(x0-rep(0.5,times=d),shape =9,rate = 2)/dgamma(x0,shape =9,rate = 2))
      
      for(i in 1:c1){
        l1[i]=sum(log(l[i:c1]))
      }
      if(max(l1)>threshold){
        break
      }
      if(max(l1)<0){
        break
      }
    }
    if(max(l1)>threshold){
      break
    }
    if(c>1000){
      break
    }
  }
  return(c)
}




opt1_cusum=function(d=2,threshold=1){
  false=0
  c=0
  repeat{
    repeat{
      c1=0
      l=numeric()
      l1=numeric()
      repeat{
        x0=rgamma(d,shape = 9,rate=2)
        c=c+1
        c1=c1+1
        print(c)
        l[c1]=prod(dgamma(x0-rep(0.5,times=d),shape =9,rate = 2)/dgamma(x0,shape =9,rate = 2))
        
        for(i in 1:c1){
          l1[i]=sum(log(l[i:c1]))
        }
        if(max(l1)>threshold){
          false=false+1
          break
        }
        if(max(l1)<0){
          break
        }
        if(c>=49){
          break
        }
        if(false>=10){
          break
        }
      }
      if(max(l1)>threshold){
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
    x0=rgamma(d,shape = 9,rate=2)+rep(0.5,times=d)
    c=c+1
    c1=length(l)+1
    print(c)
    l[c1]=prod(dgamma(x0-rep(0.5,times=d),shape =9,rate = 2)/dgamma(x0,shape =9,rate = 2))
    for(i in 1:c1){
      l1[i]=sum(log(l[i:c1]))
    }
    if(max(l1)>threshold){
      break
    }
    if(max(l1)<0){
      break
    }
    if(c>=1000){
      nondetect=1
      break
    }
  }
  if(max(l1)<=threshold){
    repeat{
      c1=0
      l=numeric()
      l1=numeric()
      repeat{
        x0=rgamma(d,shape = 9,rate=2)+rep(0.5,times=d)
        c=c+1
        c1=c1+1
        print(c)
        l[c1]=prod(dgamma(x0-rep(0.5,times=d),shape =9,rate = 2)/dgamma(x0,shape =9,rate = 2))
        
        for(i in 1:c1){
          l1[i]=sum(log(l[i:c1]))
        }
        if(max(l1)>threshold){
          break
        }
        if(max(l1)<0){
          break
        }
        if(c>=1000){
          nondetect=1
          break
        }
      }
      if(max(l1)>threshold){
        break
      }
      if(c>=1000){
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


