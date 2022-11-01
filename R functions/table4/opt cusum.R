opt_cusum=function(threshold=100,d=2){
  
  c=0
  repeat{
    l=numeric()
    l1=numeric()
    repeat{
      x0=rgamma(d,shape = 10,rate=2)
      c=c+1
      print(c)
      l[length(l)+1]=prod(dgamma(x0,shape =5,rate = 1)/dgamma(x0,shape =10,rate = 2))
      
      for(i in 1:length(l)){
        l1[i]=sum(log(l[i:length(l)]))
      }
      print(max(l1))
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




opt1_cusum=function(d=2,threshold=100){
  false=0
  c=0
  repeat{
    repeat{
      c1=0
      l=numeric()
      l1=numeric()
      repeat{
        x0=rgamma(d,shape = 10,rate=2)
        c=c+1
        c1=c1+1
        print(c)
        l[c1]=prod(dgamma(x0,shape =5,rate = 1)/dgamma(x0,shape =10,rate = 2))
        
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
        false=false+1
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
    x0=rgamma(d,shape = 5,rate=1)
    c=c+1
    c1=length(l)+1
    print(c)
    l[c1]=prod(dgamma(x0,shape =5,rate = 1)/dgamma(x0,shape =10,rate = 2))
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
      nodetect=1
      break
    }
    
  }
  if(max(l1)<=threshold){
    repeat{
      c1=0
      l=numeric()
      l1=numeric()
      repeat{
        x0=rgamma(d,shape = 5,rate=1)
        c=c+1
        c1=c1+1
        print(c)
        l[c1]=prod(dgamma(x0,shape =5,rate = 1)/dgamma(x0,shape =10,rate = 2))
        
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
          nodetect=1
          break
        }
        
      }
      if(max(l1)>threshold){
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

