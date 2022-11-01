





app_cusum=function(threshold=200,m=100,d=2){
  
  
  mu1=rep(0.5,d)
  c=0
  y0=matrix(0,m,d)
  y1=matrix(0,m,d)
  for(i in 1:m){
    y0[i,]=rgamma(n=d,shape = 9,rate = 2)
    y1[i,]=y0[i,]+mu1
    
  }
  estmu0=colMeans(y0)
  estmu1=colMeans(y1)
  estcovar0=cov(y0)
  estcovar1=cov(y1)
  repeat{
    x=rgamma(n=d,shape = 9,rate = 2)
    c=c+1
    print(c)
    x1=numeric()
    x2=numeric()
    c1=0
    repeat{
      if(c1==0){
        x2[1]=log(dmvnorm(x ,mean = estmu1,sigma = estcovar1)/dmvnorm(x ,mean = estmu0,sigma = estcovar0))
        x1[1]=x2[1]
      }
      print(max(x1))
      if(max(x1)>threshold){
        break
      }
      if(max(x1)<0){
        break
      }
      if(c1>=1){
        x2[length(x2)+1]=log(dmvnorm(x[nrow(x),],mean = estmu1,sigma = estcovar1)/dmvnorm(x[nrow(x),],mean = estmu0,sigma = estcovar0))  
        for(i in 1:length(x2)){
          x1[i]=sum(x2[i:length(x2)])
        }
      }
      print(max(x1))
      if(max(x1)>threshold){
        break
      }
      if(max(x1)<0){
        break
      }
      x=rbind(x,rgamma(n=d,shape = 9,rate = 2))
      c=c+1
      c1=c1+1
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



app1_cusum=function(threshold=2.17,m=100,d=2){
  false=0
  mu1=rep(0.5,times=d)
  covar=diag(1,d,d)
  c=0
  y0=matrix(0,m,d)
  y1=matrix(0,m,d)
  for(i in 1:m){
    y0[i,]=rgamma(n=d,shape = 9,rate = 2)
    y1[i,]=y0[i,]+mu1
    
  }
  estmu0=colMeans(y0)
  estmu1=colMeans(y1)
  estcovar0=cov(y0)
  estcovar1=cov(y1)
  repeat{
    repeat{
      x=rgamma(n=d,shape = 9,rate = 2)
      c=c+1
      print(c)
      x1=numeric()
      x2=numeric()
      c1=0
      repeat{
        if(c1==0){
          x2[1]=log(dmvnorm(x ,mean = estmu1,sigma = estcovar1)/dmvnorm(x ,mean = estmu0,sigma = estcovar0))
          x1[1]=x2[1]
        }
        print(max(x1))
        if(max(x1)>threshold){
          false=false+1
          break
        }
        if(max(x1)<0){
          break
        }
        if(c1>=1){
          x2[length(x2)+1]= log(dmvnorm(x[nrow(x),],mean = estmu1,sigma = estcovar1)/dmvnorm(x[nrow(x),],mean = estmu0,sigma = estcovar0)) 
          for(i in 1:length(x2)){
            x1[i]=sum(x2[i:length(x2)])
          }
        }
        print(max(x1))
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
  nondetect=0
  repeat{
    x=rbind(x,rgamma(n=d,shape = 9,rate = 2)+mu1)
    c=c+1
    x2[length(x2)+1]=log(dmvnorm(x[nrow(x),],mean = estmu1,sigma = estcovar1)/dmvnorm(x[nrow(x),],mean = estmu0,sigma = estcovar0)) 
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
    if(c>=1000){
      nondetect=1
      break
    }
  }
  
  if(max(x1)<=threshold){
    repeat{
      x=rgamma(n=d,shape = 9,rate = 2)+mu1
      c=c+1
      print(c)
      x1=numeric()
      x2=numeric()
      c1=0
      repeat{
        if(c1==0){
          x2[1]= log(dmvnorm(x ,mean = estmu1,sigma = estcovar1)/dmvnorm(x ,mean = estmu0,sigma = estcovar0))
          x1[1]=x2[1]
        }
        print(max(x1))
        if(max(x1)>threshold){
          break
        }
        if(max(x1)<0){
          break
        }
        if(c1>=1){
          x2[length(x2)+1]= log(dmvnorm(x[nrow(x),],mean = estmu1,sigma = estcovar1)/dmvnorm(x[nrow(x),],mean = estmu0,sigma = estcovar0))
          for(i in 1:length(x2)){
            x1[i]=sum(x2[i:length(x2)])
          }
        }
        print(max(x1))
        if(max(x1)>threshold){
          break
        }
        if(max(x1)<0){
          break
        }
        x=rbind(x,rgamma(n=d,shape = 9,rate = 2)+mu1)
        c=c+1
        c1=c1+1
        print(c) 
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
