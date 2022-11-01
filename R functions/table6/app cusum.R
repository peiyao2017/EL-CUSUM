

app_cusum=function(threshold=100,m=100,d=4,v0=list(c(-1.5,1,0.9),c(1.5,1,-0.9)),v1=list(c(-2,1,0.9),c(2,1,-0.9))){
  c=0
  y0=matrix(0,m,d)
  y1=matrix(0,m,d)
  for(i in 1:m){
    y0[i,]=rmix(n=d,pii=c(0.5,0.5),family = "Skew.normal",arg = v0)
    y1[i,]=rmix(n=d,pii=c(0.5,0.5),family = "Skew.normal",arg = v1)
    
  }
  estmu0=colMeans(y0)
  estmu1=colMeans(y1)
  estcovar0=cov(y0)
  estcovar1=cov(y1)
  repeat{
    x=rmix(n=d,pii=c(0.5,0.5),family = "Skew.normal",arg = v0)
    c=c+1
    print(c)
    x1=numeric()
    x2=numeric()
    c1=0
    repeat{
      if(c1==0){
        x2[1]=dmvnorm(x,mean = estmu1,sigma = estcovar1)/dmvnorm(x,mean = estmu0,sigma = estcovar0)
        x1[1]=log(x2[1])
      }
      print(max(x1))
      if(max(x1)>threshold){
        break
      }
      if(max(x1)<0){
        break
      }
      if(c1>=1){
        x2[length(x2)+1]=dmvnorm(x[nrow(x),],mean = estmu1,sigma = estcovar1)/dmvnorm(x[nrow(x),],mean = estmu0,sigma = estcovar0)
        print( x2[length(x2)])
        for(i in 1:length(x2)){
          x1[i]=sum(log(x2[i:length(x2)]))
        }
      }
      print(max(x1))
      if(max(x1)>threshold){
        break
      }
      if(max(x1)<0){
        break
      }
      x=rbind(x,rmix(n=d,pii=c(0.5,0.5),family = "Skew.normal",arg = v0))
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



app1_cusum=function(m=100,d=4,threshold=100,v0=list(c(-1.5,1,0.9),c(1.5,1,-0.9)),v1=list(c(-2,1,0.9),c(2,1,-0.9))){
  false=0
  c=0
  y0=matrix(0,m,d)
  y1=matrix(0,m,d)
  for(i in 1:m){
    y0[i,]=rmix(n=d,pii=c(0.5,0.5),family = "Skew.normal",arg = v0)
    y1[i,]=rmix(n=d,pii=c(0.5,0.5),family = "Skew.normal",arg = v1)
    
  }
  estmu0=colMeans(y0)
  estmu1=colMeans(y1)
  estcovar0=cov(y0)
  estcovar1=cov(y1)
  repeat{
    repeat{
      x=rmix(n=d,pii=c(0.5,0.5),family = "Skew.normal",arg = v0)
      c=c+1
      print(c)
      x1=numeric()
      x2=numeric()
      c1=0
      repeat{
        if(c1==0){
          x2[1]=dmvnorm(x,mean = estmu1,sigma = estcovar1)/dmvnorm(x,mean = estmu0,sigma = estcovar0)
          x1[1]=log(x2[1])
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
          x2[length(x2)+1]=dmvnorm(x[nrow(x),],mean = estmu1,sigma = estcovar1)/dmvnorm(x[nrow(x),],mean = estmu0,sigma = estcovar0)
          for(i in 1:length(x2)){
            x1[i]=sum(log(x2[i:length(x2)]))
          }
        }
        print(max(x1))
        if(max(x1)>threshold){
          break
        }
        if(max(x1)<0){
          break
        }
        x=rbind(x,rmix(n=d,pii=c(0.5,0.5),family = "Skew.normal",arg = v0))
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
  nodetect=0
  repeat{
    x=rbind(x,rmix(n=d,pii=c(0.5,0.5),family = "Skew.normal",arg = v1))
    c=c+1
    x2[length(x2)+1]=dmvnorm(x[nrow(x),],mean = estmu1,sigma = estcovar1)/dmvnorm(x[nrow(x),],mean = estmu0,sigma = estcovar0)
    for(i in 1:length(x2)){
      x1[i]=sum(log(x2[i:length(x2)]))
    }
    if(max(x1)>threshold){
      break
    }
    if(max(x1)<0){
      break
    }
    if(c>1000){
      nodetect=1
      break
    }
  }
  
  if(max(x1)<=threshold){
    repeat{
      x=rmix(n=d,pii=c(0.5,0.5),family = "Skew.normal",arg = v1)
      c=c+1
      print(c)
      x1=numeric()
      x2=numeric()
      c1=0
      repeat{
        if(c1==0){
          x2[1]=dmvnorm(x ,mean = estmu1,sigma = estcovar1)/dmvnorm(x ,mean = estmu0,sigma = estcovar0)
          x1[1]=log(x2[1])
        }
        print(max(x1))
        if(max(x1)>threshold){
          break
        }
        if(max(x1)<0){
          break
        }
        if(c1>=1){
          x2[length(x2)+1]=dmvnorm(x[nrow(x),],mean = estmu1,sigma = estcovar1)/dmvnorm(x[nrow(x),],mean = estmu0,sigma = estcovar0)
          for(i in 1:length(x2)){
            x1[i]=sum(log(x2[i:length(x2)]))
          }
        }
        print(max(x1))
        if(max(x1)>threshold){
          break
        }
        if(max(x1)<0){
          break
        }
        x=rbind(x,rmix(n=d,pii=c(0.5,0.5),family = "Skew.normal",arg = v1))
        c=c+1
        c1=c1+1
        print(c) 
        if(c>1000){
          nodetect=1
          break
        }
      }
      if(max(x1)>threshold){
        break
      }
      if(c>1000){
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



