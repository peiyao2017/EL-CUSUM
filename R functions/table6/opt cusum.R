opt_cusum=function(threshold=100,d=4,v0=list(c(-1.5,1,0.9),c(1.5,1,-0.9)),v1=list(c(-2,1,0.9),c(2,1,-0.9))){
  c=0
  repeat{
    l=numeric()
    l1=numeric()
    repeat{
      x0=rmix(n=d,pii=c(0.5,0.5),family = "Skew.normal",arg = v0)
      c=c+1
      print(c)
      l[length(l)+1]=prod((0.5*dsn(x0,xi=v1[[1]][1],omega=v1[[1]][2], alpha=v1[[1]][3])+0.5*dsn(x0,xi=v1[[2]][1],omega=v1[[2]][2], alpha=v1[[2]][3]))/(0.5*dsn(x0,xi=v0[[1]][1],omega=v0[[1]][2], alpha=v0[[1]][3])+0.5*dsn(x0,xi=v0[[2]][1],omega=v0[[2]][2], alpha=v0[[2]][3])))
      
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




opt1_cusum=function(d=4,threshold=100,v0=list(c(-1.5,1,0.9),c(1.5,1,-0.9)),v1=list(c(-2,1,0.9),c(2,1,-0.9))){
  false=0
  c=0
  repeat{
    repeat{
      c1=0
      l=numeric()
      l1=numeric()
      repeat{
        x0=rmix(n=d,pii=c(0.5,0.5),family = "Skew.normal",arg = v0)
        c=c+1
        c1=c1+1
        print(c)
        l[c1]=prod((0.5*dsn(x0,xi=v1[[1]][1],omega=v1[[1]][2], alpha=v1[[1]][3])+0.5*dsn(x0,xi=v1[[2]][1],omega=v1[[2]][2], alpha=v1[[2]][3]))/(0.5*dsn(x0,xi=v0[[1]][1],omega=v0[[1]][2], alpha=v0[[1]][3])+0.5*dsn(x0,xi=v0[[2]][1],omega=v0[[2]][2], alpha=v0[[2]][3])))
        
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
  nodetect=0
  repeat{
    x0=rmix(n=d,pii=c(0.5,0.5),family = "Skew.normal",arg = v1)
    c=c+1
    c1=length(l)+1
    print(c)
    l[c1]=prod((0.5*dsn(x0,xi=v1[[1]][1],omega=v1[[1]][2], alpha=v1[[1]][3])+0.5*dsn(x0,xi=v1[[2]][1],omega=v1[[2]][2], alpha=v1[[2]][3]))/(0.5*dsn(x0,xi=v0[[1]][1],omega=v0[[1]][2], alpha=v0[[1]][3])+0.5*dsn(x0,xi=v0[[2]][1],omega=v0[[2]][2], alpha=v0[[2]][3])))
    for(i in 1:c1){
      l1[i]=sum(log(l[i:c1]))
    }
    if(max(l1)>threshold){
      break
    }
    if(max(l1)<0){
      break
    }
    if(c>1000){
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
        x0=rmix(n=d,pii=c(0.5,0.5),family = "Skew.normal",arg = v1)
        c=c+1
        c1=c1+1
        print(c)
        l[c1]=prod((0.5*dsn(x0,xi=v1[[1]][1],omega=v1[[1]][2], alpha=v1[[1]][3])+0.5*dsn(x0,xi=v1[[2]][1],omega=v1[[2]][2], alpha=v1[[2]][3]))/(0.5*dsn(x0,xi=v0[[1]][1],omega=v0[[1]][2], alpha=v0[[1]][3])+0.5*dsn(x0,xi=v0[[2]][1],omega=v0[[2]][2], alpha=v0[[2]][3])))
        
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


