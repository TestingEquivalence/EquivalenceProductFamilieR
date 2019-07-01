product<-function(r,c){
  res=r %*% t(c)
  return(res)
}

l22<-function(tab1,tab2){
  t=tab1-tab2
  t=t*t
  s=sum(t)
  return(s)
}

l22_first_derivative<-function(tab1,tab2){
  res=tab1-tab2
  res=2*res
  return(res)
}

startValue<-function(tab){
  r=rowSums(tab)
  c=colSums(tab)
  return(c(r,c))
}

rc2table<-function(rc,table){
  nr=nrow(table)
  nc=ncol(table)
  r=rc[1:nr]
  c=rc[(nr+1):(nr+nc)]
  
  #normalize
  r=r/sum(r)
  c=c/sum(c)
  
  prod=product(r,c)
  return(prod)
}

fn<-function(rc, table){
  prod=rc2table(rc,table)
  dst=l22(table,prod)
  return(dst)
}

min_l22<-function(table){
   start=startValue(table)
   l=length(start)
   lowerBound=rep(0,l)
   upperBound=rep(1,l)
   opRes=optim(par=start,fn=fn, gr=NULL, method="L-BFGS-B", lower=lowerBound,
               upper=upperBound, table=table)
   return(opRes)               
}