source("distance.R")
source("asymptotic_test.R")

randomPoint<-function(tab){
  size=nrow(tab)*ncol(tab)
  x=runif(size,0,1)
  x=x/sum(x)
  m=matrix(data=x,nrow=nrow(tab), ncol=ncol(tab),byrow = TRUE)
  return(m)
}

randomExteriorPoint<-function(tab, sqEps){
  repeat{
    m = randomPoint(tab)
    res= min_l22(m)
    if (res$value>sqEps) return(m)
  }
}

linComb<-function(x,y,a){
  return((1-a)*x+a*y) 
}

linearBoundaryPoint<-function(p,q,sqEps){
  aim<-function(a){
    lc=linComb(p,q,a)
    res= min_l22(lc)
    return(res$value-sqEps)
  }
  
  aMin=uniroot(aim, c(0,1))
  return(linComb(p,q,aMin$root))
}

protoBstTest<-function(tab,n,sqEps,exteriorPoints,nSimulation){
  #calculate test statistic
  res= min_l22(tab)
  t=res$value
  
  #estimate closest boundary point
  rp=tab
  
  df<-function(x){
    r=x-tab
    sr=r*r
    sse=sum(sr)
    return(sse)
  }
  
  if (t<sqEps){
    #function(p,q,eps,distance)
    bps=lapply(X=exteriorPoints,FUN=linearBoundaryPoint,
               q=tab,sqEps=sqEps)
    dst=lapply(bps, df)
    pos=which.min(dst)
    rp=bps[[pos]]
  }
  
  
  #simulate bootstrap sample
  i=c(1:nSimulation)
  f<-function(k){
    vrp=as.vector(rp)
    v=rmultinom(n=1,size=n,prob=vrp)
    v=v/sum(v)
    m=matrix(dat=v,nrow=nrow(tab), ncol=ncol(tab))
    res= min_l22(m)
    return(res$value)
  }
  sample=lapply(i,f)
  
  #bootstrap test
  pValue=sum(sample<t)/nSimulation
  return(pValue)
}

#' The bootstrap test is based on the re-sampling method called bootstrap.
#' The bootstrap test is more precise and reliable than the asymptotic test.
#' However, it should be used carefully because the test is approximate 
#' and may be anti-conservative. 
#' In order to obtain a conservative test reducing of alpha
#' (usually halving) or slight shrinkage of the tolerance parameter epsilon
#' may be appropriate. We prefer the slight shrinkage of the tolerance parameter 
#' because it is more effective and the significance level remains unchanged.
#' \code{bootstrap_test} bootstrap test for approximate row column independence
#' in two way contingency tables. 
#' The test statistic is the minimum of the Euclidean distance 
#' between the contingency table and a product measure.
#' @param tab contingency table containing the counts of events
#' @param alpha significance level
#' @param nSimulation number of bootstrap samples, default 10000 
#' @param nExteriorPoints number of random directions to search for a boundary point,
#' default is (nrow(tab)+ncol(tab))*50
#' @return test returns the minimum tolerance parameter epsilon,
#' for which the approximate independence can be shown

bootstrap_test<-function(tab, alpha, 
                         nSimulation=10000, 
                         nExteriorPoints=0){
  #find start value for min eps
  #use for this purpose the asymptotic test with 
  #small safety margin
  eps=asymptotic_test(tab,alpha)
  sqEps=eps*eps*1.2
  
  
  n=sum(tab)
  tab=tab/n
  
  #number of search directions and seed
  if (nExteriorPoints==0) 
    nExteriorPoints=(nrow(tab)+ncol(tab))*50
  set.seed(10071977)
  
  #calculate exterior points
  f<-function(x){
    randomExteriorPoint(tab,sqEps)
  }
  
  i=c(1:nExteriorPoints)
  exteriorPoints=lapply(i, f)
  
  #calculate min epsilon
  ff<-function(x){
    set.seed(01012019)
    pval=protoBstTest(tab,n,x,exteriorPoints,nSimulation)
    return(pval-alpha)
  }
  
  #check boundary values
  #check lower bound
  lb=ff(0)
  if (lb<0) return(0)
  
  #check upper bound
  ub=ff(sqEps)
  if (ub>0) return(NA)
  
  res=uniroot(ff,c(0,sqEps))
  return(sqrt(res$root))
}
