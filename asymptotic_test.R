source("distance.R")

asympt_stdev<-function(p,derivative){
  vec = derivative
  vnsq_1  = sum(p*vec*vec)
  
  k=length(p)
  vnsq_2=0
  for (j1 in 1:k)
    for (j2 in 1:k)
      vnsq_2 = vnsq_2 + vec[j1] * vec[j2] * p[j1] * p[j2]
  
  
  vnsq  = vnsq_1 - vnsq_2
  return (sqrt(vnsq))
}


#' The asymptotic test is based on the asymptotic distribution of the test statistic. 
#' The asymptotic test needs some sufficiently large number of the observations
#' in any cell of the contingency table.
#' It should be used carefully because the test is approximate 
#' and may be anti-conservative at some points. 
#' In order to obtain a conservative test reducing of alpha  (usually halving) or
#' slight shrinkage of the tolerance parameter epsilon may be appropriate.
#' We prefer the slight shrinkage of the tolerance parameter 
#' because it is more effective and the significance level remains unchanged.
#' \code{asymptotic_test} asymptotic test for approximate row column independence
#' in two way contingency tables. 
#' The test statistic is the minimum of the Euclidean distance 
#' between the contingency table and a product measure.
#' @param tab contingency table containing the counts of events
#' @param alpha significance level
#' @return test returns the minimum tolerance parameter epsilon,
#' for which the approximate independence can be shown

asymptotic_test<-function(tab, alpha){
  #normalize tab
  n=sum(tab)
  tab=tab/n
  
  #calculate minimum distance
  res=min_l22(tab)
  q=rc2table(res$par,tab)
  
    
  vtab=as.vector(t(tab))
  der=l22_first_derivative(tab, q)
  vder=as.vector(t(der))
  
  vol = asympt_stdev(vtab,vder) / sqrt(n)
  qt=qnorm(1-alpha,0,1)
  t= res$value
  eps = t + qt*vol
  eps=sqrt(eps)
  return(eps)
}
