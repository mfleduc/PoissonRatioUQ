#Functions for the beta-prime distribution
#PDF, CDF, and RNG are implemented right now
dbetaprime <- function( x,alphaparam,betaparam, p=1,q=1 ){
  mask <- x<=0
  x[mask]=0
  xq<-x/q
  scaling <- log(p)-log(q)-lbeta(alphaparam,betaparam)
  numerator <- log(xq)*(alphaparam*p-1)
  denominator <- log(1+xq^p)*(alphaparam+betaparam)
  logpdf <- scaling + numerator - denominator
  return(exp(logpdf))
}
pbetaprime <- function(x, alphaparam,betaparam,p=1,q=1){
  mask <- x<=0
  x[mask]=0
  xqp <- (x/q)^p
  cdf <- pbeta(xqp/(1+xqp),alphaparam,betaparam)
  return(cdf)
}
rbetaprime <- function(n, alphaparam,betaparam, p=1,q=1){
  X <- rbeta(n, alphaparam,betaparam)
  Y <- X/(1-X)
  Z <- q*Y^(1/p)
  return(Z)
}
