#Functions for the beta-prime distribution
#PDF and RNG are implemented right now
dbetaprime <- function( x,alphaparam,betaparam, p=1,q=1 ){
  mask <- x<=0
  x[mask]=0
  xq<-x/q
  scaling <- p/(q*beta(alphaparam,betaparam))
  numerator <- xq^(alphaparam*p-1)
  denominator <- (1+xq^p)^(alphaparam+betaparam)
  pdf <- scaling*numerator/denominator
  return(pdf)
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
