#' Functions for doing basic calculations with the beta-prime distribution
#'
#' @title Beta-prime PDF
#' @description
#' Evaluates the Beta-Prime PDF with parameters alpha,beta,p,q (parameterization used in https://en.wikipedia.org/wiki/Beta_prime_distribution#Generalization )
#' @param x scalar or matrix. Values at which to evaluate the PDF
#' @param alphaparam Scalar value, alpha parameter of the distribution
#' @param betaparam Scalar value, beta parameter of the distribution
#' @param p Scalar value, p parameter of the distribution
#' @param q Scalar value, q parameter of the distribution
#' @returns The beta-prime pdf with input parameters evaluated at the locations in x. Check values, possible issues for very large alphaval or betaval parameters
#' @export
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
#' @title Beta-prime CDF
#' @description
#' Evaluates the Beta-Prime CDF with parameters alpha,beta,p,q (parameterization used in https://en.wikipedia.org/wiki/Beta_prime_distribution#Generalization )
#' @param x scalar or matrix. Values at which to evaluate the CDF
#' @param alphaparam Scalar value, alpha parameter of the distribution
#' @param betaparam Scalar value, beta parameter of the distribution
#' @param p Scalar value, p parameter of the distribution
#' @param q Scalar value, q parameter of the distribution
#' @returns The beta-prime cdf with input parameters evaluated at the locations in x. Check values, possible issues for very large alphaval or betaval parameters
#' @export
pbetaprime <- function(x, alphaparam,betaparam,p=1,q=1){
  mask <- x<=0
  x[mask]=0
  xqp <- (x/q)^p
  cdf <- pbeta(xqp/(1+xqp),alphaparam,betaparam)
  return(cdf)
}
#' @title Beta-prime random number generator
#' @description
#' generates random draws from the Beta-prime distribution with parameters alpha,beta,p,q (parameterization used in https://en.wikipedia.org/wiki/Beta_prime_distribution#Generalization )
#' @param n Scalar value, number of draws from the distribution
#' @param alphaparam Scalar value, alpha parameter of the distribution
#' @param betaparam Scalar value, beta parameter of the distribution
#' @param p Scalar value, p parameter of the distribution
#' @param q Scalar value, q parameter of the distribution
#' @returns A vector of length n containing draws from the desired beta-prime distribution
#' @export
rbetaprime <- function(n, alphaparam,betaparam, p=1,q=1){
  X <- rbeta(n, alphaparam,betaparam)
  Y <- X/(1-X)
  Z <- q*Y^(1/p)
  return(Z)
}
