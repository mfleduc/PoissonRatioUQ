#' Functions for doing basic calculations with the generalized beta-prime distribution
#'
#' @title generalized Beta-prime PDF
#' @description
#' Evaluates the generalized Beta-Prime PDF with parameters alpha,beta,p,q (parameterization used in https://en.wikipedia.org/wiki/Beta_prime_distribution#Generalization )
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
#' @title generalized Beta-prime CDF
#' @description
#' Evaluates the generalized Beta-Prime CDF with parameters alpha,beta,p,q (parameterization used in https://en.wikipedia.org/wiki/Beta_prime_distribution#Generalization ).This uses the fact that if X ~ BP(a,b,p,q) then (X/q)^p  ~ BP(a,b)
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
#' @title generalized Beta-prime random number generator
#' @description
#' generates random draws from the generalized Beta-prime distribution with parameters alpha,beta,p,q (parameterization used in https://en.wikipedia.org/wiki/Beta_prime_distribution#Generalization ). This uses the fact that if X ~ Beta(a,b) then q(X/(1-X))^(1/p) is BP(a,b,p,q)
#' @param n Scalar value, number of draws from the distribution
#' @param alphaparam Scalar value, alpha parameter of the distribution
#' @param betaparam Scalar value, beta parameter of the distribution
#' @param p Scalar value, p parameter of the distribution
#' @param q Scalar value, q parameter of the distribution
#' @returns A vector of length n containing draws from the desired gen. beta-prime distribution
#' @export
rbetaprime <- function(n, alphaparam,betaparam, p=1,q=1){
  X <- rbeta(n, alphaparam,betaparam)
  Y <- X/(1-X)
  Z <- q*Y^(1/p)
  return(Z)
}
#' @title Calculate the distribution of Z under a Beta-Prime assumption
#' @description given count datasets a and b, calculates the distribution of Z = lambda_a / lambda_b where lambda_i is the Poisson intensity of the channel that produced dataset i. Under a Gamma prior for each intensity, this is a Beta-Prime distribution. The Gamma pdf uses the shape/rate parameterization and the default prior is uninformative. A prior of the form x^-k_i can be obtained by setting ai=1-ki and bi=0
#' @param a Matrix. The count data for the numerator. Rows correspond to spatial locations, columns to realizations. Must have the same number of rows as b, but need not have the same number of columns. MUST be a matrix or array right now, if you want a scalar use a 1x1 array. Missing data should be replaced by NaNs.
#' @param b Matrix. The count data for the denominator. Rows correspond to spatial locations, columns to realizations. Must have the same number of rows as b, but need not have the same number of columns. MUST be a matrix or array right now, if you want a scalar use a 1x1 array. Missing data should be replaced by NaNs.
#' @param a1 Vector or scalar. The shape parameter of the Gamma prior for the upper channel. If a vector, must have the same number of entries as rows in a and b. Default is 1.
#' @param b1 Vector or scalar. The rate parameter of the Gamma prior for the upper channel. If a vector, must have the same number of entries as rows in a and b. Default is 0.
#' @param a2 Vector or scalar. The shape parameter of the Gamma prior for the lower channel. If a vector, must have the same number of entries as rows in a and b. Default is 1.
#' @param b2 Vector or scalar. The rate parameter of the Gamma prior for the lower channel. If a vector, must have the same number of entries as rows in a and b. Default is 0.
#' @returns a list containing the parameters of the beta-prime distributions at each spatial location, using the parameterization in https://en.wikipedia.org/wiki/Beta_prime_distribution#Generalization
#' @export
zbetaprime <-function(a,b,a1=1,a2=1,b1=0,b2=0){
  #Returns the parameters of the distribution given by propagating uncertainty
  #under the assumption that Z has a Beta-Prime distribution, which is derived by performing
  # a Bayesian update for the mean parameters of the Poisson distributions under the priors
  #\lambda_1 ~ Gamma(a1,b1), \lambda_2 ~ Gamma(a2,b2). This contains as a special case \lambda_i\sim x^{-k_i} by setting ai=1-ki,bi=0.
  #a,b are vectors of data, with a the counts observed in the upper channel (numerator) and b
  # the counts observed in the lower channel (denominator)
  da <- dim(a)[2]
  db <- dim(b)[2]
  if(da>1){
    suma <- rowSums(a,na.rm=TRUE)
    n1 <- rowSums(!is.nan(a))
  }else{
    n1 <-1 
    suma <- a
  }
  if(db>1){
    sumb <- rowSums(b,na.rm=TRUE)
    n2 <- rowSums(!is.nan(b))
  }else{
    n2 <- 1
    sumb <- b
  }
  alpha <- suma+a1#Calculating the parameters
  beta <- sumb+a2#
  p <- beta/beta#Will always be ones most likely
  q <- (b2+n2)/(b1+n1)
  param_list <- list('alpha'=alpha,'beta'=beta,'p'=p,'q'=q)
  return(param_list)
}