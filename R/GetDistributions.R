#' Functions for calculating distributions given data
#' @title Calculate the distribution of Z under a Gaussian assumption
#' @description given count datasets a and b, calculates the Gaussian approximation of the distribution of Z = lambda_a / lambda_b, where lambda_i is the Poisson intensity of the channel that produced dataset i, using the delta method. Assumes a Gamma prior on Poisson data to get the distribution of the intensities, and then uses the delta method to get the approximate distribution of the ratio.
#' This model assumes that realizations at different spatial locations are independent.
#' @param a Matrix. The count data for the numerator. Rows correspond to spatial locations, columns to realizations. Must have the same number of rows as b, but need not have the same number of columns. MUST be a matrix or array right now, if you want a scalar use a 1x1 array. Missing data should be replaced by NaNs.
#' @param b Matrix. The count data for the denominator. Rows correspond to spatial locations, columns to realizations. Must have the same number of rows as b, but need not have the same number of columns. MUST be a matrix or array right now, if you want a scalar use a 1x1 array. Missing data should be replaced by NaNs.
#' @param a_prior list/array. The shape and rate parameters of the Gamma prior distribution for the Poisson intensity of the numerator process. Default is (1,0), i.e. uniform prior.
#' @param b_prior list/array. The shape and rate parameters of the Gamma prior distribution for the Poisson intensity of the denominator process. Default is (1,0), i.e. uniform prior.
#' @returns a list containing the means and standard deviations of the appropriate Gaussians.
#' @export


zgaussian <- function(a,b, a_prior=c(1,0),b_prior=c(1,0)){
  # Returns the parameters of the distribution given by propagating the uncertainty in the photon counts
  #under the assumption that Z has a normal distribution
  #a = numerator, b = denominator
  da <- dim(a)[2]
  db <- dim(b)[2]
  if(da>1){
    suma <- rowSums(a,na.rm=TRUE)
    na <- rowSums(!is.nan(a))
  }else{
    na <-1 
    suma <- a
  }
  num_alpha <- suma+a_prior[1]
  num_beta <- na+a_prior[2]
  # mna <- rowMeans(a,na.rm=TRUE)
  if(db>1){
    sumb <- rowSums(b,na.rm=TRUE)
    nb <- rowSums(!is.nan(b))
  }else{
    nb <- 1
    sumb <- b
  }
  denom_alpha <- sumb+b_prior[1]
  denom_beta <- nb+b_prior[2]
  ## Calculate mean
  mu <- (num_alpha/num_beta)*(denom_beta/denom_alpha) #Mean of the Gaussian distribution for Z

  stdev <- mu*sqrt( 1/num_alpha+1/denom_alpha )
  param_list <- list('mean'=mu, 'stdev'=sigma_Z)
  return(param_list)
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
  #under the assumption that Z has a Beta-Prime distribution, equivalent to performing
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
#' @title Determine the parameters of the distribution of the temperature given the data under the assumption that Z is either a point estimate or has a Gaussian distribution, and is an affine function of T
#' @description Given count data, a prior, and the parameters of the distribution Z|T ~N(mT+z0,tausq) determines the distribution of T|data under the assumption that the counts are either Gaussian or provide a point estimate of Z.
#' @param a Matrix. The count data for the numerator. Rows correspond to spatial locations/channels, columns to realizations. Must have the same number of rows as b, but need not have the same number of columns. MUST be a matrix or array right now, if you want a scalar use a 1x1 array. Missing data should be replaced by NaNs.
#' @param b Matrix. The count data for the denominator. Rows correspond to spatial locations/channels, columns to realizations. Must have the same number of rows as b, but need not have the same number of columns. MUST be a matrix or array right now, if you want a scalar use a 1x1 array. Missing data should be replaced by NaNs.
#' @param M scalar or matrix. slope parameter of the linear regression of Z(T).
#' @param Z0 scalar or vector. Intercept parameter of the linear regression of Z(T).
#' @param TauSq scalar or matrix. Covariance matrix of the residuals of the linear regression.
#' @param priormn scalar or vector. Mean of Gaussian prior for the temperature. Default is 0.
#' @param priorvar scalar or matrix. variance of the Gaussian prior for the temperature. Default is Inf (uninformative prior).
#' @param uncertainty string. Either "Gaussian", in which case Z is assumed to have a Gaussian distribution determined by the function zgaussian, or "none", in which case mean(a)/mean(b) is used as a point estimate of Z. Default is "Gaussian"
#' @returns mean and standard deviation of the distribution T|a,b with Z=mT+z0, T|Z ~ N((a/b-z0)/m, tau^2/m), Z|a,b from zgaussian(), and T ~ N(mu0,Sigma0) where mu0 is the prior mean and Sigma0 the prior covariance matrix.
#' @export
tgivenab <- function(a,b,M,Z0,TauSq,priormn=0,priorvar=Inf,uncertainty="Gaussian"){
  # Calculate the parameters of the distribution T|a,b under the assumption that the
  # ratio Z=lambda_a/lambda_b either has a Gaussian distribution or is known exactly.
  # a and b are the data for the upper and lower channels respectively.
  # Additional assumption that Z|T \sim N(mT+z_0,tausq) and that the prior for T is
  # N(mu0,sigma0^2). These parameters are optional.
  # Setting uncertainty='none' is equivalent to saying that \lambda_a and \lambda_b are known exactly
  # and so no uncertainty should be incorporated from that. This is equivalent to the
  # algorithms of Cantrall et al and Zhang et al discussed in the associated paper.
  # Those algorithms also do not incorporate a prior distribution for T.
  #
  # By default, the assumption is that Z is Gaussian and there is no prior distribution for T.
  #
  # To obtain the distribution T|z, use the command
  # params <- tgivenab(z,1,m,z0,tausq,mu0=mu0,sigma0=sigma0,uncertainty='none')
  #
  stopifnot(tolower(uncertainty)=="gaussian"||tolower(uncertainty)=="none")#What kind of uncertainty are we going to include?
  includeprior <- all(is.finite(priorvar))
  zparams <- zgaussian(a,b)
  if(length(M)==1){
    nSpatial <- dim(a)[1] #Expanding into a spatial model
    M <- diag(nSpatial)*M
  }else if(min(dim(M))==1){
    M <- diag(M)
    nSpatial <- dim(M)[1]
  }else{
    nSpatial <- dim(M)[2]
  }
  if(length(Z0)==1){
    Z0 <- matrix(z0,nrow=dim(a)[1],ncol=1)
  }
  if(length(TauSq)==1){
    TauSq <- eye(dim(a)[1])*TauSq
  }
  TauSqInv <- solve(TauSq)
  if(length(priormn)==1){
    priormn <- matrix(priormn,nrow=nSpatial,ncol=1)
  }
  if(length(priorvar)==1&includeprior){
    priorvar <- diag(nSpatial)*priorvar
  }
  MtTsInv <- t(M)%*%TauSqInv
  if(tolower(uncertainty)=="gaussian"){
    # muT <- 1/m*(zparams$mean-z0)
    # sT <- sqrt(tausq+zparams$stdev^2)/m
    # stdev <- 1/sqrt(1/sigma0^2+1/sT^2)
    # mn <- (stdev)^2*(mu0/sigma0^2+muT/sT^2)
    # sT <- 1/sqrt(1/priorvar+m^2/tausq)
    # variance <- sT^2*(sT^2*zparams$stdev^2*m^2+tausq^2)/tausq^2
    # stdev <- sqrt(variance)
    # mn <- sT^2*(priormn/priorvar+m/tausq*(zparams$mean-z0))
    if(includeprior){
      priorvarinv <- solve(priorvar)
      SigmaT <- solve(priorvarinv+MtTsInv%*%M)
      mn <- SigmaT%*%(priorvarinv%*%priormn+(MtTsInv%*%(zparams$mean-Z0)))
      SigmaZsqrt <- diag(as.matrix(zparams$stdev))
      H <- SigmaT%*%MtTsInv%*%SigmaZsqrt
      cov <- SigmaT + H%*%t(H)
    }else{
      SigmaT <- solve(MtTsInv%*%M)
      mn <- SigmaT%*%(MtTsInv%*%(zparams$mean-Z0))
      SigmaZsqrt <- diag(as.matrix(zparams$stdev))
      H <- SigmaT%*%MtTsInv%*%SigmaZsqrt
      cov <- SigmaT + H%*%t(H)
    }
  }else if(tolower(uncertainty)=="none"){
    # stdev <- 1/sqrt(1/priorvar+m^2/tausq);
    if(includeprior){
      priorvarinv <- solve(priorvar)
      cov <- solve(priorvarinv+MtTsInv%*%M)
      mn <- cov%*%(priorvarinv%*%priormn+MtTsInv%*%(zparams$mean-Z0))
    }else{
      cov <- solve(MtTsInv%*%M)
      mn <- cov%*%(MtTsInv%*%(zparams$mean-Z0))
    }
  }
  param_list<-list('mean'=mn,'cov'=cov)
  return(param_list)

}
