#' Functions for calculating distributions given data
#' @title Calculate the distribution of Z under a Gaussian assumption
#' @description given count datasets a and b, calculates the Gaussian approximation of the distribution of Z = lambda_a / lambda_b where lambda_i is the Poisson intensity of the channel that produced dataset i. Some details are available in (Park, T., Kashyap, V. L., Siemiginowska, A., van Dyk, D. A., Zezas, A., Heinke, C., and Wargelin, B. J.: Bayesian Estimation of Hardness
#' Ratios: Modeling and Computations, In: The
#' Astrophysical Journal, 2006) and (Gehrels, N.: Confidence limits for small numbers of events in astrophysical data, The Astrophysical Journal, 1986.). This model assumes that realizations at different spatial locations are independent!
#' @param a Matrix. The count data for the numerator. Rows correspond to spatial locations, columns to realizations. Must have the same number of rows as b, but need not have the same number of columns. MUST be a matrix or array right now, if you want a scalar use a 1x1 array. Missing data should be replaced by NaNs.
#' @param b Matrix. The count data for the denominator. Rows correspond to spatial locations, columns to realizations. Must have the same number of rows as b, but need not have the same number of columns. MUST be a matrix or array right now, if you want a scalar use a 1x1 array. Missing data should be replaced by NaNs.
#' @returns a list containing the means and standard deviations of the appropriate Gaussians.
#' @export
zgaussian <- function(a,b){
  # Returns the parameters of the distribution given by propagating the uncertainty in the photon counts
  #under the assumption that Z has a normal distribution
  #a = numerator, b = denominator
  suma <- rowSums(a,na.rm=TRUE)
  mna <- rowMeans(a,na.rm=TRUE)
  sumb <- rowSums(b,na.rm=TRUE)
  mnb <- rowMeans(b,na.rm=TRUE)
  mu <- mna/mnb #Mean of the Gaussian distribution for Z
  sigma_a <- sqrt(abs(suma)+0.75)+1
  sigma_b <- sqrt(abs(sumb)+0.75)+1 # Gehrels prescription to match 1-sigma distances
  #btwn Poisson and Gaussian distributions
  #N. Gehrels. “Confidence limits for small numbers of events in astrophysical data”. In: The
  #Astrophysical Journal (1986)
  sigma_Z <- abs(mu)*sqrt( sigma_a^2/suma^2+sigma_b^2/sumb^2 )
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
  suma <- rowSums(a,na.rm=TRUE)
  n1 <- rowSums(!is.nan(a))
  sumb <- rowSums(b,na.rm=TRUE)
  n2 <- rowSums(!is.nan(b))
  alpha <- suma+a1#Calculating the parameters
  beta <- sumb+a2#
  p <- beta/beta#Will always be ones most likely
  q <- (b2+n2)/(b1+n1)
  param_list <- list('alpha'=alpha,'beta'=beta,'p'=p,'q'=q)
  return(param_list)
}
#' @title Determine the parameters of the distribution of the temperature given the data under the assumption that Z is either a point estimate or has a Gaussian distribution, and is an affine function of T
#' @param a Matrix. The count data for the numerator. Rows correspond to spatial locations, columns to realizations. Must have the same number of rows as b, but need not have the same number of columns. MUST be a matrix or array right now, if you want a scalar use a 1x1 array. Missing data should be replaced by NaNs.
#' @param b Matrix. The count data for the denominator. Rows correspond to spatial locations, columns to realizations. Must have the same number of rows as b, but need not have the same number of columns. MUST be a matrix or array right now, if you want a scalar use a 1x1 array. Missing data should be replaced by NaNs.
#' @param m scalar. slope parameter of the linear regression of Z(T). Must be a scalar for now, will be generalized to a vector in the future. Right now vectors should work but have not been tested.
#' @param z0 scalar. Intercept parameter of the linear regression of Z(T). Must be a scalar for now, will be generalized to a vector in the future. Right now vectors should work but have not been tested.
#' @param tausq scalar. Variance of the residuals of the linear regression. Must be a scalar for now, will be generalized to a vector in the future. Right now vectors should work but have not been tested.
#' @param priormn scalar or vector. Mean of Gaussian prior for the temperature. Default is 0.
#' @param priorvar scalar or matrix. variance of the Gaussian prior for the temperature. Default is Inf (uninformative prior).
#' @param uncertainty string. Either "Gaussian", in which case Z is assumed to have a Gaussian distribution determined by the function zgaussian, or "none", in which case mean(a)/mean(b) is used as a point estimate of Z. Default is "Gaussian"
#' @returns mean and standard deviation of the distribution T|a,b with Z=mT+z0, T|Z ~ N((a/b-z0)/m, tau^2/m), Z|a,b from zgaussian(), and T ~ N(mu0,Sigma0) where mu0 is the prior mean and Sigma0 the prior covariance matrix.
#' @export
tgivenab <- function(a,b,m,z0,tausq,priormn=0,priorvar=Inf,uncertainty="Gaussian"){
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
  nSpatial <- dim(a)[1] #Expanding into a spatial model
  M <- TauSq <- TauSqInv <- diag(nSpatial)
  diag(M) <- m
  Z0 <-  matrix(z0,nrow=nSpatial,ncol=1)
  diag(TauSq) <- tausq
  diag(TauSqInv) <- (1/tausq)
  if(length(priormn)==1){
    priormn <- matrix(priormn,nrow=nSpatial,ncol=1)
  }
  if(length(priorvar)==1&includeprior){
    priorvar <- diag(nSpatial)*priorvar
  }
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
      SigmaT <- solve(solve(priorvar)+M%*%(TauSqInv%*%M))
      mn <- SigmaT%*%(solve(priorvar, priormn)+(TauSqInv%*%M%*%(zparams$mean-Z0)))
      SigmaZ <- diag(nSpatial)
      diag(SigmaZ) <- zparams$stdev^2
      cov <- SigmaT + SigmaT%*%M%*%TauSqInv%*%SigmaZ%*%TauSqInv%*%M%*%SigmaT
    }else{
      SigmaT <- solve(M%*%(TauSqInv%*%M))
      mn <- SigmaT%*%((TauSqInv%*%M%*%(zparams$mean-Z0)))
      SigmaZ <- diag(nSpatial)
      diag(SigmaZ) <- zparams$stdev^2
      cov <- SigmaT + SigmaT%*%M%*%TauSqInv%*%SigmaZ%*%TauSqInv%*%M%*%SigmaT
    }
  }else if(tolower(uncertainty)=="none"){
    # stdev <- 1/sqrt(1/priorvar+m^2/tausq);
    if(includeprior){
      cov <- solve(solve(priorvar)+M%*%(TauSqInv%*%M))
      mn <- cov%*%(solve(priorvar, priormn)+(TauSqInv%*%M%*%(zparams$mean-Z0)))
    }else{
      cov <- solve(M%*%(TauSqInv%*%M))
      mn <- cov%*%((TauSqInv%*%M%*%(zparams$mean-Z0)))
    }
  }
  param_list<-list('mean'=mn,'cov'=cov)
}
