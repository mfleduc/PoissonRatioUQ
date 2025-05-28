#' Functions for calculating distributions given data
#' @title Calculate the distribution of Z under a Gaussian assumption
#' @description given count datasets a and b, calculates the Gaussian approximation of the distribution of Z = lambda_a / lambda_b where lambda_i is the Poisson intensity of the channel that produced dataset i. Some details are available in (Park, T., Kashyap, V. L., Siemiginowska, A., van Dyk, D. A., Zezas, A., Heinke, C., and Wargelin, B. J.: Bayesian Estimation of Hardness
#' Ratios: Modeling and Computations, In: The
#' Astrophysical Journal, 2006) and (Gehrels, N.: Confidence limits for small numbers of events in astrophysical data, The Astrophysical Journal, 1986.)
#' @param a vector. The count data for the numerator
#' @param b vector. The count data for the denominator
#' @returns a list containing the mean and standard deviation of the appropriate Gaussian
#' @export
zgaussian <- function(a,b){
  # Returns the parameters of the distribution given by propagating the uncertainty in the photon counts
  #under the assumption that Z has a normal distribution
  #a = numerator, b = denominator
  suma <- sum(a)
  mna <- mean(a)
  sumb <- sum(b)
  mnb <- mean(b)
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
#' @param a vector. The count data for the numerator
#' @param b vector. The count data for the denominator
#' @param a1 scalar. The shape parameter of the Gamma prior for the upper channel. Default is 1.
#' @param b1 scalar. The rate parameter of the Gamma prior for the upper channel. Default is 0.
#' @param a2 scalar. The shape parameter of the Gamma prior for the lower channel. Default is 1.
#' @param b2 scalar. The rate parameter of the Gamma prior for the lower channel. Default is 0.
#' @returns a list containing the parameters of the beta-prime distribution, using the parameterization in https://en.wikipedia.org/wiki/Beta_prime_distribution#Generalization
#' @export
zbetaprime <-function(a,b,a1=1,a2=1,b1=0,b2=0){
  #Returns the parameters of the distribution given by propagating uncertainty
  #under the assumption that Z has a Beta-Prime distribution, equivalent to performing
  # a Bayesian update for the mean parameters of the Poisson distributions under the priors
  #\lambda_1 ~ Gamma(a1,b1), \lambda_2 ~ Gamma(a2,b2). This contains as a special case \lambda_i\sim x^{-k_i} by setting ai=1-ki,bi=0.
  #a,b are vectors of data, with a the counts observed in the upper channel (numerator) and b
  # the counts observed in the lower channel (denominator)
  suma <- sum(a) #Compiling our sufficient statistics
  sumb <- sum(b)
  n1 <- length(a)
  n2 <- length(b)
  alpha <- suma+a1#Calculating the parameters
  beta <- sumb+a2#
  p<-1#Will always be this most likely
  q<- (b2+n2)/(b1+n1)
  param_list <- list('alpha'=alpha,'beta'=beta,'p'=p,'q'=q)
  return(param_list)
}
#' @title Determine the parameters of the distribution of the temperature given the data under the assumption that Z is either a point estimate or has a Gaussian distribution, and is an affine function of T
#' @param a vector. The count data for the numerator
#' @param b vector. The count data for the denominator
#' @param m scalar. slope parameter of the linear regression of Z(T)
#' @param z0 scalar. Intercept parameter of the linear regression of Z(T)
#' @param tausq scalar. Variance of the residuals of the linear regression
#' @param mu0 scalar. Mean of Gaussian prior for the temperature. Default is 0.
#' @param sigma0 scalar. Standard deviation of the Gaussian prior for the temperature. Default is Inf (uninformative prior)
#' @param uncertainty string. Either "Gaussian", in which case Z is assumed to have a Gaussian distribution determined by the function zgaussian, or "none", in which case mean(a)/mean(b) is used as a point estimate of Z. Default is "Gaussian"
#' @returns mean and standard deviation of the distribution T|a,b with Z=mT+z0 and T ~ N(mu0,sigma0^2)
#' @export
tgivenab <- function(a,b,m,z0,tausq,mu0=0,sigma0=Inf,uncertainty="Gaussian"){
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
  zparams <- zgaussian(a,b)
  if(tolower(uncertainty)=="gaussian"){
    # muT <- 1/m*(zparams$mean-z0)
    # sT <- sqrt(tausq+zparams$stdev^2)/m
    # stdev <- 1/sqrt(1/sigma0^2+1/sT^2)
    # mn <- (stdev)^2*(mu0/sigma0^2+muT/sT^2)
    sT <- 1/sqrt(1/sigma0^2+m^2/tausq)
    variance <- sT^2*(sT^2*zparams$stdev^2*m^2+tausq^2)/tausq^2
    stdev <- sqrt(variance)
    mn <- sT^2*(mu0/sigma0^2+m/tausq*(zparams$mean-z0))

  }else if(tolower(uncertainty)=="none"){
    stdev <- 1/sqrt(1/sigma0^2+m^2/tausq);
    mn <- (stdev)^2*(mu0/sigma0^2+m/tausq*(zparams$mean-z0))
  }
  param_list<-list('mean'=mn,'stdev'=stdev)
}
