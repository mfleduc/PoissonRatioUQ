#' Functions for calculating distributions given data
#' @title Calculate the distribution of Z under a Gaussian assumption
#' @description given count datasets a and b, calculates the Gaussian approximation of the distribution of Z = \lambda_a / \lambda_b where \lambda_i is the Poisson intensity of the channel that produced dataset i
#' @param a vector. The count data for the numerator
#' @param b vector. The count data for the denominator
#' @returns a list containing the mean and standard deviation of the approrpiate Gaussian
#' @export 
zgaussian <- function(a,b){
  # Returns the parameters of the distribution given by propagating the uncertainty in the photon counts
  #under the assumption that Z has a normal distribution
  #a = numerator, b = denominator
  suma <- sum(a)
  mna <- mean(a)
  sumb <- sum(b)
  mnb <- mean(b)
  mu <-mna/mnb #Mean of the Gaussian distribution for Z
  sigma_a <- sqrt(suma+0.75)+1
  sigma_b <- sqrt(sumb+0.75)+1 # Gehrels prescription to match 1-sigma distances
  #btwn Poisson and Gaussian distributions
  #N. Gehrels. “Confidence limits for small numbers of events in astrophysical data”. In: The
  #Astrophysical Journal (1986)
  sigma_Z <- mu*sqrt( sigma_a^2/suma^2+sigma_b^2/sumb^2 )
  param_list <- list('mean'=mu, 'stdev'=sigma_Z)
  return(param_list)
}
#' @title Calculate the distribution of Z under a Beta-Prime assumption
#' @description given count datasets a and b, calculates the distribution of Z = \lambda_a / \lambda_b where \lambda_i is the Poisson intensity of the channel that produced dataset i. Under a Gamma prior for each intensity, this is a Beta-Prime distribution
#' @param a vector. The count data for the numerator
#' @param b vector. The count data for the denominator
#' @param k1 scalar. The prior for the intensity of the channel generating the dataset {a}, of the form x^{-k_1}. k_1=0 is an uninformative prior, k_1=1 is scale-invariant
#' @param k2 scalar. The prior for the intensity of the channel generating the dataset {b}, of the form x^{-k_2}. k_2=0 is an uninformative prior, k_2=1 is scale-invariant
#' @returns a list containing the parameters of the beta-prime distribution, using the parameterization in https://en.wikipedia.org/wiki/Beta_prime_distribution#Generalization
#' @export 
zbetaprime <-function(a,b,k1=0,k2=0){
  #Returns the parameters of the distribution given by propagating uncertainty
  #under the assumption that Z has a Beta-Prime distribution, equivalent to performing
  # a Bayesian update for the mean parameters of the Poisson distributions under the priors
  #\lambda_1 ~ x^-k1, \lambda_2 ~ x^-k2
  #a,b are vectors of data, with a the counts observed in the upper channel (numerator) and b
  # the counts observed in the lower channel (denominator)
  suma <- sum(a) #Compiling our sufficient statistics
  sumb <- sum(b)
  n1 <- length(a)
  n2 <- length(b)
  alpha <- suma+1-k1#Calculating the parameters
  beta <- sumb+1-k2#
  p<-1#Will always be this most likely
  q<- n2/n1
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
#' @returns mean and standard deviation of the distribution T|a,b with Z=mT+z0 and T\sim N(\mu_0,\sigma_0^2)
#' @export
tgivenab <- function(a,b,m,z0,tausq,mu0=0,sigma0=Inf,uncertainty="Gaussian"){
  # Calculate the parameters of the distribution T|a,b under the assumption that the
  # ratio Z=\lambda_a/\lambda_b either has a Gaussian distribution or is known exactly.
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
  zparams <-ZGaussian(a,b)
  if(tolower(uncertainty)=="gaussian"){
    muT <- 1/m*(zparams$mean-z0)
    sT <- sqrt(tausq+zparams$stdev^2)/m
    stdev <- 1/sqrt(1/sigma0^2+1/sT^2)
    mn <- (stdev)^2*(mu0/sigma0^2+muT/sT^2)
  }else if(tolower(uncertainty)=="none"){
    stdev <- 1/sqrt(1/sigma0^2+m^2/tausq);
    mn <- (stdev)^2*(mu0/sigma0^2+m/tausq*(zparams$mean-z0))
  }
  param_list<-list('mean'=mn,'stdev'=stdev)
}
