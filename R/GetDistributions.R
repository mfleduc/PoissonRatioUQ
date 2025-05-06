ZGaussian <- function(a,b){
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
ZBetaPrime <-function(a,b,k1=0,k2=0){
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
