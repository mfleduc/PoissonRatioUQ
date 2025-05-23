#' Functions for computing the highest posterior density set of a given distribution
#'@title Highest-Posterior density
#'@description Calculates the highest-posterior density set of an abirtrary univariate distribution given a PDF fx evaluated at points x
#'@param x vector. x value of x at which the PDF is evaluated
#'@param fx vector. Values of the PDF at the points in x
#'@param alpha scalar. The desired credibility value, or the desired value of the integral of the pdf over the HPD set.
#'@returns The values of x within the alpha-highest posterior density set
#'@export
hpdset <- function(x, fx , alpha){
  #Compute the HPD alpha-credible set for an arbitrary distribution.
  # Letting $h$ be the solution to the problem
  #\[
  # \int_{f>h}f(x)dx = \alpha
  #\]
  #the HPD credible set is given by the set {x:f(x)>h}
  #Univariate functions only, probably will work best with unimodal functions
  #First: Check normalization
  fx[is.na(fx)]<- 0
  # f <- approxfun(x, fx, method="linear", yleft=0, yright=0, rule=2) #Don't think this is actually necessary
  intF <- pracma::trapz(x,fx)
  fx <- fx/intF
  #Now: Find the mode
  modefx <- max(fx, na.rm=TRUE)
  if(!is.finite(modefx)){
    return(c(NaN,NaN))
  }
  #Function to minimize
  GetCredibleSet <-function(h){
    mask <- (fx<=h)
    fxm <- fx
    fxm[mask] <- 0
    intVal <- pracma::trapz(x,fxm)
    return(intVal-alpha)
  }
  #Do rootfinding
  rootfind <-uniroot(GetCredibleSet, c(0, modefx ), tol=10^-12)
  hval<- rootfind$root
  mask <- fx>=hval
  #Return the appropriate values of x
  return(x[mask])
}
#'@title Highest-Posterior density for a Gaussian distribution
#'@description Calculates the highest-posterior density set of the univariate Gaussian with mean mu and standard deviation sigma
#'@param mu Scalar. Mean of the Gaussian distribution
#'@param sigma scalar. Standard deviation of the Gaussian distribution
#'@param alpha scalar. The desired credibility value, or the desired value of the integral of the pdf over the HPD set.
#'@returns The endpoints of the alpha-HPD interval.
#'@export
hpdintervalgaussian<-function(mu,sigma,alpha){
  #Returns the endpoints of the highest posterior density alpha-credible interval
  #for a Gaussian distribution with known mean and standard deviation
  halfwidth <- qnorm((1-alpha)/2,0,1) #Note: Halfwidth <0!!!!!
  interval <- c(  mu+sigma*halfwidth,mu-sigma*halfwidth)
  return(interval)
}
