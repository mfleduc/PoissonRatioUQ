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
  f <- approxfun(x, fx, method="linear", yleft=0, yright=0, rule=2)
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
  rootfind <-uniroot(GetCredibleSet, c(0, modefx/1.05), tol=10^-12)
  hval<- rootfind$root
  mask <- fx>=hval
  #Return the appropriate values of x
  return(x[mask])
}
hpdintervalgaussian<-function(mu,sigma,alpha){
  #Returns the endpoints of the highest posterior density alpha-credible interval
  #for a Gaussian distribution with known mean and standard deviation
  halfwidth <- qnorm((1-alpha)/2,0,1)
  interval <- c(  mu+sigma*halfwidth,mu-sigma*halfwidth)
  return(interval)
}
