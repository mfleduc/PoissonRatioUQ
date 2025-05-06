HPDSet <- function(x, fx , alpha){
  #Compute the HPD alpha-credible set for an arbitrary distribution.
  # Letting $h$ be the solution to the problem
  #\[
  # \int_{f>h}f(x)dx = \alpha
  #\]
  #the HPD credible set is given by the set {x:f(x)>h}
  #Univariate functions only, probably will work best with unimodal functions
  #First: Check normalization
  f <- approxfun(x, fx, method="linear", yleft=0, yright=0, rule=2)
  intF <- integrate(f, min(x), max(x))$value
  fx <- fx/intF
  #Now: Find the mode
  modefx <- max(fx)
  #Function to minimize
  GetCredibleSet <-function(h){
    mask <- (fx<=h)
    fxm <- fx
    fxm[mask] <- 0
    fx_approx <- approxfun(x, fxm, method="linear", yleft=0, yright=0, rule=2)
    intVal <- integrate(fx_approx,min(x),max(x))$value
    return(intVal-alpha)
  }
  #Do rootfinding
  rootfind <-uniroot(GetCredibleSet, c(0, modefx), tol=10^-12)
  hval<- rootfind$root
  mask <- fx>=hval
  #Return the appropriate values of x
  return(x[mask])
}
