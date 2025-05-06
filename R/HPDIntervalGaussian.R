HPDIntervalGaussian<-function(mu,sigma,alpha){
  #Returns the endpoints of the highest posterior density alpha-credible interval
  #for a Gaussian distribution with known mean and standard deviation
  halfwidth <- qnorm((1-alpha)/2,0,1)
  interval <- c(  mu+sigma*halfwidth,mu-sigma*halfwidth)
  return(interval)
}
