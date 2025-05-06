library(pracma)
CRPS <- function( x,pdf=NULL,cdf=NULL,xhat=NULL ){
  #Calculates the Continuous Rank Probability Score
  # for predictive distribution with input pdf/cdf and observation xhat
# Given by \int_{R}(F(y)-H(x-y))^2 dx
if(is.null(cdf)){
  stopifnot(!is.null(pdf))
  cdf <- numeric(length(x))

  f <- approxfun(x, pdf, method="linear", yleft=0, yright=0, rule=2)
  cdf[1]=0
  for(ii in 2:length(cdf)){
    cdf[ii] <- cdf[ii-1]+integrate(f,x[ii-1],x[ii])$value
  }
  cdf <- cdf/tail(cdf,1)
}
Hxmy <- ifelse( x>=xhat,1,0 )#Heaviside(xhat-y)
#approxIntegrand <- approxfun( x, (cdf-Hxmy)^2, method="linear", yleft=0, yright=0, rule=2 )
crps <- pracma::trapz(x, (cdf-Hxmy)^2 )

return(crps)
}
