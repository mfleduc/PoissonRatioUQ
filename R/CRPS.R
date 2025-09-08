#' Functions for calculating the Continuous Rank Probability Score for predictive models
#' @title Calculate the CRPS for an arbitrary univariate distribution
#' @description Calculates the CRPS of an arbitrary predictive distribution given the predictive PDF or CDF and the observed value
#' @param x vector. the values at which the PDF or CDF are evaluated
#' @param pdf vector. the predictive PDF evaluated at the values in x. Default NULL. One of PDF or CDF must not be NULL
#' @param cdf vector. the predictive CDF evaluated at the values in x. Default NULL. One of PDF or CDF must not be NULL
#' @param xhat scalar. The observed value to compare with the predictive distribution.
#' @returns The continuous rank probability score for the predictive model given the data, calculated via the integral
#' \eqn{\int_{\mathbb{R}}(F(y)-H(y-\hat{x}))^2 \mathrm{d}y}.
#' @export
CRPS <- function( x,pdf=NULL,cdf=NULL,xhat=NULL ){
  #Calculates the Continuous Rank Probability Score
  # for predictive distribution with input pdf/cdf and observation xhat
# Given by \int_{R}(F(y)-H(y-x))^2 dy
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
#' @title Calculate the CRPS for a predictive Gaussian distribution
#' @description Calculate the CRPS for a predictive Gaussian distribution with parameters mu and sigma, and observed value xhat using the closed form presented in
#' "Strictly Proper Scoring Rules, Prediction, and Estimation." by Tilmann Gneiting &Adrian E Raftery
#' @param mu scalar. the mean of the predictive distribution
#' @param sigma scalar. the standard deviation of the predictive distribution
#' @param xhat scalar. The observed value to compare with the predictive distribution.
#' @returns The continuous rank probability score for the predictive model given the data, using the closed form for the univariate Gaussian from "Strictly Proper Scoring Rules, Prediction, and Estimation." by Tilmann Gneiting &Adrian E Raftery
#' @export
CRPSgaussian <- function(mu,sigma,xhat){
  #Calculate the CRPS for a predictive Gaussian distribution with parameters mu and sigma, and observed value xhat
  # Strictly Proper Scoring Rules, Prediction, and Estimation. Tilmann Gneiting &Adrian E Raftery
  z <- (xhat-mu)/sigma# Convert to standard normal
  pdfVal <- dnorm(z, 0 ,1)#Evaluate the PDF
  cdfVal <- pnorm(z, 0, 1)#Evaluate the CDF
  crps <- -1*sigma*(1/sqrt(pi)-2*pdfVal-z*(2*cdfVal-1))
  return(crps)
}
