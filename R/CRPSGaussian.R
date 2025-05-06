CRPSGaussian <- function(mu,sigma,xhat){
  #Calculate the CRPS for a predictive Gaussian distribution with parameters mu and sigma, and observed value xhat
  # Strictly Proper Scoring Rules, Prediction, and Estimation. Tilmann Gneiting &Adrian E Raftery
  z <- (xhat-mu)/sigma# Convert to standard normal
  pdfVal <- dnorm(z, 0 ,1)#Evaluate the PDF
  cdfVal <- pnorm(z, 0, 1)#Evaluate the CDF
  crps <- -1*sigma*(1/sqrt(pi)-2*pdfVal-z*(2*cdfVal-1))
  return(crps)
}
