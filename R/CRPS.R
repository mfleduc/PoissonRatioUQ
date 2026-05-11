library(hypergeo)
library(statmod)
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
#' @title Generalized Beta prime CRPS
#' @description Compute the CRPS of a Generalized Beta Prime distribution
#' The GBP pdf is:
#'   f(z) = p/(q*B(alpha,beta)) * (z/q)^(alpha*p-1) / (1+(z/q)^p)^(alpha+beta)
#'Reference: Eq. (15) of LeDuc 2026 The Continuous Rank Probability Score of a Generalized
#'Beta-Prime Distribution and Some Special Cases
#'Due to issues with evaluating generalized hypergeometric functions at unit argument using
#'the hypergeo package, this method computes \eqn{\mathbb{E}[XF(X)]} using quadrature.
#' @param y       Observation (scalar, positive)
#' @param alpha   Shape parameter (alpha > 0)
#' @param beta    Shape parameter (beta > 1/p for finite mean)
#' @param p       Shape parameter (p > 0)
#' @param q       Scale parameter (q > 0)
#' @return CRPS value (scalar)
#' @export
CRPSgbp <- function(alpha, beta, xhat, p=1, q=1) {
  # Check finite mean condition
  y = xhat
  if (beta <= 1/p) stop("Finite mean requires beta > 1/p")
  if (any(c(y, alpha, beta, p, q) <= 0)) stop("All parameters and y must be positive")

  if(abs(beta-1)<1e-6){
    crps = crps_dagum(y,alpha,p,q)
  }else if(abs(alpha-1)<1e-6){
    crps = crps_singh_maddala(y,beta,p,q)
  }else{
    w = y^p / (q^p + y^p)
    # Mean of GBP: mu = q * B(alpha + 1/p, beta - 1/p) / B(alpha, beta)
    mu = q * exp(lbeta(alpha + 1/p, beta - 1/p) - lbeta(alpha, beta))
    # --- E[|X - y|] component ---
    # Incomplete beta B(w; alpha, beta) = betainc * B(alpha, beta)
    # pbeta returns regularized incomplete beta I_w(alpha, beta)
    inc_beta <- exp(pbeta(w, alpha, beta, log.p = TRUE) + lbeta(alpha, beta))
    # 2F1(1, alpha+beta-1; alpha+1/p; w)
    hyp2f1 <- Re(hypergeo(1, alpha + beta - 1, alpha + 1/p, w))

    E_abs <- mu - y + (2 / beta(alpha, beta)) * (
      y * inc_beta +
        y / (alpha + beta - 1) * w^(alpha - 1) * (1 - w)^beta * (1 - hyp2f1)
    )
    B_term <- beta(2*alpha + 1/p, 2*beta - 1/p)
    E_2xfx <- e_2xfx_gq(alpha, beta, p, q )
    crps <- E_abs - E_2xfx + mu
  }
  return(crps)
}
#'@title Incomplete Beta function
#'@description Evaluates the incomplete beta function at nodes w with parameters a,b
#'@param w Locations to evaluate
#'@param a the a parameter
#'@param b the b parameter
#'@returns the value of the incomplete beta function at w
#'@export
inc_beta <- function(w, a, b){pbeta(w, a, b) * beta(a, b)}
#'@title Numerical evaluation of \eqn{\mathbb{E}[2XF(X)]} for the generalized Beta-prime dsitribution
#'@description Uses a quadrature rule on [0,1] to compute the value of the expectation
#'\eqn{\mathbb{E}[2XF(X)]} for use in calculating the CRPS of the generalized Beta-prime. This is done
#'because the package hypergeo() can be poorly behaved when computing generalized
#'hypergeometric functions at argument 1
#'@param alpha scalar, the alpha parameter
#'@param beta scale, the beta parameter
#'@param p the p parameter, default 1
#'@param q the q parameter, default 1
#'@param quadrule a quadrature rule on [0,1] generated by statmod::gauss.quad. Default NULL, in which case it is generated in the function
#'@param n the number of nodes on the quadrature rule, default is \eqn{ceil(min(64, \sqrt(\alpha+\beta))}
#'@export
e_2xfx_gq <- function(alpha, beta, p=1, q=1, quadrule=NULL, n = 64) {
  # Jacobi weight exponents
  a <- alpha + 1/p   # left exponent
  b <- beta  - 1/p   # right exponent
  # Get nodes and weights for Beta(a, b) weight function
  # gauss.quad.prob returns weights that already include the 1/B(a,b) factor,
  # so the approximation is:
  #   int_0^1 u^(a-1)*(1-u)^(b-1)*g(u) du ≈ B(a,b) * sum(w_i * g(x_i))
  if(is.null(quadrule)){
    quadrule = statmod::gauss.quad.prob(round(max(n, 4*sqrt(a+b))), dist = "beta", alpha = a , beta = b )
  }
  nodes   <- quadrule$nodes
  weights <- quadrule$weights
  # Smooth function at nodes: g(u) = B(u; alpha, beta) = pbeta * beta(alpha,beta)
  g_vals <- inc_beta(nodes, alpha, beta)
  # Integral = B(a,b) * sum(w_i * g(x_i))
  integral <- exp(lbeta(a, b)) * sum(weights * g_vals)
  (2 * q / beta(alpha, beta)^2) * integral
}




