#' Functions for doing estimation of the parameters of a permanental process model using equivalent kernel representations
#' @title Estimation of point process intensities using a Permanental Process model
#' @description Given a realization of a point process, estimates the underlying intensity function by assuming that the process
#' intensity has the form \eqn{\lambda(s) = \frac{c}{2}f(s)^2}, with \eqn{f(s)} a Gaussian process. The estimation is regularized by assuming that \eqn{f(s)} lies in a RKHS given by the kernel \eqn{k(x,y)}, and the likelihood function is given by
#' \eqn{L = \sum_{i=1}^d a_i \log(\frac{c}{2}\int_{R_i}f(s)^2\mathrm{d}\mu(s))- \frac{c}{2}||f||_2^2-\frac{g}{2}||f||_{H_k}^2}. Note that if the process is not binned, all counts are equal to 1 and the integral is replaced by evaluation of \eqn{f(s)} at the location of the observation. The algorithm is a slightly generalized version of the one outlined in https://proceedings.mlr.press/v70/walder17a.html (TODO: Add link to my work when it is out)
#' @param K matrix: Kernel matrix evaluated at the data points
#' @param counts numeric: If the data represents a binned process, the number of counts in each bin. Otherwise either the total number of points observed or you can just leave it empty and the code will infer the number of points from K.
#' @param g scalar: Regularization strength. Default is 1
#' @param c scalar: Parameter such that intensity \eqn{\lambda(s) = \frac{c}{2}f(s)^2}, default is 1.
#' @returns The estimated intensity function, the parameters of the gamma approximation to its posterior at each location, and the fitted kernel coefficients of \eqn{f(s) = \sum_i\alpha_i\tilde{k}(s,s_i)}
#' @export
permproccest <- function(K,counts=NaN,g=1,c=1){
  #Gather up all the needed stuff to do the calculation: Parameters, data, etc
  if(length(counts)==1){
    #Allows counts to be the number of spatial locations, the actual counts, or left alone and inferred from K
    if(!is.nan(counts)){
      counts <- array(1,dim=c(1,counts))
    }else{
      counts <- array(1, dim=c(1, dim(K)[1]))
    }
  }
  #First: Eigendecomposition of the kernel matrix
  disp("Generating equivalent kernel matrix...")
  eigstuff <- eigen(K, symmetric=TRUE )
  eigvals_uu <- eigstuff$values
  eigvecs_uu <- eigstuff$vectors
  KTxx <- eigvecs_uu%*%diag(eigvals_uu/(c*eigvals_uu+g))%*%t(eigvecs_uu)
  # Now: Determine alpha
  disp("Determining f...")
  a0 <- array(1, dim=c(length(counts),1))
  optval <- optim( (a0), function(X)objectivefn(X,counts,KTxx,c),gr=function(X)gradientFn(X,counts,KTxx)
                   , method="CG",control=list(maxit=1000))
  alphahat <- optval$par
  #Get intensity
  f <- KTxx%*%alphahat
  intensity <- 0.5*c*f^2
  #Parameters of the posterior distributions
  Dinv <- diag(2*counts/alphahat^2)
  Sigma <- KTxx - KTxx%*%solve(Dinv+KTxx, KTxx)#Posterior covariance of f
  predvar <- diag(Sigma)
  gamma_alpha <- (f^2+predvar)^2/(2*predvar*(2*f^2+predvar))
  gamma_beta <- (f^2+predvar)/(predvar*c*(2*f^2+predvar))
  outputs <- list("lambda"=intensity, "alphaparam" = gamma_alpha, "betaparam"= gamma_beta, "kernelcoeffs" = alphahat)
  return(outputs)
}
#' @title Gradient of the RKHS-regularized Poisson point process likelihood with respect to the kernel coefficients
#' @description Gradient of the RKHS-regularized Poisson point process likelihood with respect to the kernel coefficients
#' @param alpha vector: Kernel coefficients \eqn{\vec{\alpha}} so that \eqn{f(x)= \sum \alpha_i \tilde{k}(x_i,x)}
#' @param counts vector: Count data, number observed at each observation point
#' @param K matrix: Kernel matrix for the equivalent kernel function \eqn{\tilde{k}(x,y)}
#' @returns the gradient of the objective function at alpha
#' @export
gradientFn <- function(alpha, counts, K){
  val <- zeros(length(alpha),1)
  for(ii in 1:length(counts)){
    val <- val-2*counts[ii]*K[,ii]/(K%*%alpha)
  }
  val <- val+K%*%alpha
  return(val)
}
#' @title Objective function of the RKHS-regularized Poisson point process likelihood with respect to the kernel coefficients
#' @description Objective function of the RKHS-regularized Poisson point process likelihood with respect to the kernel coefficients
#' @param alpha vector: Kernel coefficients alpha so that \eqn{f(x)= \sum \alpha_i \tilde{k}(x_i,x)}
#' @param counts vector: Count data, number observed at each observation point
#' @param K matrix: Kernel matrix for the kernel function \eqn{\tilde{k}}
#' @param c scalar: Value so that \eqn{\lambda(s)=\frac{c}{2}f(s)^2}
#' @returns the value of the objective function at alpha
#' @export
objectivefn <- function(alpha, counts, K,c){
  val <- -1*sum(counts*log((c/2)*(K%*%alpha)^2))+0.5*t(alpha)%*%K%*%(alpha)
  return(val)
}
