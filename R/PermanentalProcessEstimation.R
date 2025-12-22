#' Functions for doing estimation of the parameters of a permanental process model using equivalent kernel representations
#' @title Estimation of point process intensities using a Permanental Process model
#' @description Given Possion data, estimates the underlying intensity by assuming that it has the form \eqn{\vec{\Lambda} = \frac{c}{2}\vec{f}^2}, with \eqn{\vec{f}} a Gaussian vector. The algorithm recovers a latent vector \eqn{\vec{f}} under the prior \eqn{\vec{f}\sim\mathcal{N}(0,g^{-1}K)}, and the result is given by \eqn{\vec{f}=\tilde{K}\vec{\psi}} for some \eqn{\vec{\psi}} to be retrieved. The algorithm is adapted from the one outlined in https://proceedings.mlr.press/v70/walder17a.html and was used in https://egusphere.copernicus.org/preprints/2025/egusphere-2025-5570/ to estimate thermospheric temperature
#' @param K matrix: Kernel matrix evaluated at the spatial locations
#' @param counts numeric: The number of counts in each bin.
#' @param g scalar: Regularization strength (or prior precision). Default is 1
#' @param c scalar: Parameter such that intensity \eqn{\vec{\Lambda} = \frac{c}{2}\vec{f}^2}, default is 1.
#' @param maxiter scalar: The maximum number of iterations.
#' @returns The estimated intensities, the parameters of the gamma approximation to its posterior at each location, and the fitted kernel coefficients of \eqn{\hat{f} = \tilde{K}\vec{\psi}}
#' @export
permprocest <- function(K,counts=NaN,g=1,c=1,maxiter=300){
  #Gather up all the needed stuff to do the calculation: Parameters, data, etc
  if(length(counts)==1){
    #Allows counts to be the number of spatial locations, the actual counts, or left alone and inferred from K
    if(!is.nan(counts)){
      counts <- array(1,dim=c(counts, 1))
    }else{
      counts <- array(1, dim=c(dim(K)[1], 1))
    }
  }
  if(dim(counts)[2]>1){
    counts <- rowSums(counts)
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
                   , method="CG",control=list(maxit=maxiter,reltol=10^-6))
  alphahat <- optval$par
  #Get intensity
  f <- KTxx%*%alphahat
  intensity <- 0.5*c*f^2
  #Parameters of the posterior distributions
  # diag(Dinv) <- 0.5*(rowSums(sima20_binned[maskA,]))/optval$par^2
  Dinv <- eye(length(counts))
  diag(Dinv) <- 2*(counts/alphahat^2)
  Sigma <- KTxx - KTxx%*%solve(Dinv+KTxx, KTxx)#Posterior covariance of f
  predvar <- diag(Sigma)
  gamma_alpha <- (f^2+predvar)^2/(2*predvar*(2*f^2+predvar))
  gamma_beta <- (f^2+predvar)/(predvar*c*(2*f^2+predvar))
  outputs <- list("lambda"=intensity, "alphaparam" = gamma_alpha, "betaparam"= gamma_beta, "kernelcoeffs" = alphahat)
  return(outputs)
}
#' @title Estimation of the ratio of point process intensities using a Permanental Process model
#' @description Given a realizataion of two point processes, estimates the ratios of the underlying intensities by assuming that each intensity function is a squared Gaussian process.
#' @param K1 matrix: Kernel matrix for the numerator evaluated at the data points.
#' @param K2 matrix: Kernel matrix for the denominator evaluated at the data points. Default is K1.
#' @param counts1 numeric: Count data for the numerator, see documentation for permprocest().
#' @param counts2 numeric: Count data for the denominator, see documentation for permprocest().
#' @param g1 scalar: Regularization strength (or prior precision) for the problem in the numerator. Default is 1
#' @param c1 scalar: Parameter such that intensity \eqn{\vec{\Lambda} = \frac{c}{2}\vec{f}^2} for the problem in the numerator. Default is 1.
#' @param g2 scalar: Regularization strength (or prior precision) for the problem in the denominator. Default is g1
#' @param c2 scalar: Parameter such that intensity \eqn{\vec{\Lambda} = \frac{c}{2}\vec{f}^2} for the problem in the denominator. Default is c1.
#' @param maxiter scalar: The maximum number of iterations.
#' @returns The estimated intensity ratio (MAP) and the parameters of the Generalized Beta-Prime distribution describing the posterior
#' @export
ratioestimationpermproc <- function(K1, counts1, counts2, K2=K1,c1=1,g1=1,c2=c1,g2=g1,maxiter=300){
  print("Performing estimation for the numerator...")
  numeratorresult <- permprocest(K1, counts1, c=c1, g=g1, maxiter=maxiter)
  print("Performing estimation for the denominator...")
  denominatorresult <- permprocest(K2, counts2, c=c2, g=g2, maxiter=maxiter)
  bpalpha <- numeratorresult$alphaparam
  bpbeta <- denominatorresult$alphaparam
  bpp <- array(1, dim=c(length(bpbeta),1))
  bpq <- denominatorresult$betaparam/numeratorresult$betaparam
  pmode <- bpq*((bpalpha*bpp-1)/(bpbeta*bpp+1))^(1/bpp) #MAP estimate
  outputs <- list("ratio"=pmode, "bpalpha"=bpalpha,"bpbeta"=bpbeta,"bpp"=bpp,"bpq"=bpq)
}
#' @title Gradient of the RKHS-regularized Poisson point process likelihood with respect to the kernel coefficients
#' @description Gradient of the RKHS-regularized Poisson point process likelihood with respect to the kernel coefficients
#' @param alpha vector: Kernel coefficients \eqn{\vec{\psi}} so that \eqn{\vec{f}=\tilde{K}\vec{\psi}}
#' @param counts vector: Count data, number observed at each observation point
#' @param K matrix: Kernel matrix for the equivalent kernel function \eqn{\tilde{K}}
#' @returns the gradient of the objective function at alpha
#' @export
gradientFn <- function(alpha, counts, K){
  val <- zeros(length(alpha),1)
  KtAlpha <- (K%*%alpha)
  for(ii in 1:length(counts)){
    val <- val-2*counts[ii]*K[,ii]/KtAlpha
  }
  val <- val+KtAlpha
  return(val)
}
#' @title Objective function of the RKHS-regularized Poisson point process likelihood with respect to the kernel coefficients
#' @description Objective function of the RKHS-regularized Poisson point process likelihood with respect to the kernel coefficients
#' @param alpha vector: Kernel coefficients alpha so that \eqn{\vec{f}=\tilde{K}\vec{\psi}}
#' @param counts vector: Count data, number observed at each observation point
#' @param K matrix: Kernel matrix for the kernel function \eqn{\tilde{k}}
#' @param c scalar: Value so that \eqn{\vec{\Lambda} = \frac{c}{2}\vec{f}^2}
#' @returns the value of the objective function at alpha
#' @export
objectivefn <- function(alpha, counts, K,c){
  val <- -1*sum(counts*log((c/2)*(K%*%alpha)^2))+0.5*t(alpha)%*%K%*%(alpha)
  return(val)
}
