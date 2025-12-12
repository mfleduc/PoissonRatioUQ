#' Functions for calculating distributions given data
#' @title Calculate the distribution of Z under a Beta-Prime assumption
#' @description given count datasets a and b, calculates the distribution of Z = lambda_a / lambda_b where lambda_i is the Poisson intensity of the channel that produced dataset i. Under a Gamma prior for each intensity, this is a Beta-Prime distribution. The Gamma pdf uses the shape/rate parameterization and the default prior is uninformative. A prior of the form x^-k_i can be obtained by setting ai=1-ki and bi=0
#' @param a Matrix. The count data for the numerator. Rows correspond to spatial locations, columns to realizations. Must have the same number of rows as b, but need not have the same number of columns. MUST be a matrix or array right now, if you want a scalar use a 1x1 array. Missing data should be replaced by NaNs.
#' @param b Matrix. The count data for the denominator. Rows correspond to spatial locations, columns to realizations. Must have the same number of rows as b, but need not have the same number of columns. MUST be a matrix or array right now, if you want a scalar use a 1x1 array. Missing data should be replaced by NaNs.
#' @param a1 Vector or scalar. The shape parameter of the Gamma prior for the upper channel. If a vector, must have the same number of entries as rows in a and b. Default is 1.
#' @param b1 Vector or scalar. The rate parameter of the Gamma prior for the upper channel. If a vector, must have the same number of entries as rows in a and b. Default is 0.
#' @param a2 Vector or scalar. The shape parameter of the Gamma prior for the lower channel. If a vector, must have the same number of entries as rows in a and b. Default is 1.
#' @param b2 Vector or scalar. The rate parameter of the Gamma prior for the lower channel. If a vector, must have the same number of entries as rows in a and b. Default is 0.
#' @returns a list containing the parameters of the beta-prime distributions at each spatial location, using the parameterization in https://en.wikipedia.org/wiki/Beta_prime_distribution#Generalization
#' @export
zbetaprime <-function(a,b,a1=1,a2=1,b1=0,b2=0){
  #Returns the parameters of the distribution given by propagating uncertainty
  #under the assumption that Z has a Beta-Prime distribution, equivalent to performing
  # a Bayesian update for the mean parameters of the Poisson distributions under the priors
  #\lambda_1 ~ Gamma(a1,b1), \lambda_2 ~ Gamma(a2,b2). This contains as a special case \lambda_i\sim x^{-k_i} by setting ai=1-ki,bi=0.
  #a,b are vectors of data, with a the counts observed in the upper channel (numerator) and b
  # the counts observed in the lower channel (denominator)
  da <- dim(a)[2]
  db <- dim(b)[2]
  if(da>1){
    suma <- rowSums(a,na.rm=TRUE)
    n1 <- rowSums(!is.nan(a))
  }else{
    n1 <-1
    suma <- a
  }
  if(db>1){
    sumb <- rowSums(b,na.rm=TRUE)
    n2 <- rowSums(!is.nan(b))
  }else{
    n2 <- 1
    sumb <- b
  }
  alpha <- suma+a1#Calculating the parameters
  beta <- sumb+a2#
  p <- beta/beta#Will always be ones most likely
  q <- (b2+n2)/(b1+n1)
  param_list <- list('alpha'=alpha,'beta'=beta,'p'=p,'q'=q)
  return(param_list)
}
#' @title Determine the parameters of the distribution of the temperature given the data under the assumption that Z is either a point estimate or has a Gaussian distribution, and is an affine function of T
#' @description Given count data, a prior, and the parameters of the distribution Z|T ~N(mT+z0,tausq) determines the distribution of T|data under the assumption that the counts are either Gaussian or provide a point estimate of Z.
#' @param a Matrix. The count data for the numerator. Rows correspond to spatial locations/channels, columns to realizations. Must have the same number of rows as b, but need not have the same number of columns. MUST be a matrix or array right now, if you want a scalar use a 1x1 array. Missing data should be replaced by NaNs.
#' @param b Matrix. The count data for the denominator. Rows correspond to spatial locations/channels, columns to realizations. Must have the same number of rows as b, but need not have the same number of columns. MUST be a matrix or array right now, if you want a scalar use a 1x1 array. Missing data should be replaced by NaNs.
#' @param m scalar. Slope parameter of the model Z(T)=(mT+z_0)^p.
#' @param z0 scalar. Intercept parameter of the model Z(T)=(mT+z_0)^p.
#' @param p scalar. Exponent parameter of the model Z(T)=(mT+z_0)^p.
#' @returns Parameters of the shifted generalized Beta-Prime distribution T|a,b using the permanental process model, where the band ratio Z(T)=(mT+z_0)^p.
#' @export
tgivenab <- function(a,b,m,z0,p,spatial=TRUE){
  # Calculate the parameters of the distribution T|a,b under the assumption that the
  # ratio Z=lambda_a/lambda_b has a Beta-prime distribution.
  # a and b are the data for the upper and lower channels respectively.
  # This code uses the permanental process model under the assumption that the ratio of channel intensity parameters Z is related to the quantity of interest $T$ by $Z(T)=(MT+z_0)^p
  if(spatial){

  }else{

  }

  # zparams <- zgaussian(a,b)
  # if(length(M)==1){
  #   nSpatial <- dim(a)[1] #Expanding into a spatial model
  #   M <- diag(nSpatial)*M
  # }else if(min(dim(M))==1){
  #   M <- diag(M)
  #   nSpatial <- dim(M)[1]
  # }else{
  #   nSpatial <- dim(M)[2]
  # }
  # if(length(Z0)==1){
  #   Z0 <- matrix(Z0,nrow=dim(a)[1],ncol=1)
  # }
  # if(length(TauSq)==1){
  #   TauSq <- eye(dim(a)[1])*TauSq
  # }
  # TauSqInv <- solve(TauSq)
  # if(length(priormn)==1){
  #   priormn <- matrix(priormn,nrow=nSpatial,ncol=1)
  # }
  # if(length(priorvar)==1&includeprior){
  #   priorvar <- diag(nSpatial)*priorvar
  # }
  # MtTsInv <- t(M)%*%TauSqInv
  # if(tolower(uncertainty)=="gaussian"){
  #   # muT <- 1/m*(zparams$mean-z0)
  #   # sT <- sqrt(tausq+zparams$stdev^2)/m
  #   # stdev <- 1/sqrt(1/sigma0^2+1/sT^2)
  #   # mn <- (stdev)^2*(mu0/sigma0^2+muT/sT^2)
  #   # sT <- 1/sqrt(1/priorvar+m^2/tausq)
  #   # variance <- sT^2*(sT^2*zparams$stdev^2*m^2+tausq^2)/tausq^2
  #   # stdev <- sqrt(variance)
  #   # mn <- sT^2*(priormn/priorvar+m/tausq*(zparams$mean-z0))
  #   if(includeprior){
  #     priorvarinv <- solve(priorvar)
  #     SigmaT <- solve(priorvarinv+MtTsInv%*%M)
  #     mn <- SigmaT%*%(priorvarinv%*%priormn+(MtTsInv%*%(zparams$mean-Z0)))
  #     SigmaZsqrt <- eye(nSpatial)
  #     diag(SigmaZsqrt) <- ((zparams$stdev))
  #     H <- SigmaT%*%MtTsInv%*%SigmaZsqrt
  #     cov <- SigmaT + H%*%t(H)
  #   }else{
  #     SigmaT <- solve(MtTsInv%*%M)
  #     mn <- SigmaT%*%(MtTsInv%*%(zparams$mean-Z0))
  #     SigmaZsqrt <- eye(nSpatial)
  #     diag(SigmaZsqrt) <- ((zparams$stdev))
  #     H <- SigmaT%*%MtTsInv%*%SigmaZsqrt
  #     cov <- SigmaT + H%*%t(H)
  #   }
  # }else if(tolower(uncertainty)=="none"){
  #   # stdev <- 1/sqrt(1/priorvar+m^2/tausq);
  #   if(includeprior){
  #     priorvarinv <- solve(priorvar)
  #     cov <- solve(priorvarinv+MtTsInv%*%M)
  #     mn <- cov%*%(priorvarinv%*%priormn+MtTsInv%*%(zparams$mean-Z0))
  #   }else{
  #     cov <- solve(MtTsInv%*%M)
  #     mn <- cov%*%(MtTsInv%*%(zparams$mean-Z0))
  #   }
  # }
  # param_list<-list('mean'=mn,'cov'=cov)
  # return(param_list)

}
