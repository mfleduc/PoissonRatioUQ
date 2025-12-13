#' Functions for calculating distributions given data
#' @title Determine the parameters of the distribution of the temperature given the data under the assumption that Z is either a point estimate or has a Gaussian distribution, and is an affine function of T
#' @description Given count data, a prior, and the parameters of the distribution Z|T ~N(mT+z0,tausq) determines the distribution of T|data under the assumption that the counts are either Gaussian or provide a point estimate of Z.
#' @param a Matrix. The count data for the numerator. Rows correspond to spatial locations/channels, columns to realizations. Must have the same number of rows as b, but need not have the same number of columns. MUST be a matrix or array right now, if you want a scalar use a 1x1 array. Missing data should be replaced by NaNs.
#' @param b Matrix. The count data for the denominator. Rows correspond to spatial locations/channels, columns to realizations. Must have the same number of rows as b, but need not have the same number of columns. MUST be a matrix or array right now, if you want a scalar use a 1x1 array. Missing data should be replaced by NaNs.
#' @param m scalar or vector. Slope parameter of the model Z(T)=(mT+z_0)^p.
#' @param z0 scalar or vector. Intercept parameter of the model Z(T)=(mT+z_0)^p.
#' @param p scalar or vector. Exponent parameter of the model Z(T)=(mT+z_0)^p.
#' @param spatial boolean. If TRUE (default), run the permanental process model. If FALSE, runs the pointwise estimation scheme.
#' @param a1,a2,b1,b2 scalars: Parameters of the prior distribution for the pointwise model. See documentation for zbetaprime()
#' @param K1,K2 Matrices. Kernel matrices for the permanental process model. See ratioestimationpermproc() for details. 
#' @param g1,g2,c1,c2 scalars, other parameters for the permanental process model. See ratioestimationpermproc() for more details. This funciton assumes they are the same for both models.
#' 
#' @returns Parameters of the shifted generalized Beta-Prime distribution T|a,b using the permanental process model, where the band ratio Z(T)=(mT+z_0)^p.
#' @export
tgivenab <- function(a,b,m,z0,p,spatial=TRUE,a1=1,a2=1,b1=0,b2=0, K1=NA,K2=K1,g1=1,c1=1,g2=g1,c2=c1){
  # Calculate the parameters of the distribution T|a,b under the assumption that the
  # ratio Z=lambda_a/lambda_b has a Beta-prime distribution.
  # a and b are the data for the upper and lower channels respectively.
  # This code uses the permanental process model under the assumption that the ratio of channel intensity parameters Z is related to the quantity of interest $T$ by $Z(T)=(MT+z_0)^p
  if(spatial){
    if(any(is.na(K1))){stop("Must specify K to use the permanental process model")}
    zparams <- ratioestimationpermproc(K1, a, b, K2=K2, c1=c1,c2=c2,g1=g1,g2=g2)
  }else{
    zparams <- zbetaprime(a, b, a1=a1,a2=a2,b1=b1,b2=b2)
  }
  ## Now: Parameters for the distribution of T
  shift <- -z0/m
  # tp <- p
  tq <- zparams$bpq^(1/p)/m
  tmode <- tq*( (pmax(zparams$bpalpha*p-1,0) )/(zparams$bpbeta*p+1 ))^(1/p)-z0/m
  tparams <- list("T"=tmode,"ratio"=zparams$ratio,"shift"=shift, "bpalpha"=zparams$bpalpha,"bpbeta"=zparams$bpbeta, "bpp"=p,"bpq"=tq)
  return(tparams)
  
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
