#' Utitility functions for the package PoissonRatioUQ
#'@title Solar zenith angle
#'@description Calculates the solar zenith angle at a particular point using the latitude, longitude, time in UTC, and julian day
#'@param lat Scalar or matrix. Latitude of the points at which to calculate SZA. 0=equator, Degrees N.
#'@param lon Scalar or matrix. Longitude of the points at which to calculate SZA. Degrees E, ranging from -180 to 180
#'@param hourUTC Scalar. Time in UTC, hours since midnight
#'@param julianDay Scalar. Julian day, i.e. 1 = Jan1, 2 = Jan 2, 32 = Feb 1, etc.
#'@returns Solar zenith angle, in degrees
#'@export
getsza <- function(lat, lon, hourUTC, julianDay){
# Computes the solar zenith angle
#
  hourLocal <- hourUTC+lon/15
  lat <- lat*pi/180#Convert to radians
  delta <- -23.44*(pi/180)*cos(2*pi*(julianDay+9)/365) #Approx. declination of the sun
  h <- pi/12*(hourLocal-12) #Hour angle
  cossza <- sin(lat)*sin(delta)+cos(lat)*cos(delta)*cos(h)
  return(acosd(cossza))
}
#'@title Row variances
#'@description Calculates the variance of each row of a matrix
#'@param x matrix: The matrix
#'@param na.rm BOOL: Whether of not you want NaNs removed. Default is TRUE
#'@returns sample variance of the data in each row of the matrix
#'@export
rowVars <- function(x, na.rm=TRUE) {
  # Vectorised version of variance filter
  rowSums((x - rowMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (ncol(x) - 1)
}
#'@title Bin arrays into blocks
#'@description Takes a 2D array and divides it into blocks of constant size. Probably not the fastest possible implementation, but for now this is what we are doing.
#'@param x matrix/array: The matrix to divide into blocks
#'@param blocksize vector: dimensions of the blocks
#'@param reshapeTo2D Bool: Flag to set to TRUE if instead of wanting an output array of dimension nxmxk, you want the output to be dimension (nm)xk, i.e. you want to flatten the blocks. Default FALSE
#'@returns array blocked into chunks of the desired size
#'@export
binarray <- function(x, blocksize=c(2,2),reshapeTo2D=FALSE){
  # bin a 2x2 array into predefined blocks
  nblocks <- prod(dim(x))/prod(blocksize)
  stopifnot("Must be an integer number of blocks"=nblocks==ceiling(nblocks) )
  binnedx <- array(NaN, dim=append(blocksize, nblocks))
  counter <-0
  for(ii in seq(1, dim(x)[1]-1,by=blocksize[1])){
    for(jj in seq(1, dim(x)[2]-1,by=blocksize[2])){
      counter<-counter+1
      binnedx[,,counter] <- x[ii:(ii+blocksize[1]-1),jj:(jj+blocksize[2]-1)]
    }
  }
  if(reshapeTo2D){
    dim(binnedx)<-c(prod(blocksize),nblocks)
  }
  return(binnedx)
}
#'@title Gaussian smoother
#'@description Applies a Gaussian smoother with width sigma to the input data
#'@param sigma scalar: Width of the Gaussian smoother
#'@returns Smoothed data
#'@export
gausssmooth <- function( data, sigma ){
  n <- length(data)
  winsize <- ceiling(5*sigma)
  winhalf <-  ceiling(winsize/2)
  win <- exp( -0.5*(  1:winsize - ceiling(winsize/2) )^2/sigma^2 )
  win <- win/sum(win)
  smootheddata <- filter(  data, win, method="convolution"  )



}
