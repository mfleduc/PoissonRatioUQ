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
