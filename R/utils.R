getsza <- function(lat, hourLocal, julianDay){
# Computes the solar zenith angle
#
  lat <- lat*pi/180#Convert to radians
  delta <- -23.44*(pi/180)*cos(2*pi*(julianDay+9)/365) #Approx. declination of the sun
  h <- pi/12*(hourLocal-12) #Hour angle
  cossza <- sin(lat)*sin(delta)+cos(lat)*cos(delta)*cos(h)
  return(acos(cossza))
}
