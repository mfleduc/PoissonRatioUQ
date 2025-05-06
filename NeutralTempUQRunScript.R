
source("/home/male7736/Desktop/Research/neutral temp estimation/PoissonRatioUQ/R/GetDistributions.R")
source("/home/male7736/Desktop/Research/neutral temp estimation/PoissonRatioUQ/R/CRPSGaussian.R")
source("/home/male7736/Desktop/Research/neutral temp estimation/PoissonRatioUQ/R/CRPS.R")
source("/home/male7736/Desktop/Research/neutral temp estimation/PoissonRatioUQ/R/HPDIntervalGaussian.R")
source("/home/male7736/Desktop/Research/neutral temp estimation/PoissonRatioUQ/R/HPDSet.R")
source("/home/male7736/Desktop/Research/neutral temp estimation/PoissonRatioUQ/R/BetaPrimeFunctions.R")
# Parameters for the random data
la <- 40
lb <- 70
zhat <- la/lb
 a <- rpois(10, la)
 b <- rpois(10, lb)
# Calculate distributions for Z given the data
gaussianparams <- ZGaussian(a,b)
bpParams <- ZBetaPrime(a,b,k1=0,k2=1)

crpsGZ <- CRPSGaussian(gaussianparams$mean, gaussianparams$stdev, zhat)
hpdinterval <- HPDIntervalGaussian(gaussianparams$mean, gaussianparams$stdev, 0.95)

z <- seq(-3,3,0.01)
pdf <- dnorm(z, gaussianparams$mean, gaussianparams$stdev)
hpdSet <- HPDSet(z, pdf, 0.95)#Returns all values of Z within the HPD credible set

pdfBP <- dbetaprime(z, bpParams$alpha,bpParams$beta, bpParams$p,bpParams$q)
bprnd <- rbetaprime(40000, bpParams$alpha,bpParams$beta,bpParams$p,bpParams$q)
bpcdf <-  pbetaprime(z,bpParams$alpha,bpParams$beta, bpParams$p, bpParams$q)
#plot(z, bpcdf)
## Temperature estimation stuff: T|a,b
Tvals <- seq(0,1500,by=1)
m<- 0.0009
z0 <- 0.1088
tausq <-0.00002
# Luckily have a closed form for the temperature when Z is Gaussian
TgivenabZGaussian <- dnorm(Tvals, mean=1/m*(gaussianparams$mean-z0), sd = (m/sqrt(gaussianparams$stdev^2+m^2*tausq))^(-1))
#cdfTgivenabZGaussian <- pnorm(Tvals, mean=1/m*(gaussianparams$mean-z0), sd = (m/sqrt(gaussianparams$stdev^2+m^2*tausq))^(-1))
hpdT <- HPDIntervalGaussian(1/m*(gaussianparams$mean-z0), (m/sqrt(gaussianparams$stdev^2+m^2*tausq))^(-1), 0.95)
# Now we have to do the same for Z~BP
TgabZBP <- array(numeric(0),dim=c(length(bprnd),length(Tvals)))
for(tt in 1:length(bprnd)){
  TgabZBP[tt,]<- dnorm(Tvals, mean=1/m*( bprnd[tt]-z0), sd = sqrt(tausq/m^2))
}
TgabZBP <- colMeans(TgabZBP)
hpdTZBP <- HPDSet(Tvals, TgabZBP, 0.95)
hpdTZBP <- c(min(hpdTZBP),max(hpdTZBP))
TgabPE  <- dnorm(Tvals, mean=1/m*( mean(a)/mean(b)-z0), sd = sqrt(tausq/m^2))
#cdfTgabPE <- pnorm(Tvals, mean=1/m*( mean(a)/mean(b)-z0), sd = sqrt(tausq/m^2))
hpdPE <- HPDIntervalGaussian(1/m*( mean(a)/mean(b)-z0),sqrt(tausq/m^2),0.95)
#plot(Tvals, TgabPointEst)

crpsG <- CRPSGaussian(1/m*(gaussianparams$mean-z0), sqrt(gaussianparams$stdev^2+m^2*tausq)/m,(zhat-z0)/m)
crpsC <- CRPSGaussian(1/m*( mean(a)/mean(b)-z0),sqrt(tausq/m^2),(zhat-z0)/m)
crpsBP <- CRPS(Tvals, pdf=TgabZBP,xhat=(zhat-z0)/m )

plot(Tvals,TgivenabZGaussian,type='n',xlab='Temperature, K',ylab="PDF")
lines(Tvals, TgivenabZGaussian,type='l', col='red')
lines(Tvals, TgabZBP,type='l', col='blue')
lines(Tvals, TgabPE,type='l', col='green')
fields::xline((zhat-z0)/m)
