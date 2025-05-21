
source("../R/GetDistributions.R")
source("../R/CRPSGaussian.R")
source("../R/CRPS.R")
source("../R/HPDIntervalGaussian.R")
source("../R/HPDSet.R")
source("../R/BetaPrimeFunctions.R")
# Parameters for the random data
la <- 40
lb <- 65
zhat <- la/lb
ntrials <- 200
z <- seq(-3,3,0.005)
Tvals <- seq(0,1500,by=0.5)
m<- 0.0009
z0 <- 0.1088
tausq <-0.00002
TgabZ <- array(numeric(0),dim=c(length(z),length(Tvals)))
TgabZBP <- array(numeric(0),dim=c(1,length(Tvals)))
for(tt in 1:length(z)){
  TgabZ[tt,]<- dnorm(Tvals, mean=1/m*( z[tt]-z0), sd = sqrt(tausq/m^2))
}
crpsG <- array(numeric(0), dim=c(1,ntrials))
crpsC <- array(numeric(0), dim=c(1,ntrials))
crpsBP <- array(numeric(0), dim=c(1,ntrials))

for(nn in 1:ntrials){
  print(nn)
 a <- rpois(4, la)
 b <- rpois(4, lb)
# Calculate distributions for Z given the data
gaussianparams <- ZGaussian(a,b)
bpParams <- ZBetaPrime(a,b,k1=0,k2=1)

crpsGZ <- CRPSGaussian(gaussianparams$mean, gaussianparams$stdev, zhat)
hpdinterval <- HPDIntervalGaussian(gaussianparams$mean, gaussianparams$stdev, 0.95)
pdf <- dnorm(z, gaussianparams$mean, gaussianparams$stdev)
hpdSet <- HPDSet(z, pdf, 0.95)#Returns all values of Z within the HPD credible set
pdfBP <- dbetaprime(z, bpParams$alpha,bpParams$beta, bpParams$p,bpParams$q)

#bprnd <- rbetaprime(20000, bpParams$alpha,bpParams$beta,bpParams$p,bpParams$q)
bpcdf <-  pbetaprime(z,bpParams$alpha,bpParams$beta, bpParams$p, bpParams$q)
#plot(z, bpcdf)
## Temperature estimation stuff: T|a,b

# Luckily have a closed form for the temperature when Z is Gaussian
TgivenabZGaussian <- dnorm(Tvals, mean=1/m*(gaussianparams$mean-z0), sd = (m/sqrt(gaussianparams$stdev^2+m^2*tausq))^(-1))
#cdfTgivenabZGaussian <- pnorm(Tvals, mean=1/m*(gaussianparams$mean-z0), sd = (m/sqrt(gaussianparams$stdev^2+m^2*tausq))^(-1))
#hpdT <- HPDIntervalGaussian(1/m*(gaussianparams$mean-z0), (m/sqrt(gaussianparams$stdev^2+m^2*tausq))^(-1), 0.95)
# Now we have to do the same for Z~BP
#TgabZBP <- array(numeric(0),dim=c(length(bprnd),length(Tvals)))

for(t in 1:length(Tvals)){
TgabZBP[t] <- pracma::trapz(z, t(TgabZ[,t])*pdfBP)
}

#hpdTZBP <- HPDSet(Tvals, TgabZBP, 0.95)
#hpdTZBP <- c(min(hpdTZBP),max(hpdTZBP))
TgabPE  <- dnorm(Tvals, mean=1/m*( mean(a)/mean(b)-z0), sd = sqrt(tausq/m^2))
#cdfTgabPE <- pnorm(Tvals, mean=1/m*( mean(a)/mean(b)-z0), sd = sqrt(tausq/m^2))
#hpdPE <- HPDIntervalGaussian(1/m*( mean(a)/mean(b)-z0),sqrt(tausq/m^2),0.95)
#plot(Tvals, TgabPointEst)

crpsG[nn] <- CRPSGaussian(1/m*(gaussianparams$mean-z0), sqrt(gaussianparams$stdev^2+m^2*tausq)/m,(zhat-z0)/m)
crpsC[nn] <- CRPSGaussian(1/m*( mean(a)/mean(b)-z0),sqrt(tausq/m^2),(zhat-z0)/m)
crpsBP[nn] <- CRPS(Tvals, pdf=TgabZBP,xhat=(zhat-z0)/m )
}
plot(Tvals,TgivenabZGaussian,type='n',xlab='Temperature, K',ylab="PDF")
lines(Tvals, TgivenabZGaussian,type='l', col='red')
lines(Tvals, TgabZBP,type='l', col='blue')
lines(Tvals, TgabPE,type='l', col='green')
fields::xline((zhat-z0)/m)
