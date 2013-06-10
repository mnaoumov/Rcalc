options(warn=2, error=quote(dump.frames("BlockExp2Gamma", TRUE)))
source("PowerFunctions.R")

rexp2 = function(N) {
  rexp(N) ^ 2
}

rgamma05 = function(N) {
  rgamma(N, 0.5)
}

rgamma5 = function(N) {
  rgamma(N, 5)
}

PowerForMu (mu=c(0,0,0,0),distr=rexp2,Z=1000,V=10000)
PowerForMu (mu=c(-1/5,0,0,1/5),distr=rexp2,Z=1000,V=10000)
PowerForMu (mu=c(-2/5,0,0,2/5),distr=rexp2,Z=1000,V=10000)
PowerForMu (mu=c(-3/5,0,0,3/5),distr=rexp2,Z=1000,V=10000)
PowerForMu (mu=c(-4/5,0,0,4/5),distr=rexp2,Z=1000,V=10000)
PowerForMu (mu=c(-5/5,0,0,5/5),distr=rexp2,Z=1000,V=10000)
PowerForMu (mu=c(-6/5,0,0,6/5),distr=rexp2,Z=1000,V=10000)
PowerForMu (mu=c(-7/5,0,0,7/5),distr=rexp2,Z=1000,V=10000)

PowerForMu (mu=c(0,0,0,0),distr=rgamma05,Z=1000,V=10000)
PowerForMu (mu=c(-1/5,0,0,1/5),distr=rgamma05,Z=1000,V=10000)
PowerForMu (mu=c(-2/5,0,0,2/5),distr=rgamma05,Z=1000,V=10000)
PowerForMu (mu=c(-3/5,0,0,3/5),distr=rgamma05,Z=1000,V=10000)
PowerForMu (mu=c(-4/5,0,0,4/5),distr=rgamma05,Z=1000,V=10000)
PowerForMu (mu=c(-5/5,0,0,5/5),distr=rgamma05,Z=1000,V=10000)
PowerForMu (mu=c(-6/5,0,0,6/5),distr=rgamma05,Z=1000,V=10000)
PowerForMu (mu=c(-7/5,0,0,7/5),distr=rgamma05,Z=1000,V=10000)

PowerForMu (mu=c(0,0,0,0),distr=rgamma5,Z=1000,V=10000)
PowerForMu (mu=c(-1/5,0,0,1/5),distr=rgamma5,Z=1000,V=10000)
PowerForMu (mu=c(-2/5,0,0,2/5),distr=rgamma5,Z=1000,V=10000)
PowerForMu (mu=c(-3/5,0,0,3/5),distr=rgamma5,Z=1000,V=10000)
PowerForMu (mu=c(-4/5,0,0,4/5),distr=rgamma5,Z=1000,V=10000)
PowerForMu (mu=c(-5/5,0,0,5/5),distr=rgamma5,Z=1000,V=10000)
PowerForMu (mu=c(-6/5,0,0,6/5),distr=rgamma5,Z=1000,V=10000)
PowerForMu (mu=c(-7/5,0,0,7/5),distr=rgamma5,Z=1000,V=10000)