options(warn=2, error=quote(dump.frames("block-power-new", TRUE)))
source("PowerFunctions.R")

b=5
k=3

PowerForMu (mu=c(0,0,0),distr=rexp,Z=1000,V=10000)
PowerForMu (mu=c(-1/5,0,1/5),distr=rexp,Z=1000,V=10000)
PowerForMu (mu=c(-2/5,0,2/5),distr=rexp,Z=1000,V=10000)
PowerForMu (mu=c(-3/5,0,3/5),distr=rexp,Z=1000,V=10000)
PowerForMu (mu=c(-4/5,0,4/5),distr=rexp,Z=1000,V=10000)
PowerForMu (mu=c(-5/5,0,5/5),distr=rexp,Z=1000,V=10000)
PowerForMu (mu=c(-6/5,0,6/5),distr=rexp,Z=1000,V=10000)
PowerForMu (mu=c(-7/5,0,7/5),distr=rexp,Z=1000,V=10000)