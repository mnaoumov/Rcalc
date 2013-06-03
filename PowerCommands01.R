options(warn=2, error=quote(dump.frames("PowerFunctions01", TRUE)))
source("PowerFunctions.R")
PowerForMu (mu=c(0,0,0,0),distr=rexp,Z=2000,V=10000)
PowerForMu (mu=c(-1/5,0,0,1/5),distr=rexp,Z=2000,V=10000)