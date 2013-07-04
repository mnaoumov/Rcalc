options(warn=2, error=quote(dump.frames("accuracy-fixed", TRUE)))
source("PowerFunctions-fixed.R")

rexp2 = function(N) {
  rexp(N) ^ 2
}

rgamma5 = function(N) {
  rgamma(N, 5)
}

rgamma05 = function(N) {
  rgamma(N, 0.5)
}
 

iterations = 100000 
u = (0 : 14) / 10
M = 30

b = 5
k = 3

#accuracy exp2
A = generateRandomMatrix(rexp2, b, k)

FMCall=NULL
QMCall=NULL

for (j in 1:100){
A = generateRandomMatrix(rexp2, b, k)

Fstat = NULL

 for (i in 1 : 10000) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
SST = b * sum(xx ^ 2) 
SSE = sum((B-matrix(xx,b,k,byrow=TRUE))^2)
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))
   }

FMC = apply(outer(Fstat, b * u ^ 2/(k-1), ">"), 2, mean)
QMC = apply(outer(Fstat*(k-1), b * u ^ 2, ">"), 2, mean)

FMCall = rbind(FMCall,FMC)
QMCall = rbind(QMCall,QMC)

}

Lstat = NULL
 
  for (i in 1 : iterations) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
  Lstat = c(Lstat, Lambda(A, xx[-k])$max)
   }

FMC = apply(FMCall,2,mean)
QMC = apply(QMCall,2,mean) 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
F = 1 - pf(b * u ^ 2/(k-1) , k-1, (k-1) * (b-1))
 
tplr = NULL
tpnc = NULL
  
for (up in u) {
  tt = TailDistrMonteCarlo(A, up, M)
 tplr = c(tplr, tt$LR)
 tpnc = c(tpnc, tt$NC)
 }

rbind(FMC=round(FMC,4), F=round(F,4), QMC=round(QMC,4), Q=round(Q,4), L = round(L,4), LR=round(tplr,4), NC=round(tpnc,4))



#accuracy exp
A = generateRandomMatrix(rexp, b, k)

FMCall=NULL
QMCall=NULL

for (j in 1:100){
A = generateRandomMatrix(rexp, b, k)

Fstat = NULL

 for (i in 1 : 10000) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
SST = b * sum(xx ^ 2) 
SSE = sum((B-matrix(xx,b,k,byrow=TRUE))^2)
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))
   }

FMC = apply(outer(Fstat, b * u ^ 2/(k-1), ">"), 2, mean)
QMC = apply(outer(Fstat*(k-1), b * u ^ 2, ">"), 2, mean)

FMCall = rbind(FMCall,FMC)
QMCall = rbind(QMCall,QMC)

}

Lstat = NULL
 
  for (i in 1 : iterations) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
  Lstat = c(Lstat, Lambda(A, xx[-k])$max)
   }

FMC = apply(FMCall,2,mean)
QMC = apply(QMCall,2,mean) 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
F = 1 - pf(b * u ^ 2/(k-1) , k-1, (k-1) * (b-1))
 
tplr = NULL
tpnc = NULL
  
for (up in u) {
  tt = TailDistrMonteCarlo(A, up, M)
 tplr = c(tplr, tt$LR)
 tpnc = c(tpnc, tt$NC)
 }

rbind(FMC=round(FMC,4), F=round(F,4), QMC=round(QMC,4), Q=round(Q,4), L = round(L,4), LR=round(tplr,4), NC=round(tpnc,4))


b = 10
k = 4

#accuracy normal
A = generateRandomMatrix(rnorm, b, k)

FMCall=NULL
QMCall=NULL

for (j in 1:100){
A = generateRandomMatrix(rnorm, b, k)

Fstat = NULL

 for (i in 1 : 10000) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
SST = b * sum(xx ^ 2) 
SSE = sum((B-matrix(xx,b,k,byrow=TRUE))^2)
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))
   }

FMC = apply(outer(Fstat, b * u ^ 2/(k-1), ">"), 2, mean)
QMC = apply(outer(Fstat*(k-1), b * u ^ 2, ">"), 2, mean)

FMCall = rbind(FMCall,FMC)
QMCall = rbind(QMCall,QMC)

}

Lstat = NULL
 
  for (i in 1 : iterations) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
  Lstat = c(Lstat, Lambda(A, xx[-k])$max)
   }

FMC = apply(FMCall,2,mean)
QMC = apply(QMCall,2,mean) 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
F = 1 - pf(b * u ^ 2/(k-1) , k-1, (k-1) * (b-1))
 
tplr = NULL
tpnc = NULL
  
for (up in u) {
  tt = TailDistrMonteCarlo(A, up, M)
 tplr = c(tplr, tt$LR)
 tpnc = c(tpnc, tt$NC)
 }

rbind(FMC=round(FMC,4), F=round(F,4), QMC=round(QMC,4), Q=round(Q,4), L = round(L,4), LR=round(tplr,4), NC=round(tpnc,4))
  

#accuarcy unif
A = generateRandomMatrix(runif, b, k)

FMCall=NULL
QMCall=NULL

for (j in 1:100){
A = generateRandomMatrix(runif, b, k)

Fstat = NULL

 for (i in 1 : 10000) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
SST = b * sum(xx ^ 2) 
SSE = sum((B-matrix(xx,b,k,byrow=TRUE))^2)
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))
   }

FMC = apply(outer(Fstat, b * u ^ 2/(k-1), ">"), 2, mean)
QMC = apply(outer(Fstat*(k-1), b * u ^ 2, ">"), 2, mean)

FMCall = rbind(FMCall,FMC)
QMCall = rbind(QMCall,QMC)

}

Lstat = NULL
 
  for (i in 1 : iterations) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
  Lstat = c(Lstat, Lambda(A, xx[-k])$max)
   }

FMC = apply(FMCall,2,mean)
QMC = apply(QMCall,2,mean) 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
F = 1 - pf(b * u ^ 2/(k-1) , k-1, (k-1) * (b-1))
 
tplr = NULL
tpnc = NULL
  
for (up in u) {
  tt = TailDistrMonteCarlo(A, up, M)
 tplr = c(tplr, tt$LR)
 tpnc = c(tpnc, tt$NC)
 }

rbind(FMC=round(FMC,4), F=round(F,4), QMC=round(QMC,4), Q=round(Q,4), L = round(L,4), LR=round(tplr,4), NC=round(tpnc,4))
 

#accuarcy exp
A = generateRandomMatrix(rexp, b, k)

FMCall=NULL
QMCall=NULL

for (j in 1:100){
A = generateRandomMatrix(rexp, b, k)

Fstat = NULL

 for (i in 1 : 10000) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
SST = b * sum(xx ^ 2) 
SSE = sum((B-matrix(xx,b,k,byrow=TRUE))^2)
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))
   }

FMC = apply(outer(Fstat, b * u ^ 2/(k-1), ">"), 2, mean)
QMC = apply(outer(Fstat*(k-1), b * u ^ 2, ">"), 2, mean)

FMCall = rbind(FMCall,FMC)
QMCall = rbind(QMCall,QMC)

}

Lstat = NULL
 
  for (i in 1 : iterations) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
  Lstat = c(Lstat, Lambda(A, xx[-k])$max)
   }

FMC = apply(FMCall,2,mean)
QMC = apply(QMCall,2,mean) 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
F = 1 - pf(b * u ^ 2/(k-1) , k-1, (k-1) * (b-1))
 
tplr = NULL
tpnc = NULL
  
for (up in u) {
  tt = TailDistrMonteCarlo(A, up, M)
 tplr = c(tplr, tt$LR)
 tpnc = c(tpnc, tt$NC)
 }

rbind(FMC=round(FMC,4), F=round(F,4), QMC=round(QMC,4), Q=round(Q,4), L = round(L,4), LR=round(tplr,4), NC=round(tpnc,4))
 
#accuracy exp2 
A = generateRandomMatrix(rexp2, b, k)

FMCall=NULL
QMCall=NULL

for (j in 1:100){
A = generateRandomMatrix(rexp2, b, k)

Fstat = NULL

 for (i in 1 : 10000) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
SST = b * sum(xx ^ 2) 
SSE = sum((B-matrix(xx,b,k,byrow=TRUE))^2)
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))
   }

FMC = apply(outer(Fstat, b * u ^ 2/(k-1), ">"), 2, mean)
QMC = apply(outer(Fstat*(k-1), b * u ^ 2, ">"), 2, mean)

FMCall = rbind(FMCall,FMC)
QMCall = rbind(QMCall,QMC)

}

Lstat = NULL
 
  for (i in 1 : iterations) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
  Lstat = c(Lstat, Lambda(A, xx[-k])$max)
   }

FMC = apply(FMCall,2,mean)
QMC = apply(QMCall,2,mean) 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
F = 1 - pf(b * u ^ 2/(k-1) , k-1, (k-1) * (b-1))
 
tplr = NULL
tpnc = NULL
  
for (up in u) {
  tt = TailDistrMonteCarlo(A, up, M)
 tplr = c(tplr, tt$LR)
 tpnc = c(tpnc, tt$NC)
 }

rbind(FMC=round(FMC,4), F=round(F,4), QMC=round(QMC,4), Q=round(Q,4), L = round(L,4), LR=round(tplr,4), NC=round(tpnc,4))
 

#accuracy gamma5
A = generateRandomMatrix(rgamma5, b, k)

FMCall=NULL
QMCall=NULL

for (j in 1:100){
A = generateRandomMatrix(rgamma5, b, k)

Fstat = NULL

 for (i in 1 : 10000) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
SST = b * sum(xx ^ 2) 
SSE = sum((B-matrix(xx,b,k,byrow=TRUE))^2)
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))
   }

FMC = apply(outer(Fstat, b * u ^ 2/(k-1), ">"), 2, mean)
QMC = apply(outer(Fstat*(k-1), b * u ^ 2, ">"), 2, mean)

FMCall = rbind(FMCall,FMC)
QMCall = rbind(QMCall,QMC)

}

Lstat = NULL
 
  for (i in 1 : iterations) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
  Lstat = c(Lstat, Lambda(A, xx[-k])$max)
   }

FMC = apply(FMCall,2,mean)
QMC = apply(QMCall,2,mean) 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
F = 1 - pf(b * u ^ 2/(k-1) , k-1, (k-1) * (b-1))
 
tplr = NULL
tpnc = NULL
  
for (up in u) {
  tt = TailDistrMonteCarlo(A, up, M)
 tplr = c(tplr, tt$LR)
 tpnc = c(tpnc, tt$NC)
 }

rbind(FMC=round(FMC,4), F=round(F,4), QMC=round(QMC,4), Q=round(Q,4), L = round(L,4), LR=round(tplr,4), NC=round(tpnc,4))
 
#accuracy gamma05
A = generateRandomMatrix(rgamma05, b, k)

FMCall=NULL
QMCall=NULL

for (j in 1:100){
A = generateRandomMatrix(rgamma05, b, k)

Fstat = NULL

 for (i in 1 : 10000) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
SST = b * sum(xx ^ 2) 
SSE = sum((B-matrix(xx,b,k,byrow=TRUE))^2)
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))
   }

FMC = apply(outer(Fstat, b * u ^ 2/(k-1), ">"), 2, mean)
QMC = apply(outer(Fstat*(k-1), b * u ^ 2, ">"), 2, mean)

FMCall = rbind(FMCall,FMC)
QMCall = rbind(QMCall,QMC)

}

Lstat = NULL
 
  for (i in 1 : iterations) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
  Lstat = c(Lstat, Lambda(A, xx[-k])$max)
   }

FMC = apply(FMCall,2,mean)
QMC = apply(QMCall,2,mean) 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
F = 1 - pf(b * u ^ 2/(k-1) , k-1, (k-1) * (b-1))
 
tplr = NULL
tpnc = NULL
  
for (up in u) {
  tt = TailDistrMonteCarlo(A, up, M)
 tplr = c(tplr, tt$LR)
 tpnc = c(tpnc, tt$NC)
 }

rbind(FMC=round(FMC,4), F=round(F,4), QMC=round(QMC,4), Q=round(Q,4), L = round(L,4), LR=round(tplr,4), NC=round(tpnc,4))
 