options(warn=2, error=quote(dump.frames("new", TRUE)))
source("PowerFunctions.R")

rexp2 = function(N) {
  rexp(N) ^ 2
}

rgamma5 = function(N) {
  rgamma(N, 5)
}

rgamma05 = function(N) {
  rgamma(N, 0.5)
}

b = 10
 
k = 4
#accuracy norm

FMCall=NULL

QMCall=NULL
 
for (j in 1:100){

A = generateRandomMatrix(rnorm, b, k)
 
iterations = 10000 
Fstat = NULL
 
Lstat = NULL

 
for (i in 1 : iterations) {
   
B = t(apply(A, 1, sample, size = k))#apply(A, 1, sample, size = k)
   
xx = apply(B, 2, mean)#apply(B, 1, mean)
   
SST = b * sum(xx ^ 2) 
   
SSE = sum((B-matrix(xx,5,3,byrow=TRUE))^2)#sum((B - xx) ^ 2)
   
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))#c(Fstat, SST/ (k - 1) / (SSE / ((k - 1) * (b - 1)) ))
   
Lstat = c(Lstat, Lambda(A, xx[-k])$max)
 
}
 
 
 

u = (7 : 12) / 10
 
M = 30
 
 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
 
FMC = apply(outer(Fstat, b * u ^ 2, ">"), 2, mean)
 
QMC = apply(outer(Fstat*(k-1), b * u ^ 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
 
F = 1 - pf(b * u ^ 2 , k-1, (k-1) * (b-1))
 
FMC
 
round(F,4)
 
QMC

round(Q,4)


FMCall=rbind(FMCall,FMC)
 
QMCall=rbind(QMCall,QMC)

  } 

tplr = NULL
 
tpnc = NULL
 
 
for (up in u) {
   
tt = TailDistrMonteCarlo(A, up, M)
   
tplr = c(tplr, tt$LR)
   
tpnc = c(tpnc, tt$NC)
 
}


 

FMC
round(F,4)
QMC
round(Q,4)
L
round(tplr,4)
round(tpnc,4)

#accuarcy unif
FMCall=NULL

QMCall=NULL
 
for (j in 1:100){

A = generateRandomMatrix(runif, b, k)
 
iterations = 10000 
Fstat = NULL
 
Lstat = NULL

 
for (i in 1 : iterations) {
   
B = t(apply(A, 1, sample, size = k))#apply(A, 1, sample, size = k)
   
xx = apply(B, 2, mean)#apply(B, 1, mean)
   
SST = b * sum(xx ^ 2) 
   
SSE = sum((B-matrix(xx,5,3,byrow=TRUE))^2)#sum((B - xx) ^ 2)
   
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))#c(Fstat, SST/ (k - 1) / (SSE / ((k - 1) * (b - 1)) ))
   
Lstat = c(Lstat, Lambda(A, xx[-k])$max)
 
}
 
 
 

u = (7 : 12) / 10
 
M = 30
 
 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
 
FMC = apply(outer(Fstat, b * u ^ 2, ">"), 2, mean)
 
QMC = apply(outer(Fstat*(k-1), b * u ^ 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
 
F = 1 - pf(b * u ^ 2 , k-1, (k-1) * (b-1))
 
FMC
 
round(F,4)
 
QMC

round(Q,4)


FMCall=rbind(FMCall,FMC)
 
QMCall=rbind(QMCall,QMC)

  } 

tplr = NULL
 
tpnc = NULL
 
 
for (up in u) {
   
tt = TailDistrMonteCarlo(A, up, M)
   
tplr = c(tplr, tt$LR)
   
tpnc = c(tpnc, tt$NC)
 
}


 

FMC
round(F,4)
QMC
round(Q,4)
L
round(tplr,4)
round(tpnc,4)

#accuarcy exp
FMCall=NULL

QMCall=NULL
 
for (j in 1:100){

A = generateRandomMatrix(rexp, b, k)
 
iterations = 10000 
Fstat = NULL
 
Lstat = NULL

 
for (i in 1 : iterations) {
   
B = t(apply(A, 1, sample, size = k))#apply(A, 1, sample, size = k)
   
xx = apply(B, 2, mean)#apply(B, 1, mean)
   
SST = b * sum(xx ^ 2) 
   
SSE = sum((B-matrix(xx,5,3,byrow=TRUE))^2)#sum((B - xx) ^ 2)
   
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))#c(Fstat, SST/ (k - 1) / (SSE / ((k - 1) * (b - 1)) ))
   
Lstat = c(Lstat, Lambda(A, xx[-k])$max)
 
}
 
 
 

u = (7 : 12) / 10
 
M = 30
 
 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
 
FMC = apply(outer(Fstat, b * u ^ 2, ">"), 2, mean)
 
QMC = apply(outer(Fstat*(k-1), b * u ^ 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
 
F = 1 - pf(b * u ^ 2 , k-1, (k-1) * (b-1))
 
FMC
 
round(F,4)
 
QMC

round(Q,4)


FMCall=rbind(FMCall,FMC)
 
QMCall=rbind(QMCall,QMC)

  } 

tplr = NULL
 
tpnc = NULL
 
 
for (up in u) {
   
tt = TailDistrMonteCarlo(A, up, M)
   
tplr = c(tplr, tt$LR)
   
tpnc = c(tpnc, tt$NC)
 
}


 

FMC
round(F,4)
QMC
round(Q,4)
L
round(tplr,4)
round(tpnc,4)


#accuarcy exp2
FMCall=NULL

QMCall=NULL
 
for (j in 1:100){

A = generateRandomMatrix(rexp2, b, k)
 
iterations = 10000 
Fstat = NULL
 
Lstat = NULL

 
for (i in 1 : iterations) {
   
B = t(apply(A, 1, sample, size = k))#apply(A, 1, sample, size = k)
   
xx = apply(B, 2, mean)#apply(B, 1, mean)
   
SST = b * sum(xx ^ 2) 
   
SSE = sum((B-matrix(xx,5,3,byrow=TRUE))^2)#sum((B - xx) ^ 2)
   
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))#c(Fstat, SST/ (k - 1) / (SSE / ((k - 1) * (b - 1)) ))
   
Lstat = c(Lstat, Lambda(A, xx[-k])$max)
 
}
 
 
 

u = (7 : 12) / 10
 
M = 30
 
 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
 
FMC = apply(outer(Fstat, b * u ^ 2, ">"), 2, mean)
 
QMC = apply(outer(Fstat*(k-1), b * u ^ 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
 
F = 1 - pf(b * u ^ 2 , k-1, (k-1) * (b-1))
 
FMC
 
round(F,4)
 
QMC

round(Q,4)


FMCall=rbind(FMCall,FMC)
 
QMCall=rbind(QMCall,QMC)

  } 

tplr = NULL
 
tpnc = NULL
 
 
for (up in u) {
   
tt = TailDistrMonteCarlo(A, up, M)
   
tplr = c(tplr, tt$LR)
   
tpnc = c(tpnc, tt$NC)
 
}


 

FMC
round(F,4)
QMC
round(Q,4)
L
round(tplr,4)
round(tpnc,4)

#accuarcy gamma5
FMCall=NULL

QMCall=NULL
 
for (j in 1:100){

A = generateRandomMatrix(rgamma5, b, k)
 
iterations = 10000 
Fstat = NULL
 
Lstat = NULL

 
for (i in 1 : iterations) {
   
B = t(apply(A, 1, sample, size = k))#apply(A, 1, sample, size = k)
   
xx = apply(B, 2, mean)#apply(B, 1, mean)
   
SST = b * sum(xx ^ 2) 
   
SSE = sum((B-matrix(xx,5,3,byrow=TRUE))^2)#sum((B - xx) ^ 2)
   
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))#c(Fstat, SST/ (k - 1) / (SSE / ((k - 1) * (b - 1)) ))
   
Lstat = c(Lstat, Lambda(A, xx[-k])$max)
 
}
 
 
 

u = (7 : 12) / 10
 
M = 30
 
 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
 
FMC = apply(outer(Fstat, b * u ^ 2, ">"), 2, mean)
 
QMC = apply(outer(Fstat*(k-1), b * u ^ 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
 
F = 1 - pf(b * u ^ 2 , k-1, (k-1) * (b-1))
 
FMC
 
round(F,4)
 
QMC

round(Q,4)


FMCall=rbind(FMCall,FMC)
 
QMCall=rbind(QMCall,QMC)

  } 

tplr = NULL
 
tpnc = NULL
 
 
for (up in u) {
   
tt = TailDistrMonteCarlo(A, up, M)
   
tplr = c(tplr, tt$LR)
   
tpnc = c(tpnc, tt$NC)
 
}


 

FMC
round(F,4)
QMC
round(Q,4)
L
round(tplr,4)
round(tpnc,4)


#accuarcy gamma05
FMCall=NULL

QMCall=NULL
 
for (j in 1:100){

A = generateRandomMatrix(rgamma05, b, k)
 
iterations = 10000 
Fstat = NULL
 
Lstat = NULL

 
for (i in 1 : iterations) {
   
B = t(apply(A, 1, sample, size = k))#apply(A, 1, sample, size = k)
   
xx = apply(B, 2, mean)#apply(B, 1, mean)
   
SST = b * sum(xx ^ 2) 
   
SSE = sum((B-matrix(xx,5,3,byrow=TRUE))^2)#sum((B - xx) ^ 2)
   
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))#c(Fstat, SST/ (k - 1) / (SSE / ((k - 1) * (b - 1)) ))
   
Lstat = c(Lstat, Lambda(A, xx[-k])$max)
 
}
 
 
 

u = (7 : 12) / 10
 
M = 30
 
 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
 
FMC = apply(outer(Fstat, b * u ^ 2, ">"), 2, mean)
 
QMC = apply(outer(Fstat*(k-1), b * u ^ 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
 
F = 1 - pf(b * u ^ 2 , k-1, (k-1) * (b-1))
 
FMC
 
round(F,4)
 
QMC

round(Q,4)


FMCall=rbind(FMCall,FMC)
 
QMCall=rbind(QMCall,QMC)

  } 

tplr = NULL
 
tpnc = NULL
 
 
for (up in u) {
   
tt = TailDistrMonteCarlo(A, up, M)
   
tplr = c(tplr, tt$LR)
   
tpnc = c(tpnc, tt$NC)
 
}


 

FMC
round(F,4)
QMC
round(Q,4)
L
round(tplr,4)
round(tpnc,4)