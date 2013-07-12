options(warn=2, error=quote(dump.frames("bd-acc-with-all-distr", TRUE)))

Lambda = function(A, x, tolerance = 0.9999) {
  k = getk(A)
  if(isPointInsideDomain(A,x)==FALSE) x=.999*x
  ll = function(beta) { beta %*% x - kappa(A, beta) }
  lll = function(beta) { x - kappad(A, beta)[,1] }
  llll = function(beta) { -kappad2(A, beta) }
  lam = maxNR(ll, grad = lll, hess = llll, start = rep(0, k - 1), tol = 1 - tolerance)
  lam
}

delta = function(A, u, s, deltaIterations = 5) {
  k = getk(A)
  r = u
  V0h = V0h(A)
  
  for (i in 1 : deltaIterations) {
    tr = transformIntoSphere(A, r, s)
    while (isPointInsideDomain(A,tr)==FALSE) {
    r=.9*r
    tr = transformIntoSphere(A, r, s)}
    
    L = Lambda(A, tr)
    r = r - ((L$max - u ^ 2 / 2) / (L$est %*% V0h %*% s))[1,1]
       
  }
  
  betah = L$est
  
  (det(kappad2(A, betah)) ^ (- 1 / 2) * r ^ (k - 2) * det(V0h)) / (u ^ (k - 3) * abs(t(s) %*% V0h %*% betah))
}

TailDistrMonteCarlo = function(A, u, M, deltaIterations = 5) {
  k = getk(A)
  b = getb(A)
  
  delta2 = function(s) delta(A, u, s, deltaIterations)
  
  cb = (b ^ ((k - 1) / 2)) / (2 ^ ((k - 1) / 2 - 1) * gamma((k - 1) / 2))
  Gb = IntOnSphereMonteCarlo(delta2, k - 1, M)
  
  LugRice = 1 - pchisq(b * u ^ 2, k - 1) + cb / b * u ^ (k - 3) * exp(-b * u ^ 2 / 2) * (Gb - 1)
  ustar = u - log(Gb) / (b * u)
  NielCox = 1 - pchisq(b * ustar ^ 2, k - 1)
  list(LR=c(LugRice), NC=c(NielCox))
}

rexp2 = function(N) {
  rexp(N) ^ 2
}

rgamma5 = function(N) {
  rgamma(N, 5)
}

rgamma05 = function(N) {
  rgamma(N, 0.5)
}
 


u = (7 : 14) / 10
M = 30

b = 5
k = 3

#accuracy exp
A = generateRandomMatrix(rexp, b, k)

Fstat = NULL
Lstat = NULL
Qstat = NULL
 
  for (i in 1 : iterations) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
SST = b * sum(xx ^ 2) 
SSE = sum((B-matrix(xx,b,k,byrow=TRUE))^2)
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))
Lstat = c(Lstat, Lambda(A, xx[-k])$max)
Qstat = c(Qstat,SST)
   }
 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
FMC = apply(outer(Fstat, b * u ^ 2, ">"), 2, mean)
QMC = apply(outer(Qstat, b * u ^ 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
F = 1 - pf(b * u ^ 2 , k-1, (k-1) * (b-1))
 
tplr = NULL
tpnc = NULL
  
for (up in u) {
  tt = TailDistrMonteCarlo(A, up, M)
 tplr = c(tplr, tt$LR)
 tpnc = c(tpnc, tt$NC)
 }

rbind(FMC, F=round(F,4), QMC, Q=round(Q,4), L , LR=round(tplr,4), NC=round(tpnc,4))
 

b = 10
k = 4

#accuracy normal
A = generateRandomMatrix(rnorm, b, k)

Fstat = NULL
Lstat = NULL
Qstat = NULL
 
  for (i in 1 : iterations) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
SST = b * sum(xx ^ 2) 
SSE = sum((B-matrix(xx,b,k,byrow=TRUE))^2)
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))
Lstat = c(Lstat, Lambda(A, xx[-k])$max)
Qstat = c(Qstat,SST)
   }
 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
FMC = apply(outer(Fstat, b * u ^ 2, ">"), 2, mean)
QMC = apply(outer(Qstat, b * u ^ 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
F = 1 - pf(b * u ^ 2 , k-1, (k-1) * (b-1))
 
tplr = NULL
tpnc = NULL
  
for (up in u) {
  tt = TailDistrMonteCarlo(A, up, M)
 tplr = c(tplr, tt$LR)
 tpnc = c(tpnc, tt$NC)
 }

rbind(FMC, F=round(F,4), QMC, Q=round(Q,4), L , LR=round(tplr,4), NC=round(tpnc,4))
 

#accuarcy unif
A = generateRandomMatrix(runif, b, k)

Fstat = NULL
Lstat = NULL
Qstat = NULL
 
  for (i in 1 : iterations) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
SST = b * sum(xx ^ 2) 
SSE = sum((B-matrix(xx,b,k,byrow=TRUE))^2)
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))
Lstat = c(Lstat, Lambda(A, xx[-k])$max)
Qstat = c(Qstat,SST)
   }
 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
FMC = apply(outer(Fstat, b * u ^ 2, ">"), 2, mean)
QMC = apply(outer(Qstat, b * u ^ 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
F = 1 - pf(b * u ^ 2 , k-1, (k-1) * (b-1))
 
tplr = NULL
tpnc = NULL
  
for (up in u) {
  tt = TailDistrMonteCarlo(A, up, M)
 tplr = c(tplr, tt$LR)
 tpnc = c(tpnc, tt$NC)
 }

rbind(FMC, F=round(F,4), QMC, Q=round(Q,4), L , LR=round(tplr,4), NC=round(tpnc,4))
 

#accuarcy exp
A = generateRandomMatrix(rexp, b, k)

Fstat = NULL
Lstat = NULL
Qstat = NULL
 
  for (i in 1 : iterations) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
SST = b * sum(xx ^ 2) 
SSE = sum((B-matrix(xx,b,k,byrow=TRUE))^2)
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))
Lstat = c(Lstat, Lambda(A, xx[-k])$max)
Qstat = c(Qstat,SST)
   }
 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
FMC = apply(outer(Fstat, b * u ^ 2, ">"), 2, mean)
QMC = apply(outer(Qstat, b * u ^ 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
F = 1 - pf(b * u ^ 2 , k-1, (k-1) * (b-1))
 
tplr = NULL
tpnc = NULL
  
for (up in u) {
  tt = TailDistrMonteCarlo(A, up, M)
 tplr = c(tplr, tt$LR)
 tpnc = c(tpnc, tt$NC)
 }

rbind(FMC, F=round(F,4), QMC, Q=round(Q,4), L , LR=round(tplr,4), NC=round(tpnc,4))
  
#accuracy exp squared
A = generateRandomMatrix(rexp2, b, k)

Fstat = NULL
Lstat = NULL
Qstat = NULL
 
  for (i in 1 : iterations) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
SST = b * sum(xx ^ 2) 
SSE = sum((B-matrix(xx,b,k,byrow=TRUE))^2)
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))
Lstat = c(Lstat, Lambda(A, xx[-k])$max)
Qstat = c(Qstat,SST)
   }
 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
FMC = apply(outer(Fstat, b * u ^ 2, ">"), 2, mean)
QMC = apply(outer(Qstat, b * u ^ 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
F = 1 - pf(b * u ^ 2 , k-1, (k-1) * (b-1))
 
tplr = NULL
tpnc = NULL
  
for (up in u) {
  tt = TailDistrMonteCarlo(A, up, M)
 tplr = c(tplr, tt$LR)
 tpnc = c(tpnc, tt$NC)
 }

rbind(FMC, F=round(F,4), QMC, Q=round(Q,4), L , LR=round(tplr,4), NC=round(tpnc,4))
 
#accuracy gamma5
A = generateRandomMatrix(rgamma5, b, k)

Fstat = NULL
Lstat = NULL
Qstat = NULL
 
  for (i in 1 : iterations) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
SST = b * sum(xx ^ 2) 
SSE = sum((B-matrix(xx,b,k,byrow=TRUE))^2)
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))
Lstat = c(Lstat, Lambda(A, xx[-k])$max)
Qstat = c(Qstat,SST)
   }
 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
FMC = apply(outer(Fstat, b * u ^ 2, ">"), 2, mean)
QMC = apply(outer(Qstat, b * u ^ 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
F = 1 - pf(b * u ^ 2 , k-1, (k-1) * (b-1))
 
tplr = NULL
tpnc = NULL
  
for (up in u) {
  tt = TailDistrMonteCarlo(A, up, M)
 tplr = c(tplr, tt$LR)
 tpnc = c(tpnc, tt$NC)
 }

rbind(FMC, F=round(F,4), QMC, Q=round(Q,4), L , LR=round(tplr,4), NC=round(tpnc,4))
 
#accuracy gamma05
A = generateRandomMatrix(rgamma05, b, k)

Fstat = NULL
Lstat = NULL
Qstat = NULL
 
  for (i in 1 : iterations) {
    B = t(apply(A, 1, sample, size = k))
   xx = apply(B, 2, mean)
SST = b * sum(xx ^ 2) 
SSE = sum((B-matrix(xx,b,k,byrow=TRUE))^2)
Fstat = c(Fstat, (SST/ (k - 1)) / (SSE / ((k - 1) * (b - 1)) ))
Lstat = c(Lstat, Lambda(A, xx[-k])$max)
Qstat = c(Qstat,SST)
   }
 
L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
FMC = apply(outer(Fstat, b * u ^ 2, ">"), 2, mean)
QMC = apply(outer(Qstat, b * u ^ 2, ">"), 2, mean)
Q = 1 - pchisq(b * u ^ 2, k - 1)
F = 1 - pf(b * u ^ 2 , k-1, (k-1) * (b-1))
 
tplr = NULL
tpnc = NULL
  
for (up in u) {
  tt = TailDistrMonteCarlo(A, up, M)
 tplr = c(tplr, tt$LR)
 tpnc = c(tpnc, tt$NC)
 }

rbind(FMC, F=round(F,4), QMC, Q=round(Q,4), L , LR=round(tplr,4), NC=round(tpnc,4))
 