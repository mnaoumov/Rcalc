options(warn=2)

# used for permutations
library(gtools)

# used for maxNR
library(maxLik)

generateRandomMatrix = function(randomGenerator, b, k) {
  A = matrix(randomGenerator(b * k), b, k)
  A = A - apply(A, 1, mean)
  A = A / sqrt(mean(A ^ 2))
}

getb = function(A) dim(A)[1]
getk = function(A) dim(A)[2]

perm = function(A, i) {
  k = getk(A)
  permutations(k, k - 1, A[i,])
}

expb = function(A, i, beta) exp(perm(A, i) %*% beta)

kappa = function(A, beta) {
  b = getb(A)
  k = getk(A)
  s = 0
  for (i in 1 : b) s = s + log(sum(expb(A, i, beta)) / prod(1 : k))
  s / b
}
  
kappad = function(A, beta) {
  b = getb(A)
  k = getk(A)
  s = 0
  for (i in 1 : b) {
    expb = expb(A, i, beta)
    s = s + t(perm(A, i)) %*%  expb / sum(expb)
  }
  
  s / b
}

kappad2 = function(A, beta) {
  b = getb(A)
  k = getk(A)
  s = 0
  for (i in 1 : b) {
    perm = perm(A, i)
    expb = expb(A, i, beta)
    s = s - (t(perm) %*% expb %*% t(expb) %*% perm) / sum(expb) ^ 2 +
      (t(perm) %*% diag(c(expb)) %*% perm) / sum(expb)
  }

  s / b
}

orderedMeans = function(A) {
  b = getb(A)
  k = getk(A)
  ASorted = matrix(numeric(), b, k)
  for(i in 1 : b)
    ASorted[i,] = sort(A[i,])
  
  colMeans(ASorted)
}


getRandomMean = function(A) {
  b = getb(A)
  k = getk(A)
    B = array(dim = c(b, k))
  for(i in 1 : b)
    B[i,] = sample(A[i,], k)

  mean = apply(B, 2, mean)[1 : (k - 1)]
    mean
}

addMu = function(A, mu) {
  B = A
  b = dim(B)[1]
    for (i in 1 : b)
    B[i,] = B[i,] + mu

  B/sqrt(mean(B^2))
}

isPointInsideDomain = function(A, x, tolerance = 1) {
  k = getk(A)
  
  y = x
  if (length(y) == k - 1)
    y[length(y) + 1] = -sum(y)
  
  stopifnot(length(y) == k)
  
  for (usedIndices in 1 : (k - 1)) {
    subVectors = combinations(k, usedIndices, y, set = FALSE)
    subVectorSums = apply(subVectors, 1, sum)
    sumMins = sum(orderedMeans(A)[1 : usedIndices])
    sumMaxs = sum(orderedMeans(A)[(k - usedIndices + 1) : k])
    if (!all(subVectorSums > sumMins * tolerance))
      return(FALSE);
    if (!all(subVectorSums < sumMaxs * tolerance))
      return(FALSE);
  }
  return(TRUE)
}


putInsideDomain = function(A, x, tolerance) {
  if (isPointInsideDomain(A, x, tolerance))
  {
    return (x)
  }
  
  norm = sqrt(sum(x ^ 2))
  
  inside = 0
  outside = norm
  
  # 2^(-15) < 0.0001
  for (i in 1 : 15) {
    r = (inside + outside) / 2
    y = x / norm * r
    if (!isPointInsideDomain(A, y, tolerance))
      outside = r
    else
      inside = r
  }
  
  return (x / norm * inside)
}

Lambda = function(A, x, tolerance = 0.9999) {
  k = getk(A)
  #x = putInsideDomain(A, x, tolerance)
  
  ll = function(beta) { beta %*% x - kappa(A, beta) }
  lll = function(beta) { x - kappad(A, beta)[,1] }
  llll = function(beta) { -kappad2(A, beta) }
  lam = try(maxNR(ll, grad = lll, hess = llll, start = rep(0, k - 1), tol = 1 - tolerance), silent = TRUE)
  
  if ((length(lam) == 1) && (class(lam) == "try-error")) {
    return (list(checkNA =  1))
  }
  
  lam$checkNA = 0
  
  lam
}

Sqr = function(M) { 
z=svd(M)

S=z$u %*% diag(sqrt(z$d)) %*% t(z$v)

if (det(S) < 0) {
print("M=")
print(M)
print("S=")
print(S)
print("svd(M)$u %*% diag(sqrt(svd(M)$d)) %*% t(svd(M)$v)=")
print(svd(M)$u %*% diag(sqrt(svd(M)$d)) %*% t(svd(M)$v))
stop("det(S)<0")
}

S
}

UniformOnSphere = function(n) {
  x = rnorm(n)
  r = sqrt(sum(x ^ 2))
  x / r
}

IntOnSphereMonteCarlo = function(h, d, M) {
  sum = 0
  pointsSkipped = 0
  for(i in 1 : M) {
    funcValue = h(UniformOnSphere(d))
    if (is.na(funcValue))
      pointsSkipped = pointsSkipped + 1
    else
      sum = sum + funcValue
  }
  
  if (pointsSkipped == M) {
    return (0)
  }
  
  I = sum / (M - pointsSkipped)

  if (I < 0) {
    print("I=")
    print(I)
    print("M=")
    print(M)
    print("pointsSkipped=")
    print(pointsSkipped)

    stop("I < 0")  
}

  I
}

V0h = function(A) {
  k = getk(A)
  V = Sqr(kappad2(A, rep(0, k - 1)))  

  if (det(V) < 0) {
    print("A=")
    print(A)
    print("kappad2(A, rep(0, k - 1))=")
    print(kappad2(A, rep(0, k - 1)))
    print("V=")
    print(V)
    print("Sqr(kappad2(A, rep(0, k - 1)))=")
    print(Sqr(kappad2(A, rep(0, k - 1)))  )

    stop("det(V) < 0")
  }

V
}

transformIntoSphere = function(A, r, s) {
  (r * V0h(A) %*% s)[,1]
}

delta = function(A, u, s, deltaIterations = 5) {
  k = getk(A)
  r = u
  V0h = V0h(A)

  if(det(V0h) < 0) {
    print("V0h=")
    print(V0h)
    print("V0h(A)=")
    print(V0h(A))
    stop("det(V0h) < 0")
  }

  
  for (i in 1 : deltaIterations) {
    tr = transformIntoSphere(A, r, s)
    L = Lambda(A, tr)
    
    if (L$checkNA == 1) {
      return (NA)
    }
    r = r - ((L$max - u ^ 2 / 2) / (L$est %*% V0h %*% s))[1,1]
  }
  
  betah = L$est

  d = (det(kappad2(A, betah)) ^ (- 1 / 2) * r ^ (k - 2) * det(V0h)) / (u ^ (k - 3) * abs(t(s) %*% V0h %*% betah))

  if (d < 0) {
    print("d=")
    print(d)
    print("u=")
    print(u)
    print("r=")
    print(r)

    stop("d < 0")
  }

  d  
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


#power

PowerForMu = function(mu, distr, Z, V) {
    
  pFValues = numeric(Z)
  pLRValues = numeric(Z)
  pNCValues = numeric(Z)
  
  for (i in 1 : Z) {
A = generateRandomMatrix(distr, b, k)
A = addMu(A, mu)
A = A - apply(A, 1, mean)
A = A / sqrt(mean(A ^ 2))

Xbar = getRandomMean(A)
Fbar = c(Xbar, -sum(Xbar))
LVal = Lambda(A,Xbar)
   u = sqrt(2*LVal$max)
FVal = b*sum(Fbar ^ 2)/(mean(A^2))

xxx = NULL; for (j in 1 : V) xxx = c(xxx, apply(t(apply(A, 1, sample, size = k)), 2, mean))
xxx = matrix(xxx, k , V)

FValues = b * apply(xxx ^ 2, 2, sum)

    pFValues[i] = mean(FValues > FVal)
    tt = TailDistrMonteCarlo(A,u,M=20,deltaIterations = 5)
    pLRValues[i] = tt$LR
    pNCValues[i] = tt$NC
  }
  
  PowerF = mean(pFValues < .05)
  PowerLR = mean(pLRValues < .05)
  PowerNC = mean(pNCValues < .05)
  
  list (PowerF = PowerF, PowerLR = PowerLR, PowerNC = PowerNC)
}

options(warn=2,error = quote(dump.frames("testdump7", TRUE)))

b=10; k=4;

PowerForMu (mu=c(0,0,0,0),distr=rexp,Z=1000,V=10000)
PowerForMu (mu=c(1,0,0,1),distr=rexp,Z=1000,V=10000)
PowerForMu (mu=c(2,0,0,2),distr=rexp,Z=1000,V=10000)


