options(warn=2, error=quote(dump.frames("k-sample-accuracy-all-distr", TRUE)))

kcgf= function(tau,xd, f)
{
#This function calculates K(tau), K"(tau) and K"(tau), where
#t is a 2*(k-1) dimensional matrix, K(t) is a real number.
#Inputs:t  in K(t).      
#       xd=the population: a=1,a=2, ..., a=N.
#       f=vector of ratios of group sizes (n_2,..  .,n_{k-1}) and population size (N), 
#Outputs:cgf2$cgf:K(tau)
#        cgf2$dcgf:K'(tau)
#        cgf2$ddc gf:K"(tau)
#           
        taum=matrix(tau,length(f),2)
        nn = length(xd)
        xsd = sqrt((var(xd) * (nn - 1))/nn)
        xdb = mean(xd)
        x = (xd - xdb)/xsd
        f1=1-sum(f)
        nn1=rep(1,nn)
        x2=rbind(nn1,x)
        kt = mean(log(f1 + f%*%exp(taum%*%x2)))
    h=t(t(f * exp(taum%*%x2))/c(f1 + f %*% exp(taum%*%x2)))
        hx=t(x*t(f * exp(taum%*%x2))/c(f1 + f %*% exp(taum%*%x2)))
        hhx=rbind(h,hx)
        dk = c(h%*%t(x2)/nn)
        ddk= (cbind(rbind(diag(c(h%*%nn1)),diag(c(h%*%x))),rbind(diag(c(h%*%x)),diag(c(hx%*%x))))-hhx%*%t(hhx))/nn
        list(cgf = kt, dcgf = dk, ddcgf = ddk)
}

kLam=function(x1,t0,xd,f)
#kLam=function(x1,xd,f)

{
#This function solves  K'(t)=x for a given x=(f,x1).   
#Inputs: x1: a given k-1-dimension vector.
#        tt: an initial value for Newton's method of finding the 
#            root of  k'(t)-x for a given x.
#  xd: population.
#  f=n/N
#Outputs: Lam$lam=\Lambda(x)
#         Lam$t=t(x), 3-dimension vector, t(x)%*%x-K(t(x))=sup{t.x-K(t): t}
#         Lam$$kdd K"(t(x))
#Get tt=tt(x), tt(x)%*%x-K(t0(x))=sup{t%*%x-K(t): t}

       fx=c(f,x1)
       checkNA=0
       lf=length(f)
       tt=rep(0,2*lf)
      for (i in 1:10)
  {   
      ktt=kcgf(tt,xd,f)
      if (!is.na(det(ktt$ddcgf)) && det(ktt$ddcgf)>.00000001)
      {tt=tt+solve(ktt$ddcgf)%*%(fx-ktt$dcgf)} else {{checkNA=1}&{break}}#t0=t0+solve(ktt$ddcgf)%*%(fx-ktt$dcgf)
        }
  ktt=kcgf(tt,xd,f)
  list(lam=-ktt$cgf+fx%*%tt,t=tt,kdd=ktt$ddcgf,err=fx-ktt$dcgf,checkNA=checkNA)
}




kHu=function(xd, f, v,M)
{
#this function calculates SP tail probs for k-sample perm test based on Lam(x)
#input xd data, f group proportions for k-1, v P(Lam(X)>v),M= no of MC in integral
#first define 

        nn = length(xd)
        n = f * nn
        xsd = sqrt((var(xd) * (nn - 1))/nn)
        xdb = mean(xd)
        x = (xd - xdb)/xsd
        f1=1-sum(f)
        lf=length(f)
        tailp1=1-pchisq(nn*v^2,lf)#first order approximation
        VV=kcgf(rep(0,2*lf),x,f)$ddcgf[(lf+1):(2*lf),(lf+1):(2*lf)]#xd
        svdV=svd(VV)
        Vh=svdV$u%*%diag(svdV$d^(1/2))%*%t(svdV$v)
        Gu=0
        sdelus=matrix(0,M,lf+4)
        Mm=M                              
        for (i in 1:M){
                       s=rnorm(lf)
                       s=s/sqrt(sum(s^2))
                       s=c(Vh%*%s)
                       r=v
                       #t0=c(kLam(xb,rep(0,2*lf),x,f)$t)
                       t0=rep(0,2*lf)
                      
                       for (j in 1:6){ 
                           kl=kLam(r*s,t0,x,f)                        
                           if (kl$checkNA == 1)
                             break;
                           r=r+(v-sqrt(2*kl$lam))*sqrt(2*kl$lam)/sum(s*kl$t[(lf+1):(2*lf)])                          
                           t0=c(kl$t)
                                     }
                       kl=kLam(r*s,t0,x,f)#xd
                    
                       if (kl$checkNA==1){delus=NA}else {delus= f1*prod(f)*r^2/(v*sum(s*kl$t[(lf+1):(2*lf)])*(det(kl$kdd))^.5)}  
                       #Gu=ifelse(abs(v-sqrt(2*kl$lam))>.00001,Gu,Gu+delus)
                       if (kl$checkNA==0) Gu=Gu+delus
                       if (kl$checkNA==1) Mm=Mm-1 
sdelus[i,]=c(delus,r*s,v-sqrt(2*kl$lam),sum(s*kl$t[(lf+1):(2*lf)]),det(kcgf(kl$t,xd,f)$ddcgf))
                      #if(abs(v-sqrt(2*kl$lam))>.00001) print(c(abs(v-sqrt(2*temk)),det(kcgf(kl$t,xd,f)$ddcgf)))
                       }
        cn=(nn^(lf/2))/((2)^(lf/2-1)*gamma(lf/2))      
        tailp=ifelse(Mm==0,tailp1,tailp1+nn^(-1)*cn*v^(lf-2)*exp(-nn*v^2/2)*(Gu/Mm-1))
        tailpBN=ifelse(Mm==0,tailp1,1-pchisq(nn*(v-log(Gu/Mm)/(nn*v))^2,lf))
        list(tailpchi2=tailp1,tailp=c(tailp,tailpBN),sdelus=sdelus,Mm=Mm)
}





 pvalcal=function(xd,f,M=0,MC=0,vv=0)
                      {
                       nn=length(xd)
                       n=f*nn
                       xsd = sqrt((var(xd) * (nn - 1))/nn)
                       xdb = mean(xd)
                       x = (xd - xdb)/xsd
                       f1=1-sum(f)
                       lf=length(f)
                       #N=c(cumsum(n),nn)
                       N=c(0,cumsum(n)) 
                       xb=orderedMeans=rep(0,lf)
                       for(i in 1:lf){xb[i]=sum(x[(N[i]+1):N[i+1]])/nn}
                      
                       #xs=sort(x)
                       #for(i in 1:lf){orderedMeans[i]=sum(xs[(N[i]+1):N[i+1]])/nn}
                       #orderedMeans=c(-sum(orderedMeans),orderedMeans) 
                       vl=kLam(xb,rep(0,2*lf),xd,f)
                       v=ifelse(vv==0,sqrt(2*vl$lam),vv)
                       out=ifelse(M==0,1-pchisq(nn*v^2,lf),kHu(xd,f,v,M)$tailp)
                       #MCcal=function(MC){
                       MCpv=rep(0,MC)
                       MCsqpv=rep(0,MC)
                       MCm=MC
                       for (k in 1:MC){
                                       xs=sample(x,nn)
                                       xsb=rep(0,lf)
                                       for(i in 1:lf){xsb[i]=sum(xs[(N[i]+1):N[i+1]])/nn}
                                       #MCpv=MCpv+(kLam(xsb,rep(0,2*lf),xd,f)$lam>v^2/2)
                                       MCpv[k]=kLam(xsb,rep(0,2*lf),xd,f)$lam
                                       xsbk=c(xsb,-sum(xsb))
                                       MCsqpv[k]=nn*xsbk%*%diag(1/c(f,f1))%*%xsbk
                                      }
                       #MCpv= MCpv/(MC+1)
                                                  
                       #MCpv=ifelse(MC==0,"null",MCpv)#MCcal(MC))
                       list(SPpv=out,MCpv=MCpv,MCsqpv=MCsqpv,xb=xb)
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


u=(1:7)/10
iterations = 100000
f = c(.25,.25,.25)

#accuracy norm
xd = rnorm(40)
xsd = sqrt(var(xd)*(nn-1)/nn)
xdb = mean(xd)
xd = (xd - xdb)/xsd
lf = length(f)   #k-1
f1 = 1-sum(f)
nn = length(xd)
n = f*nn
N = c(0,cumsum(n))

Fstat = NULL
Lstat = NULL
Qstat = NULL
 
  for (i in 1 : iterations) {
   xs = sample(xd,nn)
   xsb = rep(0,lf)
   for (j in 1:lf) {xsb[j]=sum(xs[(N[j]+1):N[j+1]])/nn}
   xsbk = c(xsb,-sum(xsb))
   Fstat = c(Fstat, (nn*xsbk%*%diag(1/c(f,f1))%*%xsbk/lf)/    
((nn - nn*xsbk%*%diag(1/c(f,f1))%*%xsbk)/(nn-lf-1))  )
   Qstat = c(Qstat, nn*xsbk%*%diag(1/c(f,f1))%*%xsbk)
   Lstat = c(Lstat,kLam(xsb,rep(0,2*lf),xd,f)$lam)   
   }

L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
FMC = apply(outer(Fstat, nn * u ^ 2, ">"), 2, mean)
QMC = apply(outer(Qstat, nn * u ^ 2, ">"), 2, mean)
Q = 1 - pchisq(nn * u ^ 2, lf)
F = 1 - pf(nn * u ^ 2 , lf, (nn-lf-1))


tplr = NULL
tpnc = NULL
for (up in u) {
  tt2 = kHu(xd,f,up,100)
  tplr = c(tplr,tt2$tailp[1])
  tpnc = c(tpnc,tt2$tailp[2])
}
round(rbind(FMC, F, QMC, Q, L, LR=tplr, BN=tpnc),5)

#accuracy exp
xd = rexp(40)
xsd = sqrt(var(xd)*(nn-1)/nn)
xdb = mean(xd)
xd = (xd - xdb)/xsd
lf = length(f)   #k-1
f1 = 1-sum(f)
nn = length(xd)
n = f*nn
N = c(0,cumsum(n))

Fstat = NULL
Lstat = NULL
Qstat = NULL
 
  for (i in 1 : iterations) {
   xs = sample(xd,nn)
   xsb = rep(0,lf)
   for (j in 1:lf) {xsb[j]=sum(xs[(N[j]+1):N[j+1]])/nn}
   xsbk = c(xsb,-sum(xsb))
   Fstat = c(Fstat, (nn*xsbk%*%diag(1/c(f,f1))%*%xsbk/lf)/    
((nn - nn*xsbk%*%diag(1/c(f,f1))%*%xsbk)/(nn-lf-1))  )
   Qstat = c(Qstat, nn*xsbk%*%diag(1/c(f,f1))%*%xsbk)
   Lstat = c(Lstat,kLam(xsb,rep(0,2*lf),xd,f)$lam)   
   }

L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
FMC = apply(outer(Fstat, nn * u ^ 2, ">"), 2, mean)
QMC = apply(outer(Qstat, nn * u ^ 2, ">"), 2, mean)
Q = 1 - pchisq(nn * u ^ 2, lf)
F = 1 - pf(nn * u ^ 2 , lf, (nn-lf-1))


tplr = NULL
tpnc = NULL
for (up in u) {
  tt2 = kHu(xd,f,up,100)
  tplr = c(tplr,tt2$tailp[1])
  tpnc = c(tpnc,tt2$tailp[2])
}
round(rbind(FMC, F, QMC, Q, L, LR=tplr, BN=tpnc),5)


#accuracy gamma5
xd = rgamma(40, 5)
xsd = sqrt(var(xd)*(nn-1)/nn)
xdb = mean(xd)
xd = (xd - xdb)/xsd
lf = length(f)   #k-1
f1 = 1-sum(f)
nn = length(xd)
n = f*nn
N = c(0,cumsum(n))

Fstat = NULL
Lstat = NULL
Qstat = NULL
 
  for (i in 1 : iterations) {
   xs = sample(xd,nn)
   xsb = rep(0,lf)
   for (j in 1:lf) {xsb[j]=sum(xs[(N[j]+1):N[j+1]])/nn}
   xsbk = c(xsb,-sum(xsb))
   Fstat = c(Fstat, (nn*xsbk%*%diag(1/c(f,f1))%*%xsbk/lf)/    
((nn - nn*xsbk%*%diag(1/c(f,f1))%*%xsbk)/(nn-lf-1))  )
   Qstat = c(Qstat, nn*xsbk%*%diag(1/c(f,f1))%*%xsbk)
   Lstat = c(Lstat,kLam(xsb,rep(0,2*lf),xd,f)$lam)   
   }

L = apply(outer(Lstat, u ^ 2 / 2, ">"), 2, mean)
FMC = apply(outer(Fstat, nn * u ^ 2, ">"), 2, mean)
QMC = apply(outer(Qstat, nn * u ^ 2, ">"), 2, mean)
Q = 1 - pchisq(nn * u ^ 2, lf)
F = 1 - pf(nn * u ^ 2 , lf, (nn-lf-1))


tplr = NULL
tpnc = NULL
for (up in u) {
  tt2 = kHu(xd,f,up,100)
  tplr = c(tplr,tt2$tailp[1])
  tpnc = c(tpnc,tt2$tailp[2])
}
round(rbind(FMC, F, QMC, Q, L, LR=tplr, BN=tpnc),5)