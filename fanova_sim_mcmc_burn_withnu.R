
# Simulation and MCMC for functional ANOVA
# Multi-dimensional MH updates for GP parameters 

library(fields)
library(ggplot2)
source("func_anova_fns.R")
set.seed(133204)

### dimensions and sizes
smdm = c(10,10)
lv2 = c(2,2)
nrp = 5
ic=jc=c(-1,1)
nab = (lv2[1]-1)*(lv2[2]-1)

### set true parameters
mu0=3
#sig.mu=4; sig.a=sig.b=3; sig.ab=2; sig.eps=1
#ra.mu=.4; ra.a=ra.b=.3; ra.ab=.2; ra.eps=.1

# covariance params
rnglst = list(mean = 4, mainA = 3, mainB = 3, interac = 2, noise = 1)
siglst = list(mean = 2, mainA = 1.5, mainB = 1.5, interac = 1, noise = 0.5)
simnu = 1.5
nulst = list(mean = simnu, mainA = simnu, mainB = simnu, interac = simnu, noise = simnu)

simout = spat_anova_sim_2fact(smdm,lv2,nrp, rnglst, siglst) 

#quilt.plot(simout$locfrm[,1],simout$locfrm[,2],simout$datarr[,1,1,1])
tfrm = data.frame(XLoc=simout$locfrm[,1],YLoc=simout$locfrm[,2],ObsDat=simout$datarr[,1,1,1])
gtl = ggplot(tfrm,aes(x=XLoc,y=YLoc)) + geom_tile(aes(fill=ObsDat))


# Set up sampler
### posterior inference via Gibbs sampler
L=2000
I=J=2
nlc = nrow(simout$locfrm)
mu.samp=a.samp=b.samp=ab.samp=array(dim=c(nlc,L))
sig.samp=list(mean = rep(0,L), mainA = rep(0,L), mainB = rep(0,L), interac = rep(0,L), noise = rep(0,L))
lam.samp=list(mean = rep(0,L), mainA = rep(0,L), mainB = rep(0,L), interac = rep(0,L), noise = rep(0,L))
nu.samp=list(mean = rep(0,L), mainA = rep(0,L), mainB = rep(0,L), interac = rep(0,L), noise = rep(0,L))
mncr = array(0,c(nlc,I,J))

Cor.mu = Matern(simout$dists,range=rnglst$mean[1],smoothness=simnu)
Cor.mu.eig = eigen(Cor.mu,symmetric = TRUE,only.values = TRUE)
Cor.a = Matern(simout$dists,range=rnglst$mainA[1],smoothness=simnu)
Cor.a.eig = eigen(Cor.a,symmetric = TRUE,only.values = TRUE)
Cor.b = Matern(simout$dists,range=rnglst$mainB[1],smoothness=simnu)
Cor.b.eig = eigen(Cor.b,symmetric = TRUE,only.values = TRUE)
Cor.ab = Matern(simout$dists,range=rnglst$interac[1],smoothness=simnu)
Cor.ab.eig = eigen(Cor.ab,symmetric = TRUE,only.values = TRUE)
Cor.eps = Matern(simout$dists,range=rnglst$noise[1],smoothness=simnu)
Cor.eps.eig = eigen(Cor.eps,symmetric = TRUE,only.values = TRUE)
cordetlst = list(mean = sum(log(Cor.mu.eig$values)), mainA = sum(log(Cor.a.eig$values)), 
                 mainB = sum(log(Cor.b.eig$values)), interac = sum(log(Cor.ab.eig$values)), 
                 noise = sum(log(Cor.eps.eig$values)))

Rinv.mu = solve(Cor.mu)
Rinv.a = solve(Cor.a)
Rinv.b = solve(Cor.a)
Rinv.ab = solve(Cor.ab)
Rinv.eps = solve(Cor.eps)

C.mu=(siglst$mean[1])^2*Cor.mu
C.a=(siglst$mainA[1])^2*Cor.a
C.b=(siglst$mainB[1])^2*Cor.b
C.ab=(siglst$interac[1])^2*Cor.ab
C.eps=(siglst$noise[1])^2*Cor.eps
mu0 = 0

# Metropolis-Hastings setup
ac_mh = list(mean = 0, mainA = 0, mainB = 0, interac = 0, noise = 0)
mhcor3d = matrix(c(1,0.6,-0.4, 0.6,1,-0.8, -0.4,-0.8,1),nrow=3)
#mhcor3d = diag(rep(1,3))
mhsd3d = diag(c(0.6,1.0,0.6))
mhcv3d = mhsd3d %*% mhcor3d %*% mhsd3d
mhchl3d = t(chol(mhcv3d))


# Priors, lognormal for range
prpars = c(0.7,1.0)
prnu = c(0.25,0.25)

for(l in 1:L){
  print(l)
  
  C.mu=(siglst$mean[1])^2*Cor.mu
  C.a=(siglst$mainA[1])^2*Cor.a
  C.b=(siglst$mainB[1])^2*Cor.b
  C.ab=(siglst$interac[1])^2*Cor.ab
  C.eps=(siglst$noise[1])^2*Cor.eps
  
  ## update mu
  pprec=solve(C.mu)+(I*J*nrp)*solve(C.eps)
  y.tilde = apply(simout$datarr,1,sum) - (I*J*nrp)*mu0
  pmean=solve(pprec,solve(C.eps,y.tilde))
  mu.samp[,l] = pmean + t(chol(solve(pprec)))%*%rnorm(nlc)
  
  ## update alpha
  pprec=solve(C.a)+(I*J*nrp)*solve(C.eps)
  y.tilde = ic[1]*apply(simout$datarr[,1,,],1,sum) + ic[2]*apply(simout$datarr[,2,,],1,sum)
  pmean=solve(pprec,solve(C.eps,y.tilde))
  a.samp[,l] = pmean + t(chol(solve(pprec)))%*%rnorm(nlc)
  
  ## update beta
  pprec=solve(C.b)+(I*J*nrp)*solve(C.eps)
  y.tilde = jc[1]*apply(simout$datarr[,,1,],1,sum) + jc[2]*apply(simout$datarr[,,2,],1,sum)
  pmean=solve(pprec,solve(C.eps,y.tilde))
  b.samp[,l] = pmean + t(chol(solve(pprec)))%*%rnorm(nlc)
  
  ## update ab
  pprec=solve(C.ab)+(I*J*nrp)*solve(C.eps)
  y.tilde=rep(0,nlc)
  for(i in 1:I){ for(j in 1:J){
    y.tilde = y.tilde + ic[i]*jc[j]*apply(simout$datarr[,i,j,],1,sum)
  }}
  pmean=solve(pprec,solve(C.eps,y.tilde))
  ab.samp[,l] = pmean + t(chol(solve(pprec)))%*%rnorm(nlc)

  ## Epsilon cov params
  # mu + ic[i]*alpha + jc[j]*beta + ic[i]*jc[j]*ab
  ydv = array(0,c(nlc,nrp*I*J))
  ct0 = 0
  for(i in 1:I) { 
    for(j in 1:J) {
      mncr[,i,j] = mu.samp[,l] + ic[i]*a.samp[,l] + jc[j]*b.samp[,l] + ic[i]*jc[j]*ab.samp[,l]   
      for (k in 1:nrp) {
        ct0 = ct0 + 1
        ydv[,ct0] = simout$datarr[,i,j,k] - mncr[,i,j]
      }
    }
  }
  
  epsmh = spatial_mh3d_matern(siglst$noise,rnglst$noise,nulst$noise,Rinv.eps,Cor.eps,cordetlst$noise,ydv,
                              ct0,nlc,simout$dists,0.05,mhchl3d,2.0,prpars,prnu) 
  ac_mh$noise = ac_mh$noise + epsmh$acpt
  sig.samp$noise[l] = epsmh$sdpst
  siglst$noise[1] = epsmh$sdpst
  lam.samp$noise[l] = epsmh$lampst
  rnglst$noise[1] = epsmh$lampst
  nu.samp$noise[l] = epsmh$nupst
  nulst$noise[1] = epsmh$nupst
  Cor.eps = epsmh$CorMat
  Rinv.eps = epsmh$CorPrc
  cordetlst$noise = epsmh$lgdetcor
  
  abmh = spatial_mh3d_matern(siglst$interac,rnglst$interac,nulst$interac,Rinv.ab,Cor.ab,cordetlst$interac,ab.samp[,l],
                             nab,nlc,simout$dists,0.2,mhchl3d,4.0,prpars,prnu) 
  ac_mh$interac = ac_mh$interac + abmh$acpt
  sig.samp$interac[l] = abmh$sdpst
  siglst$interac[1] = abmh$sdpst
  lam.samp$interac[l] = abmh$lampst
  rnglst$interac[1] = abmh$lampst
  nu.samp$interac[l] = abmh$nupst
  nulst$interac[1] = abmh$nupst
  Cor.ab = abmh$CorMat
  Rinv.ab = abmh$CorPrc
  cordetlst$interac = abmh$lgdetcor
  
  amh = spatial_mh3d_matern(siglst$mainA,rnglst$mainA,nulst$mainA,Rinv.a,Cor.a,cordetlst$mainA,a.samp[,l],
                            lv2[1]-1,nlc,simout$dists,0.2,mhchl3d,4.0,prpars,prnu) 
  ac_mh$mainA = ac_mh$mainA + amh$acpt
  sig.samp$mainA[l] = amh$sdpst
  siglst$mainA[1] = amh$sdpst
  lam.samp$mainA[l] = amh$lampst
  rnglst$mainA[1] = amh$lampst
  nu.samp$mainA[l] = amh$nupst
  nulst$mainA[1] = amh$nupst
  Cor.a = amh$CorMat
  Rinv.a = amh$CorPrc
  cordetlst$mainA = amh$lgdetcor
  
  bmh = spatial_mh3d_matern(siglst$mainB,rnglst$mainB,nulst$mainB,Rinv.b,Cor.b,cordetlst$mainB,b.samp[,l],
                            lv2[2]-1,nlc,simout$dists,0.2,mhchl3d,4.0,prpars,prnu) 
  ac_mh$mainB = ac_mh$mainB + bmh$acpt
  sig.samp$mainB[l] = bmh$sdpst
  siglst$mainB[1] = bmh$sdpst
  lam.samp$mainB[l] = bmh$lampst
  rnglst$mainB[1] = bmh$lampst
  nu.samp$mainB[l] = bmh$nupst
  nulst$mainB[1] = bmh$nupst
  Cor.b = bmh$CorMat
  Rinv.b = bmh$CorPrc
  cordetlst$mainB = bmh$lgdetcor
  
  mumh = spatial_mh3d_matern(siglst$mean,rnglst$mean,nulst$mean,Rinv.mu,Cor.mu,cordetlst$mean,mu.samp[,l],
                             1,nlc,simout$dists,0.2,mhchl3d,4.0,prpars,prnu) 
  ac_mh$mean = ac_mh$mean + mumh$acpt
  sig.samp$mean[l] = mumh$sdpst
  siglst$mean[1] = mumh$sdpst
  lam.samp$mean[l] = mumh$lampst
  rnglst$mean[1] = mumh$lampst
  nu.samp$mean[l] = mumh$nupst
  nulst$mean[1] = mumh$nupst
  Cor.mu = mumh$CorMat
  Rinv.mu = mumh$CorPrc
  cordetlst$mean = mumh$lgdetcor
  
}


