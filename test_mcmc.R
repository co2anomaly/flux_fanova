
# Simulation and MCMC for functional ANOVA

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

simout = spat_anova_sim_2fact(smdm,lv2,nrp, rnglst, siglst) 

#quilt.plot(simout$locfrm[,1],simout$locfrm[,2],simout$datarr[,1,1,1])
tfrm = data.frame(XLoc=simout$locfrm[,1],YLoc=simout$locfrm[,2],ObsDat=simout$datarr[,1,1,1])
gtl = ggplot(tfrm,aes(x=XLoc,y=YLoc)) + geom_tile(aes(fill=ObsDat))


# Set up sampler
### posterior inference via Gibbs sampler
L=500
I=J=2
nlc = nrow(simout$locfrm)
mu.samp=a.samp=b.samp=ab.samp=array(dim=c(nlc,L))
sig.samp=list(mean = rep(0,L), mainA = rep(0,L), mainB = rep(0,L), interac = rep(0,L), noise = rep(0,L))
lam.samp=list(mean = rep(0,L), mainA = rep(0,L), mainB = rep(0,L), interac = rep(0,L), noise = rep(0,L))
mncr = array(0,c(nlc,I,J))

Cor.mu = Matern(simout$dists,range=rnglst$mean[1],smoothness=1.5)
Cor.mu.eig = eigen(Cor.mu,symmetric = TRUE,only.values = TRUE)
Cor.a = Matern(simout$dists,range=rnglst$mainA[1],smoothness=1.5)
Cor.a.eig = eigen(Cor.a,symmetric = TRUE,only.values = TRUE)
Cor.b = Matern(simout$dists,range=rnglst$mainB,smoothness=1.5)
Cor.b.eig = eigen(Cor.b,symmetric = TRUE,only.values = TRUE)
Cor.ab = Matern(simout$dists,range=rnglst$interac,smoothness=1.5)
Cor.ab.eig = eigen(Cor.ab,symmetric = TRUE,only.values = TRUE)
Cor.eps = Matern(simout$dists,range=rnglst$noise,smoothness=1.5)
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

#ac_sig = 0
ac_sig = list(mean = 0, mainA = 0, mainB = 0, interac = 0, noise = 0)
ac_rng = list(mean = 0, mainA = 0, mainB = 0, interac = 0, noise = 0)

# Priors, lognormal for range
prpars = c(0.7,1.0)

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
  
  sgepsmh = spatial_mh_std_dev(siglst$noise,ydv,ct0,nlc,Rinv.eps,0.1,5)
  ac_sig$noise = ac_sig$noise + sgepsmh$acpt
  sig.samp$noise[l] = sgepsmh$sdpst
  siglst$noise[1] = sgepsmh$sdpst
  
  #sgabmh = spatial_mh_std_dev(siglst$interac,ab.samp[,l],nab,nlc,Rinv.ab,0.25,5)
  #ac_sig$interac = ac_sig$interac + sgabmh$acpt
  #sig.samp$interac[l] = sgabmh$sdpst
  #siglst$interac[1] = sgabmh$sdpst

  #sgamh = spatial_mh_std_dev(siglst$mainA,a.samp[,l],lv2[1]-1,nlc,Rinv.a,0.25,5)
  #ac_sig$mainA = ac_sig$mainA + sgamh$acpt
  #sig.samp$mainA[l] = sgamh$sdpst
  #siglst$mainA[1] = sgamh$sdpst
  
  #sgbmh = spatial_mh_std_dev(siglst$mainB,b.samp[,l],lv2[2]-1,nlc,Rinv.b,0.25,5)
  #ac_sig$mainB = ac_sig$mainB + sgbmh$acpt
  #sig.samp$mainB[l] = sgbmh$sdpst
  #siglst$mainB[1] = sgbmh$sdpst

  #sgmumh = spatial_mh_std_dev(siglst$mean,mu.samp[,l],1,nlc,Rinv.mu,0.25,5)
  #ac_sig$mean = ac_sig$mean + sgmumh$acpt
  #sig.samp$mean[l] = sgmumh$sdpst
  #siglst$mean[1] = sgmumh$sdpst
  
  # Range Parameters
  lamepsmh = spatial_mh_matern_range(rnglst$noise,Rinv.eps,Cor.eps,cordetlst$noise,ydv,ct0,nlc,
                                     simout$dists,simnu,siglst$noise,0.05,prpars) 
  ac_rng$noise = ac_rng$noise + lamepsmh$acpt
  lam.samp$noise[l] = lamepsmh$lampst
  rnglst$noise[1] = lamepsmh$lampst
  Cor.eps = lamepsmh$CorMat
  Rinv.eps = lamepsmh$CorPrc
  cordetlst$noise = lamepsmh$lgdetcor
  
}


