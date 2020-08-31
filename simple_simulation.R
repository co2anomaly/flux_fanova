
####### GP ANOVA for simulated 2x2 example  #######

library(fields)


### dimensions and sizes
n=10^2
I=J=2
ic=jc=c(-1,1)
K=5


### set true parameters
mu0=3
sig.mu=4; sig.a=sig.b=3; sig.ab=2; sig.eps=1
ra.mu=.4; ra.a=ra.b=.3; ra.ab=.2; ra.eps=.1


### locations and true covariance matrices
grid.oneside=seq(0,1,length.out=sqrt(n))
locs=as.matrix(expand.grid(grid.oneside,grid.oneside))
C.mu=sig.mu^2*Matern(rdist(locs),range=ra.mu,smoothness=1.5)
C.a=sig.a^2*Matern(rdist(locs),range=ra.a,smoothness=1.5)
C.b=sig.b^2*Matern(rdist(locs),range=ra.b,smoothness=1.5)
C.ab=sig.ab^2*Matern(rdist(locs),range=ra.ab,smoothness=1.5)
C.eps=sig.eps^2*Matern(rdist(locs),range=ra.eps,smoothness=1.5)


### simulate latent fields
mu=t(chol(C.mu))%*%rnorm(n)
alpha=t(chol(C.a))%*%rnorm(n)
beta=t(chol(C.b))%*%rnorm(n)
ab=t(chol(C.ab))%*%rnorm(n)


### simulate data
y=array(dim=c(n,I,J,K))
for(i in 1:I){
  for(j in 1:J){
    for(k in 1:K){
      eps=t(chol(C.eps))%*%rnorm(n)
      y[,i,j,k] = mu0 + mu + ic[i]*alpha + jc[j]*beta + ic[i]*jc[j]*ab + eps
    }
  }
}


### plots
quilt.plot(locs[,1],locs[,2],mu,nx=sqrt(n),ny=sqrt(n))
quilt.plot(locs[,1],locs[,2],alpha,nx=sqrt(n),ny=sqrt(n))
quilt.plot(locs[,1],locs[,2],beta,nx=sqrt(n),ny=sqrt(n))
quilt.plot(locs[,1],locs[,2],ab,nx=sqrt(n),ny=sqrt(n))
quilt.plot(locs[,1],locs[,2],eps,nx=sqrt(n),ny=sqrt(n))



### posterior inference via Gibbs sampler
L=100
mu.samp=a.samp=b.samp=ab.samp=array(dim=c(n,L))

for(l in 1:L){
  print(l)
  
  ## update mu
  pprec=solve(C.mu)+(I*J*K)*solve(C.eps)
  y.tilde = apply(y,1,sum) - (I*J*K)*mu0
  pmean=solve(pprec,solve(C.eps,y.tilde))
  mu.samp[,l] = pmean + t(chol(solve(pprec)))%*%rnorm(n)
  
  ## update alpha
  pprec=solve(C.a)+(I*J*K)*solve(C.eps)
  y.tilde = ic[1]*apply(y[,1,,],1,sum) + ic[2]*apply(y[,2,,],1,sum)
  pmean=solve(pprec,solve(C.eps,y.tilde))
  a.samp[,l] = pmean + t(chol(solve(pprec)))%*%rnorm(n)
  
  ## update beta
  pprec=solve(C.b)+(I*J*K)*solve(C.eps)
  y.tilde = jc[1]*apply(y[,,1,],1,sum) + jc[2]*apply(y[,,2,],1,sum)
  pmean=solve(pprec,solve(C.eps,y.tilde))
  b.samp[,l] = pmean + t(chol(solve(pprec)))%*%rnorm(n)
  
  ## update ab
  pprec=solve(C.ab)+(I*J*K)*solve(C.eps)
  y.tilde=rep(0,n)
  for(i in 1:I){ for(j in 1:J){
    y.tilde = y.tilde + ic[i]*jc[j]*apply(y[,i,j,],1,sum)
  }}
  pmean=solve(pprec,solve(C.eps,y.tilde))
  ab.samp[,l] = pmean + t(chol(solve(pprec)))%*%rnorm(n)
  
}


### maps of true latent fields and posterior means
par(mfrow=c(1,2))

quilt.plot(locs[,1],locs[,2],mu,nx=sqrt(n),ny=sqrt(n))
quilt.plot(locs[,1],locs[,2],rowMeans(mu.samp),nx=sqrt(n),ny=sqrt(n))

quilt.plot(locs[,1],locs[,2],alpha,nx=sqrt(n),ny=sqrt(n))
quilt.plot(locs[,1],locs[,2],rowMeans(a.samp),nx=sqrt(n),ny=sqrt(n))

quilt.plot(locs[,1],locs[,2],beta,nx=sqrt(n),ny=sqrt(n))
quilt.plot(locs[,1],locs[,2],rowMeans(b.samp),nx=sqrt(n),ny=sqrt(n))

quilt.plot(locs[,1],locs[,2],ab,nx=sqrt(n),ny=sqrt(n))
quilt.plot(locs[,1],locs[,2],rowMeans(ab.samp),nx=sqrt(n),ny=sqrt(n))

par(mfrow=c(1,1))



### other posterior summaries

# posterior prob that effect of alpha larger than effect of beta
s2.a=2*a.samp^2
s2.b=2*b.samp^2
quilt.plot(locs[,1],locs[,2],rowMeans(s2.a>s2.b),nx=sqrt(n),ny=sqrt(n),col=grey.colors(100))
