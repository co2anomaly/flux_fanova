# Functional ANOVA supporting routines
library(fields)
library(reshape2)
library(ggplot2)
library(ncdf4)

spat_anova_sim_2fact = function(flddim,nlev,nrep, sp_range, sp_stddv) {
  # Simulation of random fields from ANOVA model with two factors
  #     flddim    - 2 dimensions, x and y grid sizes
  #     nlev      - 2 dimensions, number of levels for each factor
  #     nrep      - number of replicates for each level
  #     sp_range  - list of spatial range parameters for ANOVA fields
  #     sp_stddv  - list of spatial standard deviations for ANOVA fields
  
  #     Covariance parameter lists have elements:
  #     mean, mainA, mainB, interac, noise
  
  xsq = seq(1,flddim[1])
  ysq = seq(1,flddim[2])
  locs=as.matrix(expand.grid(xsq,ysq))
  n = nrow(locs)
  dstmt = rdist(locs)
  
  C.mu=(sp_stddv$mean[1])^2*Matern(dstmt,range=sp_range$mean[1],smoothness=1.5)
  C.a=(sp_stddv$mainA[1])^2*Matern(dstmt,range=sp_range$mainA[1],smoothness=1.5)
  C.b=(sp_stddv$mainB[1])^2*Matern(dstmt,range=sp_range$mainB,smoothness=1.5)
  C.ab=(sp_stddv$interac[1])^2*Matern(dstmt,range=sp_range$interac,smoothness=1.5)
  C.eps=(sp_stddv$noise[1])^2*Matern(dstmt,range=sp_range$noise,smoothness=1.5)
  
  mu=t(chol(C.mu))%*%rnorm(n)
  alpha=t(chol(C.a))%*%rnorm(n)
  beta=t(chol(C.b))%*%rnorm(n)
  ab=t(chol(C.ab))%*%rnorm(n)
  ic=jc=c(-1,1)
  
  ### simulate data
  y=array(dim=c(n,nlev[1],nlev[2],nrep))
  I = nlev[1]
  J = nlev[2]
  for(i in 1:I){
    for(j in 1:J){
      for(k in 1:nrep){
        eps=t(chol(C.eps))%*%rnorm(n)
        y[,i,j,k] = mu + ic[i]*alpha + jc[j]*beta + ic[i]*jc[j]*ab + eps
      }
    }
  }
  
  frmout = list(dists=dstmt, datarr=y, locfrm=locs, 
                mu=mu, alpha=alpha, beta=beta, ab=ab)
  return(frmout)
}

spat_anova_sim_2fact_mltlev = function(flddim,nlev,nrep, ctrlst, sp_range, sp_stddv, sp_smooth) {
  # Simulation of random fields from ANOVA model with two factors, multiple levels
  #     flddim    - 2 dimensions, x and y grid sizes
  #     nlev      - 2 dimensions, number of levels for each factor
  #     nrep      - number of replicates for each level
  #     ctrlst    - list of contrast matrices (2 elements)
  #     sp_range  - list of spatial range parameters for ANOVA fields
  #     sp_stddv  - list of spatial standard deviations for ANOVA fields
  #     sp_smooth - list of Matern smoothness parameters for ANOVA fields
  
  #     Covariance parameter lists have elements:
  #     mean, mainA, mainB, interac, noise
  
  xsq = seq(1,flddim[1])
  ysq = seq(1,flddim[2])
  locs=as.matrix(expand.grid(xsq,ysq))
  n = nrow(locs)
  dstmt = rdist(locs)
  nab = (nlev[1]-1) * (nlev[2]-1)
  I = nlev[1]
  J = nlev[2]
  
  C.mu=(sp_stddv$mean[1])^2*Matern(dstmt,range=sp_range$mean[1],smoothness=sp_smooth$mean[1])
  C.a=(sp_stddv$mainA[1])^2*Matern(dstmt,range=sp_range$mainA[1],smoothness=sp_smooth$mainA[1])
  C.b=(sp_stddv$mainB[1])^2*Matern(dstmt,range=sp_range$mainB[1],smoothness=sp_smooth$mainB[1])
  C.ab=(sp_stddv$interac[1])^2*Matern(dstmt,range=sp_range$interac[1],smoothness=sp_smooth$interac[1])
  C.eps=(sp_stddv$noise[1])^2*Matern(dstmt,range=sp_range$noise[1],smoothness=sp_smooth$noise[1])
  
  mu=t(chol(C.mu))%*%rnorm(n)
  alpha=t(chol(C.a))%*%  matrix(rnorm(n*(I-1)), nrow=n, ncol=I-1 )
  beta=t(chol(C.b))%*% matrix(rnorm(n*(J-1)), nrow=n, ncol=J-1 )
  ab=t(chol(C.ab))%*% matrix(rnorm(n*nab), nrow=n, ncol=nab )
  
  print(dim(alpha))
  ctra = ctrlst[[1]]
  ctrb = ctrlst[[2]]
  ### simulate data
  print(dim(ctra))
  y=array(dim=c(n,nlev[1],nlev[2],nrep))
  for(i in 1:I){
    for(j in 1:J){
      ictrst = ctra[i,] %x% ctrb[j,]
      for(k in 1:nrep){
        eps=t(chol(C.eps))%*%rnorm(n)
        y[,i,j,k] = mu + alpha %*% as.vector(ctra[i,]) + beta %*% as.vector(ctrb[j,]) + ab %*% as.vector(ictrst) + eps
      }
    }
  }
  
  frmout = list(dists=dstmt, datarr=y, locfrm=locs, 
                mu=mu, alpha=alpha, beta=beta, ab=ab)
  return(frmout)
}

spatial_mh_std_dev = function(cursd,datarr,nflds,nloc,corprc,mhsd,prsd) {
  # Spatial covariance standard deviation update
  # Metropolis-Hastings update of log std dev
  # Half normal prior
  #     cursd     - Current process standard deviation
  #     datarr    - Data array [nloc,nflds]
  #     nflds     - Number of fields/replicates
  #     nloc      - Number of locations
  #     corprc    - Inverse of spatial correlation matrix
  #     mhsd      - MH jumping standard deviation
  #     prsd      - Half normal prior standard deviation

  lgsdcur = log(cursd)
  lgsdprp = lgsdcur + rnorm(1,sd = mhsd)
  prpsd = exp(lgsdprp)
    
  curvr = cursd * cursd
  prpvr = prpsd * prpsd
  
  # Prior density
  curpst = -0.5 * curvr / (prsd*prsd)
  prppst = -0.5 * prpvr / (prsd*prsd)
  
  # Likelihood contributions
  if (nflds == 1) {
    dtdst = t(as.vector(datarr)) %*% corprc %*% (as.vector(datarr))
    curpst = curpst - 0.5 * dtdst / curvr - nloc * log(curvr) / 2.0
    prppst = prppst - 0.5 * dtdst / prpvr - nloc * log(prpvr) / 2.0
  }
  else {
    for (i in seq(1,nflds)) {
      dtdst = t(datarr[,i]) %*% corprc %*% datarr[,i]
      curpst = curpst - 0.5 * dtdst / curvr - nloc * log(curvr) / 2.0
      prppst = prppst - 0.5 * dtdst / prpvr - nloc * log(prpvr) / 2.0
    }
  }
  # MH Asymmetry
  curpst = curpst + lgsdcur
  prppst = prppst + lgsdprp
  
  # MH Evaluation
  logr = prppst - curpst
  logu = log(runif(1))
  #print(logr)
  if (logu[1] < logr[1]) {
    mhout = list(acpt = 1,sdpst = prpsd)
  }
  else {
    mhout = list(acpt = 0,sdpst = cursd)
  }
  
  return(mhout)
}

spatial_mh_matern_range = function(currg,currprc,currcmt,logrdet,datarr,nflds,nloc,dstmt,sptnu,sptsd,mhsd,prpars) {
  # Spatial covariance Matern range parameter update
  # Metropolis-Hastings update of log range
  # Prior?
  #     currg     - Current range parameter
  #     currprc   - Inverse of spatial correlation matrix (current)
  #     currcmt   - Correlation matrix (current)
  #     logrdet   - Log determinant for currprc, will be computed if NULL
  #     cursd     - Current process standard deviation
  #     datarr    - Data array [nloc,nflds]
  #     nflds     - Number of fields/replicates
  #     nloc      - Number of locations
  #     dstmt     - Distance matrix, assume a single set valid for all nflds
  #     sptnu     - Spatial smoothness parameter
  #     sptsd     - Marginal standard deviation
  #     mhsd      - MH jumping standard deviation
  #     prpars    - Lognormal prior mean and standard deviation
  
  lgrgcur = log(currg)
  lgsdprp = lgrgcur + rnorm(1,sd = mhsd)
  prprg = exp(lgsdprp)
  
  curvr = sptsd * sptsd
  
  # Proposal correlation structure
  prpR = Matern(dstmt,range=prprg,smoothness=sptnu)
  egprpR = eigen(prpR,only.values = TRUE)
  prprdet = sum(log(egprpR$values))
  prpRinv = solve(prpR)
  if (is.null(logrdet)) {
    egcurR = eigen(prpR,only.values = TRUE)
    currdet = sum(log(egcurR$values))
  }
  else {
    currdet = logrdet
  }
  
  # Prior density
  curpst = dnorm(lgrgcur,mean=prpars[1],sd=prpars[2],log = TRUE)
  prppst = dnorm(lgsdprp,mean=prpars[1],sd=prpars[2],log = TRUE)
  
  # Likelihood contributions
  curpst = curpst - 0.5 * nflds * logrdet
  prppst = prppst - 0.5 * nflds * prprdet
  if (nflds == 1) {
    dtdst = t(as.vector(datarr)) %*% currprc %*% (as.vector(datarr))
    dtprp = t(as.vector(datarr)) %*% prpRinv %*% (as.vector(datarr))
    curpst = curpst - 0.5 * dtdst / curvr 
    prppst = prppst - 0.5 * dtprp / curvr 
  }
  else {
    for (i in seq(1,nflds)) {
      dtdst = t(datarr[,i]) %*% currprc %*% datarr[,i]
      dtprp = t(datarr[,i]) %*% prpRinv %*% datarr[,i]
      curpst = curpst - 0.5 * dtdst / curvr 
      prppst = prppst - 0.5 * dtprp / curvr 
    }
  }
  # No MH Asymmetry due to lognormal prior
  #curpst = curpst + lgrgcur
  #prppst = prppst + lgsdprp
  
  # MH Evaluation
  logr = prppst - curpst
  logu = log(runif(1))
  print(logr)
  print(prprdet)
  print(currdet)
  if (logu[1] < logr[1]) {
    mhout = list(acpt = 1,lampst = prprg,lgdetcor=prprdet,CorMat=prpR,CorPrc=prpRinv)
  }
  else {
    mhout = list(acpt = 0,lampst = currg,lgdetcor=currdet,CorMat=currcmt,CorPrc=currprc)
  }
  # Return acceptance, parameter value, determinant, correlation and precision matrices
  return(mhout)
}


spatial_mh2d_matern = function(currg,sptsd,currprc,currcmt,logrdet,datarr,nflds,nloc,
                               dstmt,sptnu,mhsd,mhchl,prsig,prlam) {
  # Spatial covariance Matern range and variance parameter updates
  # Priors: Half normal for std dev, lognormal for range
  #     currg     - Current range parameter
  #     sptsd     - Marginal standard deviation
  #     currprc   - Inverse of spatial correlation matrix (current)
  #     currcmt   - Correlation matrix (current)
  #     logrdet   - Log determinant for currprc, will be computed if NULL
  #     cursd     - Current process standard deviation
  #     datarr    - Data array [nloc,nflds]
  #     nflds     - Number of fields/replicates
  #     nloc      - Number of locations
  #     dstmt     - Distance matrix, assume a single set valid for all nflds
  #     sptnu     - Spatial smoothness parameter
  #     mhsd      - MH jumping standard deviation
  #     mhchl     - Cholesky lower triangle for MH proposal
  #     prsig     - Half-normal prior std deviation for sigma
  #     prlam     - Lognormal prior mean, std dev for range

  lgsdcur = log(sptsd)
  lgrgcur = log(currg)
  curlg = c(lgsdcur,lgrgcur)
  
  #print(curlg)
  lgsdprp = as.vector(curlg + mhchl %*% rnorm(2,sd = mhsd))
  #print(lgsdprp)
  prpsd = exp(lgsdprp[1])
  prprg = exp(lgsdprp[2])
  
  curvr = sptsd * sptsd
  prpvr = prpsd * prpsd
  
  # Proposal correlation structure
  prpR = Matern(dstmt,range=prprg,smoothness=sptnu)
  egprpR = eigen(prpR,only.values = TRUE)
  prprdet = sum(log(egprpR$values))
  prpRinv = solve(prpR)
  if (is.null(logrdet)) {
    egcurR = eigen(prpR,only.values = TRUE)
    currdet = sum(log(egcurR$values))
  }
  else {
    currdet = logrdet
  }
  
  # Prior density
  curpst = dnorm(lgrgcur,mean=prlam[1],sd=prlam[2],log = TRUE)
  prppst = dnorm(lgsdprp[2],mean=prlam[1],sd=prlam[2],log = TRUE)
  #curpst = curpst + dnorm(lgsdcur,mean=prlam[1],sd=prlam[2],log = TRUE)
  #prppst = prppst + dnorm(lgsdprp[1],mean=prlam[1],sd=prlam[2],log = TRUE)
  curpst = curpst - 0.5 * curvr / (prsig*prsig)
  prppst = prppst - 0.5 * prpvr / (prsig*prsig)
  

  # Likelihood contributions
  curpst = curpst - 0.5 * nflds * logrdet - nloc * nflds * log(curvr) / 2.0
  prppst = prppst - 0.5 * nflds * prprdet - nloc * nflds * log(prpvr) / 2.0
  if (nflds == 1) {
    dtdst = t(as.vector(datarr)) %*% currprc %*% (as.vector(datarr))
    dtprp = t(as.vector(datarr)) %*% prpRinv %*% (as.vector(datarr))
    curpst = curpst - 0.5 * dtdst / curvr 
    prppst = prppst - 0.5 * dtprp / prpvr 
  }
  else {
    for (i in seq(1,nflds)) {
      dtdst = t(datarr[,i]) %*% currprc %*% datarr[,i]
      dtprp = t(datarr[,i]) %*% prpRinv %*% datarr[,i]
      curpst = curpst - 0.5 * dtdst / curvr 
      prppst = prppst - 0.5 * dtprp / prpvr 
    }
  }

  # MH Asymmetry
  curpst = curpst + lgsdcur
  prppst = prppst + lgsdprp[1]
  
  # No MH Asymmetry due to lognormal prior
  #curpst = curpst + lgrgcur
  #prppst = prppst + lgsdprp
  
  # MH Evaluation
  logr = prppst - curpst
  logu = log(runif(1))
  if (logu[1] < logr[1]) {
    mhout = list(acpt = 1,sdpst=prpsd,lampst = prprg,lgdetcor=prprdet,CorMat=prpR,CorPrc=prpRinv)
  }
  else {
    mhout = list(acpt = 0,sdpst=sptsd,lampst = currg,lgdetcor=currdet,CorMat=currcmt,CorPrc=currprc)
  }
  # Return acceptance, parameter value, determinant, correlation and precision matrices
  return(mhout)
}


spatial_mh3d_matern = function(sptsd,currg,sptnu,currprc,currcmt,logrdet,datarr,nflds,nloc,
                               dstmt,mhsd,mhchl,prsig,prlam,prnu) {
  # Spatial covariance Matern range and variance parameter updates
  # Priors: Half normal for std dev, lognormal for range
  #     sptsd     - Current marginal standard deviation
  #     currg     - Current range parameter
  #     sptnu     - Current spatial smoothness parameter
  #     currprc   - Inverse of spatial correlation matrix (current)
  #     currcmt   - Correlation matrix (current)
  #     logrdet   - Log determinant for currprc, will be computed if NULL
  #     cursd     - Current process standard deviation
  #     datarr    - Data array [nloc,nflds]
  #     nflds     - Number of fields/replicates
  #     nloc      - Number of locations
  #     dstmt     - Distance matrix, assume a single set valid for all nflds
  #     mhsd      - MH jumping standard deviation
  #     mhchl     - Cholesky lower triangle for MH proposal
  #     prsig     - Half-normal prior std deviation for sigma
  #     prlam     - Lognormal prior mean, std dev for range
  #     prnu      - Lognormal prior mean, std dev for smoothness
  
  lgsdcur = log(sptsd)
  lgrgcur = log(currg)
  lgnucur = log(sptnu)
  curlg = c(lgsdcur,lgrgcur,lgnucur)
  
  lgsdprp = as.vector(curlg + mhchl %*% rnorm(3,sd = mhsd))
  #print(lgsdprp)
  prpsd = exp(lgsdprp[1])
  prprg = exp(lgsdprp[2])
  prpnu = exp(lgsdprp[3])
  
  curvr = sptsd * sptsd
  prpvr = prpsd * prpsd

  # Proposal correlation structure
  cvinvld = -1
  prpR = Matern(dstmt,range=prprg,smoothness=prpnu)
  diag(prpR) = 1.0
  if (length(prpR[is.na(prpR)]) > 0) {
    cvinvld = 1
  }
  else if ((min (prpR, na.rm = TRUE) > -1.01) && (max(prpR, na.rm = TRUE) <= 1.01) ) {
    egprpR = eigen(prpR,only.values = TRUE)
    #print(min(egprpR$values))
    if (min(egprpR$values) < 1e-8) { 
      cvinvld = 1
    }
    else {
      prprdet = sum(log(egprpR$values))
      prpRinv = solve(prpR)
    }
  }
  else {
    cvinvld = 1
  }

  if (is.null(logrdet)) {
    egcurR = eigen(currcmt,only.values = TRUE)
    currdet = sum(log(egcurR$values))
  }
  else {
    currdet = logrdet
  }
 
  if (cvinvld < 0) {  
    # Prior density
    curpst = dnorm(lgrgcur,mean=prlam[1],sd=prlam[2],log = TRUE)
    prppst = dnorm(lgsdprp[2],mean=prlam[1],sd=prlam[2],log = TRUE)
    curpst = curpst + dnorm(lgnucur,mean=prnu[1],sd=prnu[2],log = TRUE)
    prppst = prppst + dnorm(lgsdprp[3],mean=prnu[1],sd=prnu[2],log = TRUE)
    curpst = curpst - 0.5 * curvr / (prsig*prsig)
    prppst = prppst - 0.5 * prpvr / (prsig*prsig)
  
    # Likelihood contributions
    curpst = curpst - 0.5 * nflds * logrdet - nloc * nflds * log(curvr) / 2.0
    prppst = prppst - 0.5 * nflds * prprdet - nloc * nflds * log(prpvr) / 2.0
    if (nflds == 1) {
      dtdst = t(as.vector(datarr)) %*% currprc %*% (as.vector(datarr))
      dtprp = t(as.vector(datarr)) %*% prpRinv %*% (as.vector(datarr))
      curpst = curpst - 0.5 * dtdst[1] / curvr 
      prppst = prppst - 0.5 * dtprp[1] / prpvr 
    }
    else {
      for (i in seq(1,nflds)) {
        dtdst = t(datarr[,i]) %*% currprc %*% datarr[,i]
        dtprp = t(datarr[,i]) %*% prpRinv %*% datarr[,i]
        curpst = curpst - 0.5 * dtdst[1] / curvr 
        prppst = prppst - 0.5 * dtprp[1] / prpvr 
      }
    }
  }
  else {
    prppst = -1.0e-6
    curpst = 0.0
  }  

  # MH Asymmetry for std dev
  curpst = curpst + lgsdcur
  prppst = prppst + lgsdprp[1]
  
  # MH Evaluation
  logr = prppst - curpst
  logu = log(runif(1))
  if ((logu[1] < logr[1]) && (cvinvld < 0)) {
    mhout = list(acpt = 1,sdpst=prpsd,lampst = prprg,nupst = prpnu,lgdetcor=prprdet,CorMat=prpR,CorPrc=prpRinv)
  }
  else {
    mhout = list(acpt = 0,sdpst=sptsd,lampst = currg,nupst = sptnu,lgdetcor=currdet,CorMat=currcmt,CorPrc=currprc)
  }
  # Return acceptance, parameter value, determinant, correlation and precision matrices
  return(mhout)
}

rand_mh_rate = function(arate,jmpsd) {
  # Calibrate MH jumping distribution
  #     arate     - Current range parameter
  #     sptsd     - Marginal standard deviation
  
  if (arate < 0.25) {
    jmpmn = -0.5 * log(2.0)
  }
  else if (arate > 0.45) {
    jmpmn = 0.5 * log(2.0)
  }
  else {
    jmpmn = 0.0
  }
  jdev = jmpmn + rnorm(1,sd=0.1)
  jmpout = jmpsd*exp(jdev[1])
  return(jmpout)
}

fit_anova_arr3d = function(datmat, contrstA, contrstB) {
    # Quick ANOVA estimate
    # Data are organized as a 3D array (factA, factB, replicate)
    dimnames(datmat)[1] = list(contrstA)
    dimnames(datmat)[2] = list(contrstB)
    mltfrm = melt(datmat,varnames=c("FactA","FactB","Rep"))
    mltfrm$Interact = mltfrm$FactA * mltfrm$FactB
    l1 = lm(mltfrm$value ~ mltfrm$FactA + mltfrm$FactB + mltfrm$Interact)
    return(l1$coefficients)
}

fit_anova_mltlev_arr3d = function(datmat, contrstA, contrstB) {
    # Quick ANOVA estimate for multi-level factors
    # Data are organized as a 3D array (factA, factB, replicate)
  
    if (is.na(ncol(contrstA))) {
        contrstA = matrix(contrstA,ncol=1)
    }
    if (is.na(ncol(contrstB))) {
        contrstB = matrix(contrstB,ncol=1)
    }
   
    amlt = melt(contrstA, varnames=c("LevA","CtrstA"), value.name="CfA")
    afrm = dcast(amlt,LevA ~ CtrstA,value.var = "CfA")
    asq = seq(2,ncol(afrm))
    names(afrm)[asq] = paste("A",names(afrm)[asq],sep="_")
    bmlt = melt(contrstB, varnames=c("LevB","CtrstB"), value.name="CfB")
    bfrm = dcast(bmlt,LevB ~ CtrstB,value.var = "CfB")
    bsq = seq(2,ncol(bfrm))
    names(bfrm)[bsq] = paste("B",names(bfrm)[bsq],sep="_")

    ctrmrg = merge(amlt,bmlt,all = TRUE)
    ctrmrg$Interact = ctrmrg$CfA * ctrmrg$CfB

    ctrcst = dcast(ctrmrg,LevA + LevB ~ CtrstA + CtrstB, value.var="Interact")
    isq = seq(3,ncol(ctrcst))
    names(ctrcst)[isq] = paste("AB",names(ctrcst)[isq],sep="_") 

    mltfrm = melt(datmat,varnames=c("LevA","LevB","Rep"))
    mltmrg = merge(mltfrm,afrm)
    mltmrg = merge(mltmrg,bfrm)
    mltmrg = merge(mltmrg,ctrcst)

    csq = seq(5,ncol(mltmrg))
    l1 = lm(mltmrg$value ~ as.matrix(mltmrg[,csq])) 
    cfout = l1$coefficients
    names(cfout) = c("Intercept",names(mltmrg)[csq])
    return(cfout)
}

# Trace plots
gp_partrace_burn = function(vrfm,cflst,pltthm,nstrt=51) {
    # Trace plot for GP parameters (stddev, range, smoothness)
    # vrfrm:  is a data frame with plot information
    #    ofile:  Plot output file name
    #    vrnm:   GP ANOVA component
    #    nrow:   Number of rows for plot
    #    width:  Plot width
    #    height: Plot height
    # cflst:  main configuration list for the analysis
    # pltthm: a ggplot theme
    nchain = as.integer(cflst$nchain)
    niter = as.integer(cfglst$nburn)
    nthin = as.integer(cfglst$nthin)
    nsamp = floor(niter / nthin)

    
    tct = nsamp - nstrt + 1
    parr = array(0,c(tct,3,nchain))
    for (i in seq(1,nchain)) {
      ncsmpl = paste0(cfglst$burn_samp_file,i,".nc")
      nc1 = nc_open(ncsmpl)
      sgnm = paste0("gp_stddev_",vrfm$vrnm)
      parr[,1,i] = ncvar_get(nc1,sgnm,start=nstrt,count=tct)
      lmnm = paste0("gp_range_",vrfm$vrnm)
      parr[,2,i] = ncvar_get(nc1,lmnm,start=nstrt,count=tct)
      nunm = paste0("gp_smoothness_",vrfm$vrnm)
      parr[,3,i] = ncvar_get(nc1,nunm,start=nstrt,count=tct)
      nc_close(nc1)
    }
    dimnames(parr)[1] = list(seq(nstrt,nsamp))
    dimnames(parr)[2] = list(c("StdDev","Range","Smoothness"))
    pmlt = melt(parr,varnames=c("Iteration","Variable","Chain"))
    pmlt$Chain = as.factor(pmlt$Chain)
    
    tstr = paste("Burn-in Trace\n",vrfm$vrnm)
    gplt = ggplot(pmlt,aes(x=Iteration,y=value,group=Chain)) + 
           facet_wrap( ~ Variable,nrow=vrfm$nrow,scales="free_y") + 
           geom_line(aes(color=Chain),size=0.3) + pltthm + ggtitle(tstr)
    print(vrfm$ofile)
    png(vrfm$ofile,width=vrfm$width,height=vrfm$height,res=150)
    print(gplt)
    dev.off()
}

gp_partrace_post = function(vrfm,cflst,pltthm,nstrt=1) {
  # Trace plot for GP parameters (stddev, range, smoothness)
  # vrfrm:  is a data frame with plot information
  #    ofile:  Plot output file name
  #    vrnm:   GP ANOVA component
  #    nrow:   Number of rows for plot
  #    width:  Plot width
  #    height: Plot height
  # cflst:  main configuration list for the analysis
  # pltthm: a ggplot theme
  nchain = as.integer(cflst$nchain)
  niter = as.integer(cfglst$npost)
  nthin = as.integer(cfglst$nthin)
  nsamp = floor(niter / nthin)
  
  tct = nsamp - nstrt + 1
  parr = array(0,c(tct,3,nchain))
  for (i in seq(1,nchain)) {
    ncsmpl = paste0(cfglst$post_samp_file,i,".nc")
    nc1 = nc_open(ncsmpl)
    sgnm = paste0("gp_stddev_",vrfm$vrnm)
    parr[,1,i] = ncvar_get(nc1,sgnm,start=nstrt,count=tct)
    lmnm = paste0("gp_range_",vrfm$vrnm)
    parr[,2,i] = ncvar_get(nc1,lmnm,start=nstrt,count=tct)
    nunm = paste0("gp_smoothness_",vrfm$vrnm)
    parr[,3,i] = ncvar_get(nc1,nunm,start=nstrt,count=tct)
    nc_close(nc1)
  }
  dimnames(parr)[1] = list(seq(nstrt,nsamp))
  dimnames(parr)[2] = list(c("StdDev","Range","Smoothness"))
  pmlt = melt(parr,varnames=c("Iteration","Variable","Chain"))
  pmlt$Chain = as.factor(pmlt$Chain)
  
  tstr = paste("Posterior Trace\n",vrfm$vrnm)
  gplt = ggplot(pmlt,aes(x=Iteration,y=value,group=Chain)) + 
    facet_wrap( ~ Variable,nrow=vrfm$nrow,scales="free_y") + 
    geom_line(aes(color=Chain),size=0.3) + pltthm + ggtitle(tstr)
  print(vrfm$ofile)
  png(vrfm$ofile,width=vrfm$width,height=vrfm$height,res=150)
  print(gplt)
  dev.off()
}

# Trace plots
gp_fldtrace_burn = function(vrfm,cflst,pltthm,nstrt=51) {
  # Trace plot for functional ANOVA component realizations
  # vrfrm:  is a data frame with plot information
  #    ofile:  Plot output file name
  #    vrnm:   GP ANOVA component
  #    nrow:   Number of rows for plot
  #    width:  Plot width
  #    height: Plot height
  #    locst:  Starting location index
  #    locfn:  Ending location index
  #    locinc: Increment for location index
  #    ndigit: Number of digits for coorrinate info
  #    nctrst: Number of contrasts for ANOVA element
  # cflst:  main configuration list for the analysis
  # pltthm: a ggplot theme
  
  nchain = as.integer(cflst$nchain)
  niter = as.integer(cfglst$nburn)
  nthin = as.integer(cfglst$nthin)
  nsamp = floor(niter / nthin)
  
  # Location sequence
  sprstr = paste("%.",vrfm$ndigit,"f",sep="")
  lcsq = seq(vrfm$locst,vrfm$locfn,by=vrfm$locinc)
  nc1 = nc_open(cflst$data_file)
  locx = ncvar_get(nc1,cfglst$location_x_name)
  locy = ncvar_get(nc1,cfglst$location_y_name)
  nc_close(nc1)
  
  xstr = sprintf(sprstr,locx[lcsq])
  ystr = sprintf(sprstr,locy[lcsq])
  lcstrs = paste0("Location (",xstr,",",ystr,")")
  lcshrt = paste0("(",xstr,",",ystr,")")
  nlc = length(lcsq)
  
  tct = nsamp - nstrt + 1
  if (vrfm$nctrst <= 1) {
      parr = array(0,c(tct,nlc,nchain))
  }
  else {
      parr = array(0,c(tct,nlc,vrfm$nctrst,nchain))
  }
  for (i in seq(1,nchain)) {
    ncsmpl = paste0(cfglst$burn_samp_file,i,".nc")
    nc1 = nc_open(ncsmpl)
    fldnm = paste0("gp_field_",vrfm$vrnm)
    for (j in seq(1,nlc)) {
        if (vrfm$nctrst == 0) {
            parr[,j,i] = ncvar_get(nc1,fldnm,start=c(nstrt,lcsq[j]),count=c(tct,1))
        }
        else if (vrfm$nctrst == 1) {
            parr[,j,i] = ncvar_get(nc1,fldnm,start=c(nstrt,lcsq[j],1),count=c(tct,1,1))
        }
        else {
            parr[,j,,i] = ncvar_get(nc1,fldnm,start=c(nstrt,lcsq[j],1),count=c(tct,1,vrfm$nctrst))
        }
    }
    nc_close(nc1)
  }
  dimnames(parr)[1] = list(seq(nstrt,nsamp))
  
  print(dim(parr))
  if (vrfm$nctrst <= 1) {
    dimnames(parr)[2] = list(lcstrs)
    pmlt = melt(parr,varnames=c("Iteration","Location","Chain"))
    pmlt$Chain = as.factor(pmlt$Chain)
  
    tstr = paste("Burn-in Trace Field Locations\n",vrfm$vrnm)
    gplt = ggplot(pmlt,aes(x=Iteration,y=value,group=Chain)) + 
           facet_wrap( ~ Location,nrow=vrfm$nrow,scales="fixed") + 
           geom_line(aes(color=Chain),size=0.3) + pltthm + ggtitle(tstr)
    print(vrfm$ofile)
    png(vrfm$ofile,width=vrfm$width,height=vrfm$height,res=150)
    print(gplt)
    dev.off()
  }
  else {
    dimnames(parr)[2] = list(lcshrt)
    pmlt = melt(parr,varnames=c("Iteration","Location","Contrast","Chain"))
    pmlt$Chain = as.factor(pmlt$Chain)
    
    tstr = paste("Burn-in Trace Field Locations\n",vrfm$vrnm)
    gplt = ggplot(pmlt,aes(x=Iteration,y=value,group=Chain)) + 
      facet_grid(Location ~ Contrast,scales="fixed") + 
      geom_line(aes(color=Chain),size=0.3) + pltthm + ggtitle(tstr)
    print(vrfm$ofile)
    png(vrfm$ofile,width=vrfm$width,height=vrfm$height,res=150)
    print(gplt)
    dev.off()
  }
}

gp_fldtrace_post = function(vrfm,cflst,pltthm,nstrt = 1) {
  # Trace plot for functional ANOVA component realizations
  # vrfrm:  is a data frame with plot information
  #    ofile:  Plot output file name
  #    vrnm:   GP ANOVA component
  #    nrow:   Number of rows for plot
  #    width:  Plot width
  #    height: Plot height
  #    locst:  Starting location index
  #    locfn:  Ending location index
  #    locinc: Increment for location index
  #    ndigit: Number of digits for coorrinate info
  #    nctrst: Number of contrasts for ANOVA element
  # cflst:  main configuration list for the analysis
  # pltthm: a ggplot theme
  
  nchain = as.integer(cflst$nchain)
  niter = as.integer(cfglst$npost)
  nthin = as.integer(cfglst$nthin)
  nsamp = floor(niter / nthin)

  # Location sequence
  sprstr = paste("%.",vrfm$ndigit,"f",sep="")
  lcsq = seq(vrfm$locst,vrfm$locfn,by=vrfm$locinc)
  nc1 = nc_open(cflst$data_file)
  locx = ncvar_get(nc1,cfglst$location_x_name)
  locy = ncvar_get(nc1,cfglst$location_y_name)
  nc_close(nc1)
  
  xstr = sprintf(sprstr,locx[lcsq])
  ystr = sprintf(sprstr,locy[lcsq])
  lcstrs = paste0("Location (",xstr,",",ystr,")")
  lcshrt = paste0("(",xstr,",",ystr,")")
  nlc = length(lcsq)
  
  tct = nsamp - nstrt + 1
  if (vrfm$nctrst <= 1) {
    parr = array(0,c(tct,nlc,nchain))
  }
  else {
    parr = array(0,c(tct,nlc,vrfm$nctrst,nchain))
  }
  for (i in seq(1,nchain)) {
    ncsmpl = paste0(cfglst$post_samp_file,i,".nc")
    nc1 = nc_open(ncsmpl)
    fldnm = paste0("gp_field_",vrfm$vrnm)
    for (j in seq(1,nlc)) {
      if (vrfm$nctrst == 0) {
        parr[,j,i] = ncvar_get(nc1,fldnm,start=c(nstrt,lcsq[j]),count=c(tct,1))
      }
      else if (vrfm$nctrst == 1) {
        parr[,j,i] = ncvar_get(nc1,fldnm,start=c(nstrt,lcsq[j],1),count=c(tct,1,1))
      }
      else {
        parr[,j,,i] = ncvar_get(nc1,fldnm,start=c(nstrt,lcsq[j],1),count=c(tct,1,vrfm$nctrst))
      }
      
    }
    nc_close(nc1)
  }

  dimnames(parr)[1] = list(seq(nstrt,nsamp))
  if (vrfm$nctrst <= 1) {
    dimnames(parr)[2] = list(lcstrs)
    pmlt = melt(parr,varnames=c("Iteration","Location","Chain"))
    pmlt$Chain = as.factor(pmlt$Chain)
 
    tstr = paste("Posterior Trace Field Locations\n",vrfm$vrnm)   
    gplt = ggplot(pmlt,aes(x=Iteration,y=value,group=Chain)) + 
      facet_wrap( ~ Location,nrow=vrfm$nrow,scales="fixed") + 
      geom_line(aes(color=Chain),size=0.3) + pltthm + ggtitle(tstr)
    print(vrfm$ofile)
    png(vrfm$ofile,width=vrfm$width,height=vrfm$height,res=150)
    print(gplt)
    dev.off()
  }
  else {
    dimnames(parr)[2] = list(lcshrt)
    pmlt = melt(parr,varnames=c("Iteration","Location","Contrast","Chain"))
    pmlt$Chain = as.factor(pmlt$Chain)
    
    tstr = paste("Posterior Trace Field Locations\n",vrfm$vrnm)
    gplt = ggplot(pmlt,aes(x=Iteration,y=value,group=Chain)) + 
      facet_grid(Location ~ Contrast,scales="fixed") + 
      geom_line(aes(color=Chain),size=0.3) + pltthm + ggtitle(tstr)
    print(vrfm$ofile)
    png(vrfm$ofile,width=vrfm$width,height=vrfm$height,res=150)
    print(gplt)
    dev.off()
  }
}

ascatter = function(vrfm,cflst,pltthm) {
  # MH Acceptance rate for scalar parameter
  # vrfrm:  is a data frame with plot information
  #    ofile:  Plot output file name
  #    vrnm:   Parameter name
  #    width:  Plot width
  #    height: Plot height
  # cflst:  main configuration list for the analysis
  # pltthm: a ggplot theme

  nchain = as.integer(cflst$nchain)
  niter = as.integer(cfglst$nburn)
  nadapt = as.integer(cfglst$burn_mh_adapt)
  nmhsv = floor(niter / nadapt)

  acparr = array(0,c(nmhsv,2,4)) 
  vsd = paste0("mh_stddev_",vrfm$vrnm)
  vrt = paste0("mh_acc_rate_",vrfm$vrnm)
  for (i in seq(1,nchain)) {
    ncsmpl = paste0(cfglst$burn_samp_file,i,".nc")
    nc1 = nc_open(ncsmpl)
    acparr[,1,i] = ncvar_get(nc1,vsd,start=c(1),
                             count=c(nmhsv))
    acparr[,2,i] = ncvar_get(nc1,vrt,start=c(1),
                             count=c(nmhsv))
    nc_close(nc1)
  }
  
  dimnames(acparr)[2] = list(c("StdDev","Rate"))
  amlt = melt(acparr,varnames=c("Iteration","Variable","Chain"))
  amlt$Chain = as.factor(amlt$Chain)
  acst = dcast(amlt,Iteration + Chain ~ Variable)
  tstr = paste("MH Acceptance\n",vrfm$vrnm,sep="")
  aplt = ggplot(acst,aes(x=StdDev,y=Rate,group=Chain)) + 
    geom_point(aes(shape=Chain)) + 
    scale_shape_manual("Chain",values = 0:3) + 
    pltthm + ggtitle(tstr)
  png(paste(vrfm$ofile),width=vrfm$width,height=vrfm$height,res=150)
  print(aplt)
  dev.off()
}

aratebox = function(vrfm,cflst,pltthm) {
  # MH Acceptance rate for scalar parameter
  # vrfrm:  is a data frame with plot information
  #    ofile:  Plot output file name
  #    vrnm:   Parameter name
  #    width:  Plot width
  #    height: Plot height
  # cflst:  main configuration list for the analysis
  # pltthm: a ggplot theme
  
  nchain = as.integer(cflst$nchain)
  niter = as.integer(cfglst$npost)
  nadapt = as.integer(cfglst$post_mh_adapt)
  nmhsv = floor(niter / nadapt)
  
  acparr = array(0,c(nmhsv,4)) 
  vrt = paste0("mh_acc_rate_",vrfm$vrnm)
  for (i in seq(1,nchain)) {
    ncsmpl = paste0(cfglst$post_samp_file,i,".nc")
    nc1 = nc_open(ncsmpl)
    acparr[,i] = ncvar_get(nc1,vrt,start=c(1),
                             count=c(nmhsv))
    nc_close(nc1)
  }
  
  amlt = melt(acparr,varnames=c("Iteration","Chain"),value.name = "Rate")
  amlt$Chain = as.factor(amlt$Chain)
  tstr = paste("MH Acceptance\n",vrfm$vrnm,sep="")
  aplt = ggplot(amlt,aes(x=Chain,y=Rate)) + 
    geom_boxplot() + 
    pltthm + ggtitle(tstr)
  png(paste(vrfm$ofile),width=vrfm$width,height=vrfm$height,res=150)
  print(aplt)
  dev.off()
}

gp_parcov_burn = function(vrnm,cflst,nstrt=51,logpar=TRUE) {
  # Empirical covariances for GP parameters (stddev, range, smoothness) - on log scale
  # vrnm:   GP ANOVA component
  # cflst:  main configuration list for the analysis
  # nstrt:  starting iteration
  # logpar: compute log transformation
  
  nchain = as.integer(cflst$nchain)
  niter = as.integer(cfglst$nburn)
  nthin = as.integer(cfglst$nthin)
  nsamp = floor(niter / nthin)
  
  tct = nsamp - nstrt + 1
  parr = array(0,c(tct,3,nchain))
  for (i in seq(1,nchain)) {
    ncsmpl = paste0(cfglst$burn_samp_file,i,".nc")
    nc1 = nc_open(ncsmpl)
    sgnm = paste0("gp_stddev_",vrnm)
    parr[,1,i] = ncvar_get(nc1,sgnm,start=nstrt,count=tct)
    lmnm = paste0("gp_range_",vrnm)
    parr[,2,i] = ncvar_get(nc1,lmnm,start=nstrt,count=tct)
    nunm = paste0("gp_smoothness_",vrnm)
    parr[,3,i] = ncvar_get(nc1,nunm,start=nstrt,count=tct)
    nc_close(nc1)
  }
  if (logpar) {
    parout = log(parr)
  }
  else {
    parout = parr
  }
  return(parout)
}

gp_parcov_post = function(vrnm,cflst,nstrt=1,logpar=TRUE) {
  # Empirical covariances for GP parameters (stddev, range, smoothness) - on log scale
  # vrnm:   GP ANOVA component
  # cflst:  main configuration list for the analysis
  # nstrt:  starting iteration
  # logpar: compute log transformation
  
  nchain = as.integer(cflst$nchain)
  niter = as.integer(cfglst$npost)
  nthin = as.integer(cfglst$nthin)
  nsamp = floor(niter / nthin)
  
  tct = nsamp - nstrt + 1
  parr = array(0,c(tct,3,nchain))
  for (i in seq(1,nchain)) {
    ncsmpl = paste0(cfglst$post_samp_file,i,".nc")
    nc1 = nc_open(ncsmpl)
    sgnm = paste0("gp_stddev_",vrnm)
    parr[,1,i] = ncvar_get(nc1,sgnm,start=nstrt,count=tct)
    lmnm = paste0("gp_range_",vrnm)
    parr[,2,i] = ncvar_get(nc1,lmnm,start=nstrt,count=tct)
    nunm = paste0("gp_smoothness_",vrnm)
    parr[,3,i] = ncvar_get(nc1,nunm,start=nstrt,count=tct)
    nc_close(nc1)
  }
  if (logpar) {
    parout = log(parr)
  }
  else {
    parout = parr
  }
  return(parout)
}

cov_to_frm = function(cvmt, vrnms = NULL) {
  if (is.null(vrnms)) {
    nmsq = paste("V",seq(1,nrow(cvmt)),sep="")
  }
  else {
    nmsq = vrnms
  }
  
  sdvc = sqrt(diag(cvmt))
  names(sdvc) = paste("SD",nmsq,sep="_")
  sfrm = as.data.frame(t(sdvc))
  
  crmt = cov2cor(cvmt)
  dimnames(crmt)[1] = list(nmsq)
  dimnames(crmt)[2] = list(nmsq)
  crmlt = melt(crmt,varnames=c("Row","Col"),value.name = "Correlation")
  
  dgmt = lower.tri(crmt)
  dimnames(dgmt)[1] = list(nmsq)
  dimnames(dgmt)[2] = list(nmsq)
  dgmlt = melt(dgmt,varnames=c("Row","Col"),value.name = "LTri")
  
  cmrg = merge(crmlt,dgmlt)
  cmsb = cmrg[cmrg$LTri,]
  
  cmsb$CorCmb = paste("Corr",cmsb$Row,cmsb$Col,sep="_")
  crvc = cmsb$Correlation
  names(crvc) = cmsb$CorCmb
  
  frmout = cbind(sfrm,t(crvc))
  return(frmout)
}

gp_scale_helmert = function(ngrp) {
  # Set up matrix of scaled Helmert contrasts, according to Kaufman and Sain
  cmthlm = contr.helmert(ngrp)
  
  crw = ngrp * (ngrp + 1) / 2
  ccl = ngrp - 1
  
  cmt = matrix(0,nrow=crw,ncol = ccl)
  rhs = rep(0,crw)
  
  tct = 0
  for (p in seq(1,ngrp)) {
    for (q in seq(p,ngrp)) {
      if (p == q) {
        rhs[tct+1] = 1.0 - 1.0 / ngrp
        cmt[tct+1,] = cmthlm[p,]^2
      }
      else {
        rhs[tct+1] = -1.0 / ngrp
        cmt[tct+1,] = cmthlm[p,] * cmthlm[q,]
      }
      tct = tct + 1
    }
  }
  sclfs = sqrt(qr.solve(cmt,rhs))
  
  # Re-scale matrix
  sclhmt = matrix(0,nrow=nrow(cmthlm),ncol=ncol(cmthlm))
  for (j in seq(1,ccl)) {
    sclhmt[,j] = sclfs[j] * cmthlm[,j]
  }
  return(sclhmt)
}

q025 = function(dat) {
  qout = quantile(dat,0.025)
}

q975 = function(dat) {
  qout = quantile(dat,0.975)
}

prob2side = function(dat) {
  p1 = 2*length(dat[dat<0]) / length(dat)
  p2 = 2*length(dat[dat>0]) / length(dat)
  pout = min(c(p1,p2))
  return(pout)
}

pstsmry = function(d1) {
  #dfrmout = data.frame(mean=mean(d1$value),sd=sd(d1$value),
  #                     q025=q025(d1$value),q975=q975(d1$value),
  #                     prob0=prob2side(d1$value))
  dfrmout = data.frame(mean=mean(d1$value),sd=sd(d1$value),
                       q025=q025(d1$value),med=median(d1$value),q975=q975(d1$value))
  return(dfrmout)
}


