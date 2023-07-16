library(fields)
library(GPvecchia)
library(Matrix)

#####  functions for sampling using vecchia approximation

## posterior sampling for multiple levels using vecchia
samp_post_vecchia_mltlev=function(y.tilde,cpars,W.eps,vecchia.specify,nrep,ntrt,nloc) {
  # Update a functional ANOVA GP via Vecchia approximation
  # y.tilde:          Vector/matrix of marginal sums of data or suitable contrasts
  # cpars:            (Matern) covariance parameters: variance, range, smoothness
  # W.eps:            ANOVA error GP precision matrix
  # vecchia.specify:  Vecchia spcification object
  # nrep:             Number of replicates per treatment
  # ntrt:             Number of treatement combinations (na x nb if complete design)
  # nloc:             Number of locations (leading dimension of y.tilde)
  
  # posterior cholesky via vecchia+IC0
  U=createU(vecchia.specify,covparms=cpars,nuggets=0)$U # prior chol
  W=Matrix::tcrossprod(U) # prior precision
  W.rev=revMat(W+(ntrt*nrep)*W.eps) # posterior precision
  V.ord=Matrix::t(ichol(W.rev)) 

  lcrv = seq(nloc,1,-1)
  nctrst = ncol(y.tilde)

  # sample from posterior, allow for more than 1 right-hand side
  #y.eps=as.numeric(W.eps%*%y.tilde[vecchia.specify$ord,])
  y.eps=W.eps%*%y.tilde[vecchia.specify$ord,]
  temp=Matrix::solve(V.ord,y.eps[lcrv,])
  zsmp = matrix(rnorm(nloc*nctrst), ncol=nctrst)
  post.samp.revord=Matrix::solve(Matrix::t(V.ord),temp+zsmp)
  ##post.samp.revord=Matrix::solve(Matrix::t(V.ord),temp)  # Post mean only
  post.samp=post.samp.revord[lcrv,]
  post.samp[vecchia.specify$ord,]=post.samp
  
  return(post.samp)
}

samp_post_vecchia_mltlev_sptm=function(y.tilde,cpars,W.eps,vecchia.specify,nrep,ntrt,nloc,tprcfct,sptmdsgn,W.epsptm) {
  # Update a functional ANOVA GP with multiple factors and time autocorrelation
  # y.tilde:          Vector/matrix of marginal sums of data or suitable contrasts
  # cpars:            (Matern) covariance parameters: variance, range, smoothness
  # W.eps:            ANOVA error GP precision matrix, spatial only, dimension (nloc,nloc)
  # vecchia.specify:  Vecchia spcification object
  # nrep:             Number of replicates per treatment (time points)
  # ntrt:             Number of treatement combinations (na x nb if complete design)
  # nloc:             Number of locations (leading dimension of y.tilde)
  # tprcfct:          Temporal precision factor (effective sample size for posterior precision)
  # sptmdsgn:         Transpose of Sparse space-time design matrix (nrep x nloc, nrep)  
  # W.epsptm          Error precision (space-time), dimension (nloc*nrep, nloc*nrep)
  
  # posterior cholesky via vecchia+IC0
  U=createU(vecchia.specify,covparms=cpars,nuggets=0)$U # prior chol
  W=Matrix::tcrossprod(U) # prior precision
  W.rev=revMat(W+(ntrt*tprcfct)*W.eps) # posterior precision
  V.ord=Matrix::t(ichol(W.rev)) 

  lcrv = seq(nloc,1,-1)
  #tlcsq = NULL
  #for (k in 1:nrep) {
  #  tlcsq = c(tlcsq,(k-1)*nloc+lcrv)
  #}
  tlcsq = seq(nrep*nloc,1,-1)
  nctrst = dim(y.tilde)[2]

  # Need to reshape
  ystk = array(0,c(nloc*nrep,nctrst))
  for (k in 1:nrep) {
    ssq = seq( (k-1)*nloc+1, k*nloc)
    ystk[ssq,] = y.tilde[vecchia.specify$ord,,k]
  }

  # sample from posterior, allow for more than 1 right-hand side
  #y.eps=as.numeric(W.eps%*%y.tilde[vecchia.specify$ord,])
  y.eps = sptmdsgn %*% W.epsptm %*% ystk
  temp=Matrix::solve(V.ord,y.eps[lcrv,])   # Need to check this
  zsmp = matrix(rnorm(nloc*nctrst), ncol=nctrst)
  post.samp.revord=Matrix::solve(Matrix::t(V.ord),temp+zsmp)
  ##post.samp.revord=Matrix::solve(Matrix::t(V.ord),temp)  # Post mean only
  post.samp=post.samp.revord[lcrv,]
  post.samp[vecchia.specify$ord,]=post.samp
  return(post.samp)
}

## function for sampling from posteriors using vecchia
samp_post_vecchia=function(y.tilde,cpars,W.eps,vecchia.specify, nrep,ntrt,nloc){
  # nrep:             Number of replicates per treatment
  # ntrt:             Number of treatement combinations (na x nb if complete design)
  # nloc:             Number of locations (leading dimension of y.tilde)
  
  # posterior cholesky via vecchia+IC0
  U=createU(vecchia.specify,covparms=cpars,nuggets=0)$U # prior chol
  W=Matrix::tcrossprod(U) # prior precision
  W.rev=revMat(W+(ntrt*nrep)*W.eps) # posterior precision
  V.ord=Matrix::t(ichol(W.rev)) 
  
  # sample from posterior
  y.eps=as.numeric(W.eps%*%y.tilde[vecchia.specify$ord])
  temp=Matrix::solve(V.ord,rev(y.eps))
  post.samp.revord=Matrix::solve(Matrix::t(V.ord),temp+rnorm(nloc))
  ##post.samp.revord=Matrix::solve(Matrix::t(V.ord),temp) # Post mean only
  post.samp=rev(post.samp.revord)
  post.samp[vecchia.specify$ord]=post.samp
  
  return(post.samp)
  
}

## function for sampling from posteriors using vecchia with time autocorrelation
samp_post_vecchia_sptm=function(y.tilde,cpars,W.eps,vecchia.specify, nrep,ntrt,nloc,tprcfct,sptmdsgn,W.epsptm){
  # y.tilde           Row/column sums from data, dimension (nloc,nrep)
  # cpars             Covariance params for ANOVA component, for Vecchia setup
  # W.eps             Error precision (spatial only), dimension (nloc,nloc)
  # nrep:             Number of replicates per treatment (time points)
  # ntrt:             Number of treatement combinations (na x nb if complete design)
  # nloc:             Number of locations (leading dimension of y.tilde)
  # tprcfct:          Temporal precision factor (effective sample size for posterior precision)
  # sptmdsgn:         Transpose of Sparse space-time design matrix (nrep x nloc, nrep)  
  # W.epsptm          Error precision (space-time), dimension (nloc*nrep, nloc*nrep)

  # posterior cholesky via vecchia+IC0
  U=createU(vecchia.specify,covparms=cpars,nuggets=0)$U # prior chol
  W=Matrix::tcrossprod(U) # prior precision
  W.rev=revMat(W+(ntrt*tprcfct)*W.eps) # posterior precision
  V.ord=Matrix::t(ichol(W.rev)) 
  
  # sample from posterior
  if (length(dim(y.tilde)) == 2) {
    vcy = y.tilde[vecchia.specify$ord,]
  } else {
    vcy = y.tilde[vecchia.specify$ord,,]
  }
  y.eps=as.numeric(sptmdsgn %*% W.epsptm %*% as.vector(vcy))
  temp=Matrix::solve(V.ord,rev(y.eps))
  post.samp.revord=Matrix::solve(Matrix::t(V.ord),temp+rnorm(nloc))
  ##post.samp.revord=Matrix::solve(Matrix::t(V.ord),temp) # Post mean only
  post.samp=rev(post.samp.revord)
  post.samp[vecchia.specify$ord]=post.samp
  
  return(post.samp)
  
}


##  incomplete cholesky decomposition (IC0)
ichol = function(M, S=NULL){
  if(!is(M, "sparseMatrix")){
    warning("Passing a dense matrix")
  }
  if(!is(M, "CsparseMatrix") || !Matrix::isTriangular(M)){
    M = as(Matrix::triu(M), "CsparseMatrix")
  }
  if(!is.null(S)){
    if(!is(S, "sparseMatrix")) S=as(Matrix::triu(S), "CsparseMatrix")
    p=S@p; i=S@i
  } else {
    p=M@p; i=M@i
  }
  vals = GPvecchia::ic0(p, i, M@x)
  Msp = Matrix::sparseMatrix(i=i, p=p, x=vals, index1=FALSE)
  Msp@x = vals
  return(as(Msp, 'dtCMatrix'))
}


## function to reverse-order a matrix
revMat=function(mat) mat[nrow(mat):1,ncol(mat):1,drop=FALSE]


spatial_mh3d_vecchia = function(sptsd,currg,sptnu, vecchia.specify, datarr,nflds,nloc,
                                mhsd,mhchl,prsig,prlam,prnu, prnumod='lognormal', sampnu=TRUE) {
  # Spatial covariance Matern range, smoothness, and variance parameter updates
  # Priors: Half normal for std dev, lognormal for range, smoothness
  # Likelihood based on Vecchia approximation
  #     sptsd     - Current marginal standard deviation
  #     currg     - Current range parameter
  #     sptnu     - Current spatial smoothness parameter
  #     vecchia.specify - Vecchia specification object
  #     datarr    - Data array [nloc,nflds]
  #     nflds     - Number of fields/replicates
  #     nloc      - Number of locations
  #     mhsd      - MH jumping standard deviation
  #     mhchl     - Cholesky lower triangle for MH proposal
  #     prsig     - Half-normal prior std deviation for sigma
  #     prlam     - Lognormal prior mean, std dev for range
  #     prnu      - Prior parameters for smoothness (2 for lognormal, 4 for beta) 
  #     prnumod   - Prior form for smoothness parameter
  #     sampnu    - Indicator for update on smoothness (FALSE = fixed)
 
  lgsdcur = log(sptsd)
  lgrgcur = log(currg)
  if (prnumod == 'lognormal') {
    lgnucur = log(sptnu)
  } else if (prnumod == 'beta') {
    nustr = (sptnu - prnu[3]) / (prnu[4] - prnu[3])
    lgnucur = log(nustr / (1.0 - nustr))
  }
  curlg = c(lgsdcur,lgrgcur,lgnucur)
 
  lgsdprp = as.vector(curlg + mhchl %*% rnorm(3,sd = mhsd))
  #print(lgsdprp)
  prpsd = exp(lgsdprp[1])
  prprg = exp(lgsdprp[2])
  if (sampnu) {
    if (prnumod == 'lognormal') {
      prpnu = exp(lgsdprp[3])
    } else if (prnumod == 'beta') {
      prpstr = exp(lgsdprp[3]) / (1.0 + exp(lgsdprp[3]))
      prpnu = prnu[3] + (prnu[4] - prnu[3]) * prpstr
    }
  }
  else {
    prpnu = sptnu
    if (prnumod == 'beta') {
      prpstr = (sptnu - prnu[3]) / (prnu[4] - prnu[3]) 
    }
  }  
  curvr = sptsd * sptsd
  prpvr = prpsd * prpsd

  # Vecchia cov params
  crcvvc = c(curvr,currg,sptnu)
  prcvvc = c(prpvr,prprg,prpnu)

  cvinvld = -1
  if (cvinvld < 0) {  
    # Prior density
    curpst = dnorm(lgrgcur,mean=prlam[1],sd=prlam[2],log = TRUE)
    prppst = dnorm(lgsdprp[2],mean=prlam[1],sd=prlam[2],log = TRUE)
    if (prnumod == 'lognormal') {
      curpst = curpst + dnorm(lgnucur,mean=prnu[1],sd=prnu[2],log = TRUE)
      prppst = prppst + dnorm(lgsdprp[3],mean=prnu[1],sd=prnu[2],log = TRUE)
    } else if (prnumod == 'beta') {
      curstr = (sptnu - prnu[3]) / (prnu[4] - prnu[3]) 
      curpst = curpst + dbeta(curstr, shape1=prnu[1], shape2=prnu[2], log=TRUE)
      prppst = prppst + dbeta(prpstr, shape1=prnu[1], shape2=prnu[2], log=TRUE)
    }
    curpst = curpst - 0.5 * curvr / (prsig*prsig)
    prppst = prppst - 0.5 * prpvr / (prsig*prsig)
  
    # Likelihood contributions
    if (nflds == 1) {
      curpst = curpst + vecchia_like_array(datarr, vecchia.specify, crcvvc, nrep=1, nuggets = 0)  
      prppst = prppst + vecchia_like_array(datarr, vecchia.specify, prcvvc, nrep=1, nuggets = 0)  
    }
    else {
      curpst = curpst + vecchia_like_array(datarr, vecchia.specify, crcvvc, nrep=nflds, nuggets=0)
      prppst = prppst + vecchia_like_array(datarr, vecchia.specify, prcvvc, nrep=nflds, nuggets=0)
    }
  }
  else {
    prppst = -1.0e-6
    curpst = 0.0
  }  

  #str1 = sprintf("  MH Vecchia Like:  Prp %.4e, Cur %.4e ",prppst,curpst)
  #print(str1) 
  # MH Asymmetry for std dev
  curpst = curpst + lgsdcur
  prppst = prppst + lgsdprp[1]
  if (prnumod == 'beta') {
    curpst = curpst + log(curstr) + log(1.0 - curstr)
    prppst = prppst + log(prpstr) + log(1.0 - prpstr)
  }

  if (is.na(prppst)) {
    prppst = -1.0e-6
    curpst = 0.0
  }

  # MH Evaluation
  logr = prppst - curpst
  logu = log(runif(1))
 
  # logr = -Inf       # Do not accept
  if ((logu[1] < logr[1]) && (cvinvld < 0)) {
    mhout = list(acpt = 1,sdpst=prpsd,lampst = prprg,nupst = prpnu)
  }
  else {
    mhout = list(acpt = 0,sdpst=sptsd,lampst = currg,nupst = sptnu)
  }
  # Return acceptance, parameter value, determinant, correlation and precision matrices
  return(mhout)
}

vecchia_like_array = function(sptarr, vecchia.specify, covpars, nrep, nuggets=0, covmodel='matern') { 

  # Evaluate Vecchia likelihood efficiently for an array of fields 
  # Only perform specification once 
  # Likelihood based on Vecchia approximation
  #     sptarr    - Data array (nloc,nrep)
  #     mhchl     - Cholesky lower triangle for MH proposal
  #     prsig     - Half-normal prior std deviation for sigma
  #     prlam     - Lognormal prior mean, std dev for range
  #     prnu      - Lognormal prior mean, std dev for smoothness
  #     sampnu    - Indicator for update on smoothness (FALSE = fixed)
 

  if(vecchia.specify$cond.yz=='zy')
    warning("cond.yz='zy' will produce a poor likelihood approximation. Use 'SGV' instead.")

  # remove NAs in data and U
  #removeNAs()
    
  # create the U matrix
  U.obj=createU(vecchia.specify,covpars,nuggets,covmodel)

  # compute the loglikelihood
  if (nrep == 1) { 
    lksum = vecchia_likelihood_U(sptarr,U.obj)
  }
  else {
    lkrep = rep(0,nrep)
    for (i in seq(1,nrep)) {
      lkrep[i] = vecchia_likelihood_U(sptarr[,i],U.obj)
    }
    lksum = sum(lkrep)
  }
  return(lksum)
}

removeNAs=function(){ # overwrites z and U.obj
    p = parent.frame()
    if(any(is.na(p$z))){

        if(length(p$nuggets)<length(p$z)) {
            new.nuggets = rep(0, length(p$z))
            new.nuggets[!is.na(p$z)] = p$nuggets
            p$nuggets = new.nuggets
        }
        
        p$nuggets[is.na(p$z)] = stats::var(p$z,na.rm=TRUE)*1e8
        p$z[is.na(p$z)] = mean(p$z,na.rm=TRUE)
    }
}

## evaluate vecchia likelihood based on U

vecchia_likelihood_U=function(z,U.obj) {

  ### output: loglikelihood (for z)
  U=U.obj$U
  latent=U.obj$latent
  zord=z[U.obj$ord.z]

  # constant
  const=sum(!latent)*log(2*pi)

  # numerator
  z1=Matrix::crossprod(U[!latent,],zord)
  quadform.num=sum(z1^2)
  logdet.num=-2*sum(log(Matrix::diag(U)))

  # denominator
  if(sum(latent)==0){ # no latents -> denominator not needed

    logdet.denom=quadform.denom=0

  } else {  # if latents, need denominator

    U.y=U[latent,]
    z2=as.numeric(U.y%*%z1)
    V.ord=U2V(U.obj)
    z3=Matrix::solve(V.ord,rev(z2),system='L')
    quadform.denom=sum(z3^2)
    logdet.denom=-2*sum(log(Matrix::diag(V.ord)))

  }

  # putting everything together
  neg2loglik=logdet.num-logdet.denom+quadform.num-quadform.denom+const
  loglik=-neg2loglik/2
  return(loglik)

}

sptm_mh4d_vecchia = function(sptsd,currg,sptnu,tmcor, vecchia.specify, datarr,nflds,ntm,nloc,
                                mhsd,mhchl,prsig,prlam,prnu, prar, prnumod='beta', prarmod='beta', sampnu=TRUE) {
  # Spatio-tempral updates for Matern range, smoothness, and variance parameters, and AR1 correlation
  # Priors: Half normal for std dev, lognormal for range, beta for smoothness, AR
  # Likelihood based on Vecchia approximation
  #     sptsd     - Current marginal standard deviation
  #     currg     - Current range parameter
  #     sptnu     - Current spatial smoothness parameter
  #     tmcor     - Current temporal autocorrelation
  #     vecchia.specify - Vecchia specification object
  #     datarr    - Data array [nloc,ntm,nflds]
  #     nflds     - Number of fields/replicates
  #     ntm       - Number of time points
  #     nloc      - Number of locations
  #     mhsd      - MH jumping standard deviation
  #     mhchl     - Cholesky lower triangle for MH proposal
  #     prsig     - Half-normal prior std deviation for sigma
  #     prlam     - Lognormal prior mean, std dev for range
  #     prnu      - Prior parameters for smoothness (2 for lognormal, 4 for beta) 
  #     prar      - Prior parameters for AR (4 for beta) 
  #     prnumod   - Prior form for smoothness parameter
  #     prarmod   - Prior form for smoothness parameter
  #     sampnu    - Indicator for update on smoothness (FALSE = fixed)

  lgsdcur = log(sptsd)
  lgrgcur = log(currg)
  if (prnumod == 'lognormal') {
    lgnucur = log(sptnu)
  } else if (prnumod == 'beta') {
    nustr = (sptnu - prnu[3]) / (prnu[4] - prnu[3])
    lgnucur = log(nustr / (1.0 - nustr))
  }
  if (prarmod == 'beta') {
    rhostr = (tmcor - prar[3]) / (prar[4] - prar[3])
    lgrhocur = log(rhostr / (1.0 - rhostr))
  }
  curlg = c(lgsdcur,lgrgcur,lgnucur,lgrhocur)
 
  lgsdprp = as.vector(curlg + mhchl %*% rnorm(4,sd = mhsd))
  #print(lgsdprp)
  prpsd = exp(lgsdprp[1])
  prprg = exp(lgsdprp[2])
  if (sampnu) {
    if (prnumod == 'lognormal') {
      prpnu = exp(lgsdprp[3])
    } else if (prnumod == 'beta') {
      prpstr = exp(lgsdprp[3]) / (1.0 + exp(lgsdprp[3]))
      prpnu = prnu[3] + (prnu[4] - prnu[3]) * prpstr
    }
  }
  else {
    prpnu = sptnu
    if (prnumod == 'beta') {
      prpstr = (sptnu - prnu[3]) / (prnu[4] - prnu[3]) 
    }
  }  
  if (prarmod == 'beta') {
    #prpast = rhostr 
    prpast = exp(lgsdprp[4]) / (1.0 + exp(lgsdprp[4]))
    prpar = prar[3] + (prar[4] - prar[3]) * prpast
  }
  curvr = sptsd * sptsd
  prpvr = prpsd * prpsd

  # Vecchia cov params
  crcvvc = c(curvr,currg,sptnu)
  prcvvc = c(prpvr,prprg,prpnu)

  cvinvld = -1
  if (cvinvld < 0) {  
    # Prior density
    curpst = dnorm(lgrgcur,mean=prlam[1],sd=prlam[2],log = TRUE)
    prppst = dnorm(lgsdprp[2],mean=prlam[1],sd=prlam[2],log = TRUE)
    if (prnumod == 'lognormal') {
      curpst = curpst + dnorm(lgnucur,mean=prnu[1],sd=prnu[2],log = TRUE)
      prppst = prppst + dnorm(lgsdprp[3],mean=prnu[1],sd=prnu[2],log = TRUE)
    } else if (prnumod == 'beta') {
      curstr = (sptnu - prnu[3]) / (prnu[4] - prnu[3]) 
      curpst = curpst + dbeta(curstr, shape1=prnu[1], shape2=prnu[2], log=TRUE)
      prppst = prppst + dbeta(prpstr, shape1=prnu[1], shape2=prnu[2], log=TRUE)
    }
    if (prnumod == 'beta') {
      curast = (tmcor - prar[3]) / (prar[4] - prar[3]) 
      curpst = curpst + dbeta(curast, shape1=prar[1], shape2=prar[2], log=TRUE)
      prppst = prppst + dbeta(prpast, shape1=prar[1], shape2=prar[2], log=TRUE)
    }
    curpst = curpst - 0.5 * curvr / (prsig*prsig)
    prppst = prppst - 0.5 * prpvr / (prsig*prsig)
  
    # Likelihood contributions
    curpst = curpst + vecchia_like_array_sptm(datarr, vecchia.specify, crcvvc, tmcor, ntm, nflds, nuggets=0)
    prppst = prppst + vecchia_like_array_sptm(datarr, vecchia.specify, prcvvc, prpar, ntm, nflds, nuggets=0)
  }
  else {
    prppst = -1.0e-6
    curpst = 0.0
  }  

  #str1 = sprintf("  MH Vecchia Like:  Prp %.4e, Cur %.4e ",prppst,curpst)
  #print(str1) 
  # MH Asymmetry for std dev
  curpst = curpst + lgsdcur
  prppst = prppst + lgsdprp[1]
  if (prnumod == 'beta') {
    curpst = curpst + log(curstr) + log(1.0 - curstr)
    prppst = prppst + log(prpstr) + log(1.0 - prpstr)
  }
  if (prarmod == 'beta') {
    curpst = curpst + log(curast) + log(1.0 - curast)
    prppst = prppst + log(prpast) + log(1.0 - prpast)
  }

  if (is.na(prppst)) {
    prppst = -1.0e-6
    curpst = 0.0
  }

  # MH Evaluation
  logr = prppst - curpst
  logu = log(runif(1))
 
  # logr = -Inf       # Do not accept
  if ((logu[1] < logr[1]) && (cvinvld < 0)) {
    mhout = list(acpt = 1,sdpst=prpsd,lampst = prprg,nupst = prpnu,rhopst = prpar)
  }
  else {
    mhout = list(acpt = 0,sdpst=sptsd,lampst = currg,nupst = sptnu,rhopst = tmcor)
  }
  # Return acceptance, parameter value, determinant, correlation and precision matrices
  return(mhout)
}

tmqf = function(zvc, prcmt) {
    # Quick calculation of temporal quadratic form
    qfrm = Matrix::t(zvc) %*% prcmt %*% zvc
    return(as.vector(qfrm))
}

vecchia_like_array_sptm = function(sptarr, vecchia.specify, covpars, tmrho, ntm, nfld, nuggets=0, covmodel='matern') { 

  # Evaluate spatio-temporal Vecchia likelihood efficiently for an array of fields 
  # Only perform specification once 
  # Likelihood based on Vecchia approximation
  #     sptarr    - Data array (nloc,ntm,nfld)
  #     vecchia.specify - GPvecchia specification object 
  #     covpars   - Matern covariance parameter vector (variance, range, smoothness) 
  #     tmrho     - Temporal autocorrelation
  #     ntm       - Number of time points
  #     nfld      - Number of field 'replicates' 

  if(vecchia.specify$cond.yz=='zy')
    warning("cond.yz='zy' will produce a poor likelihood approximation. Use 'SGV' instead.")

  # remove NAs in data and U
  #removeNAs()
  nloc = dim(sptarr)[1]    

  # create the U matrix
  U.obj=createU(vecchia.specify,covpars,nuggets,covmodel)
  lgu = log(Matrix::diag(U.obj$U))

  # Set up temporal precision 
  tmdg = c(0,rep(tmrho*tmrho,ntm-2),0)
  tpltz = toeplitz(c(1,-1.0*tmrho,rep(0,ntm-2)))
  prctm = (1.0 / (1-tmrho*tmrho)) * Matrix( diag(tmdg) + tpltz,sparse=TRUE)
  d2 = Matrix::determinant(prctm,logarithm=TRUE)

  # Full determinant
  d1fl = -2.0 * (ntm * sum(lgu) - nloc * d2$modulus)

  lkcmp = rep(0,nfld)

  # compute the loglikelihood
  for (k in 1:nfld) {
    # Quad form stuff
    ztrns = array(0,c(nloc,ntm))
    for (t in 1:ntm) {
      zord = sptarr[U.obj$ord.z,t,k]
      zchk = Matrix::crossprod(U.obj$U,zord)
      ztrns[,t] = zchk@x
    }
    qftms = apply(ztrns,1,tmqf, prcmt = prctm)
    lkcmp[k] = -0.5 * (d1fl + sum(qftms) + nloc*nrp*log(2*pi))
    
  }

  lksum = sum(lkcmp)
  return(lksum)
}

