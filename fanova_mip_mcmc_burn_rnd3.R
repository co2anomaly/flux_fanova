
# Read flux MIP data for 2-way factorial 
# Multi-dimensional MH updates for GP parameters 
# Expected command line arguments
#    config: CSV file with MCMC info
#    chain:  Chain number
# Round 2: Incorporate MH updates for GP parameters
# Alternate prior for smoothness

library(ncdf4)
library(jsonlite)
library(GPvecchia)
source("func_anova_fns.R")
source("func_anova_vecchia.R")

# If running interactively 
#config = "config/fanova_mip_namer_prdev.csv"
#chain = 2

args=(commandArgs(TRUE))

config = args[1]
chain = as.integer(args[2])
ptxt = sprintf("Config: %s, Chain %d",config,chain)
print(ptxt)

# Read configuration
exmpcfg = read.csv(config,header = TRUE,stringsAsFactors = FALSE)
cfglst = as.list(exmpcfg$Value)
names(cfglst) = exmpcfg$Field

prinfo = fromJSON(txt = cfglst$json_file)

nchain = as.integer(cfglst$nchain)
niter = as.integer(cfglst$nburn)
nthin = as.integer(cfglst$nthin)
nsamp = floor(niter / nthin)
nadapt = as.integer(cfglst$burn_mh_adapt)
nmhsv = floor(niter / nadapt)
sadapt = as.integer(cfglst$save_adapt)

insd = as.integer(cfglst$burn_seed) + chain*as.integer(cfglst$burn_inc)
set.seed(insd)
nloc = as.integer(cfglst$nloc)
nmnA = as.integer(cfglst$nleva)
nmnB = as.integer(cfglst$nlevb)
nab = (nmnA-1) * (nmnB-1)
nrp = as.integer(cfglst$nrep)

ctxt = sprintf("Configuration parsed. Factor A has %d levels. Factor B has %d levels. %d replicates\n",nmnA,nmnB,nrp)
print(ctxt)

lgfl = paste0(cfglst$burn_log_file,chain,".txt")
write(ptxt,file=lgfl)
write(ctxt,file=lgfl,append=TRUE)

# Read data and locations
nc1 = nc_open(cfglst$data_file)
dtarr = ncvar_get(nc1,cfglst$data_variable)
ctrstA = ncvar_get(nc1,"contrast_A")
ctrstB = ncvar_get(nc1,"contrast_B")
locx = ncvar_get(nc1,cfglst$location_x_name)
locy = ncvar_get(nc1,cfglst$location_y_name)
nc_close(nc1)

ctrlon = -110
if ((cfglst$great_circ == "yes") || (cfglst$great_circ == "Yes")) {
    locs = as.matrix(cbind( (locx - ctrlon) * cos(pi*locy/180),locy))
    dstmt = rdist.earth(locs,miles = FALSE)
} else {
    locs=as.matrix(cbind(locx,locy))
    dstmt = rdist(locs)
}

ctxt = "Data and locations read"
print(ctxt)
write(ctxt,file=lgfl,append=TRUE)

if (is.na(ncol(ctrstA))) {
    ctrstA = matrix(ctrstA,ncol=1)
}
if (is.na(ncol(ctrstB))) {
    ctrstB = matrix(ctrstB,ncol=1)
}


# Covariance params, setup default values
# Use list structures convenient for GPvecchia (variance,range,smoothness)
cvprlst = list( mean = c(200,2000,1.5), mainA = c(1000,150,1.5), mainB = c(1000,150,1.5), 
                interact = c(500,100,1.5), noise = c(300,100,1.5) )
gpgrps = c("mean","mainA","mainB","interact","noise")
gpvnms = c("gp_stddev","gp_range","gp_smoothness")

# Read in inital values - all elements stored in init files
ncinit = paste0(cfglst$burn_init_file,chain,".nc")
nc1 = nc_open(ncinit)
for (g1 in seq(1,length(gpgrps))) {
  sgnm = paste0("gp_stddev_",gpgrps[g1])
  sdtmp = ncvar_get(nc1,sgnm)
  cvprlst[[gpgrps[g1]]][1] = sdtmp^2
  lmnm = paste0("gp_range_",gpgrps[g1])
  cvprlst[[gpgrps[g1]]][2] = ncvar_get(nc1,lmnm)
  nunm = paste0("gp_smoothness_",gpgrps[g1])
  cvprlst[[gpgrps[g1]]][3] = ncvar_get(nc1,nunm)
}
mu.cur = as.vector(ncvar_get(nc1,"gp_field_mean"))
a.cur = ncvar_get(nc1,"gp_field_mainA")
b.cur = ncvar_get(nc1,"gp_field_mainB")
ab.cur = ncvar_get(nc1,"gp_field_interact")
nc_close(nc1)

if (nmnA == 2) {
    dim(a.cur) = c(nloc,nmnA-1)
}
if (nmnB == 2) {
    dim(b.cur) = c(nloc,nmnB-1)
}
if (nab == 1) {
    dim(ab.cur) = c(nloc,nab)
}

print(cvprlst)

# For reverting
mu.prv = mu.cur
a.prv = a.cur 
b.prv = b.cur 
ab.prv = ab.cur 
cvprprv = cvprlst

# Output arrays
mu.samp=array(dim=c(sadapt,nloc))
a.samp=array(dim=c(sadapt,nloc,nmnA-1))
b.samp=array(dim=c(sadapt,nloc,nmnB-1))
ab.samp=array(dim=c(sadapt,nloc,nab))
sig.samp=list(mean = rep(0,sadapt), mainA = rep(0,sadapt), mainB = rep(0,sadapt), interact = rep(0,sadapt), noise = rep(0,sadapt))
lam.samp=list(mean = rep(0,sadapt), mainA = rep(0,sadapt), mainB = rep(0,sadapt), interact = rep(0,sadapt), noise = rep(0,sadapt))
nu.samp=list(mean = rep(0,sadapt), mainA = rep(0,sadapt), mainB = rep(0,sadapt), interact = rep(0,sadapt), noise = rep(0,sadapt))
mncr = array(0,c(nloc,nmnA,nmnB))

# Vecchia setup
m=20
vec.spec=vecchia_specify(locs,m,cond.yz='y')


# Metropolis-Hastings setup
ac_mh = list(mean = 0, mainA = 0, mainB = 0, interact = 0, noise = 0)
acrt_mh = list(mean = 0, mainA = 0, mainB = 0, interact = 0, noise = 0)
jpsd_mh = list(mean = 1.0, mainA = 1.0, mainB = 1.0, interact = 1.0, noise = 1.0)
for (g1 in seq(1,length(gpgrps))) {
    jpsd_mh[[gpgrps[g1] ]] = prinfo$mh[[ gpgrps[g1] ]][["stddev"]]
}

mhcor3d = matrix(c(1,0.6,-0.4, 0.6,1,-0.8, -0.4,-0.8,1),nrow=3)
mhsd3d = diag(c(0.6,1.0,0.6))
mhcv3d = mhsd3d %*% mhcor3d %*% mhsd3d
mhchl3d = t(chol(mhcv3d))

# Optional refinement of MH proposal
chol_mh = list(mean = mhchl3d, mainA = mhchl3d, mainB = mhchl3d,
               interact = mhchl3d, noise = mhchl3d) 
for (g1 in seq(1,length(gpgrps))) {
    if ( (!is.null(prinfo$mh[[ gpgrps[g1] ]][["corvec"]])) &&
         ( !is.null(prinfo$mh[[ gpgrps[g1] ]][["sdvec"]])) ) {
        mhcor3d = matrix(prinfo$mh[[ gpgrps[g1] ]][["corvec"]],nrow=3)
        mhsd3d = diag( prinfo$mh[[ gpgrps[g1] ]][["sdvec"]]  )
        mhcv3d = mhsd3d %*% mhcor3d %*% mhsd3d
        chol_mh[[gpgrps[g1] ]] =  t(chol(mhcv3d))
    }
}

# Sampling file
ncsmpl = paste0(cfglst$burn_samp_file,chain,".nc")

# Priors, from JSON
mu0 = 0
prpars = c(0.7,1.0)
prnu = c(0.25,0.25)
pr_gp_sig = list(mean = 1.0, mainA = 1.0, mainB = 1.0, interact = 1.0, noise = 1.0)
pr_gp_lam = list(mean = prpars, mainA = prpars, mainB = prpars, interact = prpars, noise = prpars)
pr_gp_nu = list(mean = prnu, mainA = prnu, mainB = prnu, interact = prnu, noise = prnu)
for (g1 in seq(1,length(gpgrps))) {
    pr_gp_sig[[gpgrps[g1] ]] = prinfo$prior[[ gpgrps[g1]]][["gp_stddev"]][["stddev"]]
    pr_gp_lam[[gpgrps[g1] ]] = c(prinfo$prior[[ gpgrps[g1]]][["gp_range"]][["mean"]],
                                 prinfo$prior[[ gpgrps[g1]]][["gp_range"]][["stddev"]])
    pr_gp_nu[[gpgrps[g1] ]] = c(prinfo$prior[[ gpgrps[g1]]][["gp_smoothness"]][["shape1"]],
                                prinfo$prior[[ gpgrps[g1]]][["gp_smoothness"]][["shape2"]],
                                prinfo$prior[[ gpgrps[g1]]][["gp_smoothness"]][["lim1"]],
                                prinfo$prior[[ gpgrps[g1]]][["gp_smoothness"]][["lim2"]])
}

# Bookkeeping
saveidx = 0
aidx = 0

# Initialize error (epsilon) ANOVA cov
U.eps=createU(vec.spec,covparms=cvprlst$noise,nuggets=0)$U
W.eps=Matrix::tcrossprod(U.eps) # precision matrix


for(l in 1:niter){
  
  str1 = sprintf("Iteration %04d",l)
  #print(str1)
 
  # Thinning/output
  if ( (l %% nthin) == 0) {
    print(str1)
      save = 1
      sct = (saveidx %% sadapt) + 1
      saveidx = saveidx + 1
  }
  else {
      save = 0
  }
  if ( (l %% nadapt) == 0 ) {
      # MH Update
      aidx = aidx + 1
      gpgrps = c("mean","mainA","mainB","interact","noise")
      for (g1 in seq(1,length(gpgrps))) {
          acrt_mh[[gpgrps[g1] ]] = 1.0 * ac_mh[[gpgrps[g1] ]] / nadapt

          # Output result
          ncmh = nc_open(ncsmpl,write = TRUE)
          sdnm = paste0("mh_stddev_",gpgrps[g1])
          ncvar_put(ncmh,sdnm,jpsd_mh[[gpgrps[g1] ]],start = aidx,count=1)
          rtnm = paste0("mh_acc_rate_",gpgrps[g1])
          ncvar_put(ncmh,rtnm,acrt_mh[[gpgrps[g1] ]],start = aidx,count=1)
          nc_close(ncmh)
          
          # Update
          jpsd_mh[[gpgrps[g1] ]] = rand_mh_rate(acrt_mh[[gpgrps[g1] ]], jpsd_mh[[gpgrps[g1] ]])
          ac_mh[[gpgrps[g1] ]] = 0
      }  
        
      # Logging
      ctxt = sprintf("\nIteration %d\n  Acc Rate Mean: %.4f, Noise: %.4f\n  %s",
                     l,acrt_mh$mean,acrt_mh$noise, Sys.time())
      write(ctxt,file=lgfl,append=TRUE)
      
  }
  
  ### State/parameter updates
  ### Use previously computed inverse correlation matrices where possible
  ### Posterior mean, Cholesky of posterior covariance still needed 
 
  ## update mu
  ysum = apply(dtarr,1,sum) - (nmnA*nmnB*nrp)*mu0
  mu.cur = samp_post_vecchia(ysum,cvprlst$mean,W.eps,vec.spec,nrp,nmnA*nmnB,nloc)
  
  ## update alpha - multiple columns
  print("Update alpha")
  ysma = apply(dtarr,1:2,sum) %*% ctrstA
  if (nmnA > 2) {
    a.cur = samp_post_vecchia_mltlev(ysma,cvprlst$mainA,W.eps,vec.spec,nrp,nmnA*nmnB,nloc)
  }
  else {
    a.cur = samp_post_vecchia(ysma,cvprlst$mainA,W.eps,vec.spec,nrp,nmnA*nmnB,nloc)
    dim(a.cur) = c(nloc,nmnA-1)
  }
  
  ## update beta - multiple columns 
  ysmb = apply(dtarr,c(1,3),sum) %*% ctrstB
  if (nmnB > 2) { 
    b.cur = samp_post_vecchia_mltlev(ysmb,cvprlst$mainB,W.eps,vec.spec,nrp,nmnA*nmnB,nloc)
  }
  else {
    b.cur = samp_post_vecchia(ysmb,cvprlst$mainB,W.eps,vec.spec,nrp,nmnA*nmnB,nloc)
    dim(b.cur) = c(nloc,nmnB-1)
  }
  
  ## update ab - loop
  ysmab = apply(dtarr,1:3,sum)
  y.tilde=array(0,c(nloc,nab))
  for(i in 1:nmnA) { 
    for(j in 1:nmnB) {
      ictrst = ctrstA[i,] %x% ctrstB[j,]
      y.tilde = y.tilde + ysmab[,i,j] %x% t(ictrst)
    }
  }
  if ( (nmnA > 2) | (nmnB > 2) ) {
    ab.cur = samp_post_vecchia_mltlev(y.tilde,cvprlst$interact,W.eps,vec.spec,nrp,nmnA*nmnB,nloc)
  }
  else {
    ab.cur = samp_post_vecchia(y.tilde,cvprlst$interact,W.eps,vec.spec,nrp,nmnA*nmnB,nloc)
    dim(ab.cur) = c(nloc,nab)
  }
  
  ## Epsilon cov params
  # mu + ic[i]*alpha + jc[j]*beta + ic[i]*jc[j]*ab
  ydv = array(0,c(nloc,nrp*nmnA*nmnB))
  ct0 = 0
  for(i in 1:nmnA) { 
    for(j in 1:nmnB) {
      ictrst = ctrstA[i,] %x% ctrstB[j,]
      tst1 = mu.cur + a.cur %*% as.vector(ctrstA[i,]) + b.cur %*% as.vector(ctrstB[j,]) + ab.cur %*% as.vector(ictrst)
      mncr[,i,j] = as.vector(tst1)
      for (k in 1:nrp) {
        ct0 = ct0 + 1
        ydv[,ct0] = dtarr[,i,j,k] - mncr[,i,j]
      }
    }
  }
  
  tmpsd = sqrt(cvprlst$noise[1])
  epsmh = spatial_mh3d_vecchia(tmpsd,cvprlst$noise[2],cvprlst$noise[3], vec.spec,ydv,
                               ct0,nloc, jpsd_mh$noise,chol_mh$noise,pr_gp_sig$noise[1],pr_gp_lam$noise,pr_gp_nu$noise,prnumod='beta')
  ac_mh$noise = ac_mh$noise + epsmh$acpt
  cvprlst$noise = c(epsmh$sdpst^2, epsmh$lampst, epsmh$nupst)

  # If new params, update Vecchia
  if (epsmh$acpt == 1) {  
    U.eps=createU(vec.spec,covparms=cvprlst$noise,nuggets=0)$U
    W.eps=Matrix::tcrossprod(U.eps) # precision matrix
  }

  # Other GP parameter updates
  tmpsd = sqrt(cvprlst$interact[1])
  abmh = spatial_mh3d_vecchia(tmpsd,cvprlst$interact[2],cvprlst$interact[3], vec.spec, ab.cur,
                              nab,nloc, jpsd_mh$interact,chol_mh$interact, 
                              pr_gp_sig$interact[1],pr_gp_lam$interact,pr_gp_nu$interact,prnumod='beta')
  ac_mh$interact = ac_mh$interact + abmh$acpt
  cvprlst$interact = c(abmh$sdpst^2, abmh$lampst, abmh$nupst)

  tmpsd = sqrt(cvprlst$mainA[1])
  amh = spatial_mh3d_vecchia(tmpsd,cvprlst$mainA[2],cvprlst$mainA[3], vec.spec, a.cur,
                              nmnA-1,nloc, jpsd_mh$mainA,chol_mh$mainA, 
                              pr_gp_sig$mainA[1],pr_gp_lam$mainA,pr_gp_nu$mainA,prnumod='beta')
  ac_mh$mainA = ac_mh$mainA + amh$acpt
  cvprlst$mainA = c(amh$sdpst^2, amh$lampst, amh$nupst)

  tmpsd = sqrt(cvprlst$mainB[1])
  bmh = spatial_mh3d_vecchia(tmpsd,cvprlst$mainB[2],cvprlst$mainB[3], vec.spec, b.cur,
                              nmnB-1,nloc, jpsd_mh$mainB,chol_mh$mainB, 
                              pr_gp_sig$mainB[1],pr_gp_lam$mainB,pr_gp_nu$mainB,prnumod='beta')
  ac_mh$mainB = ac_mh$mainB + bmh$acpt
  cvprlst$mainB = c(bmh$sdpst^2, bmh$lampst, bmh$nupst)

  tmpsd = sqrt(cvprlst$mean[1])
  mumh = spatial_mh3d_vecchia(tmpsd,cvprlst$mean[2],cvprlst$mean[3], vec.spec, mu.cur,
                              1,nloc, jpsd_mh$mean,chol_mh$mean, 
                              pr_gp_sig$mean[1],pr_gp_lam$mean,pr_gp_nu$mean,prnumod='beta')
  ac_mh$mean = ac_mh$mean + mumh$acpt
  cvprlst$mean = c(mumh$sdpst^2, mumh$lampst, mumh$nupst)

  # Final check for NAs, etc.
  natot = length(mu.cur[is.na(mu.cur)]) 
  if (nmnA > 2) { 
    natot = natot + length(a.cur@x[is.na(a.cur@x)])
  } else {
    natot = natot + length(a.cur[is.na(a.cur)])
  }
  if (nmnB > 2) { 
    natot = natot + length(b.cur@x[is.na(b.cur@x)])
  } else {
    natot = natot + length(b.cur[is.na(b.cur)])
  }
  if ( (nmnA > 2) | (nmnB > 2) ) {
    natot = natot + length(ab.cur@x[is.na(ab.cur@x)])
  } else {
    natot = natot + length(ab.cur[is.na(ab.cur)])
  }

  if (natot > 0) {
    mu.cur = mu.prv
    a.cur = a.prv 
    b.cur = b.prv 
    ab.cur = ab.prv 
    cvprlst = cvprprv
    ac_mh$mean = ac_mh$mean - mumh$acpt
    ac_mh$mainA = ac_mh$mainA - amh$acpt
    ac_mh$mainB = ac_mh$mainB - bmh$acpt
    ac_mh$interact = ac_mh$interact - abmh$acpt
    ac_mh$noise = ac_mh$noise - epsmh$acpt
    if (epsmh$acpt == 1) {  
      U.eps=createU(vec.spec,covparms=cvprlst$noise,nuggets=0)$U
      W.eps=Matrix::tcrossprod(U.eps) # precision matrix
    }
    ctxt = sprintf("\n  %d Missing at Iteration %d\n", 
                    natot, l)
    write(ctxt,file=lgfl,append=TRUE)
  } else {            
    mu.prv = mu.cur
    a.prv = a.cur 
    b.prv = b.cur 
    ab.prv = ab.cur 
    cvprprv = cvprlst
  }

  # Save/store output
  if (save == 1) {
    for (g1 in seq(1,length(gpgrps))) {
        sig.samp[[gpgrps[g1] ]][sct] = sqrt(cvprlst[[gpgrps[g1] ]][1])
        lam.samp[[gpgrps[g1] ]][sct] = cvprlst[[gpgrps[g1] ]][2]
        nu.samp[[gpgrps[g1] ]][sct] = cvprlst[[gpgrps[g1] ]][3]
    }
    mu.samp[sct,] = mu.cur
    if (nmnA > 2) { 
      a.samp[sct,,] = as.matrix(a.cur)
    }
    else {
      a.samp[sct,,] = a.cur
    }
    if (nmnB > 2) {  
      b.samp[sct,,] = as.matrix(b.cur)
    }
    else {
      b.samp[sct,,] = b.cur
    }
    if ( (nmnA > 2) | (nmnB > 2) ) {
      ab.samp[sct,,] = as.matrix(ab.cur)
    }
    else {
      ab.samp[sct,,] = ab.cur
    }
    if ( (saveidx %% sadapt) == 0) {
        ftxt = sprintf('Saving: saveidx %d, sct %d',saveidx,sct) 
        print(ftxt)
        # Write output
        ast = saveidx + 1 - sadapt
        nc1 = nc_open(ncsmpl,write=TRUE)
        for (g1 in seq(1,length(gpgrps))) {
          sgnm = paste0("gp_stddev_",gpgrps[g1])
          ncvar_put(nc1,sgnm,sig.samp[[gpgrps[g1]]],start = c(ast),count=c(sadapt))
          lmnm = paste0("gp_range_",gpgrps[g1])
          ncvar_put(nc1,lmnm,lam.samp[[gpgrps[g1]]],start = c(ast),count=c(sadapt))
          nunm = paste0("gp_smoothness_",gpgrps[g1])
          ncvar_put(nc1,nunm,nu.samp[[gpgrps[g1]]],start = c(ast),count=c(sadapt))
        }
        ncvar_put(nc1,"gp_field_mean",mu.samp,start = c(ast,1),count=c(sadapt,nloc))
        ncvar_put(nc1,"gp_field_mainA",a.samp,start = c(ast,1,1),count=c(sadapt,nloc,nmnA-1))
        ncvar_put(nc1,"gp_field_mainB",b.samp,start = c(ast,1,1),count=c(sadapt,nloc,nmnB-1))
        ncvar_put(nc1,"gp_field_interact",ab.samp,start = c(ast,1,1),count=c(sadapt,nloc,nab))
        nc_close(nc1)
    }
  }
}

# Save final states
ctxt = "Saving final states"
print(ctxt)
write(ctxt,file=lgfl,append=TRUE)

ncfnl = paste0(cfglst$burn_final_file,chain,".nc")
nc1 = nc_open(ncfnl,write=TRUE)
for (g1 in seq(1,length(gpgrps))) {
  sgnm = paste0("gp_stddev_",gpgrps[g1])
  sdtmp = sqrt(cvprlst[[gpgrps[g1]]][1])
  ncvar_put(nc1,sgnm,sdtmp)
  lmnm = paste0("gp_range_",gpgrps[g1])
  ncvar_put(nc1,lmnm,cvprlst[[gpgrps[g1]]][2])
  nunm = paste0("gp_smoothness_",gpgrps[g1])
  ncvar_put(nc1,nunm,cvprlst[[gpgrps[g1]]][3])
}
ncvar_put(nc1,"gp_field_mean",mu.cur)
ncvar_put(nc1,"gp_field_mainA",a.cur)
ncvar_put(nc1,"gp_field_mainB",b.cur)
ncvar_put(nc1,"gp_field_interact",ab.cur)
nc_close(nc1)


