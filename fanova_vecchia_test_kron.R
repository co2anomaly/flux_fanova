# Test Vecchia setup with Kronecker product for time 

library(fields)
library(ncdf4)
library(jsonlite)
library(GPvecchia)
library(Matrix)
source("func_anova_tidy.R")
source("func_anova_vecchia.R")

# Temporary
config = "config/fanova_mip_namer_prdev.csv"
chain = 1

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

# Read in inital values
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

# Vecchia setup
m=20
vec.spec=vecchia_specify(locs,m,cond.yz='y')


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
    pr_gp_nu[[gpgrps[g1] ]] = c(prinfo$prior[[ gpgrps[g1]]][["gp_smoothness"]][["mean"]],
                                prinfo$prior[[ gpgrps[g1]]][["gp_smoothness"]][["stddev"]])
}

# Bookkeeping
saveidx = 0
aidx = 0

# Initialize error (epsilon) ANOVA cov
U.eps=createU(vec.spec,covparms=cvprlst$noise,nuggets=0)$U
W.eps=Matrix::tcrossprod(U.eps) # precision matrix

# Set up standard matrix, set sparse = TRUE
ntm = 3
rhotm = 0.3
tmdg = c(0,rep(rhotm*rhotm,ntm-2),0)
tpltz = toeplitz(c(1,-1.0*rhotm,rep(0,ntm-2)))
prctm = (1.0 / (1-rhotm*rhotm)) * Matrix( diag(tmdg) + tpltz,sparse=TRUE)

W.spctm = prctm %x% W.eps
# Force symmetric?
W.spctmu = forceSymmetric(W.spctm)
  
mncr = array(0,c(nloc,nmnA,nmnB))
## Epsilon cov params
# mu + ic[i]*alpha + jc[j]*beta + ic[i]*jc[j]*ab
ydv = array(0,c(nloc,nrp*nmnA*nmnB))
lkhld = rep(0,nrp*nmnA*nmnB)
ct0 = 0
print("Starting likelihood")
for(i in 1:nmnA) { 
  for(j in 1:nmnB) {
    ictrst = ctrstA[i,] %x% ctrstB[j,]
    tst1 = mu.cur + a.cur %*% as.vector(ctrstA[i,]) + b.cur %*% as.vector(ctrstB[j,]) + ab.cur %*% as.vector(ictrst)
    mncr[,i,j] = as.vector(tst1)
    for (k in 1:nrp) {
      ct0 = ct0 + 1
      ydv[,ct0] = dtarr[,i,j,k] - mncr[,i,j]
      lkhld = vecchia_likelihood(ydv[,ct0],vec.spec,cvprlst$noise,nuggets=0)
    }
  }
}
l1txt = sprintf("Likelihood 1 Evaluation: %.6e",sum(lkhld))
print(l1txt)
 
lk2 = vecchia_like_array(ydv[,22:24], vec.spec, cvprlst$noise, nrep=nrp)
l2txt = sprintf("Likelihood 2 Evaluation: %.6e",lk2)
print(l2txt)

# Possible posterior precision shortcuts   
#W.spctmu = forceSymmetric(W.spctm)

sptdsgn = Matrix(diag(nrow=nloc), sparse=TRUE) 
tmdsgn = Matrix(rep(1,3), nrow=3)
fldsgn = tmdsgn %x% sptdsgn

dtprc = (Matrix::t(fldsgn)) %*% W.spctmu %*% fldsgn
dtprcu = forceSymmetric(dtprc)

tmprc = (Matrix::t(tmdsgn)) %*% prctm %*% tmdsgn

print(tmprc)
print('Summary of sptm prec / spatial only')
print(summary(dtprcu@x  / W.eps@x))

# Determinant of full precision
d1 = Matrix::determinant(W.spctmu,logarithm=TRUE)

U.obj = createU(vec.spec,cvprlst$noise,nuggets=0)
lgu = log(Matrix::diag(U.obj$U))

d2 = Matrix::determinant(prctm,logarithm=TRUE)
d1fl = -2.0 * (nrp * sum(lgu) - nloc*d2$modulus)

d1txt = sprintf("Full Det Evaluation: %.4e",d1$modulus)
print(d1txt)
 
d2txt = sprintf("Det Comp Evaluation: %.4e",d1fl)
print(d2txt)

# Quad form stuff
ztrns = array(0,c(nloc,nrp))
for (k in 1:nrp) {
  zord = ydv[U.obj$ord.z,k+21]
  zchk = Matrix::crossprod(U.obj$U,zord)
  ztrns[,k]=zchk@x
}

tmqf = function(zvc, prcmt) {
    # Quick calculation of temporal quadratic form
    qfrm = Matrix::t(zvc) %*% prcmt %*% zvc
    return(as.vector(qfrm))
}
qftms = apply(ztrns,1,tmqf, prcmt = prctm)

lk3 = vecchia_like_array(ydv[,22:24], vec.spec, cvprlst$noise, nrep=nrp)
l3txt = sprintf("Likelihood 3 Evaluation: %.6e",lk3)
print(l3txt)

l3b = -0.5 * (d1fl + sum(qftms) + nloc*nrp*log(2*pi)) 
l3btxt = sprintf("Likelihood 3B Evaluation: %.6e",l3b)
print(l3btxt)

