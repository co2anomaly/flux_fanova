# Generate initial guesses for MIP functional ANOVA states and parameters

library(fields)
library(ncdf4)
library(jsonlite)
library(plyr)
library(colorspace)
library(GPvecchia)
source("func_anova_fns.R")

cfgfile = "config/fanova_mip_namer_prdev.csv"
exmpcfg = read.csv(cfgfile,header = TRUE,stringsAsFactors = FALSE)
cfglst = as.list(exmpcfg$Value)
names(cfglst) = exmpcfg$Field

prinfo = fromJSON(txt = cfglst$json_file)

insd = as.integer(cfglst$init_seed)
set.seed(insd)
nchain = as.integer(cfglst$nchain)
nloc = as.integer(cfglst$nloc)

nmnA = as.integer(cfglst$nleva)
nmnB = as.integer(cfglst$nlevb)
nab = (nmnA-1) * (nmnB-1)
nrp = as.integer(cfglst$nrep)

# Read in Dataset
nc1 = nc_open(cfglst$data_file)
dtarr = ncvar_get(nc1,cfglst$data_variable)
ctrstA = ncvar_get(nc1,"contrast_A")
ctrstB = ncvar_get(nc1,"contrast_B")
locx = ncvar_get(nc1,cfglst$location_x_name)
locy = ncvar_get(nc1,cfglst$location_y_name)
latmip = as.vector(ncvar_get(nc1,"latitude"))
lonmip = as.vector(ncvar_get(nc1,"longitude"))
nc_close(nc1)

nloc = length(latmip)
lcfrm = data.frame(LocID = seq(1,nloc),Longitude = lonmip,Latitude = latmip)

theme_mat = theme_bw() 
theme_mat$axis.title.x$size = 14
theme_mat$axis.text.x$size = 12
theme_mat$axis.title.y$size = 14
theme_mat$axis.text.y$size = 12
theme_mat$plot.title$size = 16
theme_mat$plot.title$hjust = 0.5
theme_mat$legend.text$size = 12
theme_mat$strip.text$size = 14


locs = cbind(locx,locy)
dsts = rdist(locs)

# Coef estimates to initialize
cfbsln = aaply(dtarr,1,.fun = fit_anova_mltlev_arr3d, .progress="text", contrstA=ctrstA, contrstB=ctrstB)
sdbsln = apply(dtarr,1:3,sd)
gpsd = median(sdbsln) * 0.5
cfsd = sd(cfbsln)

# groupings of GPs
vct = 0
varlst = list()
gpgrps = c("mean","mainA","mainB","interact","noise")
gpvnms = c("gp_stddev","gp_range","gp_smoothness")
gplng = c("overall mean","main effect factor A","main effect factor B","interaction","replicate error")


# Parameter sequences
m=20
vec.spec=vecchia_specify(locs,m,cond.yz='y')

nusq = seq(0.5,3.0,by=0.5)
nsmth = length(nusq)

nrng = 10
mxdst = ceiling(0.5 * max(dsts) / 1000) * 1000    # Original was 0.2
dstsq = seq(1.0 * mxdst/nrng, mxdst, length.out = nrng)

nsd = 10
sdfld = apply(cfbsln,2,sd)
sdmx = ceiling(1.5 * max(sdfld) / 100) * 100    # Default was 2
sdsq = seq(1.*sdmx/nsd, sdmx, length.out = nsd)

mulks = array(0,c(10,10,length(nusq)))
for (i in seq(1,nsd)) { 
    for (j in seq(1,nrng)) {
        for (k in seq(1,nsmth)) {
            cvvc = c(sdsq[i] * sdsq[i], dstsq[j], nusq[k])
            tmplk = vecchia_likelihood(cfbsln[,"Intercept"], vec.spec, cvvc, nuggets = 0)
            strrp = sprintf("Sigma: %.1f, Range: %.1f, Nu: %.1f    LogLik: %.4e",sdsq[i], dstsq[j], nusq[k], tmplk)
            print(strrp)
            mulks[i,j,k] = tmplk
        }
    }
}

dimnames(mulks)[1] = list(sdsq)
dimnames(mulks)[2] = list(dstsq)
dimnames(mulks)[3] = list(nusq)
mumlt = melt(mulks,varnames = c("StdDev","Range","Smooth"))
max_mean = mumlt[mumlt$value == max(mumlt$value,na.rm = TRUE),]
mumlt$LkDev = mumlt$value - max_mean$value[1] 

#r7 = diverge_hcl(7,h=c(240,0),c=80,l=c(50,100),power=3)
gmu = ggplot(mumlt,aes(x=Range,y=LkDev)) + geom_point(aes(color=StdDev)) + facet_wrap(~ Smooth,nrow=2) + 
      xlab("Range [km]") + 
      scale_y_continuous("Log Like Diff",limits=c(-10000,0)) + 
      ggtitle("Mean Likelihood Grid Search") + theme_mat 
pnm = paste0(cfglst$output_dir,"/MIP_Like_Mean.pdf")
pdf(pnm,width = 9,height=6)
print(gmu)
dev.off()

# Main A 
alks = array(0,c(10,10,length(nusq)))
for (i in seq(1,nsd)) { 
  for (j in seq(1,nrng)) {
    for (k in seq(1,nsmth)) {
      cvvc = c(sdsq[i] * sdsq[i], dstsq[j], nusq[k])
      tmplk = 0
      for (q in seq(1,nmnA-1)) {
          clnm = paste("A",q,sep="_")
          tmplk = tmplk + vecchia_likelihood(cfbsln[,clnm], vec.spec, cvvc, nuggets = 0)
      }
      strrp = sprintf("Sigma: %.1f, Range: %.1f, Nu: %.1f    LogLik: %.4e",sdsq[i], dstsq[j], nusq[k], tmplk)
      print(strrp)
      alks[i,j,k] = tmplk
    }
  }
}

dimnames(alks)[1] = list(sdsq)
dimnames(alks)[2] = list(dstsq)
dimnames(alks)[3] = list(nusq)
amlt = melt(alks,varnames = c("StdDev","Range","Smooth"))
max_mainA = amlt[amlt$value == max(amlt$value,na.rm = TRUE),]
amlt$LkDev = amlt$value - max_mainA$value[1] 

gmna = ggplot(amlt,aes(x=Range,y=LkDev)) + geom_point(aes(color=StdDev)) + facet_wrap(~ Smooth,nrow=2) + 
       xlab("Range [km]") + 
       scale_y_continuous("Log Like Diff",limits=c(-10000,0)) + 
       ggtitle("Model Effect Likelihood Grid Search") + theme_mat 
pnm = paste0(cfglst$output_dir,"/MIP_Like_MainA.pdf")
pdf(pnm,width = 9,height=6)
print(gmna)
dev.off()


# Main B 
blks = array(0,c(10,10,length(nusq)))
for (i in seq(1,nsd)) { 
  for (j in seq(1,nrng)) {
    for (k in seq(1,nsmth)) {
      cvvc = c(sdsq[i] * sdsq[i], dstsq[j], nusq[k])
      tmplk = 0
      for (q in seq(1,nmnB-1)) {
        clnm = paste("B",q,sep="_")
        tmplk = tmplk + vecchia_likelihood(cfbsln[,clnm], vec.spec, cvvc, nuggets = 0)
      }
      strrp = sprintf("Sigma: %.1f, Range: %.1f, Nu: %.1f    LogLik: %.4e",sdsq[i], dstsq[j], nusq[k], tmplk)
      print(strrp)
      blks[i,j,k] = tmplk
    }
  }
}

dimnames(blks)[1] = list(sdsq)
dimnames(blks)[2] = list(dstsq)
dimnames(blks)[3] = list(nusq)
bmlt = melt(blks,varnames = c("StdDev","Range","Smooth"))
max_mainB = bmlt[bmlt$value == max(bmlt$value,na.rm = TRUE),]
bmlt$LkDev = bmlt$value - max_mainB$value[1] 

gmnb = ggplot(bmlt,aes(x=Range,y=LkDev)) + geom_point(aes(color=StdDev)) + facet_wrap(~ Smooth,nrow=2) + 
       xlab("Range [km]") + 
       scale_y_continuous("Log Like Diff",limits=c(-10000,0)) + 
       ggtitle("Data Effect Likelihood Grid Search") + theme_mat 
pnm = paste0(cfglst$output_dir,"/MIP_Like_MainB.pdf")
pdf(pnm,width = 9,height=6)
print(gmnb)
dev.off()



# Interact
ilks = array(0,c(10,10,length(nusq)))
for (i in seq(1,nsd)) { 
  for (j in seq(1,nrng)) {
    for (k in seq(1,nsmth)) {
      cvvc = c(sdsq[i] * sdsq[i], dstsq[j], nusq[k])
      tmplk = 0
      for (p in seq(1,nmnA-1)) {
        for (q in seq(1,nmnB-1)) {
          clnm = paste("AB",p,q,sep="_")
          tmplk = tmplk + vecchia_likelihood(cfbsln[,clnm], vec.spec, cvvc, nuggets = 0)
        }
      }
      strrp = sprintf("Sigma: %.1f, Range: %.1f, Nu: %.1f    LogLik: %.4e",sdsq[i], dstsq[j], nusq[k], tmplk)
      print(strrp)
      ilks[i,j,k] = tmplk
    }
  }
}

dimnames(ilks)[1] = list(sdsq)
dimnames(ilks)[2] = list(dstsq)
dimnames(ilks)[3] = list(nusq)
imlt = melt(ilks,varnames = c("StdDev","Range","Smooth"))
max_interact = imlt[imlt$value == max(imlt$value,na.rm = TRUE),]
imlt$LkDev = imlt$value - max_interact$value[1] 

gmni = ggplot(imlt,aes(x=Range,y=LkDev)) + geom_point(aes(color=StdDev)) + facet_wrap(~ Smooth,nrow=2) + 
       xlab("Range [km]") + 
       scale_y_continuous("Log Like Diff",limits=c(-10000,0)) + 
       ggtitle("Interaction Effect Likelihood Grid Search") + theme_mat 
pnm = paste0(cfglst$output_dir,"/MIP_Like_Interact.pdf")
pdf(pnm,width = 9,height=6)
print(gmni)
dev.off()

# Output
cvprgrd = list( mean = c(max_mean$StdDev[1]^2, max_mean$Range[1], max_mean$Smooth[1]),
                mainA = c(max_mainA$StdDev[1]^2, max_mainA$Range[1], max_mainA$Smooth[1]),
                mainB = c(max_mainB$StdDev[1]^2, max_mainB$Range[1], max_mainB$Smooth[1]),
                interact = c(max_interact$StdDev[1]^2, max_interact$Range[1], max_interact$Smooth[1]),
                noise = c(300^2,150,1.0))
save(cvprgrd,file = cfglst$like_est)
