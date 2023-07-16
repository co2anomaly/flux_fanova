# Plot prior densities

library(tidyverse)
library(colorspace)
library(jsonlite)
library(scales)

theme_mat = theme_bw() 
theme_mat$axis.title.x$size = 14
theme_mat$axis.text.x$size = 12
theme_mat$axis.title.y$size = 14
theme_mat$axis.text.y$size = 12
theme_mat$plot.title$size = 16
theme_mat$plot.title$hjust = 0.5
theme_mat$legend.text$size = 12
theme_mat$strip.text$size = 14


# Parameter Ranges
nusq <- seq(0.25,3.75,by=0.05)
lmsq <- c(seq(3,27,by=3),seq(30,3000,by=30))
sgsq <- c(2,seq(20,2000,by=20))
arsq <- seq(-1,1,by=0.05)

# Read JSON
prinfo <- fromJSON(txt = "../config/MIP_Africa_SpTm_Prior_Beta.json")
gpgrps <- c("mean","mainA","mainB","interact","noise")

for (j in seq(1,length(gpgrps))) {
  nustr <- (nusq - prinfo$prior[[gpgrps[j] ]][["gp_smoothness"]][["lim1"]]) /
           (prinfo$prior[[gpgrps[j] ]][["gp_smoothness"]][["lim2"]] - prinfo$prior[[gpgrps[j] ]][["gp_smoothness"]][["lim1"]])
  nudns <- dbeta(nustr, shape1=prinfo$prior[[gpgrps[j] ]][["gp_smoothness"]][["shape1"]], 
                        shape2=prinfo$prior[[gpgrps[j] ]][["gp_smoothness"]][["shape2"]], log=FALSE) /
           (prinfo$prior[[gpgrps[j] ]][["gp_smoothness"]][["lim2"]] - prinfo$prior[[gpgrps[j] ]][["gp_smoothness"]][["lim1"]])
  nutbl <- tibble(ParVal = nusq,Density = nudns,
                  ANOVAComp = rep(gpgrps[j],length(nusq)))
  if (j == 1) {
    nufnl <- nutbl
  }
  else {
    nufnl <- rbind(nufnl,nutbl)
  }
  
  
  lmdns <- dnorm(log(lmsq),mean = prinfo$prior[[ gpgrps[j]]][["gp_range"]][["mean"]],
                 sd = prinfo$prior[[ gpgrps[j]]][["gp_range"]][["stddev"]]) / lmsq
  lmtbl <- tibble(ParVal = lmsq,Density = lmdns,ANOVAComp = rep(gpgrps[j],length(lmsq)) )
  if (j == 1) {
    lmfnl <- lmtbl
  }
  else {
    lmfnl <- rbind(lmfnl,lmtbl)
  }
  
  sgdns <- 2.0 * dnorm(sgsq, mean = 0, 
                       sd = prinfo$prior[[ gpgrps[j]]][["gp_stddev"]][["stddev"]])
  sgtbl <- tibble(ParVal = sgsq,Density = sgdns,ANOVAComp = rep(gpgrps[j],length(sgsq)) )
  if (j == 1) {
    sgfnl <- sgtbl
  }
  else {
    sgfnl <- rbind(sgfnl,sgtbl)
  }
  
  if (gpgrps[j] == "noise") {
    arstr <- (arsq - prinfo$prior[[gpgrps[j] ]][["ar_correlation"]][["lim1"]]) /
      (prinfo$prior[[gpgrps[j] ]][["ar_correlation"]][["lim2"]] - prinfo$prior[[gpgrps[j] ]][["ar_correlation"]][["lim1"]])
    ardns <- dbeta(arstr, shape1=prinfo$prior[[gpgrps[j] ]][["ar_correlation"]][["shape1"]], 
                   shape2=prinfo$prior[[gpgrps[j] ]][["ar_correlation"]][["shape2"]], log=FALSE) /
          (prinfo$prior[[gpgrps[j] ]][["ar_correlation"]][["lim2"]] - prinfo$prior[[gpgrps[j] ]][["ar_correlation"]][["lim1"]])
    artbl <- tibble(ParVal = arsq,Density = ardns,
                  ANOVAComp = rep(gpgrps[j],length(arsq)))
  }
}

c5 = qualitative_hcl(5,palette = "dark3")
sgfnl <- sgfnl %>% 
  mutate(CompFact = recode(ANOVAComp, mean = 'mu',mainA = 'alpha * "," * beta', mainB = 'alpha * "," * beta',
                           interact = '(alpha * beta) * "," * epsilon',
                           noise = '(alpha * beta) * "," * epsilon'),
         CompOrd = recode(ANOVAComp,mean = 1, mainA = 2, mainB = 2,
                          interact = 3, noise = 3) )
sgfnl <- sgfnl %>% mutate(CompFact = fct_reorder(CompFact,CompOrd))
g1 <- ggplot(sgfnl,aes(x=ParVal,y=Density,group=ANOVAComp)) + 
      geom_line(aes(color=CompFact),size=1.2) + 
      scale_x_continuous("GP Standard Deviation") + ylab("Density") + 
      ggtitle("MIP GP Prior Standard Deviation") + 
      scale_colour_manual("Component",values=c5[c(1,2,4)],labels = parse_format()) + 
      theme_mat
pdf("MIP_Prior_Sigma.pdf",width = 8,height=6)
print(g1)
dev.off()

# Range
lmfnl <- lmfnl %>% 
  mutate(CompFact = recode(ANOVAComp, mean = 'mu * "," * alpha',mainA = 'mu * "," * alpha', 
                           mainB = 'beta * "," * (alpha * beta) * "," * epsilon',
                           interact = 'beta * "," * (alpha * beta) * "," * epsilon',
                           noise = 'beta * "," * (alpha * beta) * "," * epsilon'),
         CompOrd = recode(ANOVAComp,mean = 1, mainA = 1, mainB = 2,
                          interact = 2, noise = 2) )
lmfnl <- lmfnl %>% mutate(CompFact = fct_reorder(CompFact,CompOrd))
g2 <- ggplot(lmfnl,aes(x=ParVal,y=Density,group=ANOVAComp)) + 
  geom_line(aes(color=CompFact),size=1.2) + 
  scale_x_continuous("GP Range") + ylab("Density") + 
  ggtitle("MIP GP Prior Range") + 
  scale_colour_manual("Component",values=c5[c(1,3)],labels = parse_format()) + 
  theme_mat
pdf("MIP_Prior_Lambda.pdf",width = 8,height=6)
print(g2)
dev.off()

# Smoothness
g3 <- ggplot(nutbl,aes(x=ParVal,y=Density)) + 
  geom_line(color=c5[1],size=1.2) + 
  scale_x_continuous("GP Smoothness",limits=c(0,4)) + ylab("Density") + 
  ggtitle("MIP GP Prior Smoothness") + 
  theme_mat
pdf("MIP_Prior_Nu.pdf",width = 8,height=6)
print(g3)
dev.off()

# AR Correlation
artbl <- artbl %>% 
  mutate(CompFact = recode(ANOVAComp, 
                           noise = 'epsilon'),
         CompOrd = recode(ANOVAComp,noise = 1) )
artbl <- artbl %>% mutate(CompFact = fct_reorder(CompFact,CompOrd))
g4 <- ggplot(artbl,aes(x=ParVal,y=Density,group=ANOVAComp)) + 
  geom_line(aes(color=CompFact),size=1.2) + 
  scale_x_continuous("AR Correlation") + ylab("Density") + 
  ggtitle("MIP Error AR Correlation") + 
  scale_colour_manual("Component",values=c5[3],labels = parse_format()) + 
  theme_mat
pdf("MIP_Prior_Rho.pdf",width = 8,height=6)
print(g4)
dev.off()
