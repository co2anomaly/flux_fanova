# Correlation functions

library(GPvecchia)
library(tidyr)
library(dplyr)
library(ggplot2)
library(colorspace)
xsq = seq(1,500,by=10)
xdst = abs(outer(xsq,xsq,FUN = "-"))

theme_mat = theme_bw() 
theme_mat$axis.title.x$size = 14
theme_mat$axis.text.x$size = 12
theme_mat$axis.title.y$size = 14
theme_mat$axis.text.y$size = 12
theme_mat$plot.title$size = 16
theme_mat$plot.title$hjust = 0.5
theme_mat$legend.text$size = 12
theme_mat$strip.text$size = 14


lmsq = seq(0.5,3.0,by=0.5)
stnrng = 100

for (j in seq(1,length(lmsq))) {
    crgp = MaternFun(xdst,covparms = c(1,stnrng,lmsq[j]))
    gptbl = tibble(Dist = xdst[,1], MtrnCor = crgp[,1])
    gptbl = gptbl %>% mutate(Smooth = lmsq[j])
    if (j == 1) {
        gpfnl = gptbl      
    } else {
        gpfnl = rbind(gpfnl,gptbl)
    }
}
gpfnl = gpfnl %>% mutate(SmthFact = factor(Smooth))

c6 = qualitative_hcl(6,palette = "dark3")
g1 <- ggplot(gpfnl,aes(x=Dist,y=MtrnCor,group=SmthFact)) + 
  geom_line(aes(color=SmthFact),linewidth=1.2) + 
  scale_x_continuous("Distance") + ylab("Correlation") + 
  ggtitle("Matern Correlation, Range = 100") + 
  scale_colour_manual("Smoothness",values=c6) + 
  geom_hline(yintercept = exp(-1),linetype=2, linewidth=0.8) + 
  theme_mat
pdf("MaternCorr_Stein_Rng100.pdf",width = 8,height=6)
print(g1)
dev.off()

# Use "effective range" parameterization
effrng = 100

for (j in seq(1,length(lmsq))) {
  stnrng = effrng / sqrt(2.0 * lmsq[j])
  crgp = MaternFun(xdst,covparms = c(1,stnrng,lmsq[j]))
  gptbl = tibble(Dist = xdst[,1], MtrnCor = crgp[,1])
  gptbl = gptbl %>% mutate(Smooth = lmsq[j])
  if (j == 1) {
    gpfnl = gptbl      
  } else {
    gpfnl = rbind(gpfnl,gptbl)
  }
}
gpfnl = gpfnl %>% mutate(SmthFact = factor(Smooth))

c6 = qualitative_hcl(6,palette = "dark3")
g2 <- ggplot(gpfnl,aes(x=Dist,y=MtrnCor,group=SmthFact)) + 
  geom_line(aes(color=SmthFact),linewidth=1.2) + 
  scale_x_continuous("Distance") + ylab("Correlation") + 
  ggtitle("Matern Correlation, Effective Range = 100") + 
  scale_colour_manual("Smoothness",values=c6) + 
  geom_hline(yintercept = exp(-1),linetype=2, linewidth=0.8) + 
  theme_mat
pdf("MaternCorr_EffRng_Rng100.pdf",width = 8,height=6)
print(g2)
dev.off()
