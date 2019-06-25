## -- C-STAD: plotting results ------------ ##
## Date: 24/06/2019
## Author: Ugofilippo Basellini
## Comments: 
##       - FIGURE A2: Residuals 
##
## ------------------------------------------ ##

rm(list = ls())
library(lattice)

## read data
cou <- "SWE"
sex <- "F"
n.sim <- 1600
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Results")
name <- paste0(cou,"_",sex,"_",n.sim,"simul_results.Rdata")
load(name)
ResCSTAD_SWE <- DevRes

## load SWE results
cou <- "DNK"
name <- paste0(cou,"_",sex,"_",n.sim,"simul_results.Rdata")
load(name)
ResCSTAD_DNK <- DevRes

setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Functions")
source("CohortToPeriod.R")

ResCSTAD_SWE <- PeriodMat_FromCohort(ages = xo,
                                     cohorts = as.numeric(colnames(ResCSTAD_SWE)),
                                     CohortMat = ResCSTAD_SWE,t.end = 2015)
ResCSTAD_DNK <- PeriodMat_FromCohort(ages = xo,
                                     cohorts = as.numeric(colnames(ResCSTAD_DNK)),
                                     CohortMat = ResCSTAD_DNK,t.end = 2015)

## plot
n.breaks <- 12
res.min <- max(abs(min(ResCSTAD_DNK,na.rm = T)),abs(min(ResCSTAD_SWE,na.rm = T)))
res.max <- max(abs(max(ResCSTAD_DNK,na.rm = T)),abs(max(ResCSTAD_SWE,na.rm = T)))
if(res.min >= res.max){
  my.breaks <- seq(-res.min,res.min,
                   length.out = n.breaks)
}else{
  my.breaks <- seq(-res.max,res.max,
                   length.out = n.breaks)
}
my.breaks <- seq(-5.8,5.8,length.out = n.breaks)
my.cols <- brewer.pal((n.breaks-1), "RdBu")

y.cstad <- as.numeric(colnames(ResCSTAD_SWE))
list.plot <- list(x = xo, y = y.cstad, type = c("SWE","DNK"))
grid.plot <- expand.grid(list.plot)
grid.plot$Z <- c(c(ResCSTAD_SWE),c(ResCSTAD_DNK))
cex.lab <- 1.5
(A3 <- levelplot(Z ~ y * x | type, grid.plot, layout = c(2, 1), at = my.breaks, 
          col.regions = my.cols, colorkey = list(col = my.cols),
          xlab=list("Year",cex=cex.lab),
          ylab=list("Age",cex=cex.lab),
          par.settings=list(panel.background=list(col="grey")),
          panel = function(...){
            panel.levelplot(...)
            panel.text(1930,60,labels=expression(c[1]),cex=1.3)
            panel.text(1975,60,labels=expression(c[2]),cex=1.3)
            panel.text(1995,60,labels=expression(c[3]),cex=1.3)
            panel.abline(a=c(-c1[1]+1),b=1,lty=2,lwd=2)
            panel.abline(a=c(-c2[1]+1),b=1,lty=2,lwd=2)
            panel.abline(a=c(-c3[1]+1),b=1,lty=2,lwd=2)
          })
)

setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/Paper/Figures")
pdf("FA2.pdf",width = 10,height=5.5)
print(A3)
dev.off()
