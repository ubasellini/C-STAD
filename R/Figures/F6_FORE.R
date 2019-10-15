## -- C-STAD: plotting results ------------ ##
## Date: 24/06/2019
## Author: Ugofilippo Basellini
## Comments: 
##       - FIGURE 5: Forecasts e40, g40
##
## ------------------------------------------ ##

rm(list = ls())

## read data
cou <- "SWE"
sex <- "F"
n.sim <- 1600
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Results")
name <- paste0(cou,"_",sex,"_",n.sim,"simul_results.Rdata")
load(name)

## CI
lev.p <- 0.8
e40.actSWE <- e40.act
g40.actSWE <- g40.act
e40_2DSWE <- e40_2D
g40_2DSWE <- g40_2D
E40_cstad_MEAN_SWE <- apply(E40sim,1,median)
E40_cstad_UP_SWE <- apply(E40sim,1,quantile,prob=1- (1-lev.p)/2)
E40_cstad_LOW_SWE <- apply(E40sim,1,quantile,prob=(1-lev.p)/2)
G40_cstad_MEAN_SWE <- apply(G40sim,1,median)
G40_cstad_UP_SWE <- apply(G40sim,1,quantile,prob=1- (1-lev.p)/2)
G40_cstad_LOW_SWE <- apply(G40sim,1,quantile,prob=(1-lev.p)/2)
c1SWE <- c1
c2SWE <- c2
c3SWE <- c3
couSWE <- "Sweden"

## load DNK results
cou <- "DNK"
name <- paste0(cou,"_",sex,"_",n.sim,"simul_results.Rdata")
load(name)
## CI
lev.p <- 0.8
e40.actDNK <- e40.act
g40.actDNK <- g40.act
e40_2DDNK <- e40_2D
g40_2DDNK <- g40_2D
E40_cstad_MEAN_DNK <- apply(E40sim,1,median)
E40_cstad_UP_DNK <- apply(E40sim,1,quantile,prob=1- (1-lev.p)/2)
E40_cstad_LOW_DNK <- apply(E40sim,1,quantile,prob=(1-lev.p)/2)
G40_cstad_MEAN_DNK <- apply(G40sim,1,median)
G40_cstad_UP_DNK <- apply(G40sim,1,quantile,prob=1- (1-lev.p)/2)
G40_cstad_LOW_DNK <- apply(G40sim,1,quantile,prob=(1-lev.p)/2)
c1DNK <- c1
c2DNK <- c2
c3DNK <- c3
couDNK <- "Denmark"

## plotting
require(RColorBrewer)
my.cols <- brewer.pal(8, "Set1")
my.colsT <- adjustcolor(my.cols, alpha=0.3)
cex.axis <- 0.75
cex4 <- 1.25
cex16 <- 0.75
cex18 <- cex16+0.2
lwd.x <- 2.5
cex.coh <- 1.5
cex.cou <- 1.3
cex.leg <- 1.2
cex.x.axis <- 1.1
cex.y.axis <- 1.1
cex.y.lab <-  1.75
cex.main <- 1.75
my.pspline <- brewer.pal(n=9, name = 'Oranges')[7]

c1_mid <- mean(c(c1DNK[length(c1DNK)],c1SWE[length(c1SWE)]))
c2_mid <- mean(c(c2DNK[length(c2DNK)],c2SWE[length(c2SWE)]))

setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/Paper/Figures")
pdf("F6.pdf",width = 10,height=5.5)
## set overall margins
par(mfrow = c(1,2),
    oma = c(1.25,0.7,0.1,0.25),
    mar = c(1.75,2.35,1.5,0.9))
## bottom, left, top, right
xlim <- range(y)
ylimE40fit <- range(e40.actDNK,e40.actSWE,E40_cstad_UP_DNK,E40_cstad_UP_SWE)
plot(c1,e40.actDNK,t="n",ylim=ylimE40fit,
     xlim=xlim,cex.main=1,axes=F)
title(main = expression(paste("Observed and completed ",e[40])),
      cex.main=cex.main)
axis(1,cex.axis=cex.x.axis,padj = -0.5);
axis(2,las=2,cex.axis=cex.y.axis,hadj = 0.75)
mtext(text=expression(e[40]), side=2, las=3,
      line=1.5, cex=cex.y.lab)
grid();box()
abline(v=c1_mid,lty=2)
abline(v=c2_mid,lty=2)
points(c1,e40.actDNK,t="p",pch=16,col=my.cols[2],cex=cex16)
text(1890,32,couDNK,col=my.cols[2],cex=cex.cou)
points(c1SWE,e40.actSWE,t="p",pch=16,col=my.cols[3],cex=cex16)
text(1860,35,couSWE,col=my.cols[3],cex=cex.cou)
lines(y[!y%in%c1SWE],e40_2DSWE[!y%in%c1SWE],col=my.pspline,lwd=1.5,lty=1)
lines(y[!y%in%c1DNK],e40_2DDNK[!y%in%c1DNK],col=my.pspline,lwd=1.5,lty=1)
xx <- c(y[!y%in%c1DNK],rev(y[!y%in%c1DNK]))
yy <- c(E40_cstad_LOW_DNK[!y%in%c1],rev(E40_cstad_UP_DNK[!y%in%c1]))
polygon(xx, yy, border = my.cols[2], col=my.colsT[2])
xx <- c(y[!y%in%c1SWE],rev(y[!y%in%c1SWE]))
yy <- c(E40_cstad_LOW_SWE[!y%in%c1SWE],rev(E40_cstad_UP_SWE[!y%in%c1SWE]))
polygon(xx, yy, border = my.cols[3], col=my.colsT[3])
lines(y[!y%in%c1DNK],E40_cstad_MEAN_DNK[!y%in%c1DNK],col=my.cols[2],lwd=2.5)
lines(y[!y%in%c1SWE],E40_cstad_MEAN_SWE[!y%in%c1SWE],col=my.cols[3],lwd=2.5)
text(1865,45,expression(c[1]),col=1,cex=cex.coh)
text(1915,45,expression(c[2]),col=1,cex=cex.coh)
text(1955,40,expression(c[3]),col=1,cex=cex.coh)
legend("topleft",c("observed","C-STAD","2D P-spline"),pch=c(16,NA,NA),col=c(1,1,my.pspline),pt.cex=0.8,
       bty="n",lwd=c(1,2,1.5),cex=cex.leg,lty=c(NA,1,1),inset = 0.01)
legend("topleft",c("","",""),pch=c(NA,15,NA),col=adjustcolor(1,alpha.f = 0.3),
       bty="n",lwd=c(NA,NA,NA),cex=cex.leg,lty=c(NA,1,NA),inset = 0.01)

## right: G40
ylimG40fit <- range(g40.actDNK,g40.actSWE,G40_cstad_LOW_DNK,G40_cstad_LOW_SWE)
plot(c1,g40.actDNK,t="n",ylim=ylimG40fit,
     xlim=xlim,cex.main=1,axes=F)
title(main = expression(paste("Observed and completed ",G[40])),
      cex.main=cex.main)
axis(1,cex.axis=cex.x.axis,padj = -0.5);
axis(2,las=2,cex.axis=cex.y.axis,hadj = 0.75,
     at=seq(0.12,0.24,0.02),labels = c("12","14","16","18","20","22","24"))
mtext(text=expression(G[40]), side=2, las=3,
      line=1.5, cex=cex.y.lab)
grid();box()
abline(v=c1_mid,lty=2)
abline(v=c2_mid,lty=2)
points(c1,g40.actDNK,t="p",pch=16,col=my.cols[2],cex=cex16)
text(1890,0.23,couDNK,col=my.cols[2],cex=cex.cou)
points(c1SWE,g40.actSWE,t="p",pch=16,col=my.cols[3],cex=cex16)
text(1868,0.20,couSWE,col=my.cols[3],cex=cex.cou)
lines(y[!y%in%c1DNK],g40_2DDNK[!y%in%c1DNK],col=my.pspline,lwd=1.5,lty=1)
lines(y[!y%in%c1SWE],g40_2DSWE[!y%in%c1SWE],col=my.pspline,lwd=1.5,lty=1)
xx <- c(y[!y%in%c1DNK],rev(y[!y%in%c1DNK]))
yy <- c(G40_cstad_LOW_DNK[!y%in%c1],rev(G40_cstad_UP_DNK[!y%in%c1]))
polygon(xx, yy, border = my.cols[2], col=my.colsT[2])
xx <- c(y[!y%in%c1SWE],rev(y[!y%in%c1SWE]))
yy <- c(G40_cstad_LOW_SWE[!y%in%c1SWE],rev(G40_cstad_UP_SWE[!y%in%c1SWE]))
polygon(xx, yy, border = my.cols[3], col=my.colsT[3])

lines(y[!y%in%c1DNK],G40_cstad_MEAN_DNK[!y%in%c1],col=my.cols[2],lwd=2.5)
lines(y[!y%in%c1SWE],G40_cstad_MEAN_SWE[!y%in%c1SWE],col=my.cols[3],lwd=2.5)
text(1865,0.145,expression(c[1]),col=1,cex=cex.coh)
text(1915,0.145,expression(c[2]),col=1,cex=cex.coh)
text(1955,0.2,expression(c[3]),col=1,cex=cex.coh)

cex.x.lab <- cex.y.lab
title(xlab = "Cohort",cex.lab=cex.x.lab,
      outer = TRUE, line = 0.3)

dev.off()
