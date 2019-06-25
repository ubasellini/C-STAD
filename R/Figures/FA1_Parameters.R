## -- C-STAD: plotting results ------------ ##
## Date: 24/06/2019
## Author: Ugofilippo Basellini
## Comments: 
##       - FIGURE A1: Parameters 
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

cor(diff(s_cstad[1:(n1+n2)]),diff(bU_cstad[1:(n1+n2)]))

couSWE <- "Sweden"
c1SWE <- c1
c2SWE <- c2
c3SWE <- c3
n1SWE <- n1
n2SWE <- n2
n3SWE <- n3
lev.p <- 0.8
S_cstad_MEAN_SWE <- apply(Ssim,1,median)
S_cstad_UP_SWE <- apply(Ssim,1,quantile,prob=1- (1-lev.p)/2)
S_cstad_LOW_SWE <- apply(Ssim,1,quantile,prob=(1-lev.p)/2)
BL_cstad_MEAN_SWE <- apply(BLsim,1,median)
BL_cstad_UP_SWE <- apply(BLsim,1,quantile,prob=1- (1-lev.p)/2)
BL_cstad_LOW_SWE <- apply(BLsim,1,quantile,prob=(1-lev.p)/2)
BU_cstad_MEAN_SWE <- apply(BUsim,1,median)
BU_cstad_UP_SWE <- apply(BUsim,1,quantile,prob=1- (1-lev.p)/2)
BU_cstad_LOW_SWE <- apply(BUsim,1,quantile,prob=(1-lev.p)/2)
CL_cstad_MEAN_SWE <- apply(CLsim,1,median)
CL_cstad_UP_SWE <- apply(CLsim,1,quantile,prob=1- (1-lev.p)/2)
CL_cstad_LOW_SWE <- apply(CLsim,1,quantile,prob=(1-lev.p)/2)
DL_cstad_MEAN_SWE <- apply(DLsim,1,median)
DL_cstad_UP_SWE <- apply(DLsim,1,quantile,prob=1- (1-lev.p)/2)
DL_cstad_LOW_SWE <- apply(DLsim,1,quantile,prob=(1-lev.p)/2)

## load DNK results
cou <- "DNK"
name <- paste0(cou,"_",sex,"_",n.sim,"simul_results.Rdata")
load(name)

cor(diff(s_cstad[1:(n1+n2)]),diff(bU_cstad[1:(n1+n2)]))

couDNK <- "Denmark"
c1DNK <- c1
c2DNK <- c2
c3DNK <- c3
n1DNK <- n1
n2DNK <- n2
n3DNK <- n3
lev.p <- 0.8
S_cstad_MEAN_DNK <- apply(Ssim,1,median)
S_cstad_UP_DNK <- apply(Ssim,1,quantile,prob=1- (1-lev.p)/2)
S_cstad_LOW_DNK <- apply(Ssim,1,quantile,prob=(1-lev.p)/2)
BL_cstad_MEAN_DNK <- apply(BLsim,1,median)
BL_cstad_UP_DNK <- apply(BLsim,1,quantile,prob=1- (1-lev.p)/2)
BL_cstad_LOW_DNK <- apply(BLsim,1,quantile,prob=(1-lev.p)/2)
BU_cstad_MEAN_DNK <- apply(BUsim,1,median)
BU_cstad_UP_DNK <- apply(BUsim,1,quantile,prob=1- (1-lev.p)/2)
BU_cstad_LOW_DNK <- apply(BUsim,1,quantile,prob=(1-lev.p)/2)
CL_cstad_MEAN_DNK <- apply(CLsim,1,median)
CL_cstad_UP_DNK <- apply(CLsim,1,quantile,prob=1- (1-lev.p)/2)
CL_cstad_LOW_DNK <- apply(CLsim,1,quantile,prob=(1-lev.p)/2)
DL_cstad_MEAN_DNK <- apply(DLsim,1,median)
DL_cstad_UP_DNK <- apply(DLsim,1,quantile,prob=1- (1-lev.p)/2)
DL_cstad_LOW_DNK <- apply(DLsim,1,quantile,prob=(1-lev.p)/2)

## plotting
require(RColorBrewer)
my.cols <- brewer.pal(8, "Set1")
my.colsT <- adjustcolor(my.cols, alpha=0.3)
cex.axis <- 0.75
cex4 <- 1.25
cex16 <- 0.8
cex18 <- cex16+0.2
lwd.x <- 2
cex.coh <- 1.5
cex.leg <- 1
lwd.ln <- 2.5

setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/Paper/Figures")
pdf("FA1.pdf",width = 10,height=5.5)
## set overall margins
par(mfrow = c(2,3),
    oma = c(2.25,0.4,0.1,0.25),
    mar = c(1.25,2.3,1.2,0.3))
## bottom, left, top, right

## S
xlim <- range(y)
ylimS <- range(S_cstad_UP_DNK,S_cstad_UP_SWE,S_cstad_LOW_DNK,S_cstad_LOW_SWE)
plot(c1DNK,S_cstad_UP_DNK[1:n1DNK],t="n",ylim=ylimS,xlim=xlim,
     cex.main=1,axes=F)
title(main = expression(s))
axis(1,cex.axis=cex.axis,padj = -1);
axis(2,las=2,cex.axis=cex.axis,hadj = 0.75)
grid();box()
abline(v=c1[length(c1)],lty=2)
abline(v=mean(c3SWE[1],c3DNK[1]),lty=2)
xx <- c(y, rev(y))
yy <- c(S_cstad_LOW_DNK, rev(S_cstad_UP_DNK))
polygon(xx, yy, border = my.cols[2], col=my.colsT[2])
yy <- c(S_cstad_LOW_SWE, rev(S_cstad_UP_SWE))
polygon(xx, yy, border = my.cols[3], col=my.colsT[3])
lines(y,S_cstad_MEAN_DNK,col=my.cols[2],lwd=lwd.ln)
lines(y,S_cstad_MEAN_SWE,col=my.cols[3],lwd=lwd.ln)
text(1865,12,expression(c[1]),col=1,cex=cex.coh)
text(1915,3,expression(c[2]),col=1,cex=cex.coh)
text(1955,3,expression(c[3]),col=1,cex=cex.coh)
legend("topleft",c("DNK","SWE"),inset=0.01,
       col = my.cols[2:3],
       bty="n",lwd=2,lty=c(1,1),cex = cex.leg)
legend("topleft", inset=0.01,
       c("",""),
       col=my.colsT[2:3], lwd=3,
       lty = c(NA,NA),
       pch=c(15,15), pt.cex=2.5,
       bty="n",cex = cex.leg)

## bU
ylimBU <- range(BU_cstad_UP_DNK,BU_cstad_UP_SWE,BU_cstad_LOW_DNK,BU_cstad_LOW_SWE)
plot(c1DNK,BU_cstad_UP_DNK[1:n1DNK],t="n",ylim=ylimBU,xlim=xlim,
     cex.main=1,axes=F)
title(main = expression(b[U]))
axis(1,cex.axis=cex.axis,padj = -1);
axis(2,las=2,cex.axis=cex.axis,hadj = 0.75)
grid();box()
abline(v=c1[length(c1)],lty=2)
abline(v=mean(c3SWE[1],c3DNK[1]),lty=2)
yy <- c(BU_cstad_LOW_DNK, rev(BU_cstad_UP_DNK))
polygon(xx, yy, border = my.cols[2], col=my.colsT[2])
yy <- c(BU_cstad_LOW_SWE, rev(BU_cstad_UP_SWE))
polygon(xx, yy, border = my.cols[3], col=my.colsT[3])
lines(y,BU_cstad_MEAN_DNK,col=my.cols[2],lwd=lwd.ln)
lines(y,BU_cstad_MEAN_SWE,col=my.cols[3],lwd=lwd.ln)
text(1865,1.2,expression(c[1]),col=1,cex=cex.coh)
text(1915,0.85,expression(c[2]),col=1,cex=cex.coh)
text(1955,0.85,expression(c[3]),col=1,cex=cex.coh)

## VOID
ylimBU <- range(BU_cstad_UP_DNK,BU_cstad_UP_SWE,BU_cstad_LOW_DNK,BU_cstad_LOW_SWE)
plot(c1DNK,BU_cstad_UP_DNK[1:n1DNK],t="n",ylim=ylimBU,xlim=xlim,
     cex.main=1,axes=F)

## bL
ylimBL <- range(BL_cstad_UP_DNK,BL_cstad_UP_SWE,BL_cstad_LOW_DNK,BL_cstad_LOW_SWE)
plot(c1,BL_cstad_UP_DNK[1:n1DNK],t="n",ylim=ylimBL,xlim=xlim,
     cex.main=1,axes=F)
title(main = expression(b[L]))
axis(1,cex.axis=cex.axis,padj = -1);
axis(2,las=2,cex.axis=cex.axis,hadj = 0.75)
# mtext(text=expression(b[L]), side=2, las=3,
#       line=1.25, cex=1.25)
# mtext(text="Cohort", side=1, las=1,
#       line=1.4, cex=1.25)
grid();box()
abline(v=c1[length(c1)],lty=2)
abline(v=mean(c3SWE[1],c3DNK[1]),lty=2)
yy <- c(BL_cstad_LOW_DNK, rev(BL_cstad_UP_DNK))
polygon(xx, yy, border = my.cols[2], col=my.colsT[2])
yy <- c(BL_cstad_LOW_SWE, rev(BL_cstad_UP_SWE))
polygon(xx, yy, border = my.cols[3], col=my.colsT[3])
lines(y,BL_cstad_MEAN_DNK,col=my.cols[2],lwd=lwd.ln)
lines(y,BL_cstad_MEAN_SWE,col=my.cols[3],lwd=lwd.ln)
text(1865,1.8,expression(c[1]),col=1,cex=cex.coh)
text(1915,0.3,expression(c[2]),col=1,cex=cex.coh)
text(1955,0.3,expression(c[3]),col=1,cex=cex.coh)

## cL
ylimCL <- range(CL_cstad_UP_DNK,CL_cstad_UP_SWE,CL_cstad_LOW_DNK,CL_cstad_LOW_SWE)
plot(c1,CL_cstad_UP_DNK[1:n1DNK],t="n",ylim=ylimCL,xlim=xlim,
     cex.main=1,axes=F)
title(main = expression(c[L]))
axis(1,cex.axis=cex.axis,padj = -1);
axis(2,las=2,cex.axis=cex.axis,hadj = 0.75)
grid();box()
abline(v=c1[length(c1)],lty=2)
abline(v=mean(c3SWE[1],c3DNK[1]),lty=2)
yy <- c(CL_cstad_LOW_DNK, rev(CL_cstad_UP_DNK))
polygon(xx, yy, border = my.cols[2], col=my.colsT[2])
yy <- c(CL_cstad_LOW_SWE, rev(CL_cstad_UP_SWE))
polygon(xx, yy, border = my.cols[3], col=my.colsT[3])
lines(y,CL_cstad_MEAN_DNK,col=my.cols[2],lwd=lwd.ln)
lines(y,CL_cstad_MEAN_SWE,col=my.cols[3],lwd=lwd.ln)
text(1865,0.015,expression(c[1]),col=1,cex=cex.coh)
text(1915,-0.03,expression(c[2]),col=1,cex=cex.coh)
text(1955,-0.03,expression(c[3]),col=1,cex=cex.coh)

## dL
ylimDL <- range(DL_cstad_UP_DNK,DL_cstad_UP_SWE,DL_cstad_LOW_DNK,DL_cstad_LOW_SWE)
plot(c1,DL_cstad_UP_DNK[1:n1DNK],t="n",ylim=ylimDL,xlim=xlim,
     cex.main=1,axes=F)
title(main = expression(d[L]))
axis(1,cex.axis=cex.axis,padj = -1);
axis(2,las=2,cex.axis=cex.axis,hadj = 0.75)
grid();box()
abline(v=c1[length(c1)],lty=2)
abline(v=mean(c3SWE[1],c3DNK[1]),lty=2)
yy <- c(DL_cstad_LOW_DNK, rev(DL_cstad_UP_DNK))
polygon(xx, yy, border = my.cols[2], col=my.colsT[2])
yy <- c(DL_cstad_LOW_SWE, rev(DL_cstad_UP_SWE))
polygon(xx, yy, border = my.cols[3], col=my.colsT[3])
lines(y,DL_cstad_MEAN_DNK,col=my.cols[2],lwd=lwd.ln)
lines(y,DL_cstad_MEAN_SWE,col=my.cols[3],lwd=lwd.ln)
text(1865,2e-4,expression(c[1]),col=1,cex=cex.coh)
text(1915,-5e-4,expression(c[2]),col=1,cex=cex.coh)
text(1955,-5e-4,expression(c[3]),col=1,cex=cex.coh)

cex.x.lab <- 2
title(xlab = "Cohort",cex.lab=cex.x.lab,
      outer = TRUE, line = 1)

dev.off()
