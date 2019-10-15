## -- C-STAD: plotting results ------------ ##
## Date: 24/06/2019
## Author: Ugofilippo Basellini
## Comments: 
##       - FIGURE 6a: Rates over time
##
## ------------------------------------------ ##

## cleaning the workspace
rm(list=ls(all=TRUE))
library(RColorBrewer)
library(fields)
library(colorspace)

##---- FIGURE: PARAMETERS -------

## read data
cou <- "SWE"
sex <- "F"
n.sim <- 1600
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Results")
name <- paste0(cou,"_",sex,"_",n.sim,"simul_results.Rdata")
load(name)
MXactSWE <- exp(LMX.act)
c1SWE <- c1
c2SWE <- c2
c3SWE <- c3

## CI for bootstrap
lev.p <- 0.8
ETA.MEANswe <- ETA.LOWswe <- ETA.UPswe <- matrix(NA,dim(ETAS)[1],dim(ETAS)[2])
i <- j <- 1
for (i in 1:dim(ETAS)[1]){
  for (j in 1:dim(ETAS)[2]){
    ETA.MEANswe[i,j] <- median(ETAS[i,j,,],na.rm=T)
    ETA.UPswe[i,j] <- quantile(ETAS[i,j,,],prob=1- (1-lev.p)/2,na.rm=T)
    ETA.LOWswe[i,j] <- quantile(ETAS[i,j,,],prob=(1-lev.p)/2,na.rm=T)
  }
}
# matplot(xo,ETA.MEANswe[1:length(xo),],t="l",col = rainbow_hcl(dim(ETAS)[2]),lty=1)
# matplot(40:120,ETA.LOWswe,t="l",col = rainbow_hcl(dim(ETAS)[2]),lty=1)

## read data
cou <- "DNK"
name <- paste0(cou,"_",sex,"_",n.sim,"simul_results.Rdata")
load(name)
MXactDNK <- exp(LMX.act)
c1DNK <- c1
c2DNK <- c2
c3DNK <- c3
## CI for bootstrap
lev.p <- 0.8
ETA.MEANdnk <- ETA.LOWdnk <- ETA.UPdnk <- matrix(NA,dim(ETAS)[1],dim(ETAS)[2])
i <- j <- 1
for (i in 1:dim(ETAS)[1]){
  for (j in 1:dim(ETAS)[2]){
    ETA.MEANdnk[i,j] <- median(ETAS[i,j,,],na.rm=T)
    ETA.UPdnk[i,j] <- quantile(ETAS[i,j,,],prob=1- (1-lev.p)/2,na.rm=T)
    ETA.LOWdnk[i,j] <- quantile(ETAS[i,j,,],prob=(1-lev.p)/2,na.rm=T)
  }
}

## colors
col.stad <- brewer.pal(n = 8, name = 'Blues')[8]
col.stadT <- adjustcolor(col.stad, alpha=0.3)
my.orange <- brewer.pal(n=9, name = 'Oranges')[6]
my.green <- brewer.pal(n=9, name = 'Greens')[6]
my.purple <- brewer.pal(n=9, name = 'Purples')[7]
my.cols <- c(col.stad,my.orange,my.green,my.purple)
my.colsT <- adjustcolor(my.cols, alpha=0.3)
col.obs <- "grey70"
cex.x.axis <- 1.1
cex.y.axis <-  1.1
cex.x.lab <- 2
cex.y.lab <-  1.75
cex.leg <- 1.15
cex.coh <- 1.5
cex.title <- 1.7 
cex.obs <- 0.85
cex.age <- 1
lwd.pt <- 0.9
lwd.mean <- 2.5

##--- RATES ------
ylim <- range(t(MXactSWE),finite=T)
ylim[1] <- ylim[1]+0.0005
ylim[2] <- ylim[2]-1.9
ylim

## SAVE
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/Paper/Figures")
pdf("F7a.pdf",width = 10,height = 5.5)
par(mfrow = c(1,2),
    oma = c(1.25,1.2,0.1,0.25),
    mar = c(1.75,2.3,1.2,0.3))

## SWE
## bottom, left, top, right
matplot(y,t(MXactSWE),pch=16,log="y",
     ylim=ylim,xlab="",ylab="",t="n",axes=F)
axis(1,cex.axis=cex.x.axis,padj = -0.5)
axis(2,las=3,cex.axis=cex.y.axis,at=c(1e-3,1e-2,1e-1,1e0),
     labels = c("0.001","0.01","0.1","1"))
grid();box()
# mtext("Cohort",side=1,cex=cex.x.lab,line=2.4)
mtext("Log-mortality",side=2,cex=cex.y.lab,line=2,las=3)
title(main="Sweden", cex.main=cex.title)
j <- 40
whi.age <- which(xo==j)
points(y,t(MXactSWE[whi.age,]),pch=1,col=my.cols[1],cex=cex.obs,lwd=lwd.pt)
xx <- c(y,rev(y))
yy <- c(t(exp(ETA.LOWswe[whi.age,])),rev(t(exp(ETA.UPswe[whi.age,]))))
polygon(xx, yy, border = my.cols[1], col=my.colsT[1])
lines(y,t(exp(ETA.MEANswe[whi.age,])),t="l",lwd=lwd.mean,col=my.cols[1])
text(1848,0.005,paste("age =",j),lwd=1,cex=cex.age,font=2,col=my.cols[1])

j <- 60
whi.age <- which(xo==j)
points(y,t(MXactSWE[whi.age,]),pch=1,col=my.cols[2],cex=cex.obs,lwd=lwd.pt)
yy <- c(t(exp(ETA.LOWswe[whi.age,])),rev(t(exp(ETA.UPswe[whi.age,]))))
polygon(xx, yy, border = my.cols[2], col=my.colsT[2])
lines(y,t(exp(ETA.MEANswe[whi.age,])),t="l",lwd=lwd.mean,col=my.cols[2])
text(1868,0.029,paste("age =",j),lwd=1,cex=cex.age,font=2,col = my.cols[2])

j <- 80
whi.age <- which(xo==j)
points(y,t(MXactSWE[whi.age,]),pch=1,col=my.cols[3],cex=cex.obs,lwd=lwd.pt)
yy <- c(t(exp(ETA.LOWswe[whi.age,])),rev(t(exp(ETA.UPswe[whi.age,]))))
polygon(xx, yy, border = my.cols[3], col=my.colsT[3])
lines(y,t(exp(ETA.MEANswe[whi.age,])),t="l",lwd=lwd.mean,col=my.cols[3])
text(1848,0.07,paste("age =",j),lwd=1,cex=cex.age,font=2,col=my.cols[3])

j <- 100
whi.age <- which(xo==j)
points(y,t(MXactSWE[whi.age,]),pch=1,col=my.cols[4],cex=cex.obs,lwd=lwd.pt)
yy <- c(t(exp(ETA.LOWswe[whi.age,])),rev(t(exp(ETA.UPswe[whi.age,]))))
polygon(xx, yy, border = my.cols[4], col=my.colsT[4])
lines(y,t(exp(ETA.MEANswe[whi.age,])),t="l",lwd=lwd.mean,col=my.cols[4])
text(1868,0.23,paste("age =",j),lwd=1,cex=cex.age,font=2,col=my.cols[4])
# lines(cohorts,t(exp(LMX.smo2D[whi.smo,])),t="l",lwd=2,col=my.orange[6])

abline(v=c1SWE[length(c1SWE)],lty=2)
abline(v=c2SWE[length(c2SWE)],lty=2)
text(1868,0.002,expression(c[1]),lwd=1,cex=cex.coh,font=2)
text(1917,0.0006,expression(c[2]),lwd=1,cex=cex.coh,font=2)
text(1942,0.0006,expression(c[3]),lwd=1,cex=cex.coh,font=2)

legend("bottomleft",c("observed","C-STAD"),pch=c(1,NA),col=1,pt.cex=cex.obs,
       bty="n",lwd=c(lwd.pt,lwd.mean),cex=cex.leg,lty=c(NA,1),inset = 0.01)
legend("bottomleft",c("",""),pch=c(NA,15),col=adjustcolor(1,alpha.f = 0.3),
       bty="n",lwd=c(NA,NA),cex=cex.leg,lty=c(NA,1),inset = 0.01)


## DNK
## bottom, left, top, right
matplot(y,t(MXactDNK),pch=16,log="y",
        ylim=ylim,xlab="",ylab="",t="n",axes=F)
axis(1,cex.axis=cex.x.axis,padj = -0.5)
axis(2,las=3,cex.axis=cex.y.axis,at=c(1e-3,1e-2,1e-1,1e0),
     labels = c("0.001","0.01","0.1","1"))
grid();box()
title(main="Denmark", cex.main=cex.title)
j <- 40
whi.age <- which(xo==j)
points(y,t(MXactDNK[whi.age,]),pch=1,col=my.cols[1],cex=cex.obs,lwd=lwd.pt)
xx <- c(y,rev(y))
yy <- c(t(exp(ETA.LOWdnk[whi.age,])),rev(t(exp(ETA.UPdnk[whi.age,]))))
polygon(xx, yy, border = my.cols[1], col=my.colsT[1])
lines(y,t(exp(ETA.MEANdnk[whi.age,])),t="l",lwd=lwd.mean,col=my.cols[1])
text(1848,0.005,paste("age =",j),lwd=1,cex=cex.age,font=2,col=my.cols[1])

j <- 60
whi.age <- which(xo==j)
points(y,t(MXactDNK[whi.age,]),pch=1,col=my.cols[2],cex=cex.obs,lwd=lwd.pt)
yy <- c(t(exp(ETA.LOWdnk[whi.age,])),rev(t(exp(ETA.UPdnk[whi.age,]))))
polygon(xx, yy, border = my.cols[2], col=my.colsT[2])
lines(y,t(exp(ETA.MEANdnk[whi.age,])),t="l",lwd=lwd.mean,col=my.cols[2])
text(1868,0.029,paste("age =",j),lwd=1,cex=cex.age,font=2,col = my.cols[2])

j <- 80
whi.age <- which(xo==j)
points(y,t(MXactDNK[whi.age,]),pch=1,col=my.cols[3],cex=cex.obs,lwd=lwd.pt)
yy <- c(t(exp(ETA.LOWdnk[whi.age,])),rev(t(exp(ETA.UPdnk[whi.age,]))))
polygon(xx, yy, border = my.cols[3], col=my.colsT[3])
lines(y,t(exp(ETA.MEANdnk[whi.age,])),t="l",lwd=lwd.mean,col=my.cols[3])
text(1848,0.07,paste("age =",j),lwd=1,cex=cex.age,font=2,col=my.cols[3])

j <- 100
whi.age <- which(xo==j)
points(y,t(MXactDNK[whi.age,]),pch=1,col=my.cols[4],cex=cex.obs,lwd=lwd.pt)
yy <- c(t(exp(ETA.LOWdnk[whi.age,])),rev(t(exp(ETA.UPdnk[whi.age,]))))
polygon(xx, yy, border = my.cols[4], col=my.colsT[4])
lines(y,t(exp(ETA.MEANdnk[whi.age,])),t="l",lwd=lwd.mean,col=my.cols[4])
text(1868,0.2,paste("age =",j),lwd=1,cex=cex.age,font=2,col=my.cols[4])

abline(v=c1DNK[length(c1DNK)],lty=2)
abline(v=c2DNK[length(c2DNK)],lty=2)
text(1868,0.002,expression(c[1]),lwd=1,cex=cex.coh,font=2)
text(1917,0.0006,expression(c[2]),lwd=1,cex=cex.coh,font=2)
text(1942,0.0006,expression(c[3]),lwd=1,cex=cex.coh,font=2)

cex.x.lab <- cex.y.lab
title(xlab = "Cohort",cex.lab=cex.x.lab,
      outer = TRUE, line = 0.3)

dev.off()

