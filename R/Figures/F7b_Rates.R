## -- C-STAD: plotting results ------------ ##
## Date: 24/06/2019
## Author: Ugofilippo Basellini
## Comments: 
##       - FIGURE 6b: Rates over age
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
LMXcstadSWE <- LMXcstad

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

## read data
cou <- "DNK"
name <- paste0(cou,"_",sex,"_",n.sim,"simul_results.Rdata")
load(name)
MXactDNK <- exp(LMX.act)

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
my.colsT <- adjustcolor(my.cols, alpha=0.2)
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
pdf("F7b.pdf",width = 10,height = 5.5)
par(mfrow = c(1,2),
    oma = c(1.25,1.2,0.1,0.25),
    mar = c(1.75,2.3,1.2,0.3))

## SWE
## bottom, left, top, right
matplot(xo,MXactSWE,pch=16,log="y",
     ylim=ylim,xlab="",ylab="",t="n",axes=F)
axis(1,cex.axis=cex.x.axis,padj = -0.5)
axis(2,las=3,cex.axis=cex.y.axis,at=c(1e-3,1e-2,1e-1,1e0),
     labels = c("0.001","0.01","0.1","1"))
grid();box()
# mtext("Cohort",side=1,cex=cex.x.lab,line=2.4)
mtext("Log-mortality",side=2,cex=cex.y.lab,line=2,las=3)
title(main="Sweden", cex.main=cex.title)

coh <- y[length(y)]
points(xo,MXactSWE[,which(y==coh)],pch=1,col=my.cols[3],cex=cex.obs,lwd=lwd.pt)
xx <- c(xo,rev(xo))
yy <- c(exp(ETA.LOWswe[1:length(xo),which(y==coh)]),rev(exp(ETA.UPswe[1:length(xo),which(y==coh)])))
polygon(xx, yy, border = my.cols[3], col=my.colsT[3])
lines(xo,exp(ETA.MEANswe[1:length(xo),which(y==coh)]),t="l",lwd=lwd.mean,col=my.cols[3])
text(85,0.01,paste(coh),lwd=1,cex=cex.age,font=2,col=my.cols[3])

coh <- 1925
points(xo,MXactSWE[,which(y==coh)],pch=1,col=my.cols[2],cex=cex.obs,lwd=lwd.pt)
yy <- c(exp(ETA.LOWswe[1:length(xo),which(y==coh)]),rev(exp(ETA.UPswe[1:length(xo),which(y==coh)])))
polygon(xx, yy, border = my.cols[2], col=my.colsT[2])
lines(xo,exp(ETA.MEANswe[1:length(xo),which(y==coh)]),t="l",lwd=lwd.mean,col=my.cols[2])
text(45,0.005,paste(coh),lwd=1,cex=cex.age,font=2,col = my.cols[2])

coh <- y[1]
points(xo,MXactSWE[,which(y==coh)],pch=1,col=my.cols[1],cex=cex.obs,lwd=lwd.pt)
yy <- c(exp(ETA.LOWswe[1:length(xo),which(y==coh)]),rev(exp(ETA.UPswe[1:length(xo),which(y==coh)])))
polygon(xx, yy, border = my.cols[1], col=my.colsT[1])
lines(xo,exp(ETA.MEANswe[1:length(xo),which(y==coh)]),t="l",lwd=lwd.mean,col=my.cols[1])
text(70,0.1,paste(coh),lwd=1,cex=cex.age,font=2,col=my.cols[1])

legend("topleft",c("observed","C-STAD"),pch=c(1,NA),col=1,
       bty="n",lwd=c(lwd.pt,lwd.mean),cex=cex.leg,lty=c(NA,1),inset = 0.01)
legend("topleft",c("",""),pch=c(NA,15),col=adjustcolor(1,alpha.f = 0.3),
       bty="n",lwd=c(NA,NA),cex=cex.leg,lty=c(NA,1),inset = 0.01)


## DNK
## bottom, left, top, right
matplot(xo,MXactDNK,pch=16,log="y",
        ylim=ylim,xlab="",ylab="",t="n",axes=F)
axis(1,cex.axis=cex.x.axis,padj = -0.5)
axis(2,las=3,cex.axis=cex.y.axis,at=c(1e-3,1e-2,1e-1,1e0),
     labels = c("0.001","0.01","0.1","1"))
grid();box()
title(main="Denmark", cex.main=cex.title)

coh <- y[length(y)]
points(xo,MXactDNK[,which(y==coh)],pch=1,col=my.cols[3],cex=cex.obs,lwd=lwd.pt)
xx <- c(xo,rev(xo))
yy <- c(exp(ETA.LOWdnk[1:length(xo),which(y==coh)]),rev(exp(ETA.UPdnk[1:length(xo),which(y==coh)])))
polygon(xx, yy, border = my.cols[3], col=my.colsT[3])
lines(xo,exp(ETA.MEANdnk[1:length(xo),which(y==coh)]),t="l",lwd=lwd.mean,col=my.cols[3])
text(85,0.01,paste(coh),lwd=1,cex=cex.age,font=2,col=my.cols[3])

coh <- 1925
points(xo,MXactDNK[,which(y==coh)],pch=1,col=my.cols[2],cex=cex.obs,lwd=lwd.pt)
yy <- c(exp(ETA.LOWdnk[1:length(xo),which(y==coh)]),rev(exp(ETA.UPdnk[1:length(xo),which(y==coh)])))
polygon(xx, yy, border = my.cols[2], col=my.colsT[2])
lines(xo,exp(ETA.MEANdnk[1:length(xo),which(y==coh)]),t="l",lwd=lwd.mean,col=my.cols[2])
text(45,0.005,paste(coh),lwd=1,cex=cex.age,font=2,col = my.cols[2])


coh <- y[1]
points(xo,MXactDNK[,which(y==coh)],pch=1,col=my.cols[1],cex=cex.obs,lwd=lwd.pt)
yy <- c(exp(ETA.LOWdnk[1:length(xo),which(y==coh)]),rev(exp(ETA.UPdnk[1:length(xo),which(y==coh)])))
polygon(xx, yy, border = my.cols[1], col=my.colsT[1])
lines(xo,exp(ETA.MEANdnk[1:length(xo),which(y==coh)]),t="l",lwd=lwd.mean,col=my.cols[1])
text(70,0.1,paste(coh),lwd=1,cex=cex.age,font=2,col=my.cols[1])

cex.x.lab <- cex.y.lab
title(xlab = "Ages",cex.lab=cex.x.lab,
      outer = TRUE, line = 0.1)

dev.off()

