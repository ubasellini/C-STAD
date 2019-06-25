## -- C-STAD: plotting results ------------ ##
## Date: 24/06/2019
## Author: Ugofilippo Basellini
## Comments: 
##       - FIGURE 7: Age-at-death distributions
##
## ------------------------------------------ ##

## cleaning the workspace
rm(list=ls(all=TRUE))
library(RColorBrewer)
library(fields)
library(colorspace)

## read data
cou <- "SWE"
sex <- "F"
n.sim <- 1600
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Results")
name <- paste0(cou,"_",sex,"_",n.sim,"simul_results.Rdata")
load(name)
ADDobsSWE <- ADDobs
ADDcstadSWE <- ADDcstad

## CI for bootstrap
lev.p <- 0.8
ADD.MEANswe <- ADD.LOWswe <- ADD.UPswe <- matrix(NA,dim(ADD)[1],dim(ADD)[2])
i <- j <- 1
for (i in 1:dim(ADD)[1]){
  for (j in 1:dim(ADD)[2]){
    ADD.MEANswe[i,j] <- median(ADD[i,j,,],na.rm=T)
    ADD.UPswe[i,j] <- quantile(ADD[i,j,,],prob=1- (1-lev.p)/2,na.rm=T)
    ADD.LOWswe[i,j] <- quantile(ADD[i,j,,],prob=(1-lev.p)/2,na.rm=T)
  }
}

## read data
cou <- "DNK"
name <- paste0(cou,"_",sex,"_",n.sim,"simul_results.Rdata")
load(name)
ADDobsDNK <- ADDobs
ADDcstadDNK <- ADDcstad

## CI for bootstrap
lev.p <- 0.8
ADD.MEANdnk <- ADD.LOWdnk <- ADD.UPdnk <- matrix(NA,dim(ADD)[1],dim(ADD)[2])
i <- j <- 1
for (i in 1:dim(ADD)[1]){
  for (j in 1:dim(ADD)[2]){
    ADD.MEANdnk[i,j] <- median(ADD[i,j,,],na.rm=T)
    ADD.UPdnk[i,j] <- quantile(ADD[i,j,,],prob=1- (1-lev.p)/2,na.rm=T)
    ADD.LOWdnk[i,j] <- quantile(ADD[i,j,,],prob=(1-lev.p)/2,na.rm=T)
  }
}

## colors
display.brewer.pal(n=9, name = 'Purples')
col.stad <- brewer.pal(n = 8, name = 'Blues')[8]
col.stadT <- adjustcolor(col.stad, alpha=0.3)
my.orange <- brewer.pal(n=9, name = 'Oranges')[6]
my.green <- brewer.pal(n=9, name = 'Greens')[6]
my.purple <- brewer.pal(n=9, name = 'Purples')[7]
my.cols <- c(col.stad,my.orange,my.green,my.purple)
my.colsT <- adjustcolor(my.cols, alpha=0.2)
cex.x.axis <- 1.1
cex.y.axis <-  1.1
cex.x.lab <- 2
cex.y.lab <-  1.75
cex.leg <- 1.15
cex.coh <- 1.5
cex.title <- 1.7 
cex.obs <- 0.85
cex.age <- 0.9
lwd.pt <- 0.9
lwd.mean <- 2.5

##--- RATES ------
ylim <- c(0,max(ADDobsSWE,na.rm=T))
ylim[2] <- 5000
x <- 40:120

## SAVE
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/Paper/Figures")
pdf("F7.pdf",width = 10,height = 5.5)
par(mfrow = c(1,2),
    oma = c(1.25,1.2,0.1,0.25),
    mar = c(1.75,2.3,1.2,0.3))

## SWE
## bottom, left, top, right
matplot(xo,ADDobsSWE,pch=16,xlim=range(x),
     ylim=ylim,xlab="",ylab="",t="n",axes=F)
axis(1,cex.axis=cex.x.axis,padj = -0.5)
axis(2,las=2,cex.axis=cex.y.axis,at=seq(0,5000,1000),
     labels = c("0","1","2","3","4","5"))
grid();box()
mtext("x 10^3", 3, cex=0.6, at=40)
mtext("Age-at-death distributions",side=2,cex=cex.y.lab,line=2,las=3)
title(main="Sweden", cex.main=cex.title)

coh <- y[length(y)]
points(xo,ADDobsSWE[,which(y==coh)],pch=1,col=my.cols[3],cex=cex.obs,lwd=lwd.pt)
xx <- c(x,rev(x))
yy <- c(ADD.LOWswe[1:length(x),which(y==coh)],rev(ADD.UPswe[1:length(x),which(y==coh)]))
polygon(xx, yy, border = my.cols[3], col=my.colsT[3])
lines(x,ADD.MEANswe[1:length(x),which(y==coh)],t="l",lwd=lwd.mean,col=my.cols[3])
# lines(xo,exp(LMXcstadSWE[1:length(xo),which(y==coh)]),t="l",lwd=lwd.mean,col=my.cols[3],lty=1)
text(108,3000,paste(coh),lwd=1,cex=cex.age,font=2,col=my.cols[3])
# lines(cohorts,t(exp(LMX.smo2D[whi.smo,])),t="l",lwd=2,col=my.orange[6])

coh <- 1925
points(xo,ADDobsSWE[,which(y==coh)],pch=1,col=my.cols[2],cex=cex.obs,lwd=lwd.pt)
yy <- c(ADD.LOWswe[1:length(x),which(y==coh)],rev(ADD.UPswe[1:length(x),which(y==coh)]))
polygon(xx, yy, border = my.cols[2], col=my.colsT[2])
lines(x,ADD.MEANswe[1:length(x),which(y==coh)],t="l",lwd=lwd.mean,col=my.cols[2])
text(92,2100,paste(coh),lwd=1,cex=cex.age,font=2,col = my.cols[2])

coh <- y[1]
points(xo,ADDobsSWE[,which(y==coh)],pch=1,col=my.cols[1],cex=cex.obs,lwd=lwd.pt)
yy <- c(ADD.LOWswe[1:length(x),which(y==coh)],rev(ADD.UPswe[1:length(x),which(y==coh)]))
polygon(xx, yy, border = my.cols[1], col=my.colsT[1])
lines(x,ADD.MEANswe[1:length(x),which(y==coh)],t="l",lwd=lwd.mean,col=my.cols[1])
text(65,3000,paste(coh),lwd=1,cex=cex.age,font=2,col=my.cols[1])

legend("topleft",c("observed","C-STAD"),pch=c(1,NA),col=1,
       bty="n",lwd=c(lwd.pt,lwd.mean),cex=cex.leg,lty=c(NA,1),inset = 0.01)
legend("topleft",c("",""),pch=c(NA,15),col=adjustcolor(1,alpha.f = 0.3),
       bty="n",lwd=c(NA,NA),cex=cex.leg,lty=c(NA,1),inset = 0.01)


## DNK
## bottom, left, top, right
matplot(xo,ADDobsDNK,pch=16,xlim=range(x),
        ylim=ylim,xlab="",ylab="",t="n",axes=F)
axis(1,cex.axis=cex.x.axis,padj = -0.5)
axis(2,las=2,cex.axis=cex.y.axis,at=seq(0,5000,1000),
     labels = c("0","1","2","3","4","5"))
mtext("x 10^3", 3, cex=0.6, at=40)
grid();box()
title(main="Denmark", cex.main=cex.title)

coh <- y[length(y)]
points(xo,ADDobsDNK[,which(y==coh)],pch=1,col=my.cols[3],cex=cex.obs,lwd=lwd.pt)
xx <- c(x,rev(x))
yy <- c(ADD.LOWdnk[1:length(x),which(y==coh)],rev(ADD.UPdnk[1:length(x),which(y==coh)]))
polygon(xx, yy, border = my.cols[3], col=my.colsT[3])
lines(x,ADD.MEANdnk[1:length(x),which(y==coh)],t="l",lwd=lwd.mean,col=my.cols[3])
text(107,3000,paste(coh),lwd=1,cex=cex.age,font=2,col=my.cols[3])

coh <- 1925
points(xo,ADDobsDNK[,which(y==coh)],pch=1,col=my.cols[2],cex=cex.obs,lwd=lwd.pt)
yy <- c(ADD.LOWdnk[1:length(x),which(y==coh)],rev(ADD.UPdnk[1:length(x),which(y==coh)]))
polygon(xx, yy, border = my.cols[2], col=my.colsT[2])
lines(x,ADD.MEANdnk[1:length(x),which(y==coh)],t="l",lwd=lwd.mean,col=my.cols[2])
text(91,2100,paste(coh),lwd=1,cex=cex.age,font=2,col = my.cols[2])

coh <- y[1]
points(xo,ADDobsDNK[,which(y==coh)],pch=1,col=my.cols[1],cex=cex.obs,lwd=lwd.pt)
yy <- c(ADD.LOWdnk[1:length(x),which(y==coh)],rev(ADD.UPdnk[1:length(x),which(y==coh)]))
polygon(xx, yy, border = my.cols[1], col=my.colsT[1])
lines(x,ADD.MEANdnk[1:length(x),which(y==coh)],t="l",lwd=lwd.mean,col=my.cols[1])
text(65,3000,paste(coh),lwd=1,cex=cex.age,font=2,col=my.cols[1])

cex.x.lab <- cex.y.lab
title(xlab = "Ages",cex.lab=cex.x.lab,
      outer = TRUE, line = 0.1)

dev.off()

