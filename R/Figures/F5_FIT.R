## -- C-STAD: plotting results ------------ ##
## Date: 24/06/2019
## Author: Ugofilippo Basellini
## Comments: 
##       - FIGURE 4: Estimated e40, g40
##
## ------------------------------------------ ##

rm(list = ls())

## load SWE results
cou <- "SWE"
sex <- "F"
n.sim <- 1600
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Results")
name <- paste0(cou,"_",sex,"_",n.sim,"simul_results.Rdata")
load(name)

couSWE <- "Sweden"
c1SWE <- c1
e40.actSWE <- e40.act
g40.actSWE <- g40.act
e40.hatSWE <- e40_c1
g40.hatSWE <- g40_c1

## load DNK results
cou <- "DNK"
name <- paste0(cou,"_",sex,"_",n.sim,"simul_results.Rdata")
load(name)

couDNK <- "Denmark"
c1DNK <- c1
e40.actDNK <- e40.act
g40.actDNK <- g40.act
e40.hatDNK <- e40_c1
g40.hatDNK <- g40_c1

## plotting
require(RColorBrewer)
my.cols <- brewer.pal(8, "Set1")
cex.axis <- 0.75
cex4 <- 0.8
cex16 <- 1.1
lwd.x <- 2.5
cex.coh <- 1.5
cex.cou <- 1.3
cex.leg <- 1.15
cex.x.axis <- 1.1
cex.y.axis <- 1.1
cex.y.lab <-  1.75
cex.main <- 1.75

my.cols[3] <- brewer.pal(9, "Greens")[6]
my.cols[4] <- brewer.pal(9, "BuPu")[8]
display.brewer.pal(n = 9, name = 'BuPu')

setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/Paper/Figures")
pdf("F5.pdf",width = 10,height=5.5)
## set overall margins
par(mfrow = c(1,2),
    oma = c(1.25,0.7,0.1,0.25),
    mar = c(1.75,2.35,1.5,0.9))
## bottom, left, top, right

## left: E40
ylimE40fit <- range(e40.actDNK,e40.actSWE,e40.hatDNK,e40.hatSWE)
plot(c1,e40.actDNK,t="n",pch=16,col=my.cols[2],ylim=ylimE40fit,
     cex.main=1,axes=F,xlim=range(c1DNK,c1SWE))
title(main = expression(paste("Observed and fitted ",e[40])),
      cex.main=cex.main)
axis(1,cex.axis=cex.x.axis,padj = -0.5);
axis(2,las=2,cex.axis=cex.y.axis,hadj = 0.75)
mtext(text=expression(e[40]), side=2, las=3,
      line=1.25, cex=cex.y.lab)
grid();box()
points(c1,e40.actDNK,pch=16,col=my.cols[2],cex=cex16)
text(1895,33,couDNK,col=my.cols[2],cex=cex.cou)
points(c1SWE,e40.actSWE,pch=16,col=my.cols[3],cex=cex16)
text(1878,37,couSWE,col=my.cols[3],cex=cex.cou)
points(c1,e40.hatDNK,pch=4,col=my.cols[5],lwd=lwd.x,cex=cex4)
points(c1SWE,e40.hatSWE,pch=4,col=my.cols[4],lwd=lwd.x,cex=cex4)
legend("topleft",c("Observed","C-STAD"),pch = c(16,4),
       bty="n",lwd=2,lty=c(NA,NA),cex = cex.leg)

## right: G40
ylimG40fit <- range(g40.actDNK,g40.actSWE,g40.hatDNK,g40.hatSWE)
plot(c1,g40.actDNK,t="n",pch=16,col=my.cols[2],ylim=ylimG40fit,
     cex.main=1,axes=F,xlim=range(c1DNK,c1SWE))
title(main = expression(paste("Observed and fitted ",G[40])),
      cex.main=cex.main)
axis(1,cex.axis=cex.x.axis,padj = -0.5);
axis(2,las=2,cex.axis=cex.y.axis,hadj = 0.75,
     at=seq(0.18,0.24,0.02),labels = c("18","20","22","24"))
mtext(text=expression(G[40]), side=2, las=3,
      line=1.25, cex=cex.y.lab)
grid();box()
points(c1,g40.actDNK,pch=16,col=my.cols[2],cex=cex16)
text(1895,0.22,couDNK,col=my.cols[2],cex=cex.cou)
points(c1SWE,g40.actSWE,pch=16,col=my.cols[3],cex=cex16)
text(1880,0.19,couSWE,col=my.cols[3],cex=cex.cou)
points(c1,g40.hatDNK,pch=4,col=my.cols[5],lwd=lwd.x,cex=cex4)
points(c1SWE,g40.hatSWE,pch=4,col=my.cols[4],lwd=lwd.x,cex=cex4)

cex.x.lab <- cex.y.lab
title(xlab = "Cohort",cex.lab=cex.x.lab,
      outer = TRUE, line = 0.3)
par(mfrow = c(1,1))
dev.off()
