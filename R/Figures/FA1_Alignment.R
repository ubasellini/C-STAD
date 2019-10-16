## -- C-STAD: plotting results ------------ ##
## Date: 24/06/2019
## Author: Ugofilippo Basellini
## Comments: 
##       - FIGURE A1: Aligned distributions 
##
## ------------------------------------------ ##

## clean the workspace
rm(list = ls())

## load useful packages
library(MortalitySmooth)
library(colorspace)
library(viridis)
library(fields)

## -- DATA & FUNCTIONS ----------

## load C-STAD functions
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Functions")
source("C-STAD_Functions.R")
source("BootDxFUN.R")

## load data 
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Data")
cou <- "SWE"        ## SWE or DNK
sex <- "F"          ## only F 
name <- paste0(cou,"coh",sex,".Rdata")  ## Females
load(name)    

## age dimensions
ages <- as.numeric(rownames(cE))
cohorts <- as.numeric(colnames(cE))
age.start <- 40
xo <- age.start:110     ## HMD ages
x <- age.start:120      ## expanded ages
mo <- length(xo)
m <- length(x)
delta <- 0.1
xs <- seq(min(x), max(x), delta)  ## ages at fine grid
ms <- length(xs)

## cohort dimensions
year.start <- 1835
year.end <- 2015 - age.start - 5
y <- year.start:year.end  ## use for both SWE and DNK (to have same time-range of analysis + reliable data)
n <- length(y)
coly <- rainbow_hcl(n)  
## (1835 = first cohort with data observed at all ages in DNK)
## (1970 = 5 years before last cohort with data observed at age 40)

## cohorts first Lexis parallelogram (c1)
c_breve <- cohorts[min(which(is.na(cE[nrow(cE),])))] - 1
c1 <- y[1]:c_breve  ## 1905 = last cohort with fully observed data
n1 <-  length(c1)

## starting data
E <- cE[ages%in%xo,cohorts%in%y]
MX.act <- cMx[ages%in%xo,cohorts%in%y]
Z <- Zna <- cZ[ages%in%xo,cohorts%in%y]
Z[E==0] <- 0
W <- matrix(1,mo,n)
W[is.na(E)] <- 0  ## zero weight where data is missing

## expand data for extrapolation 
EA <- rbind(E,matrix(100,(m-mo),n))
ZA <- rbind(Z,matrix(100,(m-mo),n))
Wup <- cbind(matrix(1,(m-mo),n1),       
             matrix(0,(m-mo),(n-n1)))   ## 1 to consider 110-120 of c1 for standard
WA <- rbind(W,Wup)

## weights augmented (repeat each weight 10 times) for the smooth standard
One <- matrix(rep(1,10),1/delta,1)
Ws <- kronecker(WA,One)       ## expanded matrix of weights
Ws <- Ws[c(1:ms),]            ## remove last 9 weights
WA[EA==0] <- 0                ## zero weight where data equal to zero
WA[which(x>xo[mo]),] <- 0     ## zero weight for ages above 110

## log death rates
LMX.act <- log(MX.act)
matplot(xo,LMX.act,t="l",lty=1,col = coly,xlab="Age",
        main=paste("Observed Cohort Mortality,",cou,sex,y[1],"-",y[n]))

## B-splines parameters
xl <- min(x)
xr <- max(x)
xmin <- round(xl - 0.01 * (xr - xl),3)
xmax <- round(xr + 0.01 * (xr - xl),3)
ndx <- floor(m/5)
yl <- min(y)
yr <- max(y)
ymin <- round(yl - 0.01 * (yr - yl),3)
ymax <- round(yr + 0.01 * (yr - yl),3)
ndy <- floor(n/5)
deg <- 3

## B-splines bases
Bx <- MortSmooth_bbase(x, xmin, xmax, ndx, deg)
Bxs <- MortSmooth_bbase(xs, xmin, xmax, ndx, deg)
By <- MortSmooth_bbase(y, ymin, ymax, ndy, deg)
B <- kronecker(By,Bx)
Bs <- kronecker(By,Bxs)

## -- STANDARD ----------

## 2D smooth (optimal parameters)
if(cou=="DNK"){
  lambdaX.hat <- 10^2.5
  lambdaY.hat <- 10^3 
}else if(cou=="SWE"){
  lambdaX.hat <- 10^2.5
  lambdaY.hat <- 10^3.5
}
smooth2D <- Mort2Dsmooth(x=x,y=y,Z=ZA,offset=log(EA),W=WA,
                         ndx=c(ndx,ndy),method = 3,
                         lambdas = c(lambdaX.hat,lambdaY.hat))
plot(smooth2D, palette = "terrain.colors")
LMX.smo2D <- matrix(Bs %*% c(smooth2D$coefficients),ms,n)
LMX.smo2DW <- LMX.smo2D*Ws
LMX.smo2DW[LMX.smo2DW==0] <- NA
par(mfrow=c(1,3))
matplot(xo,LMX.act,t="l",lty=1,col = coly,xlab="Age",ylim = range(LMX.act,LMX.smo2D,finite=T),
        main=paste("Observed Mortality,",cou,y[1],"-",y[n]))
matplot(xs,LMX.smo2D,t="l",lty=1,col = coly,xlab="Age",ylim = range(LMX.act,LMX.smo2D,finite=T),
        main=paste("Smooth lmx extrapolated"))
matplot(xs,LMX.smo2DW,t="l",lty=1,col = coly,xlab="Age",ylim = range(LMX.act,LMX.smo2D,finite=T),
        main=paste("Smooth lmx observed"))
par(mfrow=c(1,1))

## actual and 2Dsmo e40, g40
e40.act <- apply(exp(LMX.act[,y%in%c1]),2,lifetable.ex,x=xo,sex=sex)
g40.act <- apply(exp(LMX.act[,y%in%c1]),2,GINI_func,ages=xo,sex=sex)
e40_2D <- apply(exp(LMX.smo2D[xs%in%xo,]),2,lifetable.ex,x=xo,sex=sex)
g40_2D <- apply(exp(LMX.smo2D[xs%in%xo,]),2,GINI_func,ages=xo,sex=sex)

par(mfrow=c(1,2))
plot(c1,e40.act,ylim=range(e40.act,e40_2D),pch=16,
     main=paste0("E",xo[1]),ylab="",cex.lab=1.25,xlim=range(y),xlab="Cohort")
lines(y,e40_2D,col=5,lwd=2)
plot(c1,g40.act,ylim=range(g40.act,g40_2D),pch=16,
     main=paste0("G",xo[1]),ylab="",cex.lab=1.25,xlim=range(y),xlab="Cohort")
lines(y,g40_2D,col=5,lwd=2)
par(mfrow=c(1,1))

## FX
FX.smo2D <- apply(exp(LMX.smo2D),2,dx_from_mx,age=xs)
FX.smo2DW <- FX.smo2D*Ws
FX.smo2DW[FX.smo2DW==0] <- NA
par(mfrow=c(1,2))
matplot(xs, FX.smo2D, lty=1, t="l", col=coly,xlab="Age",ylab="fx",
        main="Smooth Dx, extrapolated",cex.lab=1.25)
matplot(xs, FX.smo2DW, lty=1, t="l", col=coly,xlab="Age",ylab="fx",
        main="Smooth Dx, observed",cex.lab=1.25)
par(mfrow=c(1,1))

## compute modal age at death (it's pretty smooth here due to 2D smoothing)
M_2D <- xs[apply(FX.smo2D, 2, which.max)] + (delta/2)
plot(y, M_2D, t="o", lwd=2, pch=16, main="Modal Age at Death (smooth)")

## STANDARD DISTRIBUTION
s_2D <- M_2D - M_2D[1]

## derive aligned distributions
FX.align <- matrix(0, nrow=ms, ncol=n)
for(i in 1:n){
  FX.align[,i] <- fx_shift(age=xs,fx=FX.smo2D[,i],shift=-s_2D[i],ndx = ndx,deg = deg)
}
FX.alignW <- FX.align*Ws
FX.alignW[FX.alignW==0] <- NA
par(mfrow=c(1,2))
matplot(xs, FX.align, lty=1, t="l", col=coly,xlab="Age",ylab="fx",
        main="Smooth Dx, extrapolated",cex.lab=1.25)
matplot(xs, FX.alignW, lty=1, t="l", col=coly,xlab="Age",ylab="fx",
        main="Smooth Dx, observed",cex.lab=1.25)
par(mfrow=c(1,1))

## Standard = mean of the aligned densities
FXallmeanW <- exp(apply(log(FX.alignW), 1, mean, na.rm=T))  
FXstand <- FXallmeanW
matplot(xs, FX.alignW, lty=1, t="l", col=coly,xlab="Age",ylab="fx",
        main="Smooth Dx, extrapolated",cex.lab=1.25)
lines(xs, FXallmeanW, lwd=3,col=1)
legend("bottomleft", c("Standard"),col=1, lwd=3,lty = 1,
       bg="white", pt.cex=1.2,bty="n",cex = 1.5)


ylim <- range(FX.smo2DW,na.rm=T)

## colors
# display.brewer.pal(n=9, name = 'Purples')
# col.stad <- brewer.pal(n = 8, name = 'Blues')[8]
# col.stadT <- adjustcolor(col.stad, alpha=0.3)
# my.orange <- brewer.pal(n=9, name = 'Oranges')[6]
# my.green <- brewer.pal(n=9, name = 'Greens')[6]
# my.purple <- brewer.pal(n=9, name = 'Purples')[7]
# my.cols <- c(col.stad,my.orange,my.green,my.purple)
# my.colsT <- adjustcolor(my.cols, alpha=0.2)
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
my.col <- colorRampPalette(c("purple3",
                             "blue2", "cyan2","green", "yellow","goldenrod2"))(n)
my.col <- rainbow_hcl(n)

## SAVE
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/Paper/Figures")
pdf("FA0.pdf",width = 10,height = 5.5)
par(mfrow = c(1,2),
    oma = c(1.25,1.2,0.1,0.25),
    mar = c(1.75,1.3,1.2,0.1))
## bottom, left, top, right

## Smooth dx, observed
matplot(xs,FX.smo2DW,xlim=range(x),
        ylim=ylim,xlab="",ylab="",t="n",axes=F)
axis(1,cex.axis=cex.x.axis,padj = -0.5)
axis(2,cex.axis=cex.y.axis,at=seq(0,0.04,0.01),
     labels = c("0","0.01","0.02","0.03","0.04"))
grid();box()
title(main="Smooth distributions", cex.main=cex.title)
matlines(xs,FX.smo2DW,t="l",lty=1,col=my.col)
# lines(xs,apply(FX.smo2DW,1,mean,na.rm=T),lwd=3)
image.plot(smallplot=c(.15,.45, .84,.88),axis.args=list(cex.axis=0.8),
           legend.only=TRUE, zlim=range(y),
           col=my.col, nlevel=n,
           horizontal = TRUE)

## Aligned dx
matplot(xs,FX.smo2DW,xlim=range(x),
        ylim=ylim,xlab="",ylab="",t="n",axes=F)
axis(1,cex.axis=cex.x.axis,padj = -0.5)
axis(2,las=2,cex.axis=cex.y.axis,at=seq(0,0.04,0.01),
     labels = rep("",5))
grid();box()
title(main="Aligned distributions", cex.main=cex.title)
abline(v=xs[which.max(apply(FX.alignW,1,mean,na.rm=T))],lty=2)

matlines(xs,FX.alignW,t="l",lty=1,col=my.col)
lines(xs,apply(FX.alignW,1,mean,na.rm=T),lwd=3)
legend("topright", c("Standard"),col=1, lwd=3,lty = 1,
       bg="white", pt.cex=1.2,bty="n",cex = 1.5)

cex.x.lab <- cex.y.lab
title(xlab = "Ages",cex.lab=cex.x.lab,
      outer = TRUE, line = 0.1)

dev.off()




