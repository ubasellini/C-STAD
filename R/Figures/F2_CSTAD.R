## -- C-STAD: MODEL OVERVIEW  --------------- ##
## Date: 16/07/2019
## Author: Ugofilippo Basellini
## Comments: 
##       - FIGURE 1: Schematic C-STAD overview
##
## ------------------------------------------ ##

## cleaning the workspace
rm(list=ls(all=TRUE))

## packages
library(MortalitySmooth)
library(pBrackets)

## load C-STAD functions
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Functions")
source("C-STAD_Functions.R")

## STANDARD: gompertz
dgom <- function(x, theta){
  ## Density Function
  alpha <- theta[1]
  beta <- theta[2]
  haza <- alpha * exp(beta*x)
  surv <- exp(alpha/beta * (1 - exp(beta*x)))
  d <- haza * surv
  return(d)
}
a <- 0.000018
b <- 0.12

## ages
delta <- 0.1
x <- 40:110
m <- length(x)
xs <- seq(min(x),max(x),by=delta)
ms <- length(xs)

## B-splines parameters
xl <- min(x)
xr <- max(x)
xmin <- round(xl - 0.01 * (xr - xl),3)
xmax <- round(xr + 0.01 * (xr - xl),3)
ndx <- floor(m/3)
deg <- 3
nbx <- ndx+deg

## B-splines bases
B <- MortSmooth_bbase(x, xmin, xmax, ndx, deg)
Bs <- MortSmooth_bbase(xs, xmin, xmax, ndx, deg)

## standard
FXstand <- dgom(xs,theta = c(a,b))  ## STANDARD
Mstand <- xs[which.max(FXstand)]

## coeff standard
ages.add.l <- 50
ages.add.r <- 40
delta1 <- 0.1

## define new augmented age-axis 
xA <- c(rev(seq(from=x[1]-delta1, by=-delta1,length=ages.add.l/delta1)), 
        xs, 
        seq(from=x[m]+delta1, by=delta1, length=ages.add.r/delta1))
mA <- length(xA)

## new B-splines parameters on augmented axis
xlA <- min(xA)
xrA <- max(xA)
xminA <- round(xlA - 0.01 * (xrA - xlA),3)
xmaxA <- round(xrA + 0.01 * (xrA - xlA),3)
ndxA <- floor((xA[mA]-xA[1])/3)

## augmented B-splines
BA <- MortSmooth_bbase(xA, xminA, xmaxA, ndx=ndxA, deg=deg)
nbxA <- ncol(BA)

## Derive coefficients of the standard
Standard <- coeff_stand(age=xs,fx=FXstand,ndx=ndxA,deg=deg,
                        ages.add.l=ages.add.l,ages.add.r=ages.add.r)
coeff_Stand <- Standard$betasA

## PURE SHIFTING 
s <- 7
dx.shift <- fx_shift(age = xs,fx = FXstand,shift = s,ndx = ndx,deg = deg)
ws <- xs - s

## C-STAD 

## segment axis
plo <- which(xA<=floor(Mstand+s))
pup <- which(xA>floor(Mstand+s))
xlo <- xA[plo]
xup <- xA[pup]

## segmenting a linear transformation function
## below the mode of year i
bL <- 1.2
cL <- -0.002
dL <- 0.0002
x.low <- xlo-s-Mstand
wL <- Mstand + bL*x.low + cL*(x.low^2) + dL*(x.low^3)
## above the mode of year i
bU <- 0.8
x.up <- xup-s-Mstand
wU <- Mstand + bU*x.up
## unique transformation function
wb <- c(wL, wU)

## B-splines on transformed ages
Bwb <- MortSmooth_bbase(x=c(wb),
                        xminA,xmaxA,ndx=ndxA,deg=deg)
## transformed density
fwb <- as.vector(exp(Bwb%*%coeff_Stand))
fwb <- fwb[xA%in%xs]
fwb <- 10*fwb/sum(fwb)
## hazard
eta <- log(mx_from_dx(fwb))
wb <- wb[xA%in%xs]

## plot
par(mfrow=c(1,2))
plot(xs,FXstand,t="l",ylim=range(FXstand,fwb))
lines(xs,dx.shift,col=2,lwd=2,lty=1)
lines(xs,fwb,col=4)
plot(xs, xs, t="l",ylim=range(xs,wb))
lines(xs, ws, col=2, lwd=3)
lines(xs, wb, col=4, lwd=3)
par(mfrow=c(1,1))

## PLOTTING
cex.axis <- 1
cex.pars <- 1.25
cex.leg <- 1.1
my.cols <- c("grey30","royalblue2","sienna2")

setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/Paper/Figures")
pdf("F2.pdf",width = 10,height=5.5)
## set overall margins
par(mfrow = c(1,2),
    oma = c(1.25,0.4,0,0.25),
    mar = c(1.25,3,1,0.1))
## bottom, left, top, right
plot(xs, FXstand, t="n",ylim=range(FXstand,dx.shift,(fwb+0.0008)),
     xlab="",ylab="",
     axes=FALSE,cex.main=1.25,
     main="Standard and transformed distributions")
axis(1,at=seq(20,120,20),cex.axis=cex.axis,padj = -1);
axis(2,las=2,at=seq(-0.01,0.06,0.01),cex.axis=cex.axis,hadj = 0.75)
mtext(text="f(x), f(t(x))", side=2, las=3,
      line=2.25, cex=1.35)
mtext(text="Ages", side=1, las=1,
      line=1.4, cex=1.35)
grid();box()
segments(x0=xs[which.max(dx.shift)], x1=xs[which.max(dx.shift)],
         y0=-10, y1=max(fwb), col=1,
         lwd=1, lty=2)
text(x=xs[which.max(dx.shift)]+0.75, y=0.04810, labels = expression(M^g),
     col=1, cex=1.2)
segments(x0=(Mstand), x1=(Mstand),
         y0=-10, y1=max(FXstand), col=my.cols[1],
         lwd=1, lty=2)
text(x=xs[which.max(FXstand)]-0.1, 
     y=max(FXstand)+0.002, labels = expression(M^f),
     col=my.cols[1], cex=1.2)

lines(xs, FXstand, col=my.cols[1], t="l", lwd=3,lty=1)
lines(xs, dx.shift, col=my.cols[2], t="l",lwd=3)
lines(xs, fwb, col=my.cols[3], lwd=3)

h.s <- 0.01020

h.arr <- 0.0020
fwb2 <- round(fwb,4)
x.left <- xs[which(fwb2==h.arr)][2]
x.right <- xs[which(fwb2==h.arr)][3]
delta <- 0.43
arrows(xs[which.max(dx.shift)]+delta, h.arr, x.right-delta, h.arr, code=2, col=my.cols[3], lwd=2, length=0.1)
arrows(xs[which.max(dx.shift)]-delta, h.arr, x.left+delta, h.arr, code=2, col=my.cols[3], lwd=2, length=0.1)

## write pars

brackets(Mstand, h.s, xs[which.max(dx.shift)], h.s, lwd=2, col=my.cols[2], h=0.00080)
text((Mstand+s/2), h.s+0.002, expression(s),
     col=my.cols[2], font=3, cex=cex.pars)
text(90, h.arr + 0.002, expression(b[U]),
     col=my.cols[3], font=3, cex=cex.pars)
text(71, h.arr + 0.002, expression(paste(b[L],",",c[L],",",d[L])),
     col=my.cols[3], font=3, cex=cex.pars)

legend("topleft",c("Standard f(x)",
                   expression(paste(g[1](x),"=",f(x-s))),
                   expression(paste(g[2](x),"=",f,"[",t(x,theta),"]"))),
       pch = NA,lty = c(1,1,1),col = my.cols,bty="n",cex = cex.leg,lwd = c(3,3,3))


#### TRANSFORMATION
plot(xs, xs, t="n",
     xlab="",
     ylab="",
     xlim=c(39, 112),
     ylim=range(xs,ws,wb),
     axes=FALSE,cex.main=1.25,
     main="Transformation functions")
axis(1,at=seq(20,120,20),cex.axis=cex.axis,padj = -1)
axis(2,at=seq(0,120,20),las=2,cex.axis=cex.axis,hadj = 0.75)
mtext(text="Transformed ages, t(x)", side=2, las=3,
      line=1.8, cex=1.35)
mtext(text="Actual ages, x", side=1, las=1,
      line=1.4, cex=1.35)
grid();box()

lines(xs, xs, col=my.cols[1], lwd=3, lty=1)
lines(xs, ws, col=my.cols[2], lwd=3)
lines(xs, wb, col=my.cols[3], lwd=3)
text(102.9, 86.1, expression(b[U]),
     col=my.cols[3], font=3, cex=cex.pars)
text(58, 25, expression(paste(b[L],",",c[L],",",d[L])),
     col=my.cols[3], font=3, cex=cex.pars)
brackets(110.7, 110.3, 110.7, 110.3-s, lwd=2, col=my.cols[2], h=1)
text(113.6, 107, expression(s),
     col=my.cols[2], font=3, cex=cex.pars)
legend("topleft",
       c(expression(paste("No transf., ",t(x,theta),"=",x)),
         expression(paste("Shifting transf., ",t(x,theta),"=",x-s)),
         expression(paste("Segmented transf., ",t(x,theta)))),
       pch = NA,lty = 1,col = my.cols,bty="n",
       cex = cex.leg,lwd = c(3,3,3))

segments(x0=(Mstand+s), x1=(Mstand+s),
         y0=Mstand, y1=-10, col=1,
         lwd=1, lty=2)
text(x=(Mstand+s)+4, y=36.9, labels = expression(M^g),
     col=1, cex=1.2)

par(mfrow=c(1,1))
dev.off()


