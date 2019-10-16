## ------------------------------------------------------------- ##
##  C-STAD model: fitting and forecasting cohort mortality
##  Last update: 24/06/2019
##  Authors: Ugofilippo Basellini 
##           Soren Kjaergaard
##           Giancarlo Camarda
##
##  sessionInfo() details:
##  
##  R version 3.5.1 (2018-07-02)
##  Platform: x86_64-apple-darwin15.6.0 (64-bit)
##  Running under: macOS  10.14.5
##  
##  locale: en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
##
##  attached base packages:
##  splines  grid  stats  graphics  grDevices  utils  datasets 
##  methods  base     
## 
##  other attached packages:
##  rootSolve_1.7  vars_1.5-3  lmtest_0.9-36  urca_1.3-0  strucchange_1.5-1  sandwich_2.5-0       
##  zoo_1.8-3  MASS_7.3-50  MortalitySmooth_2.3.4  svcm_0.1.2  Matrix_1.2-14  shape_1.4.4          
##  plotrix_3.7-4  lattice_0.20-35  colorspace_1.3-2  fields_9.6  maps_3.3.0  spam_2.2-0  
##  dotCall64_1.0-0  RColorBrewer_1.1-2
## ------------------------------------------------------------- ##

## clean the workspace
rm(list = ls())

## load useful packages
library(MortalitySmooth)
library(vars)
library(colorspace)

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

## Derive coefficients of the standard (need to augment x-axis)
ages.add.l <- 90
ages.add.r <- 80
if (cou == "DNK"){
  ages.add.l <- ages.add.l+10
  ages.add.r <- ages.add.r+10
}
delta1 <- 1

## define new augmented age-axis 
xA <- c(rev(seq(from=x[1]-delta1, by=-delta1,length=ages.add.l/delta1)), 
        x, 
        seq(from=x[m]+delta1, by=delta1, length=ages.add.r/delta1))
mA <- length(xA)
xAs <- seq(xA[1],xA[mA],delta)

## new B-splines parameters on augmented axis
xlA <- min(xA)
xrA <- max(xA)
xminA <- round(xlA - 0.01 * (xrA - xlA),3)
xmaxA <- round(xrA + 0.01 * (xrA - xlA),3)
ndxA <- floor((xA[mA]-xA[1])/4)

## augmented B-splines
BA <- MortSmooth_bbase(xA, xminA, xmaxA, ndx=ndxA, deg=deg)
BAs <- MortSmooth_bbase(xAs, xminA, xmaxA, ndx=ndxA, deg=deg)
nbxA <- ncol(BA)

## coefficients & mode
lambdaSTAND <- 1
Standard <- coeff_stand(age=xs,fx=FXstand,ndx=ndxA,deg=deg,
                        ages.add.l=ages.add.l,ages.add.r=ages.add.r,
                        lambda=lambdaSTAND)
coeff_Stand <- Standard$betasA
(Mstand <-  xAs[which.max(BAs%*%coeff_Stand)])

## -- C1 ----------

## C-STAD ESTIMATION on c1
E1 <- cE[ages%in%xo,cohorts%in%c1]
MX1.act <- cMx[ages%in%xo,cohorts%in%c1]
Z1 <- cZ[ages%in%xo,cohorts%in%c1]
Z1[E1==0] <- 0
LMX1.act <- log(MX1.act)
W1 <- matrix(1,mo,n1)
W1[E1==0] <- 0      ## zero weight where no exposure

## expand data for 1D smoothing
E1A <- rbind(E1,matrix(100,m-mo,n1))
Z1A <- rbind(Z1,matrix(100,m-mo,n1))
W1A <- rbind(W1,matrix(0,m-mo,n1))

## find unsmooth M and s
LMX.smo1D <- FX.smo1D <-  matrix(NA,ms,n1)
i <- 1
for (i in 1:n1){
  lmx.smo <- lmx_smooth(age = x,y = Z1A[,i],e = E1A[,i],w =  W1A[,i],
                          ndx = ndx,deg = deg)
  LMX.smo1D[,i] <- Bxs %*% lmx.smo$coef
  FX.smo1D[,i] <- dx_from_mx(age=xs,mx=exp(LMX.smo1D[,i]))
}

## compute modal age at death 
M_c1 <- xs[apply(FX.smo1D, 2, which.max)] + (delta/2)
s_c1 <- M_c1 - Mstand
plot(c1, M_c1, t="o", lwd=2, pch=16, main="Modal Age at Death")
points(c1, M_2D[y%in%c1], t="o", lwd=2, pch=16, col=2)

## Break point of the age axis for each year 
PLO1 <- PUP1 <-  list()
XLO1 <- XUP1 <- list()
for(i in 1:n1){
  PLO1[[i]] <- which(xA<=M_c1[i])
  PUP1[[i]] <- which(xA>M_c1[i])
  XLO1[[i]] <- xA[PLO1[[i]]]
  XUP1[[i]] <- xA[PUP1[[i]]]
}

## empty vectors and matrices to store results 
PLOT=F
bL_c1 <- bU_c1 <- cL_c1 <- dL_c1 <- numeric(n1)
DXcstad_c1 <- LMXcstad_c1 <-  matrix(NA, m, n1)

## Estimation
i <- 1
for(i in 1:n1){
  cat("fitting cohort", c1[i], "\n")
  conv.stad <- FALSE
  ## starting values
  if (i==1){
  ## start from 1 if first year
  start.value <- c(1,1,0,0)
  }else{
  ## start from previously estimated pars
  start.value <- c(bL_c1[i-1],bU_c1[i-1],cL_c1[i-1],dL_c1[i-1])
  }
  start.value <- c(1,1,0,0)
  ## MLE
  opt <- optim(par=start.value, fn=MLE_obj_FUN_Cohort,x=x,xA=xA,
               Mstand=Mstand, shat=s_c1[i],
               xlo=XLO1[[i]], xup=XUP1[[i]],
               coeff.stand=coeff_Stand, 
               Dx=Z1[,i], Ex=E1[,i], wei=W1[,i],
               xmin=xminA, xmax=xmaxA,
               ndx=ndxA, deg=deg)
  if (opt$convergence != 0) break
  ## assign
  shat <- s_c1[i]  
  bLhat <- opt$par[1]
  bUhat <- opt$par[2]
  cLhat <- opt$par[3]
  dLhat <- opt$par[4]

  ## compute dx:
  ## segment a linear transformation function
  ## below the mode
  x.low <- XLO1[[i]]-shat-Mstand
  wL <- Mstand + bLhat*x.low + cLhat*(x.low^2) + dLhat*(x.low^3)
  ## above the mode
  x.up <- XUP1[[i]]-shat-Mstand
  wU <- Mstand + bUhat*x.up
  ## unique transformation function
  wb <- c(wL, wU)
  ## B-splines on transformed ages
  Bwb <- MortSmooth_bbase(x=c(wb),
                          xminA,xmaxA,ndx=ndxA,deg=deg)
  ## transformed density
  fwb <- as.vector(exp(Bwb%*%coeff_Stand))
  fwb <- fwb[xA%in%x]
  fwb <- fwb/sum(fwb)
  ## hazard
  eta <- log(mx_from_dx(fwb))
  ## save parameters, density and mx
  bL_c1[i] <- bLhat
  bU_c1[i] <- bUhat
  cL_c1[i] <- cLhat
  dL_c1[i] <- dLhat
  DXcstad_c1[,i] <- fwb
  LMXcstad_c1[,i] <- eta
  ## plotting
  if(PLOT){
    plot(xo,LMX.act[,i],pch=16,main=c1[i],ylim = range(LMX.act[,y%in%c1],finite=T),
         xlim=range(xo),xlab="Age")
    lines(xs,LMX.smo1D[,i],col=2,lwd=2,lty=2)
    lines(xs,LMX.smo2D[,i],col=3,lwd=2,lty=2)
    lines(x,LMXcstad_c1[,i],col=4,lwd=2,lty=1)
    abline(v=xo[mo],lty=3);abline(v=M_c1[i],lty=3)
    legend("topleft",c("Observed","Smooth 1D","Smooth 2D","C-STAD"),pch = c(16,NA,NA,NA),
           lty = c(0,2,2,1),col = c(1,2,3,4),bty="n",cex = 1.3,lwd = 3)
    # locator(1)
  }
}

## plot estimated parameters
par(mfrow=c(2,3))
plot(c1,s_c1,t="o",pch=16,main="s",cex.main=2,ylab="")
plot(c1,bL_c1,t="o",pch=16,cex.main=2,main=expression(b[L]),ylab="")
plot(c1,bU_c1,t="o",pch=16,cex.main=2,main=expression(b[U]),ylab="")
plot(c1,cL_c1,t="o",pch=16,cex.main=2,main=expression(c[L]),ylab="")
plot(c1,dL_c1,t="o",pch=16,cex.main=2,main=expression(d[L]),ylab="")
par(mfrow=c(1,1))

## actual and fitted life expectancy & gini
e40_c1 <- apply(exp(LMXcstad_c1[1:mo,]),2,lifetable.ex,x=xo,sex=sex)
g40_c1 <- apply(exp(LMXcstad_c1[1:mo,]),2,GINI_func,ages=xo,sex=sex)

par(mfrow = c(1,2))
plot(c1,e40.act,ylim=range(e40.act,e40_c1),pch=16,xlab="cohort",
     main=paste0("E",xo[1]),ylab="",cex.lab=1.25)
points(c1,e40_2D[y%in%c1],col=5,pch=4,lwd=2)
points(c1,e40_c1,col=2,pch=4,lwd=2)
legend("topleft",c("Observed","2D","C-STAD"),pch = c(16,4,4),
       lty = c(0,0,0),col = c(1,5,2),bty="n",cex = 1.3,lwd = 3)
plot(c1,g40.act,ylim=range(g40.act,g40_c1),pch=16,xlab="cohort",
     main=paste0("G",xo[1]),ylab="",cex.lab=1.25)
points(c1,g40_2D[y%in%c1],col=5,pch=4,lwd=2)
points(c1,g40_c1,col=2,pch=4,lwd=2)
par(mfrow = c(1,1))

## -- C2 ----------

## select wide cohort range
c.fore <- (c1[n1]+1):y[n]   ## 1970 = last cohort with 7 observed ages (need to estimate three parameters)
n.fore <- length(c.fore)
coly.fore <- rainbow(n.fore)

## get observed Exp, Deaths and Mx (incomplete cohorts)
E.fore <- cE[ages%in%xo,cohorts%in%c.fore]
MX.fore <- cMx[ages%in%xo,cohorts%in%c.fore]
LMX.fore <- log(MX.fore)
Z.fore <- cZ[ages%in%xo,cohorts%in%c.fore]
Z.fore[E.fore==0] <- 0
W.fore <- matrix(1,mo,n.fore)
W.fore[E.fore==0] <- 0
W.fore[is.na(E.fore)] <- 0  ## zero weight where data is missing

## plot
matplot(xo,LMX.act,t="l",lty=1,col = coly,ylim = range(LMX.act,LMX.fore,finite=T),
        main="Log-Mortality, Observed, full and partial cohorts")
matlines(xo, LMX.fore, lty=1, t="l", col=coly.fore)

## derive LT densities (truncated)
DX.fore <- matrix(NA,mo,n.fore)
i <- 1
for (i in 1:n.fore){
  mx <- MX.fore[,i]
  if (any(is.na(mx))){
    whi <- which(is.na(mx))
  }else{
    whi <- length(mx)
  }
  x.run <- xo[-whi]
  mx.run <- mx[-whi]
  dx.run <- dx_from_mx(age=x.run,mx=mx.run)
  DX.fore[-whi,i] <- dx.run
}
matplot(xo, DX.fore, lty=1, t="l", col=coly.fore,xlab="Age",cex.lab=1.25,
        main="LT death distributions, partial cohorts")

## compute raw mode of LT dx 
## (to be compared afterwards with Mode from smooth rates)
Mraw <- xo[apply(DX.fore,2,which.max)]
plot(c.fore,Mraw)

## smooth 1D
## expand data for 1D smoothing
E.foreA <- rbind(E.fore,matrix(100,m-mo,n.fore))
Z.foreA <- rbind(Z.fore,matrix(100,m-mo,n.fore))
W.foreA <- rbind(W.fore,matrix(0,m-mo,n.fore))

## find unsmooth M and s
LMX.foreSMO1D <- FX.foreSMO1D <-  matrix(NA,ms,n.fore)
for (i in 1:n.fore){
  lmx.smo <- lmx_smooth(age = x,y = Z.foreA[,i],e = E.foreA[,i],w =  W.foreA[,i],
                          ndx = ndx,deg = deg)
  LMX.foreSMO1D[,i] <- Bxs %*% lmx.smo$coef
  FX.foreSMO1D[,i] <- dx_from_mx(age=xs,mx=exp(LMX.foreSMO1D[,i]))
}

## compute modal age at death 
M2.fore <- xs[apply(FX.foreSMO1D, 2, which.max)] + (delta/2)

## find cut-off cohort, last cohort where we observe the mode
## e.g. the observed LT distribution decreases after the peak 
## (c tilde in the paper)
if(cou == "DNK"){
  c.cutoff <- 1927  ## 1927 for DNK 
}else if(cou == "SWE"){
  c.cutoff <- 1925  ## 1925 
}

## check if cutoff is fine 
## (we want at least two points abobe M to estimate bU)
plot(xo,DX.fore[,which(c.fore==c.cutoff)],pch=16,ylab="fx",
     ylim=range(0,DX.fore[,which(c.fore==c.cutoff)],finite=T))
lines(xs,FX.foreSMO1D[,which(c.fore==c.cutoff)],col=2,lwd=2)
abline(v=M2.fore[which(c.fore==c.cutoff)])

## divide forecast period in c2 and c3
c2 <- c.fore[c.fore<=c.cutoff]
c3 <- c.fore[!c.fore%in%c2]
n2 <- length(c2)
n3 <- length(c3)

## data in c2
E2 <- cE[ages%in%xo,cohorts%in%c2]
MX2 <- cMx[ages%in%xo,cohorts%in%c2]
LMX2 <- log(MX2)
Z2 <- cZ[ages%in%xo,cohorts%in%c2]
Z2[E2==0] <- 0
W2 <- matrix(1,mo,n2)
W2[E2==0] <- 0
W2[is.na(E2)] <- 0  ## zero weight where data is missing

## observed modes in c2 (from 1D smooth)
M_c2 <- M2.fore[c.fore%in%c2]
plot(c1, M_c1, t="o", lwd=2, pch=16, main="Modal Age at Death",xlab="cohort",
     xlim=range(c1,c2),ylim=range(M_c1,M_c2))
points(c2,M_c2,col=2,pch=16,t="o")

## Shifting parameter in c2
s_c2 <- M_c2 - Mstand

## MLE ESTIMATION 
## Break point of the age axis for each year 
PLO2 <- PUP2 <- XLO2 <- XUP2 <- list()
for(i in 1:n2){
  PLO2[[i]] <- which(xA<=M_c2[i])
  PUP2[[i]] <- which(xA>M_c2[i])
  XLO2[[i]] <- xA[PLO2[[i]]]
  XUP2[[i]] <- xA[PUP2[[i]]]
}

## empty vectors and matrices to store results 
PLOT=F
bL_c2 <- bU_c2 <- cL_c2 <- dL_c2 <- numeric(n2)
DXcstad_c2 <- LMXcstad_c2 <- matrix(NA, m, n2)

## Estimation
i <- 1
for(i in 1:n2){
  cat("fitting cohort", c2[i], "\n")
  conv.stad <- FALSE
  ## starting values
  if (i==1){
    ## start from last value of c1 in first year 
    start.value <- c(bL_c1[n1],bU_c1[n1],cL_c1[n1],dL_c1[n1])
  }else{
    ## start from previously estimated pars
    start.value <- c(bL_c2[i-1],bU_c2[i-1],cL_c2[i-1],dL_c2[i-1])
  }
  start.value <- c(1,1,0,0)
  ## MLE
  opt <- optim(par=start.value, fn=MLE_obj_FUN_Cohort,x=x,xA=xA,
               Mstand=Mstand, shat=s_c2[i],
               xlo=XLO2[[i]], xup=XUP2[[i]],
               coeff.stand=coeff_Stand, 
               Dx=Z2[,i], Ex=E2[,i], wei=W2[,i],
               xmin=xminA, xmax=xmaxA,
               ndx=ndxA, deg=deg,control = list(maxit=1000))
  if (opt$convergence != 0) break
  
  ## assign 
  shat <- s_c2[i]
  bLhat <- opt$par[1]
  bUhat <- opt$par[2]
  cLhat <- opt$par[3]
  dLhat <- opt$par[4]
  
  ## compute dx:
  ## segment a linear transformation function
  ## below the mode
  x.low <- XLO2[[i]]-shat-Mstand
  wL <- Mstand + bLhat*x.low + cLhat*(x.low^2) + dLhat*(x.low^3)
  ## above the mode
  x.up <- XUP2[[i]]-shat-Mstand
  wU <- Mstand + bUhat*x.up
  ## unique transformation function
  wb <- c(wL, wU)
  ## B-splines on transformed ages
  Bwb <- MortSmooth_bbase(x=c(wb),
                          xminA,xmaxA,ndx=ndxA,deg=deg)
  ## transformed density
  fwb <- as.vector(exp(Bwb%*%coeff_Stand))
  fwb <- fwb[xA%in%x]
  fwb <- fwb/sum(fwb)

  ## hazard
  eta <- log(mx_from_dx(fwb))
  ## save parameters, density and mx
  bL_c2[i] <- bLhat
  bU_c2[i] <- bUhat
  cL_c2[i] <- cLhat
  dL_c2[i] <- dLhat
  DXcstad_c2[,i] <- fwb
  LMXcstad_c2[,i] <- eta
  ## plotting
  if(PLOT){
    plot(xo,LMX2[,i],pch=16,main=c2[i],ylim = range(LMX2,LMXcstad_c2,finite=T),
         xlim=range(xo),xlab="Age")
    lines(xs,LMX.foreSMO1D[,i],col=2,lwd=2,lty=2)
    lines(xs,LMX.smo2D[,which(y==c2[i])],col=3,lwd=2,lty=2)
    lines(x,LMXcstad_c2[,i],col=4,lwd=2,lty=1)
    abline(v=M_c2[i])
    abline(v=xo[mo],lty=3)
    legend("topleft",c("Observed","Smooth 1D","Smooth 2D","C-STAD"),pch = c(16,NA,NA,NA),
           lty = c(0,2,2,1),col = c(1,2,3,4),bty="n",cex = 1.3,lwd = 3)
    # locator(1)
  }
}

## plot estimated parameters
if (PLOT){
  par(mfrow=c(2,3))
  plot(c1,s_c1,t="o",pch=16,main="s",cex.main=2,ylab="",
       ylim=range(s_c1,s_c2),xlim=range(c1,c2))
  points(c2,s_c2,col=2,pch=16,t="o")
  plot(c1,bL_c1,t="o",pch=16,main=expression(b[L]),cex.main=2,ylab="",
       ylim=range(bL_c1,bL_c2),xlim=range(c1,c2))
  points(c2,bL_c2,col=2,pch=16,t="o")
  plot(c1,bU_c1,t="o",pch=16,main=expression(b[U]),cex.main=2,ylab="",
       ylim=range(bU_c1,bU_c2),xlim=range(c1,c2))
  points(c2,bU_c2,col=2,pch=16,t="o")
  plot(c1,cL_c1,t="o",pch=16,main=expression(c[L]),cex.main=2,ylab="",
       ylim=range(cL_c1,cL_c2),xlim=range(c1,c2))
  points(c2,cL_c2,col=2,pch=16,t="o")
  plot(c1,dL_c1,t="o",pch=16,main=expression(d[L]),cex.main=2,ylab="",
       ylim=range(dL_c1,dL_c2),xlim=range(c1,c2))
  points(c2,dL_c2,col=2,pch=16,t="o")
  par(mfrow=c(1,1))
}

## forecast life expectancy & gini
e40_c2 <- apply(exp(LMXcstad_c2[1:mo,]),2,lifetable.ex,x=xo,sex=sex)
g40_c2 <- apply(exp(LMXcstad_c2[1:mo,]),2,GINI_func,ages=xo,sex=sex)

if (PLOT){
  par(mfrow = c(1,2),mar=c(4.5,3,1,1))
  ## bottom, left, top, right
  plot(c1,e40.act,ylim=range(e40.act,e40_c1,e40_c2),pch=16,
       main=paste0("E",xo[1]),xlab="Cohorts",ylab="",cex.lab=1.25,xlim=range(c1,c2))
  lines(y[y%in%c(c1,c2)],e40_2D[y%in%c(c1,c2)],col=5,lwd=2)
  points(c1,e40_c1,col=2,pch=4,lwd=2)
  lines(c2,e40_c2,col=2,pch=4,lwd=2)
  abline(v=(c1[n1]+0.5),lty=3)
  legend("topleft",c("Observed","2D smo", "C-STAD"),pch = c(16,4,4),
         lty = c(0,0,0),col = c(1,5,2),bty="n",cex = 0.9,lwd = 2)
  plot(c1,g40.act,ylim=range(g40.act,g40_c1,g40_c2),pch=16,
       main=paste0("G",xo[1]),xlab="Cohorts",ylab="",cex.lab=1.25,xlim=range(c1,c2))
  lines(y[y%in%c(c1,c2)],g40_2D[y%in%c(c1,c2)],col=5,lwd=2)
  points(c1,g40_c1,col=2,pch=4,lwd=2)
  abline(v=(c1[n1]+0.5),lty=3)
  lines(c2,g40_c2,col=2,pch=4,lwd=2)
  par(mfrow = c(1,1))
}

## -- C3 ----------

## extend ts of parameters
c12 <- c(c1,c2)
n12 <- length(c12)
s_c12 <- c(s_c1,s_c2)
bU_c12 <- c(bU_c1,bU_c2)
bL_c12 <- c(bL_c1,bL_c2)
cL_c12 <- c(cL_c1,cL_c2)
dL_c12 <- c(dL_c1,dL_c2)

## VAR for (S, BU) in first differences
s_diff <- diff(s_c12)
bU_diff <- diff(bU_c12)
(COR_s_bU_diff <- cor(s_diff,bU_diff))  ## extremely high!!

df_var <- data.frame(s1=s_diff,bU1=bU_diff)
df_var.ts <- ts(df_var, start = y[1])
var_m1 <- VAR(y=df_var.ts,p = 1,type = "const")
var_resid <- residuals(var_m1)
var_fitted <- fitted(var_m1)
pred.var <- predict(var_m1,n.ahead = n3)

s_diff_fore <- pred.var$fcst$s1[,1]
s_c3 <- s_c12[n12]+cumsum(s_diff_fore)
bU_diff_fore <- pred.var$fcst$bU1[,1]
bU_c3 <- bU_c12[n12]+cumsum(bU_diff_fore)
M_c3 <- Mstand+s_c3

if (PLOT){
  par(mfrow=c(1,3))
  plot(y,M_2D,main="Mode",
       ylim=range(M_2D,M_c3))
  points(c1,M_c1 , t="o", lwd=2, pch=16)
  points(c2,M_c2,col=2,pch=16,t="o")
  points(c3,M_c3,col=4,pch=16,t="o")
  
  plot(c(c1,c2),s_c12,main="S",
       ylim=range(s_c12,s_c3),xlim=range(y))
  lines(c3,s_c3,col=4,pch=16,t="o",lwd=2)
  
  plot(c(c1,c2),bU_c12,main="bU",
       ylim=range(bU_c12,bU_c3),xlim=range(y))
  lines(c3,bU_c3,col=4,pch=16,t="o")
  par(mfrow=c(1,1))
}

## data in c3
E3 <- cE[ages%in%xo,cohorts%in%c3]
MX3 <- cMx[ages%in%xo,cohorts%in%c3]
LMX3 <- log(MX3)
Z3 <- cZ[ages%in%xo,cohorts%in%c3]
Z3[E3==0] <- 0
W3 <- matrix(1,mo,n3)
W3[E3==0] <- 0
W3[is.na(E3)] <- 0  ## zero weight where data is missing

## MLE ESTIMATION 
## Break point of the age axis for each year 
PLO3 <- PUP3 <- XLO3 <- XUP3 <- list()
for(i in 1:n3){
  PLO3[[i]] <- which(xA<=M_c3[i])
  PUP3[[i]] <- which(xA>M_c3[i])
  XLO3[[i]] <- xA[PLO3[[i]]]
  XUP3[[i]] <- xA[PUP3[[i]]]
}

## acceptable range for LOW parameters
pc.change <- 0.75
rangeB <- diff(range(bL_c12))
rangeC <- diff(range(cL_c12))
rangeD <- diff(range(dL_c12))
bL_min <- min(bL_c12)-pc.change*rangeB
bL_max <- max(bL_c12)+pc.change*rangeB
cL_min <- min(cL_c12)-pc.change*rangeC
cL_max <- max(cL_c12)+pc.change*rangeC
dL_min <- min(dL_c12)-pc.change*rangeD
dL_max <- max(dL_c12)+pc.change*rangeD

par(mfrow=c(1,3))
plot(bL_c12,ylim=range(bL_min,bL_max),t="l");abline(h=bL_min);abline(h=bL_max)
plot(cL_c12,ylim=range(cL_min,cL_max),t="l");abline(h=cL_min);abline(h=cL_max)
plot(dL_c12,ylim=range(dL_min,dL_max),t="l");abline(h=dL_min);abline(h=dL_max)
par(mfrow=c(1,1))

## empty vectors and matrices to store results 
PLOT=F
bL_c3 <- cL_c3 <- dL_c3 <- rep(NA,n3)
DXcstad_c3 <- LMXcstad_c3 <- matrix(NA, m, n3)

## Estimation 
i <- 1
n.grid.fit <- 500
for(i in 1:n3){
  cat("fitting cohort", c3[i], "\n")
  ## declare explicitely inputs of cleversearch function
  xlo_2 <- XLO3[[i]]
  shat <- s_c3[i]
  bUhat <- bU_c3[i]
  xup_2 <- XUP3[[i]]
  coeff.stand <- coeff_Stand 
  Dx_2 <- Z3[,i] 
  Ex_2 <- E3[,i]
  wei_2 <- W3[,i]
  MLE_obj_FUN_Cohort_LOW_CLEVER <- function(pars){
    ## starting b, c, d LOW
    bL <- pars[1]
    cL <- pars[2]
    dL <- pars[3]
    ## segment a linear transformation function
    ## below the mode
    x.low <- xlo_2-shat-Mstand
    wL <- Mstand + bL*x.low + cL*(x.low^2) + dL*(x.low^3)
    ## above the mode
    x.up <- xup_2-shat-Mstand
    wU <- Mstand + bUhat*x.up
    ## unique transformation function
    wb <- c(wL, wU)
    ## B-splines on transformed ages
    Bwb <- MortSmooth_bbase(x=c(wb),
                            xminA,xmaxA,ndx=ndxA,deg=deg)
    ## transformed density
    dwb <- as.vector(exp(Bwb%*%coeff.stand))
    dwb <- dwb[xA%in%x]
    dwb <- dwb/sum(dwb)
    ## hazard
    eta <- log(mx_from_dx(dx=dwb))
    mu <- exp(eta)
    ## minimise minus the Log-Likelihood (maximise the LL)
    Lij <- -sum(wei_2*(Dx_2 * eta[1:length(Dx_2)]-Ex_2*mu[1:length(Ex_2)]), na.rm = T)
    return(Lij)
  }
  ## start values
  if (i==1){
    start.value <- c(bL_c2[n2],cL_c2[n2],dL_c2[n2])
  }else{
    start.value <- c(bL_c3[i-1],cL_c3[i-1],dL_c3[i-1])
  }
  ## optimize
  opt <- cleversearch(fn = MLE_obj_FUN_Cohort_LOW_CLEVER, 
                      lower = c(bL_min,cL_min,dL_min),
                      upper = c(bL_max,cL_max,dL_max), 
                      startvalue = start.value, 
                      ngrid=n.grid.fit, logscale=FALSE)
  
  ## assign 
  shat <- s_c3[i]
  bUhat <- bU_c3[i]
  bLhat <- opt$par[1]
  cLhat <- opt$par[2]
  dLhat <- opt$par[3]
  
  ## C-STAD model
  ## compute dx:
  ## segment a linear transformation function
  ## below the mode
  x.low <- XLO3[[i]]-shat-Mstand
  wL <- Mstand + bLhat*x.low + cLhat*(x.low^2) + dLhat*(x.low^3)
  ## above the mode
  x.up <- XUP3[[i]]-shat-Mstand
  wU <- Mstand + bUhat*x.up
  ## unique transformation function
  wb <- c(wL, wU)
  ## B-splines on transformed ages
  Bwb <- MortSmooth_bbase(x=c(wb),
                          xminA,xmaxA,ndx=ndxA,deg=deg)
  ## transformed density
  fwb <- as.vector(exp(Bwb%*%coeff_Stand))
  fwb <- fwb[xA%in%x]
  fwb <- fwb/sum(fwb)
  
  ## hazard
  eta <- log(mx_from_dx(fwb))
  ## save parameters, density and mx
  bL_c3[i] <- bLhat
  cL_c3[i] <- cLhat
  dL_c3[i] <- dLhat
  DXcstad_c3[,i] <- fwb
  LMXcstad_c3[,i] <- eta
  ## plotting
  if (PLOT==TRUE){
    plot(xo,LMX3[,i],pch=16,main=c3[i],ylim = range(LMX3,LMXcstad_c3,finite=T),
         xlim=range(xo),xlab="Age")
    lines(xs,LMX.smo2D[,which(y==c3[i])],col=3,lwd=2,lty=2)
    lines(x,LMXcstad_c3[,i],col=4,lwd=2,lty=1)
    abline(v=M_c3[i])
    abline(v=xo[mo],lty=3)
    legend("topleft",c("Observed","Smo 2D","C-STAD"),pch = c(16,NA,NA),
           lty = c(0,2,1),col = c(1,3,4),bty="n",cex = 1.3,lwd = 3)
    # locator(1)
  }
}

## plot forecast/estimated parameters
if (PLOT){
  par(mfrow=c(2,3))
  plot(c1,s_c1,t="o",pch=16,main="s",cex.main=2,ylab="",
       ylim=range(s_c1,s_c2,s_c3),xlim=range(y))
  points(c2,s_c2,col=2,pch=16,t="o")
  points(c3,s_c3,col=4,pch=16,t="o")
  plot(c1,bU_c1,t="o",pch=16,main=expression(b[U]),cex.main=2,ylab="",
       ylim=range(bU_c1,bU_c2,bU_c3),xlim=range(y))
  points(c2,bU_c2,col=2,pch=16,t="o")
  points(c3,bU_c3,col=4,pch=16,t="o")
  plot(c1,bL_c1,t="n",pch=16,cex.main=2,ylab="",xlab="",
       ylim=range(bL_c1,bL_c2),xlim=range(c1,c2),axes=F)
  plot(c1,bL_c1,t="o",pch=16,main=expression(b[L]),cex.main=2,ylab="",
       ylim=range(bL_c1,bL_c2,bL_c3,bL_max,bL_min,finite=T),xlim=range(y))
  abline(h=bL_min);abline(h=bL_max)
  points(c2,bL_c2,col=2,pch=16,t="o")
  points(c3,bL_c3,col=3,pch=16,t="l")
  plot(c1,cL_c1,t="o",pch=16,main=expression(c[L]),cex.main=2,ylab="",
       ylim=range(cL_c1,cL_c2,cL_c3,cL_max,cL_min,finite=T),xlim=range(y))
  abline(h=cL_min);abline(h=cL_max)
  points(c2,cL_c2,col=2,pch=16,t="o")
  points(c3,cL_c3,col=3,pch=16,t="l")
  plot(c1,dL_c1,t="o",pch=16,main=expression(d[L]),cex.main=2,ylab="",
       ylim=range(dL_c1,dL_c2,dL_c3,dL_min,dL_max,finite=T),xlim=range(y))
  abline(h=dL_min);abline(h=dL_max)
  points(c2,dL_c2,col=2,pch=16,t="o")
  points(c3,dL_c3,col=3,pch=16,t="l")
  par(mfrow=c(1,1))
}

## forecast life expectancy & gini
e40_c3 <- apply(exp(LMXcstad_c3[1:mo,]),2,lifetable.ex,x=xo,sex=sex)
g40_c3 <- apply(exp(LMXcstad_c3[1:mo,]),2,GINI_func,ages=xo,sex=sex)

## PLOT
if (PLOT){
  par(mfrow = c(1,2))
  ## bottom, left, top, right
  plot(c1,e40.act,ylim=range(e40.act,e40_c1,e40_c3,e40_2D),pch=16,
       main=paste0("E",xo[1]),ylab="",cex.lab=1.25,xlim=range(y),xlab="Cohort")
  lines(y,e40_2D,col=5,lwd=2)
  points(c1,e40_c1,col=2,pch=4,lwd=2)
  lines(c2,e40_c2,col=2,pch=4,lwd=2)
  lines(c3,e40_c3,col=3,pch=4,lwd=2)
  abline(v=c2[1]);abline(v=c3[1])
  plot(c1,g40.act,ylim=range(g40.act,g40_c1,g40_c3),pch=16,
       main=paste0("G",xo[1]),ylab="",cex.lab=1.25,xlim=range(y),xlab="Cohort")
  lines(y,g40_2D,col=5,lwd=2)
  points(c1,g40_c1,col=2,pch=4,lwd=2)
  lines(c2,g40_c2,col=2,pch=4,lwd=2)
  lines(c3,g40_c3,col=3,pch=4,lwd=2)
  par(mfrow = c(1,1))
}

## overall CSTAD results
all_cohorts <- c(c1,c2,c3)
FXcstad <- cbind(DXcstad_c1,DXcstad_c2,DXcstad_c3)
LMXcstad <- cbind(LMXcstad_c1,LMXcstad_c2,LMXcstad_c3)
s_cstad <- c(s_c12,s_c3)
bU_cstad <- c(bU_c12,bU_c3)
bL_cstad <- c(bL_c12,bL_c3)
cL_cstad <- c(cL_c12,cL_c3)
dL_cstad <- c(dL_c12,dL_c3)
e40_cstad <- c(e40_c1,e40_c2,e40_c3)
g40_cstad <- c(g40_c1,g40_c2,g40_c3)
ADDobs <- apply(exp(LMX.act[1:mo,]), 2, lifetable.dx, x=xo, sex="F")
ADDcstad <- apply(exp(LMXcstad), 2, lifetable.dx, x=x, sex="F")

## -- CONSTRUCT CI via BOOTSTRAP ----------
set.seed(2019)
n.boot <- 10
n.simul <- 10
col.sim <- rainbow_hcl(n.simul)
col.boot <- rainbow_hcl(n.boot)

## FITTED DEATHS & Poisson residuals
DXcstad <- E*exp(LMXcstad[1:mo,])
res.t1 <- sign(Zna-DXcstad)
res.t2 <- Zna * log(ifelse(Zna==0,1e-8,Zna)/ifelse(DXcstad==0,1e-8,DXcstad))
res.t3 <- Zna - DXcstad
DevRes <- res.t1 * sqrt(2) * (res.t2 - res.t3)^0.5

## STORING bootstrap RESULTS
M <- S <- BU <- BL <- CL <- DL <- E40 <- G40 <- array(NA,dim=c(n,n.boot,n.simul))
ETAS <- ADD <-  array(NA,dim=c(m,n,n.boot,n.simul))

t_start <- Sys.time()
boot.num <- 1
sim.num <- 1
for (boot.num in 1:n.boot){
  cat("bootstrap number",boot.num,"\n")
  
  ## residual bootstrap 
  DXboot <- matrix(NA,mo,n)
  j <- 1
  for (j in 1:n){
    DXboot[,j] <- genPoissonResBoot(Zobs=matrix(Zna[,j],mo,1),
                                    Zhat=matrix(DXcstad[,j],mo,1),
                                    rescale.dths = F, centre.res = F)
  }

  ## refit CSTAD ## refit CSTAD
  ## refit CSTAD ## refit CSTAD
  
  ## NEW STANDARD
  ZA_boot <- rbind(DXboot,matrix(100,(m-mo),n))
  smooth2D_boot <- Mort2Dsmooth(x=x,y=y,Z=ZA_boot,offset=log(EA),W=WA,
                                ndx=c(ndx,ndy),method = 3,
                                lambdas = c(lambdaX.hat,lambdaY.hat))
  LMX.smo2D_boot <- matrix(Bs %*% c(smooth2D_boot$coefficients),ms,n)
  LMX.smo2DW_boot <- LMX.smo2D_boot*Ws
  LMX.smo2DW_boot[LMX.smo2DW_boot==0] <- NA
  
  ## Fx and alignment
  FX.smo2D_boot <- apply(exp(LMX.smo2D_boot),2,dx_from_mx,age=xs)
  FX.smo2DW_boot <- FX.smo2D_boot*Ws
  FX.smo2DW_boot[FX.smo2DW_boot==0] <- NA
  
  ## compute modal age at death (it's pretty smooth here due to 2D smoothing)
  M_2D_boot <- xs[apply(FX.smo2D_boot, 2, which.max)] + (delta/2)
  Mstand_boot <- M_2D_boot[1]
  
  ## STANDARD DISTRIBUTION
  s_2D_boot <- M_2D_boot - Mstand_boot
  
  ## derive aligned distributions
  FX.align_boot <- matrix(0, nrow=ms, ncol=n)
  for(i in 1:n){
    FX.align_boot[,i] <- fx_shift(age=xs,fx=FX.smo2D_boot[,i],shift=-s_2D_boot[i],ndx = ndx,deg = deg)
  }
  FX.alignW_boot <- FX.align_boot*Ws
  FX.alignW_boot[FX.alignW_boot==0] <- NA
  
  ## Standard = mean of the aligned densities
  FXallmeanW_boot <- exp(apply(log(FX.alignW_boot), 1, mean, na.rm=T))
  FXstand_boot <- FXallmeanW_boot
  
  ## coefficients & mode
  Standard_boot <- coeff_stand(age=xs,fx=FXstand_boot,ndx=ndxA,deg=deg,
                               ages.add.l=ages.add.l,ages.add.r=ages.add.r,
                               lambda = lambdaSTAND)
  coeff_Stand_boot <- Standard_boot$betasA
  Mstand_boot <- xAs[which.max(BAs%*%coeff_Stand_boot)]  ## reduced bias for CHE
  ## C-STAD ESTIMATION on c1
  Z1_boot <- DXboot[,1:n1]
  
  ## first way to compute mode: BIC
  M_c1_boot <- rep(NA,n1)
  
  ## find unsmooth M and s
  LMX.smo1D_boot <- FX.smo1D_boot <-  matrix(NA,ms,n1)
  i <- 1
  for (i in 1:n1){
    lmx.smo <- lmx_smooth(age = x,y = Z1A[,i],e = E1A[,i],w =  W1A[,i],
                            ndx = ndx,deg = deg)
    LMX.smo1D_boot[,i] <- Bxs %*% lmx.smo$coef
    FX.smo1D_boot[,i] <- dx_from_mx(age=xs,mx=exp(LMX.smo1D_boot[,i]))
  }
  
  ## compute modal age at death 
  M_c1_boot <- xs[apply(FX.smo1D_boot, 2, which.max)] + (delta/2)
  s_c1_boot <- M_c1_boot - Mstand_boot

  ## Break point of the age axis for each year 
  PLO1_boot <- PUP1_boot <-  list()
  XLO1_boot <- XUP1_boot <- list()
  for(i in 1:n1){
    PLO1_boot[[i]] <- which(xA<=M_c1_boot[i])
    PUP1_boot[[i]] <- which(xA>M_c1_boot[i])
    XLO1_boot[[i]] <- xA[PLO1_boot[[i]]]
    XUP1_boot[[i]] <- xA[PUP1_boot[[i]]]
  }
  
  ## empty vectors and matrices to store results 
  PLOT=F
  bL_c1_boot <- bU_c1_boot <- cL_c1_boot <- dL_c1_boot <- numeric(n1)
  DXcstad_c1_boot <- LMXcstad_c1_boot <-  matrix(NA, m, n1)
  
  ## Estimation
  cat("fitting cohorts c1", "\n")
  for(i in 1:n1){
    conv.stad <- FALSE
    ## starting values
    if (i==1){
      ## start from 1 if first year
      start.value <- c(1,1,0,0)
    }else{
      ## start from previously estimated pars
      start.value <- c(bL_c1_boot[i-1],bU_c1_boot[i-1],
                       cL_c1_boot[i-1],dL_c1_boot[i-1])
    }
    start.value <- c(1,1,0,0)
    ## MLE
    opt <- optim(par=start.value, fn=MLE_obj_FUN_Cohort,x=x,xA=xA,
                 Mstand=Mstand_boot, shat=s_c1_boot[i],
                 xlo=XLO1_boot[[i]], xup=XUP1_boot[[i]],
                 coeff.stand=coeff_Stand_boot, 
                 Dx=Z1_boot[,i], Ex=E1[,i], wei=W1[,i],
                 xmin=xminA, xmax=xmaxA,
                 ndx=ndxA, deg=deg)
    # if (opt$convergence != 0) break
    ## assign
    shat <- s_c1_boot[i]  
    bLhat <- opt$par[1]
    bUhat <- opt$par[2]
    cLhat <- opt$par[3]
    dLhat <- opt$par[4]
    
    ## compute dx:
    ## segment a linear transformation function
    ## below the mode
    x.low <- XLO1_boot[[i]]-shat-Mstand_boot
    wL <- Mstand_boot + bLhat*x.low + cLhat*(x.low^2) + dLhat*(x.low^3)
    ## above the mode
    x.up <- XUP1_boot[[i]]-shat-Mstand_boot
    wU <- Mstand_boot + bUhat*x.up
    ## unique transformation function
    wb <- c(wL, wU)
    ## B-splines on transformed ages
    Bwb <- MortSmooth_bbase(x=c(wb),
                            xminA,xmaxA,ndx=ndxA,deg=deg)
    ## transformed density
    fwb <- as.vector(exp(Bwb%*%coeff_Stand_boot))
    fwb <- fwb[xA%in%x]
    fwb <- fwb/sum(fwb)
    
    ## hazard
    eta <- log(mx_from_dx(fwb))
    ## save parameters, density and mx
    bL_c1_boot[i] <- bLhat
    bU_c1_boot[i] <- bUhat
    cL_c1_boot[i] <- cLhat
    dL_c1_boot[i] <- dLhat
    DXcstad_c1_boot[,i] <- fwb
    LMXcstad_c1_boot[,i] <- eta
    ## plotting
    if(PLOT){
      plot(xo,LMX.act[,i],pch=16,main=c1[i],ylim = range(LMXcstad_c1,LMX.act,finite=T),
           xlim=range(xo),xlab="Age")
      lines(x,LMXcstad_c1[,i],col=4,lwd=2,lty=1)
      lines(x,LMXcstad_c1_boot[,i],col=5,lwd=2,lty=1)
      abline(v=xo[mo],lty=3)
      legend("topleft",c("Obs","C-STAD","Boot"),pch = c(16,NA,NA),
             lty = c(0,1,1),col = c(1,4,5),bty="n",cex = 1.3,lwd = 3)
      # locator(1)
    }
  }
  
  ## refit CSTAD in c2 ## refit CSTAD in c2
  ## refit CSTAD in c2 ## refit CSTAD in c2
  
  ## smooth 1D c2
  ## expand data for 1D smoothing
  Z2_boot <- DXboot[,1:n2+n1]
  
  ## find unsmooth M and s
  LMX.foreSMO1D_boot <- FX.foreSMO1D_boot <-  matrix(NA,ms,n2)
  i <- 1
  for (i in 1:n2){
    z2boot <- c(Z2_boot[,i],rep(99,10))
    lmx.smo <- lmx_smooth(age = x,y = Z.foreA[,i],e = E.foreA[,i],w =  W.foreA[,i],
                            ndx = ndx,deg = deg)
    LMX.foreSMO1D_boot[,i] <- Bxs %*% lmx.smo$coef
    FX.foreSMO1D_boot[,i] <- dx_from_mx(age=xs,mx=exp(LMX.foreSMO1D_boot[,i]))
  }
  
  ## compute modal age at death 
  M_c2_boot <- xs[apply(FX.foreSMO1D_boot, 2, which.max)] + (delta/2)
  s_c2_boot <- M_c2_boot - Mstand_boot
  
  ## MLE ESTIMATION 
  ## Break point of the age axis for each year 
  PLO2_boot <- PUP2_boot <- XLO2_boot <- XUP2_boot <- list()
  for(i in 1:n2){
    PLO2_boot[[i]] <- which(xA<=M_c2_boot[i])
    PUP2_boot[[i]] <- which(xA>M_c2_boot[i])
    XLO2_boot[[i]] <- xA[PLO2_boot[[i]]]
    XUP2_boot[[i]] <- xA[PUP2_boot[[i]]]
  }
  
  ## empty vectors and matrices to store results 
  PLOT=F
  bL_c2_boot <- bU_c2_boot <- cL_c2_boot <- dL_c2_boot <- numeric(n2)
  DXcstad_c2_boot <- LMXcstad_c2_boot <- matrix(NA, m, n2)
  
  ## Estimation
  cat("fitting cohorts c2", "\n")
  for(i in 1:n2){
    conv.stad <- FALSE
    ## starting values
    if (i==1){
      ## start from last value of c1 in first year 
      start.value <- c(bL_c1_boot[n1],bU_c1_boot[n1],cL_c1_boot[n1],dL_c1_boot[n1])
    }else{
      ## start from previously estimated pars
      start.value <- c(bL_c2_boot[i-1],bU_c2_boot[i-1],cL_c2_boot[i-1],dL_c2_boot[i-1])
    }
    start.value <- c(1,1,0,0)
    ## MLE
    opt <- optim(par=start.value, fn=MLE_obj_FUN_Cohort,x=x,xA=xA,
                 Mstand=Mstand_boot, shat=s_c2_boot[i],
                 xlo=XLO2_boot[[i]], xup=XUP2_boot[[i]],
                 coeff.stand=coeff_Stand_boot, 
                 Dx=Z2_boot[,i], Ex=E2[,i], wei=W2[,i],
                 xmin=xminA, xmax=xmaxA,
                 ndx=ndxA, deg=deg)
    # if (opt$convergence != 0) break
    
    ## assign 
    shat <- s_c2_boot[i]
    bLhat <- opt$par[1]
    bUhat <- opt$par[2]
    cLhat <- opt$par[3]
    dLhat <- opt$par[4]
    
    ## compute dx:
    ## segment a linear transformation function
    ## below the mode
    x.low <- XLO2_boot[[i]]-shat-Mstand_boot
    wL <- Mstand_boot + bLhat*x.low + cLhat*(x.low^2) + dLhat*(x.low^3)
    ## above the mode
    x.up <- XUP2_boot[[i]]-shat-Mstand_boot
    wU <- Mstand_boot + bUhat*x.up
    ## unique transformation function
    wb <- c(wL, wU)
    ## B-splines on transformed ages
    Bwb <- MortSmooth_bbase(x=c(wb),
                            xminA,xmaxA,ndx=ndxA,deg=deg)
    ## transformed density
    fwb <- as.vector(exp(Bwb%*%coeff_Stand_boot))
    fwb <- fwb[xA%in%x]
    fwb <- fwb/sum(fwb)
    
    ## hazard
    eta <- log(mx_from_dx(fwb))
    ## save parameters, density and mx
    bL_c2_boot[i] <- bLhat
    bU_c2_boot[i] <- bUhat
    cL_c2_boot[i] <- cLhat
    dL_c2_boot[i] <- dLhat
    DXcstad_c2_boot[,i] <- fwb
    LMXcstad_c2_boot[,i] <- eta
    ## plotting
    if(PLOT){
      plot(xo,LMX2[,i],pch=16,main=c2[i],ylim = range(LMX2,LMXcstad_c2,finite=T),
           xlim=range(xo),xlab="Age")
      lines(x,LMXcstad_c2[,i],col=4,lwd=2,lty=1)
      lines(x,LMXcstad_c2_boot[,i],col=5,lwd=2,lty=1)
      abline(v=M_c2[i])
      abline(v=xo[mo],lty=3)
      # locator(1)
    }
  }
  
  ## refit CSTAD in c3 ## refit CSTAD in c3
  ## refit CSTAD in c3 ## refit CSTAD in c3
  Z3_boot <- DXboot[,1:n3+n1+n2]
  
  ## extend ts of parameters
  s_c12_boot <- c(s_c1_boot,s_c2_boot)
  bU_c12_boot <- c(bU_c1_boot,bU_c2_boot)
  bL_c12_boot <- c(bL_c1_boot,bL_c2_boot)
  cL_c12_boot <- c(cL_c1_boot,cL_c2_boot)
  dL_c12_boot <- c(dL_c1_boot,dL_c2_boot)
  
  ## VAR (S, BU)
  s_diff_boot <- diff(s_c12_boot)
  bU_diff_boot <- diff(bU_c12_boot)
  df_var_boot <- data.frame(s1=s_diff_boot,bU1=bU_diff_boot)
  df_var_boot.ts <- ts(df_var_boot, start = y[1])
  var_m1_boot <- VAR(y=df_var_boot.ts,p = 1,type = "const")
  var_resid_boot <- residuals(var_m1_boot)
  var_fitted_boot <- fitted(var_m1_boot)

  n.grid.boot <- 250
  cat("fitting cohorts c3", "\n")
  for (sim.num in 1:n.simul){
    cat("simulation", sim.num, "\n")
    
    sample.res1 <- sample(c(var_resid_boot[,1]),size = c(n12-2),replace = T)
    sample.res2 <- sample(c(var_resid_boot[,2]),size = c(n12-2),replace = T)
    df.sim.ts <- var_fitted_boot + cbind(sample.res1,sample.res2)
    var_sim <- VAR(y=df.sim.ts,p = 1,type = "const")
    boot_fore <- predict(var_sim,n.ahead = n3)
    s1f1.boot <- boot_fore$fcst$s1[,1]
    bU1f1.boot <- boot_fore$fcst$bU1[,1]
    s_c3_boot <- s_c12_boot[n12]+cumsum(s1f1.boot)
    bU_c3_boot <- bU_c12_boot[n12]+cumsum(bU1f1.boot)
    M_c3_boot <- s_c3_boot + Mstand_boot

    ## for each simulation, estimate BL, CL, DL
    
    ## MLE ESTIMATION 
    ## Break point of the age axis for each year 
    PLO3_boot <- PUP3_boot <- XLO3_boot <- XUP3_boot <- list()
    for(i in 1:n3){
      PLO3_boot[[i]] <- which(xA<=M_c3_boot[i])
      PUP3_boot[[i]] <- which(xA>M_c3_boot[i])
      XLO3_boot[[i]] <- xA[PLO3_boot[[i]]]
      XUP3_boot[[i]] <- xA[PUP3_boot[[i]]]
    }
    
    ## acceptable range for LOW parameters
    rangeB_boot <- diff(range(bL_c12_boot))
    rangeC_boot <- diff(range(cL_c12_boot))
    rangeD_boot <- diff(range(dL_c12_boot))
    bL_min_boot <- min(bL_c12_boot)-pc.change*rangeB_boot
    bL_max_boot <- max(bL_c12_boot)+pc.change*rangeB_boot
    cL_min_boot <- min(cL_c12_boot)-pc.change*rangeC_boot
    cL_max_boot <- max(cL_c12_boot)+pc.change*rangeC_boot
    dL_min_boot <- min(dL_c12_boot)-pc.change*rangeD_boot
    dL_max_boot <- max(dL_c12_boot)+pc.change*rangeD_boot
    
    ## empty vectors and matrices to store results 
    PLOT=F
    bL_c3_boot <- cL_c3_boot <- dL_c3_boot <- rep(NA,n3)
    DXcstad_c3_boot <- LMXcstad_c3_boot <- matrix(NA, m, n3)
    
    ## Estimation 
    i <- 1
    for(i in 1:n3){
      ## declare explicitely inputs of cleversearch function
      xlo_2 <- XLO3_boot[[i]]
      shat <- s_c3_boot[i]
      bUhat <- bU_c3_boot[i]
      xup_2 <- XUP3_boot[[i]]
      coeff.stand <- coeff_Stand_boot 
      Dx_2 <- Z3_boot[,i] 
      Ex_2 <- E3[,i]
      wei_2 <- W3[,i]
      MLE_obj_FUN_Cohort_LOW_CLEVER <- function(pars){
        ## starting b, c, d LOW
        bL <- pars[1]
        cL <- pars[2]
        dL <- pars[3]
        ## segment a linear transformation function
        ## below the mode
        x.low <- xlo_2-shat-Mstand
        wL <- Mstand + bL*x.low + cL*(x.low^2) + dL*(x.low^3)
        ## above the mode
        x.up <- xup_2-shat-Mstand
        wU <- Mstand + bUhat*x.up
        ## unique transformation function
        wb <- c(wL, wU)
        ## B-splines on transformed ages
        Bwb <- MortSmooth_bbase(x=c(wb),
                                xminA,xmaxA,ndx=ndxA,deg=deg)
        ## transformed density
        dwb <- as.vector(exp(Bwb%*%coeff.stand))
        dwb <- dwb[xA%in%x]
        dwb <- dwb/sum(dwb)
        ## hazard
        eta <- log(mx_from_dx(dx=dwb))
        mu <- exp(eta)
        ## minimise minus the Log-Likelihood (maximise the LL)
        Lij <- -sum(wei_2*(Dx_2 * eta[1:length(Dx_2)]-Ex_2*mu[1:length(Ex_2)]), na.rm = T)
        return(Lij)
      }
      ## start values
      if (i==1){
        start.value <- c(bL_c2_boot[n2],cL_c2_boot[n2],dL_c2_boot[n2])
      }else{
        start.value <- c(bL_c3_boot[i-1],cL_c3_boot[i-1],dL_c3_boot[i-1])
      }
      ## optimize
      opt <- cleversearch(fn = MLE_obj_FUN_Cohort_LOW_CLEVER, 
                          lower = c(bL_min_boot,cL_min_boot,dL_min_boot),
                          upper = c(bL_max_boot,cL_max_boot,dL_max_boot), 
                          startvalue = start.value, 
                          ngrid=n.grid.boot, logscale=FALSE)
      
      ## assign 
      shat <- s_c3_boot[i]
      bUhat <- bU_c3_boot[i]
      bLhat <- opt$par[1]
      cLhat <- opt$par[2]
      dLhat <- opt$par[3]
      
      ## C-STAD model
      ## compute dx:
      ## segment a linear transformation function
      ## below the mode
      x.low <- XLO3_boot[[i]]-shat-Mstand_boot
      wL <- Mstand_boot + bLhat*x.low + cLhat*(x.low^2) + dLhat*(x.low^3)
      ## above the mode
      x.up <- XUP3_boot[[i]]-shat-Mstand_boot
      wU <- Mstand_boot + bUhat*x.up
      ## unique transformation function
      wb <- c(wL, wU)
      ## B-splines on transformed ages
      Bwb <- MortSmooth_bbase(x=c(wb),
                              xminA,xmaxA,ndx=ndxA,deg=deg)
      ## transformed density
      fwb <- as.vector(exp(Bwb%*%coeff_Stand_boot))
      fwb <- fwb[xA%in%x]
      fwb <- fwb/sum(fwb)
      
      ## hazard
      eta <- log(mx_from_dx(fwb))
      ## save parameters, density and mx
      bL_c3_boot[i] <- bLhat
      cL_c3_boot[i] <- cLhat
      dL_c3_boot[i] <- dLhat
      DXcstad_c3_boot[,i] <- fwb
      LMXcstad_c3_boot[,i] <- eta
      ## plotting
      if (PLOT==TRUE){
        plot(xo,LMX3[,i],pch=16,main=c3[i],ylim = range(LMX3,LMXcstad_c3,finite=T),
             xlim=range(xo),xlab="Age")
        lines(xs,LMX.smo2D[,which(y==c3[i])],col=3,lwd=2,lty=2)
        lines(x,LMXcstad_c3[,i],col=4,lwd=2,lty=1)
        lines(x,LMXcstad_c3_boot[,i],col=5,lwd=2,lty=1)
        abline(v=M_c3[i],col=4);abline(v=M_c3_boot[i],col=5)
        abline(v=xo[mo],lty=3)
        legend("topleft",c("Observed","Smo 2D","C-STAD"),pch = c(16,NA,NA),
               lty = c(0,2,1),col = c(1,3,4),bty="n",cex = 1.3,lwd = 3)
        # locator(1)
      }
    }
    ## save results
    M[,boot.num,sim.num] <- c(M_c1_boot,M_c2_boot,M_c3_boot)
    S[,boot.num,sim.num] <- c(s_c1_boot,s_c2_boot,s_c3_boot)
    BL[,boot.num,sim.num] <- c(bL_c1_boot,bL_c2_boot,bL_c3_boot)
    BU[,boot.num,sim.num] <- c(bU_c1_boot,bU_c2_boot,bU_c3_boot)
    CL[,boot.num,sim.num] <- c(cL_c1_boot,cL_c2_boot,cL_c3_boot)
    DL[,boot.num,sim.num] <- c(dL_c1_boot,dL_c2_boot,dL_c3_boot)
    ETAS[,,boot.num,sim.num] <- cbind(LMXcstad_c1_boot,LMXcstad_c2_boot,LMXcstad_c3_boot)
    ADD[,,boot.num,sim.num] <- apply(exp(ETAS[,,boot.num,sim.num]), 2, lifetable.dx, x=x, sex="F")
    E40[,boot.num,sim.num] <- apply(exp(ETAS[1:mo,,boot.num,sim.num]),2,lifetable.ex,x=xo,sex=sex)
    G40[,boot.num,sim.num] <- apply(exp(ETAS[1:mo,,boot.num,sim.num]),2,GINI_func,ages=xo,sex=sex)
  }
}
t_end <- Sys.time()
(t_boot <- round(t_end - t_start))
matplot(x,ADD[,,boot.num,sim.num],t="l",lty=1,lwd=0.8,col = coly)
col.sim <- rainbow_hcl(n.boot*n.simul)
Msim <- matrix(c(M),nrow=n,ncol=n.boot*n.simul)
Ssim <- matrix(c(S),nrow=n,ncol=n.boot*n.simul)
BLsim <- matrix(c(BL),nrow=n,ncol=n.boot*n.simul)
CLsim <- matrix(c(CL),nrow=n,ncol=n.boot*n.simul)
DLsim <- matrix(c(DL),nrow=n,ncol=n.boot*n.simul)
BUsim <- matrix(c(BU),nrow=n,ncol=n.boot*n.simul)
E40sim <- matrix(c(E40),nrow=n,ncol=n.boot*n.simul)
G40sim <- matrix(c(G40),nrow=n,ncol=n.boot*n.simul)

##### PLOTS -------------
lev.p <- 0.8
plot(all_cohorts,s_cstad,t="o",pch=16,lwd=2,ylim=range(s_cstad,S,finite=T))
matlines(all_cohorts,Ssim,lwd=0.8,lty=1,col = col.sim)
plot(all_cohorts,bU_cstad,t="o",pch=16,lwd=2,ylim=range(bU_cstad,BU,finite=T))
matlines(all_cohorts,BUsim,lwd=0.8,lty=1,col = col.sim)
plot(all_cohorts,bL_cstad,t="o",pch=16,lwd=2,ylim=range(bL_cstad,BL,finite=T))
matlines(all_cohorts,BLsim,lwd=0.8,lty=1,col = col.sim)
plot(all_cohorts,cL_cstad,t="o",pch=16,lwd=2,ylim=range(cL_cstad,CL,finite=T))
matlines(all_cohorts,CLsim,lwd=0.8,lty=1,col = col.sim)
plot(all_cohorts,dL_cstad,t="o",pch=16,lwd=2,ylim=range(dL_cstad,DL,finite=T))
matlines(all_cohorts,DLsim,lwd=0.8,lty=1,col = col.sim)
plot(all_cohorts,e40_cstad,t="o",pch=16,lwd=2,ylim=range(e40_cstad,E40,finite=T))
matlines(all_cohorts,E40sim,lwd=0.8,lty=1,col = col.sim)
lines(all_cohorts,apply(E40sim,1,median),col=2,lwd=2)
lines(all_cohorts,apply(E40sim,1,quantile,prob=1- (1-lev.p)/2),col=2,lwd=2,lty=2)
lines(all_cohorts,apply(E40sim,1,quantile,prob=(1-lev.p)/2),col=2,lwd=2,lty=2)
plot(all_cohorts,g40_cstad,t="o",pch=16,lwd=2,ylim=range(g40_cstad,G40,finite=T))
matlines(all_cohorts,G40sim,lwd=0.8,lty=1,col = col.sim)
lines(all_cohorts,apply(G40sim,1,median),col=2,lwd=2)
lines(all_cohorts,apply(G40sim,1,quantile,prob=1- (1-lev.p)/2),col=2,lwd=2,lty=2)
lines(all_cohorts,apply(G40sim,1,quantile,prob=(1-lev.p)/2),col=2,lwd=2,lty=2)

## keep variables to save results
var.keep <- c("cou","c1","c2","c3","n1","n2","n3","xo",
              "e40.act","g40.act","y","n.boot","n.simul",
              "e40_c1","e40_c2","e40_c3","e40_cstad","E40sim","e40_2D",
              "g40_c1","g40_c2","g40_c3","g40_cstad","G40sim","g40_2D",
              "ETAS","t_boot","Msim","Ssim",
              "s_cstad","S","bU_cstad","BUsim","bL_cstad","BLsim","CLsim","DLsim",
              "cL_cstad","CL","dL_cstad","DL","e0act","g0act",
              "Mstand","M_c1","M_c2","M_c3","sex","e0hat","g0hat",
              "LMX.act","LMXcstad","LMX.smo2D","DevRes",
              "cE","cMx","cZ","DXcstad","MX.act_low","c1b","c2b",
              "ADDobs","ADDcstad","ADD")

rm(list=setdiff(ls(),var.keep))
name <- paste0(cou,"_",sex,"_",n.boot*n.simul,"simul_results.Rdata")

# ## save results
# setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Results")
# save.image(name)
