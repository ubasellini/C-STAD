## ------------------------------------------------------------- ##
##  C-STAD model: fitting and forecasting cohort mortality
##  Last update: 17/06/2019
##  Author: UB with inputs from SK and GC
##  Comments:
##  Out-of-sample validation of C-STAD vs 2D P-spline
## ------------------------------------------------------------- ##

## clean the workspace
rm(list = ls())

## load useful packages
library(MortalitySmooth)
library(forecast)
library(colorspace)
library(vars)
library(DescTools)

## load C-STAD functions
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Functions")
source("C-STAD_Functions.R")

## load data
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Data")
cou <- "DNK"        ## SWE or DNK
sex <- "F"          ## only F 
name <- paste0(cou,"coh",sex,".Rdata")  ## Females
load(name)    

## STARTING DATA 
## age dimensions
age.start <- 40
xo <- age.start:110     ## HMD ages
x <- age.start:120      ## expanded ages
mo <- length(xo)
m <- length(x)
delta <- 0.1
xs <- seq(min(x), max(x), delta)  ## ages at fine grid
ms <- length(xs)

## cohort dimensions
year.start <- 1835  ## first cohort with data observed at all ages in DNK
year.end <- 1905    ## for all exercises in the out-of-sample 
y <- year.start:year.end  ## use for both SWE and DNK (to have same time-range of analysis + reliable data)
n <- length(y)
coly <- rainbow_hcl(n)  

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

## augmented age-axis (for standard)
ages.add.l <- 90
ages.add.r <- 80
delta1 <- 1
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

## observed data
ages <- as.numeric(row.names(cE))
cohorts <- as.numeric(colnames(cE))
E.obs <- cE[ages%in%xo,cohorts%in%y]
MX.obs <- cMx[ages%in%xo,cohorts%in%y]
Z.obs <- cZ[ages%in%xo,cohorts%in%y]

## actual and fitted life expectancy & gini
e40.obs <- apply(MX.obs,2,lifetable.ex,x=xo,sex=sex)
g40.obs <- apply(MX.obs,2,GINI_func,ages=xo,sex=sex)

## set out-of-sample exercises
ih.opt <- seq(10,35,5)
n.opt <- length(ih.opt)
MAE.err <- MAPE.err <- RMSE.err <- matrix(NA,n.opt,4)
t <- 6
## loop over scenarios
for(t in 1:n.opt){
  cat("out-of-sample exercise",t,"\n")
  cat("forecast horizon",ih.opt[t],"\n")
  
  ## select out-of-sample duration
  ih <- ih.opt[t]
  
  ## re-load data at each exercise
  setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Data")
  load(name)    
  
  ## NA years for the out-of-sample
  years.na <- c((1905-ih+1):1905)
  
  # insert NA in the years we do not observe (out-of-sample)
  for(i in 1:length(years.na)){
    cE[((nrow(cE)+1)-i):nrow(cE),cohorts%in%years.na[i]] <- NA
    cZ[((nrow(cE)+1)-i):nrow(cE),cohorts%in%years.na[i]] <- NA
    cMx[((nrow(cE)+1)-i):nrow(cE),cohorts%in%years.na[i]] <- NA
  }
  
  ## cohorts first Lexis parallelogram (c1)
  c1 <- y[1]:(years.na[1]-1) 
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

  ## -- STANDARD ----------
  
  ## find the optimal smoothing parameters of 2D P-spline model ##
  smooth2D <- Mort2Dsmooth(x=x,y=y,Z=ZA,offset=log(EA),W=WA,
                           ndx=c(ndx,ndy),method = 1)
  
  ## subset with observed data
  LMX.smo2D <- matrix(Bs %*% c(smooth2D$coefficients),ms,n)
  LMX.smo2DW <- LMX.smo2D*Ws
  LMX.smo2DW[LMX.smo2DW==0] <- NA
  
  ## actual and 2Dsmo e40, g40
  e40_2D <- apply(exp(LMX.smo2D[xs%in%xo,]),2,lifetable.ex,x=xo,sex=sex)
  g40_2D <- apply(exp(LMX.smo2D[xs%in%xo,]),2,GINI_func,ages=xo,sex=sex)
  
  ## FX
  FX.smo2D <- apply(exp(LMX.smo2D),2,dx_from_mx,age=xs)
  FX.smo2DW <- FX.smo2D*Ws
  FX.smo2DW[FX.smo2DW==0] <- NA

  ## compute modal age at death (it's pretty smooth here due to 2D smoothing)
  M_2D <- xs[apply(FX.smo2D, 2, which.max)] + (delta/2)

  ## STANDARD DISTRIBUTION
  s_2D <- M_2D - M_2D[1]
  
  ## derive aligned distributions
  FX.align <- matrix(0, nrow=ms, ncol=n)
  for(i in 1:n){
    FX.align[,i] <- fx_shift(age=xs,fx=FX.smo2D[,i],shift=-s_2D[i],ndx = ndx,deg = deg)
  }
  FX.alignW <- FX.align*Ws
  FX.alignW[FX.alignW==0] <- NA

  ## Standard = mean of the aligned densities
  FXallmeanW <- exp(apply(log(FX.alignW), 1, mean, na.rm=T))  
  FXstand <- FXallmeanW

  ## coefficients & mode
  lambdaSTAND <- 1
  Standard <- coeff_stand(age=xs,fx=FXstand,ndx=ndxA,deg=deg,
                          ages.add.l=ages.add.l,ages.add.r=ages.add.r,
                          lambda=lambdaSTAND)
  coeff_Stand <- Standard$betasA
  Mstand <-  xAs[which.max(BAs%*%coeff_Stand)]
  
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
    # lmx.smo <- try(Mort1Dsmooth(x=x,y=Z1A[,i],offset = log(E1A[,i]),
    #                             w = W1A[,i],ndx = ndx,deg = deg),silent = T)
    # if (class(lmx.smo) == "Mort1Dsmooth"){
    #   LMX.smo1D[,i] <- Bxs %*% lmx.smo$coefficients
    # }else{
      lmx.smo <- lmx_smooth(age = x,y = Z1A[,i],e = E1A[,i],w =  W1A[,i],
                            ndx = ndx,deg = deg)
      LMX.smo1D[,i] <- Bxs %*% lmx.smo$coef
    # }
    FX.smo1D[,i] <- dx_from_mx(age=xs,mx=exp(LMX.smo1D[,i]))
  }
  
  ## compute modal age at death 
  M_c1 <- xs[apply(FX.smo1D, 2, which.max)] + (delta/2)
  s_c1 <- M_c1 - Mstand
  
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
    # cat("fitting cohort", c1[i], "\n")
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
           xlim=range(x),xlab="Age")
      lines(xs,LMX.smo1D[,i],col=2,lwd=2,lty=2)
      lines(xs,LMX.smo2D[,i],col=3,lwd=2,lty=2)
      lines(x,LMXcstad_c1[,i],col=4,lwd=2,lty=1)
      abline(v=xo[mo],lty=3);abline(v=M_c1[i],lty=3)
      legend("topleft",c("Observed","Smooth 1D","Smooth 2D","C-STAD"),pch = c(16,NA,NA,NA),
             lty = c(0,2,2,1),col = c(1,2,3,4),bty="n",cex = 1.3,lwd = 3)
      # locator(1)
    }
  }
  
  ## actual and fitted life expectancy & gini
  e40_c1 <- apply(exp(LMXcstad_c1[1:mo,]),2,lifetable.ex,x=xo,sex=sex)
  g40_c1 <- apply(exp(LMXcstad_c1[1:mo,]),2,GINI_func,ages=xo,sex=sex)
  
  ## -- C2 ----------
  
  ## select wide cohort range
  c.fore <- (c1[n1]+1):1905   
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
  
  ## compute raw mode of LT dx 
  ## (to be compared afterwards with Mode from smooth rates)
  Mraw <- xo[apply(DX.fore,2,which.max)]
  
  
  ## smooth 1D
  ## expand data for 1D smoothing
  E.foreA <- rbind(E.fore,matrix(100,m-mo,n.fore))
  Z.foreA <- rbind(Z.fore,matrix(100,m-mo,n.fore))
  W.foreA <- rbind(W.fore,matrix(0,m-mo,n.fore))
  
  ## find unsmooth M and s
  LMX.foreSMO1D <- FX.foreSMO1D <-  matrix(NA,ms,n.fore)
  for (i in 1:n.fore){
    # lmx.smo <- try(Mort1Dsmooth(x=x,y=Z.foreA[,i],offset = log(E.foreA[,i]),
    #                             w = W.foreA[,i],ndx = ndx,deg = deg),silent = T)
    # if (class(lmx.smo) == "Mort1Dsmooth"){
    #   LMX.foreSMO1D[,i] <- Bxs %*% lmx.smo$coefficients
    # }else{
      lmx.smo <- lmx_smooth(age = x,y = Z.foreA[,i],e = E.foreA[,i],w =  W.foreA[,i],
                            ndx = ndx,deg = deg)
      LMX.foreSMO1D[,i] <- Bxs %*% lmx.smo$coef
    # }
    FX.foreSMO1D[,i] <- dx_from_mx(age=xs,mx=exp(LMX.foreSMO1D[,i]))
  }
  
  ## compute modal age at death 
  M2.fore <- xs[apply(FX.foreSMO1D, 2, which.max)] + (delta/2)
  
  if(cou == "DNK" & ih == 5 ){
    c.cutoff <- 1905    
  }else if(cou == "DNK" & ih == 10){
    c.cutoff <- 1905
  }else if(cou == "DNK" & ih == 15){
    c.cutoff <- 1905
  }else if(cou == "DNK" & ih == 20){
    c.cutoff <- 1905  
  }else if(cou == "DNK" & ih == 25){
    c.cutoff <- 1904   
  }else if(cou == "DNK" & ih == 30){
    c.cutoff <- 1899
  }else if(cou == "DNK" & ih == 35){
    c.cutoff <- 1895
  }else if(cou == "SWE" & ih == 10){
    c.cutoff <- 1905
  }else if(cou == "SWE" & ih == 15){
    c.cutoff <- 1905
  }else if(cou == "SWE" & ih == 20){
    c.cutoff <- 1905  
  }else if(cou == "SWE" & ih == 25){
    c.cutoff <- 1903   
  }else if(cou == "SWE" & ih == 30){
    c.cutoff <- 1899 
  }else if(cou == "SWE" & ih == 35){
    c.cutoff <- 1894}
  
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
  
  ## Shifting parameter for constrained forecast
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
  bL_c2 <- bU_c2 <- cL_c2 <- dL_c2 <- numeric(n2)
  DXcstad_c2 <- LMXcstad_c2 <- matrix(NA, m, n2)
  
  ## Estimation
  i <- 1
  for(i in 1:n2){
    # cat("fitting cohort", c2[i], "\n")
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
           xlim=range(x),xlab="Age")
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
  
  ## forecast life expectancy & gini
  e40_c2 <- apply(exp(LMXcstad_c2[1:mo,]),2,lifetable.ex,x=xo,sex=sex)
  g40_c2 <- apply(exp(LMXcstad_c2[1:mo,]),2,GINI_func,ages=xo,sex=sex)
  
  
  ## -- C3 ----------
  if (c.cutoff<1905){
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
    
    
    ## empty vectors and matrices to store results 
    bL_c3 <- cL_c3 <- dL_c3 <- rep(NA,n3)
    DXcstad_c3 <- LMXcstad_c3 <- matrix(NA, m, n3)
    
    ## Estimation 
    i <- 1
    n.grid.fit <- 500
    for(i in 1:n3){
      # cat("fitting cohort", c3[i], "\n")
      ## declare explicitely inputs of cleversearch function
      xlo_2 <- XLO3[[i]]
      shat <- s_c3[i]
      bUhat <- bU_c3[i]
      xup_2 <- XUP3[[i]]
      coeff.stand <- coeff_Stand 
      if (is.vector(Z3)){
        Dx_2 <- Z3
        Ex_2 <- E3
        wei_2 <- W3
      }else{
        Dx_2 <- Z3[,i]  
        Ex_2 <- E3[,i]
        wei_2 <- W3[,i]
      }
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
      if (PLOT){
        plot(xo,LMX3[,i],pch=16,main=c3[i],ylim = range(LMX3,LMXcstad_c3,finite=T),
             xlim=range(x),xlab="Age")
        lines(xs,LMX.smo2D[,which(y==c3[i])],col=3,lwd=2,lty=2)
        lines(x,LMXcstad_c3[,i],col=4,lwd=2,lty=1)
        abline(v=M_c3[i])
        abline(v=xo[mo],lty=3)
        legend("topleft",c("Observed","Smo 2D","C-STAD"),pch = c(16,NA,NA),
               lty = c(0,2,1),col = c(1,3,4),bty="n",cex = 1.3,lwd = 3)
        # locator(1)
      }
    }
    
    ## forecast life expectancy & gini
    if (is.vector(Z3)){
      e40_c3 <- lifetable.ex(x = xo,mx = exp(LMXcstad_c3[1:mo]),sex=sex)
      g40_c3 <- GINI_func(ages=xo,mx = exp(LMXcstad_c3[1:mo]),sex=sex)
    }else{
      e40_c3 <- apply(exp(LMXcstad_c3[1:mo,]),2,lifetable.ex,x=xo,sex=sex)
      g40_c3 <- apply(exp(LMXcstad_c3[1:mo,]),2,GINI_func,ages=xo,sex=sex)
    }
    e40_cstad <- c(e40_c1,e40_c2,e40_c3)
    g40_cstad <- c(g40_c1,g40_c2,g40_c3)
  }else{
    e40_cstad <- c(e40_c1,e40_c2)
    g40_cstad <- c(g40_c1,g40_c2) 
  }
  
  
  #--------- calcualate error measures -----------# 
  rmse <- function(actual,predicted){return(sqrt(mean((actual-predicted)^2)))}
  
  e40_cstad_out <- e40_cstad[y %in% (y[n]-ih+1):y[n]]
  g40_cstad_out <- g40_cstad[y %in% (y[n]-ih+1):y[n]]
  e40_2D_out <- e40_2D[y %in% (y[n]-ih+1):y[n]]
  g40_2D_out <- g40_2D[y %in% (y[n]-ih+1):y[n]]
  e40_obs_out <- e40.obs[y %in% (y[n]-ih+1):y[n]]
  g40_obs_out <- g40.obs[y %in% (y[n]-ih+1):y[n]]
  
  RMSE.err[t,1] <- rmse(actual=e40_obs_out, predicted=e40_cstad_out)
  RMSE.err[t,2] <- rmse(actual=e40_obs_out, predicted=e40_2D_out)
  RMSE.err[t,3] <- rmse(actual=g40_obs_out*100, predicted=g40_cstad_out*100)
  RMSE.err[t,4] <- rmse(actual=g40_obs_out*100, predicted=g40_2D_out*100)
  
  MAPE.err[t,1] <- MAPE(x=e40_cstad_out, ref=e40_obs_out, na.rm = T)
  MAPE.err[t,2] <- MAPE(x=e40_2D_out, ref=e40_obs_out, na.rm = T)
  MAPE.err[t,3] <- MAPE(x=g40_cstad_out*100, ref=g40_obs_out*100, na.rm = T)
  MAPE.err[t,4] <- MAPE(x=g40_2D_out*100, ref=g40_obs_out*100, na.rm = T)

  MAE.err[t,1] <- MAE(x=e40_cstad_out, ref=e40_obs_out, na.rm = T)
  MAE.err[t,2] <- MAE(x=e40_2D_out, ref=e40_obs_out, na.rm = T)
  MAE.err[t,3] <- MAE(x=g40_cstad_out*100, ref=g40_obs_out*100, na.rm = T)
  MAE.err[t,4] <- MAE(x=g40_2D_out*100, ref=g40_obs_out*100, na.rm = T)
}

colnames(RMSE.err) <- c("e CSTAD","e 2Dsmo","g CSTAD","g 2Dsmo")
colnames(MAPE.err) <- c("e CSTAD","e 2Dsmo","g CSTAD","g 2Dsmo")
colnames(MAE.err) <- c("e CSTAD","e 2Dsmo","g CSTAD","g 2Dsmo")
rownames(RMSE.err) <- rownames(MAPE.err) <- rownames(MAE.err) <- ih.opt

round(RMSE.err,digits = 3)
# round(MAPE.err,digits = 6)
# round(MAE.err,digits = 4)

var.keep <- c("RMSE.err","MAPE.err","MAE.err","ih.opt","cou")
rm(list=setdiff(ls(),var.keep))

## save results
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Out-of-sample/Results")
name <- paste0(cou,".Rdata")
save.image(name)
