## ------------------------------------------------------------- ##
##  C-STAD model: fitting and forecasting cohort mortality
##  Last update: 17/10/2019
##  Author: Ugofilippo Basellini
##  Comments:
##  Out-of-sample validation of LC period model
## ------------------------------------------------------------- ##

## clean the workspace
rm(list = ls())

## load useful packages
library(demography)
library(colorspace)
library(DescTools)

## -- DATA & FUNCTIONS ----------

## load C-STAD functions
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Functions")
source("LCfun.R")

## load data 
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Data")
cou <- "SWE"        ## SWE or DNK
sex <- "F"          ## only F 
name <- paste0(cou,"coh",sex,".Rdata")  
load(name)    
PLOT <- F

## age dimensions
all.ages <- as.numeric(rownames(cE))
all.cohorts <- as.numeric(colnames(cE))
age.start <- 40
ages <- x <- age.start:110     ## HMD ages
m <- length(x)

## set out-of-sample exercises
ih.opt <- seq(10,35,5)
n.opt <- length(ih.opt)
MAE.errLC <- MAPE.errLC <- RMSE.errLC <- matrix(NA,n.opt,2)
FittingPeriod <- matrix(NA,n.opt,2)
t <- 1

## loop over scenarios
for (t in 1:n.opt){
  cat("out-of-sample exercise",t,"\n")
  cat("forecast horizon",ih.opt[t],"\n")
  
  ## OUT-OF-SAMPLE SCENARIO
  H <- ih.opt[t]

  ## cohort dimensions
  year.start <- 1835 - H
  year.end <- max(all.cohorts)
  cohorts <- y <- year.start:year.end  
  n <- length(y)
  
  ## starting data
  E <- cE[all.ages%in%x,all.cohorts%in%y]
  MX <- cMx[all.ages%in%x,all.cohorts%in%y]
  Z <- cZ[all.ages%in%x,all.cohorts%in%y]
  
  ## ages, cohorts and periods matrices
  AGES <- matrix(x,nrow=m,ncol=n)  
  COHORTS <- matrix(rep(y,each=m),nrow=m)
  PERIODS <- COHORTS + AGES
  
  ## cohorts first Lexis parallelogram (c1)
  c_breve <- 1905     ## for all exercises in the out-of-sample 
  c1 <- y[1]:c_breve  
  n1 <-  length(c1)
  cBT <- tail(c1,H)
  
  ## compute e40 and g40
  e40.act <- apply(MX[,y%in%c1],2,lifetable.ex,x=x,sex=sex)
  g40.act <- apply(MX[,y%in%c1],2,GINI_func,ages=x,sex=sex)
  e40.BT <- tail(e40.act,H)
  g40.BT <- tail(g40.act,H)
  if (PLOT){
    par(mfrow=c(1,2))
    plot(c1,e40.act);points(cBT,e40.BT,pch=16)
    plot(c1,g40.act);points(cBT,g40.BT,pch=16)
    par(mfrow=c(1,1))
  }
  
  ## periods for LC PERIOD
  periods <- (cBT[1]+age.start):PERIODS[m,which(y==c_breve)]
  np <- length(periods)
  
  ## PERIOD DATA
  MXp <- Zp <- Ep <- matrix(NA,m,np)
  j <- 1
  for (j in 1:np){
    MXp[,j] <- rev(MX[PERIODS == periods[j]])
    Zp[,j] <- rev(Z[PERIODS == periods[j]])
    Ep[,j] <- rev(E[PERIODS == periods[j]])
  }
  colnames(MXp) <- colnames(Zp) <- colnames(Ep) <- periods
  if (PLOT){
    par(mfrow=c(1,3))
    matplot(x,MXp,t="l",log="y",col=rainbow_hcl(np))
    matplot(x,Zp,t="l",col=rainbow_hcl(np))
    matplot(x,Ep,t="l",col=rainbow_hcl(np))
    par(mfrow=c(1,1))
  }
  
  ## fitting and forecast horizon
  periods.fit <- periods[1]:(periods[np]-H)
  n.fit <- length(periods.fit)
  periods.fore <- (periods.fit[n.fit]+1):periods[np]
  n.fore <- length(periods.fore)
  Ep.fit <- Ep[,periods%in%periods.fit]
  Zp.fit <- Zp[,periods%in%periods.fit]
  MXp.fit <- MXp[,periods%in%periods.fit]
  Ep.fore <- Ep[,periods%in%periods.fore]
  Zp.fore <- Zp[,periods%in%periods.fore]
  MXp.fore <- MXp[,periods%in%periods.fore]
  
  ## fitting data for LC
  LCfitData <- demogdata(data=MXp.fit, pop=Ep.fit, ages=ages, years=periods.fit,
                         type="mortality", label=cou, name="female", lambda=0)
  
  ## fit standard LC model
  iter <- 0
  for (iter in 0:m){
    ages.lc <- x[1]:x[m-iter]
    m.lc <- length(ages.lc)
    ## tro to fit
    LCfit <- try(lca(LCfitData, adjust="dt",years = periods.fit,
                     ages = ages.lc,interpolate = T,series = "female",
                     max.age = ages.lc[length(ages.lc)]),silent=T)
    if (class(LCfit)!="try-error") break
  }
  
  ## LC pars
  AlphaLC <- LCfit$ax
  BetaLC <- LCfit$bx
  KappaLC <- LCfit$kt
  if (PLOT){
    par(mfrow=c(1,3))
    plot(ages.lc,AlphaLC,main = expression(alpha[x]));
    plot(ages.lc,BetaLC,main = expression(beta[x]));
    plot(periods.fit,KappaLC,main = expression(kappa[t]))
    par(mfrow=c(1,1))
  }
  
  ## fitted rates and deaths
  ETAhatLC <- LCfit$fitted$y
  DXhatLC <- exp(ETAhatLC)*Ep.fit[x%in%ages.lc,]
  
  ## central forecast
  Kappats <- ts(c(KappaLC), start = periods.fit[1])
  modKappa <- Arima(Kappats, order=c(0,1,0), include.drift=TRUE)
  predK <- forecast(modKappa,h=H)
  ETAforeLC <- LCetaFUN(AlphaLC,BetaLC,predK$mean)
  
  if (PLOT){
    matplot(ages.lc,ETAhatLC,t="l",lty=1,lwd = 0.8,col=rainbow_hcl(n.fit))
    matlines(ages.lc,ETAforeLC,t="l",lty=1,lwd = 0.8,col=rainbow(H))
  }
  
  ## summary measures
  e40obs <- apply(MXp.fit,2,lifetable.ex,x=x,sex=sex)
  g40obs <- apply(MXp.fit,2,GINI_func,ages=x,sex=sex)
  e40BT <- apply(MXp.fore,2,lifetable.ex,x=x,sex=sex)
  g40BT <- apply(MXp.fore,2,GINI_func,ages=x,sex=sex)
  e40foreLC <- apply(exp(ETAforeLC),2,lifetable.ex,x=ages.lc,sex=sex)
  g40foreLC <- apply(exp(ETAforeLC),2,GINI_func,ages=ages.lc,sex=sex)
  
  if (PLOT){
    par(mfrow=c(1,2))
    plot(periods.fit,e40obs,ylim=range(e40obs,e40BT),xlim=range(periods),pch=19)
    points(periods.fore,e40BT)
    lines(periods.fore,e40foreLC,lwd=2,col=2)
    plot(periods.fit,g40obs,ylim=range(g40obs,g40BT),xlim=range(periods),pch=19)
    points(periods.fore,g40BT)
    lines(periods.fore,g40foreLC,lwd=2,col=2)
    par(mfrow=c(1,1))
  }
  
  ## extract LC cohorts from fitted and forecast rates
  ETA_LC <- cbind(ETAhatLC,ETAforeLC)
  AGES_LC <- matrix(ages.lc,nrow=length(ages.lc),ncol=np)  
  PERIODS_LC <- matrix(rep(periods,each=length(ages.lc)),
                       nrow=length(ages.lc))
  COHORTS_LC <- PERIODS_LC - AGES_LC
  
  LMXc.LC <- matrix(NA,length(ages.lc),H)
  j <- 1
  for (j in 1:H){
    LMXc.LC[,j] <- ETA_LC[COHORTS_LC==cBT[j]]
  }
  if (PLOT){
    matplot(ages,log(MX[,y%in%cBT]),t="l",lty=1,lwd = 0.8,col=rainbow_hcl(H))
    matlines(ages.lc,LMXc.LC,t="l",lty=1,lwd = 0.8,col=rainbow(H))
  }
  
  e40LCbt <- apply(exp(LMXc.LC),2,lifetable.ex,x=ages.lc,sex=sex)
  g40LCbt <- apply(exp(LMXc.LC),2,GINI_func,ages=ages.lc,sex=sex)
  if (PLOT){
    par(mfrow=c(1,2))
    plot(c1,e40.act);points(cBT,e40.BT,pch=16);lines(cBT,e40LCbt,col=2,lwd=2)
    plot(c1,g40.act);points(cBT,g40.BT,pch=16);lines(cBT,g40LCbt,col=2,lwd=2)
    par(mfrow=c(1,1))
  }
  
  RMSE.errLC[t,1] <- rmse(actual=e40.BT, predicted=e40LCbt)
  RMSE.errLC[t,2] <- rmse(actual=g40.BT*100, predicted=g40LCbt*100)
  
  MAPE.errLC[t,1] <- MAPE(x=e40LCbt, ref=e40.BT, na.rm = T)
  MAPE.errLC[t,2] <- MAPE(x=g40LCbt*100, ref=g40.BT*100, na.rm = T)
  
  MAE.errLC[t,1] <- MAE(x=e40LCbt, ref=e40.BT, na.rm = T)
  MAE.errLC[t,2] <- MAE(x=g40LCbt*100, ref=g40.BT*100, na.rm = T)
  
  FittingPeriod[t,1] <- periods[1]
  FittingPeriod[t,2] <- periods.fit[n.fit]
}

colnames(RMSE.errLC) <- colnames(MAPE.errLC) <- 
  colnames(MAE.errLC) <- c("e LC","g LC")
rownames(RMSE.errLC) <- rownames(MAPE.errLC) <- rownames(MAE.errLC) <- ih.opt

var.keep <- c("RMSE.errLC","MAPE.errLC","MAE.errLC","ih.opt","cou","FittingPeriod")
rm(list=setdiff(ls(),var.keep))

## save results
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Out-of-sample/Results")
name <- paste0(cou,"_LC.Rdata")
save.image(name)

