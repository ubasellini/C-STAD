##-- C-STAD MODEL #### C-STAD MODEL #### C-STAD MODEL --##
##-- C-STAD MODEL #### C-STAD MODEL #### C-STAD MODEL --##
##
##  R code to compare the forecast accuracy of the C-STAD
##  vs the 2D P-spline model and the LC model
##  
##  Authors: Ugofilippo Basellini
##           Soren Kjaegaard
##           Giancarlo Camarda
##  Last update: 14/10/2019
##
##------------------------------------------------------##

## cleaning the workspace
rm(list=ls(all=TRUE))

## load Results
setwd("~/Documents/Demography/Work/STADcohorts/99_Github/C-STAD/R/Out-of-sample/Results")
cou <- "SWE" 
name <- paste0(cou,".Rdata")
load(name)
load(paste0(cou,"_LC.Rdata"))  ## add LC results

## reshape table
RMSE_swe <- MAE_swe <- MAPE_swe <- matrix(NA,12,3) 
i <- 1
for (i in 1:6){
  RMSE_swe[2*i-1,1:2] <- RMSE.err[i,1:2]
  RMSE_swe[2*i-1,3] <- RMSE.errLC[i,1]    ## include LC
  RMSE_swe[2*i,1:2] <- RMSE.err[i,3:4]
  RMSE_swe[2*i,3] <- RMSE.errLC[i,2]      ## include LC
  MAE_swe[2*i-1,1:2] <- MAE.err[i,1:2]
  MAE_swe[2*i-1,3] <- MAE.errLC[i,1]      ## include LC
  MAE_swe[2*i,1:2] <- MAE.err[i,3:4]
  MAE_swe[2*i,3] <- MAE.errLC[i,2]        ## include LC
  MAPE_swe[2*i-1,1:2] <- MAPE.err[i,1:2]
  MAPE_swe[2*i-1,3] <- MAPE.errLC[i,1]    ## include LC
  MAPE_swe[2*i,1:2] <- MAPE.err[i,3:4]
  MAPE_swe[2*i,3] <- MAPE.errLC[i,2]      ## include LC
}

## load Denmark
cou <- "DNK"     
name <- paste0(cou,".Rdata")
load(name)
load(paste0(cou,"_LC.Rdata"))  ## add LC results

RMSE_dnk <- MAE_dnk <- MAPE_dnk <- matrix(NA,12,3) 
for (i in 1:6){
  RMSE_dnk[2*i-1,1:2] <- RMSE.err[i,1:2]
  RMSE_dnk[2*i-1,3] <- RMSE.errLC[i,1]    ## include LC
  RMSE_dnk[2*i,1:2] <- RMSE.err[i,3:4]
  RMSE_dnk[2*i,3] <- RMSE.errLC[i,2]      ## include LC
  MAE_dnk[2*i-1,1:2] <- MAE.err[i,1:2]
  MAE_dnk[2*i-1,3] <- MAE.errLC[i,1]      ## include LC
  MAE_dnk[2*i,1:2] <- MAE.err[i,3:4]
  MAE_dnk[2*i,3] <- MAE.errLC[i,2]        ## include LC
  MAPE_dnk[2*i-1,1:2] <- MAPE.err[i,1:2]
  MAPE_dnk[2*i-1,3] <- MAPE.errLC[i,1]    ## include LC
  MAPE_dnk[2*i,1:2] <- MAPE.err[i,3:4]
  MAPE_dnk[2*i,3] <- MAPE.errLC[i,2]      ## include LC
}

## join
RMSE_table <- cbind(RMSE_swe,RMSE_dnk,rep(ih.opt,each=2))
MAE_table <- cbind(MAE_swe,MAE_dnk,rep(ih.opt,each=2))
MAPE_table <- cbind(MAPE_swe,MAPE_dnk,rep(ih.opt,each=2))
rownames(RMSE_table) <-  rownames(MAE_table) <- 
  rownames(MAPE_table) <- rep(c("e40","g40"),6)
colnames(RMSE_table) <- colnames(MAE_table) <- 
  colnames(MAPE_table) <- c(rep(c("CSTAD","2Dsmo","LC"),2),"horiz")

## table 1 in paper
round(RMSE_table,digits = 2)
round(RMSE_table,digits = 3)

## table B1 in paper
round(MAE_table,digits = 2)
round(MAE_table,digits = 3)

## table B2 in paper
100*round(MAPE_table,digits = 4)
100*round(MAPE_table,digits = 5)

