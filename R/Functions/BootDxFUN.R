genPoissonResBoot <- function(Zobs,Zhat,rescale.dths=F,centre.res=F){
  require(rootSolve)
  ## dimensions
  m <- nrow(Zobs)
  n <- ncol(Zobs)
  mn <- m*n
  ## compute Poisson Deviance Residuals 
  res.t1 <- sign(Zobs-Zhat)
  res.t2 <- Zobs * log(ifelse(Zobs==0,1e-8,Zobs)/ifelse(Zhat==0,1e-8,Zhat))
  res.t3 <- Zobs - Zhat
  PoiDevRes <- res.t1 * sqrt(2) * (res.t2 - res.t3)^0.5
  ## resample residuals
  poidevres <- c(PoiDevRes)
  poidevres <- poidevres[!is.na(poidevres)]  ## remove NA
  if (centre.res){
    poidevres <- poidevres - mean(poidevres)
  }
  PoiDevResResample <- matrix(sample(x=poidevres,replace = T,size = mn),m,n)
  ## find DxHat via residual bootstap (Renshaw & Haberman 2008)
  ## univariate Newton-Raphson
  ZBoot <- matrix(NA,m,n)
  for (j in 1:n){
    for (i in 1:m){
      if (is.na(Zobs[i,j])){
        ## assign NA if observed death is NA
        ZBoot[i,j] <- NA
      }else{
        d_hat <- Zhat[i,j]
        a_hat <- log(d_hat)
        r_star <- PoiDevResResample[i,j]
        r_sign <- sign(r_star)
        c_hat_star <- (r_star^2)/2 - d_hat
        g <- function(d){
          d*log(d) - d*(1+a_hat) - c_hat_star
        }
        dg <- function(d){
          log(d) - a  
        }
        if (r_star < 0 & c_hat_star > 0 | 
            r_star==-sqrt(2*d_hat) & c_hat_star==0){
          ZBoot[i,j] <- 0
        }else{
          start <- max(1e-6, d_hat + r_sign * 0.5 * d_hat)
          ZBoot[i,j] <- multiroot(f = g, start = start, 
                          jacfunc = dg, positive = TRUE)$root
        }
      }
    }
  }
  ## need to rescale to number of observed deaths?
  if (rescale.dths) ZBoot <- ZBoot*sum(Zobs,na.rm = T)/sum(ZBoot,na.rm = T)
  return(ZBoot)
}