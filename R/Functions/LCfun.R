




## functions log(hazard) LC given parameters
LCetaFUN <- function(alpha,beta,kappa){
  One <- matrix(1, nrow=length(kappa), ncol=1)
  p1 <- alpha %*% t(One)
  p2 <- beta %*% t(kappa)
  eta <- p1+p2
  return(eta)
}

## function to compute rmse 
rmse <- function(actual,predicted){
  return(sqrt(mean((actual-predicted)^2)))
}

## return LT e(x) at first age of life table
lifetable.ex <- function(x, mx, sex="M", ax=NULL){
  m <- length(x)
  n <- c(diff(x), NA)
  if(is.null(ax)){
    ax <- rep(0,m)
    if(x[1]!=0 | x[2]!=1){
      ax <- n/2
      ax[m] <- 1 / mx[m]
    }else{    
      if(sex=="F"){
        if(mx[1]>=0.107){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.053 + 2.800*mx[1]
        }
      }
      if(sex=="M"){
        if(mx[1]>=0.107){
          ax[1] <- 0.330
        }else{
          ax[1] <- 0.045 + 2.684*mx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- 1 / mx[m]
    }
  }
  qx  <- n*mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]/mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  return(ex[1])
}

## function to compute Gini starting from mx
GINI_func <- function(ages,mx,ax=NULL,sex="M"){
  ## adjust mx column (remove zeros and NaN)
  mx.adj <- mx 
  whi.na <- which(is.na(mx.adj))
  if (length(whi.na) > 0) mx.adj <- mx.adj[-whi.na]
  whi.zero <- which(mx.adj==0)
  if (length(whi.zero) > 0) mx.adj <- mx.adj[-whi.zero]
  mx <- mx.adj
  ## adjust ages
  m <- length(mx)
  ages <- ages[1:m]
  n <- c(diff(ages), NA)
  ## build ax
  if(is.null(ax)){
    ax <- rep(0,m)
    if(ages[1]!=0 | ages[2]!=1){
      ax <- n/2
      ax[m] <- 1 / mx[m]
    }else{    
      if(sex=="F"){
        if(mx[1]>=0.107){
          ax[1] <- 0.350
        }else{
          ax[1] <- 0.053 + 2.800*mx[1]
        }
      }
      if(sex=="M"){
        if(mx[1]>=0.107){
          ax[1] <- 0.330
        }else{
          ax[1] <- 0.045 + 2.684*mx[1]
        }
      }
      ax[-1] <- n[-1]/2
      ax[m] <- 1 / mx[m]
    }
  }
  ## build cx
  cx <- contr <- rep(NA,m)
  cx <- ax - 1/2
  ## compute lx and qx
  qx  <- n*mx / (1 + (n - ax) * mx)
  qx[m] <- 1
  px  <- 1-qx
  lx  <- cumprod(c(1,px))*100000
  dx  <- -diff(lx)
  Lx  <- n*lx[-1] + ax*dx
  lx <- lx[-(m+1)]
  Lx[m] <- lx[m]/mx[m]
  Lx[is.na(Lx)] <- 0 ## in case of NA values
  Lx[is.infinite(Lx)] <- 0 ## in case of Inf values
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx
  ## build ax.hat
  ax.hat <- (1-(2/3)*qx+cx*(2-qx-(6/5)*cx))/(2-qx)
  if (ages[1]==0){
    ax.hat[1] <- ax[1]*(1-qx[1]*(3+0.831*ax[1])/(2+qx[1]))
  }
  lx.f <- c(lx[-1],0)
  contr <- lx.f^2+ax.hat*(lx^2-lx.f^2)
  gini.v <- 1- sum(contr)/
    (ex[1]*(lx[1]^2))
  return(gini.v)
}
