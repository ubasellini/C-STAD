PeriodMat_FromCohort <- function(ages,cohorts,t.end=NULL,CohortMat){
  ## derive all periods
  m <- length(ages)
  n <- length(cohorts)
  mn <- m*n
  t.start <- min(cohorts) + min(ages)
  if (is.null(t.end)) t.end <- max(cohorts) + min(ages)
  t <- t.start:t.end
  nt <- length(t)
  ## matrix with period indicators
  Period_Indic <- matrix(rep(cohorts,each=m) + rep(ages,n),
                         nrow=m,ncol = n)
  ## matrix to store results
  PeriodMat <- matrix(NA,nrow = m,ncol = nt)
  rownames(PeriodMat) <- rownames(Period_Indic) <- ages
  colnames(Period_Indic) <- cohorts
  colnames(PeriodMat) <- t
  for (i in 1:m){
    for (j in 1:n){
      which.age <- ages[i]
      which.per <- Period_Indic[i,j]
      PeriodMat[i,which(t==which.per)] <- CohortMat[i,j]
    }
  }
  return(PeriodMat)
}
