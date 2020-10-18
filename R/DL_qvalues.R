DLqvalues<-function(DL_data,K,estimation="auto")
{
  require(robust)
  require(qvalue)
  loadings<-DL_data# [,1:as.numeric(K)]
  resscale <- apply(loadings, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim=estimation)$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_DL<-qval$qvalues
  padj <- p.adjust(reschi2test,method="bonferroni")
  return(data.frame(p.values=reschi2test, q.values=q.values_DL,padj=padj))
}
###The covRob function selects a robust covariance estimator that is likely to provide a good estimate in a reasonable amount of time. Presently this selection is based on the problem size. The Donoho-Stahel estimator is used if there are less than 1000 observations and less than 10 variables or less than 5000 observations and less than 5 variables. If there are less than 50000 observations and less than 20 variables then the MCD is used. For larger problems, the Orthogonalized Quadrant Correlation estimator is used.
# Simqvalue=DLqvalues(DLsim1,10)