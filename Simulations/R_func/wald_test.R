#' Wald test
#' C  : linear combination of testing
#' fit: model fitting
#' sum: summary(fit)
#' df : degree of freedom (depends on C)
wald = function(est,cov,C){
  df = nrow(C)
  Ccov = C%*%cov%*%t(C)
  X.w = t(C%*%est)%*%ginv(Ccov)%*%(C%*%est) #test.stat
  res = 1-pchisq(X.w,df)
  # res = c(X.w,1-pchisq(X.w,df))
  # names(res) = c("Test.Stat","P-value")
  return(res)
}