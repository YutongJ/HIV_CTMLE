# fluctuating the initial Qn_bar
# input:
# Y   : The outcome vector
# A   : A The treatment
# W   : The covariates
# Qn  : A list of outcome regression estimates evaluated on observed data
# gn  : A list of propensity regression estimates evaluated on observed data
# a_0 : A list of fixed treatment values
# output:
# Qnstar : fluctuated Qn

# Y = dat_AA[,1]; A=AA_k; Qn=Qn10; gn=gn10; a_0=nams
# Q=Qn10[[1]]; g=gn10[[1]];a=nams[1]
fluctuate_Q <- function(Y, A, W, Qn, gn, a_0, Qnstar_only = F) {
  fluctuate_Q_apply <- function(a, Q, g) {
    offset0 <- qlogis(Q)

    Ha <- as.numeric(A == a) / g
    # fit the logistic regression to get epsilons with starting value = 0
    fit <- stats::glm(
      Y ~ -1 + offset(off) + Ha,
      start = 0,
      data = data.frame(Y = Y, off = offset0, Ha = Ha), family = "binomial")
    if (!fit$converged){
      # if not converge, re-fit without starting values
      fit <- stats::glm(
        Y ~ -1 + offset(off) + Ha,
        data = data.frame(Y = Y, off = offset0, Ha = Ha), family = "binomial")
    }
    # fluctuating the initial estimator
    Qnstar <- stats::predict(
      fit,
      type = "response",
      newdata = data.frame(off = offset0, Ha = 1 / g)
    ) 
    list(est = Qnstar, eps = fit$coef)
  }
  
  QnStar <- mapply(a = a_0, Q = Qn, g = gn, FUN = fluctuate_Q_apply, SIMPLIFY = FALSE) # it will loop four times for each trt level
  
  
  if (Qnstar_only==F){
    return(QnStar)
  }
  if (Qnstar_only==T){
    return(do.call("list",lapply(QnStar, function(x){return(x$est)})))
  }
}
