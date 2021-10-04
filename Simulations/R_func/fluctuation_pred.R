# fluctuating prediction Qn_bar
# input:
# Y   : The outcome vector
# A   : A The treatment
# W   : The covariates
# Qn  : A list of outcome regression estimates evaluated on observed data
# gn  : A list of propensity regression estimates evaluated on observed data
# a_0 : A list of fixed treatment values
# alpha:given the coefficent of fluctuation parameter and calculate observed Qnstar
# output:
# Qnstar : fluctuated Qn
pred_fluctuate_Q <- function(Qn, gn, a_0, alpha){
  pred_fluctuate_Q_apply <- function(a, Q, g, alp){
    plogis(qlogis(Q) + 1/g*alp)
  }
  Qnstar <- mapply(a = a_0, Q = Qn, g = gn, alp=alpha, FUN = pred_fluctuate_Q_apply, SIMPLIFY = FALSE)
  return(Qnstar)
}

