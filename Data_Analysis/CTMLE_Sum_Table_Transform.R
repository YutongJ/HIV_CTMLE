options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)
arguments = matrix(unlist(strsplit(args,"=")),ncol=2,byrow = T)

for (args_i in 1:length(args)) {
  if(nchar(arguments[args_i,2])>5){
    assign(arguments[args_i,1],arguments[args_i,2])
  }
  else{
    assign(arguments[args_i,1],as.numeric(arguments[args_i,2]))
  }
}

ls()

# anti = 1



# Antibodies
antibody <- c("vrc26", "1074", "PGT121", "PGT145")
s_max <- c(229, 303, 313, 294)


# Summary Table prep
# dir = "H:/Research/RF_TMLE/Results/20201012_HIV_CTMLE/AA/"
dir = paste0(path_save, antibody[anti],"/AA/")
dir.save = paste0(path_save, antibody[anti],"/")


# i_pos <- 19:20
i_pos <- 19:(18+s_max[anti])

# sig_level
level <- 0.95

est_list = vector(mode = "list", length = s_max[anti])
for (i in i_pos){
  AA <- get(load(paste0(dir, "AA_", i, ".rda")))
  est <- unlist(AA$estimates$CTMLE$est_CTMLE_cv)
  cov_cv <- AA$estimates$CTMLE$cov_CTMLE_cv
  
  # # naive wald
  # lower <- est - qnorm(0.975) * sqrt(diag(cov_cv))
  # lower <- round(lower, 3)
  # upper <- est + qnorm(0.975) * sqrt(diag(cov_cv))
  # upper <- round(upper,3)
  
  # logit transformed
  logitMean <- list(f       = function(eff){ qlogis(eff) }, # qlogis: logit
                    f_inv   = function(eff){ plogis(eff) }, # plogis: expit
                    h       = function(est0){ est0 },
                    fh_grad = function(est0){ 1/(est0 - est0^2) })
  
  out <- matrix(NA, nrow = length(est), ncol = 3)
  
  thisC <- do.call(logitMean$h, args = list(est0 = est))
  
  f_thisC <- do.call(logitMean$f, args = list(eff = thisC))
  
  # gradient for each AA level 
  grad <- mapply(logitMean$fh_grad, est)
  
  # sd.err for each AA level
  thisSe <- sqrt(grad * diag(cov_cv) * grad)
    
  for (j in seq_along(est)) {
    transformCI <- rep(f_thisC[j], 3) +
      stats::qnorm(c(0.5, (1 - level) / 2, (1 + level) / 2)) *
      rep(thisSe[j], 3)
    out[j, ] <- do.call(logitMean$f_inv, args = list(eff = transformCI))
  }
  row.names(out) <- names(est)
  colnames(out) <- c("est", "cil", "ciu")
  
  # result format presented in the table: "est (cil, ciu)"
  val <- paste0(round(out[,1],5)," (",round(out[,2],5), ", ", round(out[,3],5),")")
  ord <- order(est,decreasing = T)
  
  est_list[[i-18]] <- c(rbind(names(est)[ord],val[ord]))
}


save(est_list, file = paste0(dir.save, "est_list_transform.rda"))









