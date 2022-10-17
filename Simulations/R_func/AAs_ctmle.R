# datt1=dat;w_i=87;s=c(5,10,20,30,40);l=lel;n_fold=4
OR_AAk_tmle_ctmle <- function(datt1,w_i,s = c(5,10,50,100,200), l = lel, n_fold=4){
  # identify the column for AA_k
  ind <- which(colnames(datt1)==paste0("W",w_i))
  
  # results from cross validation
  cv_top <- cv_triplets(fold=n_fold, datt=datt1, wi=w_i,levels=l, selections=s, mnds=10, tolg = g)
  
  # results from full data
  full_top <- AAk_step2(datt = datt1, wi = w_i, mnds=10, ntrees=1000, levels=l, selections=s, tolg = g)
  
  
  #############   TMLE
  # get TMLE point estimates for tops
  est_tops_TMLE <- full_top$est_tops_TMLE
  est_tops_TMLE_cv <- cv_top$cv_est_TMLE
  
  # get TMLE covariance matrix using full data
  cov_TMLE_full <- full_top$cov_tops_TMLE
  # get cross-validation TMLE covariance matrix for tops (mean of covariance matrices for different folds)
  cov_TMLE_cv <- mapply(FUN=function(x,i){Reduce("+", do.call("c", lapply(x, function(x){return(x[i])})))/n_fold}, i = seq(s), MoreArgs = list(x = cv_top$cv_cov_TMLE), SIMPLIFY = F)
  names(cov_TMLE_cv) <- s
  
  #############   CTMLE
  # get the top k from ctmle
  cv_top1 <- as.numeric(names(cv_top[[1]])[1])
  # selected triplets based on ctmle results
  res_cv <- full_top$triplets[[which(s==cv_top1)]]
  
  # QnStar = list length 4 with Qbar(a, W_i), i = 1,..,n
  # gn = list length 4 with g(a | W_i), i = 1,..,n
  est_CTMLE <- apply(do.call("cbind",res_cv[[3]]), MARGIN = 2, mean)
  est_CTMLE_cv <- cv_top$cv_est_CTMLE[[as.character(cv_top1)]]
  
  # get CTMLE covariance matrix using full data
  cov_CTMLE_full <- full_top$cov_tops_CTMLE[[which(s==cv_top1)]]
  # get cross-validation CTMLE covariance matrix for tops (mean of covariance matrices for different folds)
  cov_CTMLE_cvs <- mapply(FUN=function(x,i){Reduce("+", do.call("c", lapply(x, function(x){return(x[i])})))/n_fold}, i = seq(s), MoreArgs = list(x = cv_top$cv_cov_CTMLE), SIMPLIFY = F)
  cov_CTMLE_cv <- cov_CTMLE_cvs[[which(s==cv_top1)]]
  
  # contrast to test whether all estimates are equal to each other
  Con <- rep(0,times=l*(l-1))
  Con[seq(1,l*(l-1),l)+seq(l-1)-1] <- 1
  Con[seq(1,l*(l-1),l)+seq(l-1)] <- -1
  C <- matrix(Con,ncol = l,byrow = T)
  
  
  res_list <- list()
  # save the results for TMLE (different tops)
  for (ii in seq(s)) {
    res_list[[ii]] = list(est_TMLE = est_tops_TMLE[[ii]],
                          cov_TMLE_full = cov_TMLE_full[[ii]],
                          cov_TMLE_cv   = cov_TMLE_cv[[ii]],
                          p_full = wald(est = est_tops_TMLE[[ii]], cov = cov_TMLE_full[[ii]],C),
                          p_cv = wald(est = est_tops_TMLE_cv[[ii]], cov = cov_TMLE_cv[[ii]],C),
                          est_TMLE_cv = est_tops_TMLE_cv[[ii]])
  }
  
  # save the results for CTMLE
  res_list[[length(s)+1]] <- list(est_CTMLE = est_CTMLE,
                                  cov_CTMLE_full = cov_CTMLE_full,
                                  cov_CTMLE_cv   = cov_CTMLE_cv,
                                  p_full = wald(est = est_CTMLE, cov = cov_CTMLE_full,C),
                                  p_cv = wald(est = est_CTMLE_cv, cov = cov_CTMLE_cv,C),
                                  est_CTMLE_cv = est_CTMLE_cv)
  names(res_list) = c(paste0("TMLE_",s), "CTMLE")
  
  # hypotheses test for all pairs in res_list (TMLE & CTMLE)
  hypo_test_p_full = lapply(res_list, function(x){return(wald(est = x[[1]], cov = x[[2]],C = C))})
  hypo_test_p_cv = lapply(res_list, function(x){return(wald(est = x[[6]], cov = x[[3]],C = C))})
  
  return(list(estimates = res_list,    # results list of 3: (1)est;  (2)cov;  (3)p-value (pairwise comparisons)
              p_value_full   = hypo_test_p_full,  # linear combination of all pairwise comparison with covariance_full
              p_value_cv   = hypo_test_p_cv  # linear combination of all pairwise comparison with covariance_cv
  ))
}

# test1 <- OR_AAk_tmle_ctmle(datt1=dat, w_i = 10, s = c(5,10,50,100,200), l = lel, n_fold=4)

