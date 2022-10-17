# datt=dat;mnds=10;ntrees=1000;levels=4;wi=87;selections=c(5,10,20,50,70,100);tolg=0.01#;top=tops
AAk_step2 <- function(datt,wi,mnds=10,ntrees=1000,levels,selections,tolg=0.01){
  # ic_cov for tops
  ic_cov_func <- function(Qn_v, gn_v){
    ic_matrix <- mapply(a = seq(levels), gn_a = gn_v, Qn_a = Qn_v,
                        function(a, Qn_a, gn_a){as.numeric(datt[,ind] == a)/gn_a * ((datt[,1]) - Qn_a) + Qn_a - mean(Qn_a)}, 
                        SIMPLIFY = TRUE)
    # covariance matrix like the one that comes out of drtmle
    ic_cov <- cov(ic_matrix) / samples.size
    return(ic_cov)
  }
  
  # identify the column for AA_k
  ind <- which(colnames(datt)==paste0("W",wi))
  
  newdat=list()
  for(i in seq(levels)){
    temp = datt
    temp[,ind]=rep(as.character(i),nrow(datt))#levels(temp[,ind])=as.factor(seq(levels))
    newdat[[i]] = temp
  }
  rm(temp,i)
  
  # fit RF for OR
  # set.seed(sim)
  fit1.RF <- randomForest::randomForest(factor(Y)~.,data = datt, maxnodes=mnds, ntree=ntrees) #
  # get the importance score for each covariate
  imp1 = round(randomForest::importance(fit1.RF), 2)
  # get the names of ranked covariates (based on importance)
  imp_name1 = names(sort(imp1[,1],decreasing = T))
  
  Qn10 = list()
  # generate Qn list
  for (j in seq(levels)) {
    Qn10[[j]] = predict(fit1.RF, newdata = newdat[[j]], type = "prob")[,2]
  }
  
  
  AA_k = datt[,ind]
  datt_k = datt[,-c(1,ind)]
  
  
  # the triplets contains three list in it: (1)gn; (2)Qn; (3)QnStar.
  triplets = list()
  # set the initial neg log-lik to Infinity
  loglik <- Inf
  # get Qnstar_full for tops 
  Qnstar_full_tops <- list()
  cov_tops_TMLE <- list()
  cov_tops_CTMLE <- list()
  
  for(i in seq(selections)){
    # select top-ranked AA's for PS regression
    datt_k1 = datt_k[,colnames(datt_k)%in%imp_name1[1:selections[i]]]
    
    # fit a RF for propensity score (PS) outside of drtmle
    fit_PS <- randomForest::randomForest(factor(AA_k) ~ ., data=datt_k1, maxnodes=mnds, ntree=ntrees)
    
    # get P(A = xxx | W) by calling predict
    gn_W <- predict(fit_PS, newdata = datt_k1, type = "prob")
    # truncation of gn at tolg (all gn estimates less than tolg(e.g. 0.01) will be replaced by tolg)
    gn_W[gn_W < tolg] <- tolg
    # generate gn list
    gn10 <- split(gn_W, col(gn_W))
    
    # set-up the initial triplet (e.x. with selections=5 covariates in it)
    if(i==1){
      # fluctuate the initial Qn based on gn
      QnStar <- fluctuate_Q(Y = datt[,1], A = AA_k, W = datt_k, Qn = Qn10, gn = gn10, a_0 = seq(levels))
    }else{
      # fluctuate the initial Qn based on gn
      QnStar <- fluctuate_Q(Y = datt[,1], A = AA_k, W = datt_k, Qn = triplets[[1]][[2]], gn = gn10, a_0 = seq(levels))
    }
    
    # obtain the observed Qn to calculate empirical loss
    Qn_obs <- do.call("cbind",lapply(QnStar,FUN=function(x){return(x[[1]])}))
    Qn_obs <- Qn_obs[cbind(1:nrow(datt),as.numeric(AA_k))]
    
    # negative log-likelihood ---empirical loss
    loglik1 <- -1/nrow(datt)*sum(datt[,1]*log(Qn_obs)+(1-datt[,1])*log(1-Qn_obs))
    
    # if the neg-loglik is smaller than that of Qnstar_bar_(k-1), accept the candidate triplet
    if(loglik1<loglik){
      loglik = loglik1
      triplets[[i]] <- list(gn10, Qn10, do.call("list",lapply(QnStar,function(x){x[[1]]})))
    }
    
    # if the neg-loglik is larger than that of Qnstar_bar_(k-1), then set Qn_bar = Qnstar_bar_(k-1) and go back to steps above
    if(loglik1>loglik){
      tmp <- triplets[[i-1]] # reset the Qn_bar to the updated one and re-estimate Qnstar
      QnStar1 <- fluctuate_Q(Y = datt[,1], A = AA_k, W = datt_k, Qn = tmp[[3]], gn = gn10, a_0 = seq(levels))
      # obtain the observed Qn to calculate empirical loss
      Qn_obs <- do.call("cbind",lapply(QnStar1,FUN=function(x){return(x[[1]])}))
      Qn_obs <- Qn_obs[cbind(1:nrow(datt),as.numeric(AA_k))]
      # negative log-likelihood ---empirical loss
      loglik1 <- -1/nrow(datt)*sum(datt[,1]*log(Qn_obs)+(1-datt[,1])*log(1-Qn_obs))
      loglik = loglik1
      
      triplets[[i]] <- list(gn10, Qn10, do.call("list",lapply(QnStar1,function(x){x[[1]]})))
    }
    # print(loglik1)
    
    
    # save Qnstar_full list for each top in a list
    Qnstar_full_tops[[i]] <- fluctuate_Q(Y = datt[,1], A = AA_k, W = datt_k, Qn = Qn10, gn = gn10, a_0 = seq(levels), Qnstar_only = T)
    
    # calculate covariance matrix for each top
    cov_tops_TMLE[[i]]  <- ic_cov_func(gn_v = gn10, Qn_v = Qnstar_full_tops[[i]])
    
    # calculate covariance matrix for each top
    cov_tops_CTMLE[[i]] <- ic_cov_func(gn_v = triplets[[i]][[1]], Qn_v = triplets[[i]][[3]])
    
  }
  names(cov_tops_TMLE) = selections
  names(cov_tops_CTMLE) = selections
  
  # get TMLE estimates for tops
  est_tops_TMLE <- do.call("list", lapply(Qnstar_full_tops, function(x){return(unlist(lapply(x, mean)))}))
  # rename for better correspondence
  names(est_tops_TMLE) <- paste0("tops",selections)
  
  return(list(triplets = triplets,
              est_tops_TMLE = est_tops_TMLE,
              cov_tops_TMLE = cov_tops_TMLE,
              cov_tops_CTMLE = cov_tops_CTMLE))
}
