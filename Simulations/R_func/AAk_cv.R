# cross validation
# fold=4;datt=dat;mnds=10;ntrees=1000;levels=4;wi=37;selections=c(5,10,20,100,200);tolg=0.01#;top=tops
# output:
# sortrisk    : sorted mean empirical risks from CTMLE
# cv_cov_TMLE : cross-validation TMLE variance matrix
#               list of 4(folds), estimated covariance matrix for four folds, within each fold, several matrix corresponding to different tops
# cv_cov_CTMLE : cross-validation CTMLE variance matrix
#               list of 4(folds), estimated covariance matrix for four folds, within each fold, several matrix corresponding to different tops

cv_triplets <- function(fold=4,datt=dat,wi,mnds=10,ntrees=1000,levels,selections,tolg=0.01){
  ###### funciton1: given Qnstar, get empirical risk
  neg_loglik <- function(Qn_obs, Y){
    # negative log-likelihood ---empirical loss
    loglik1 <- -1/length(Y)*sum(Y*log(Qn_obs)+(1-Y)*log(1-Qn_obs))
    return(loglik1)
  }
  
  ###### funciton2: ic_cov for tops
  ic_cov_func <- function(Qn_v, gn_v,dat_v){
    ic_matrix <- mapply(a = seq(levels), gn_a = gn_v, Qn_a = Qn_v,
                        function(a, Qn_a, gn_a){as.numeric(dat_v[,ind] == a)/gn_a * ((dat_v[,1]) - Qn_a) + Qn_a - mean(Qn_a)}, 
                        SIMPLIFY = TRUE)
    # covariance matrix like the one that comes out of drtmle
    ic_cov <- cov(ic_matrix) / samples.size
    return(ic_cov)
  }
  
  
  # dat list
  datlist_train <- list()
  datlist_validate <- list()
  
  num = samples.size/fold
  for(i in seq(fold)){
    datlist_train[[i]] <- datt[-((i*num-num+1):(i*num)),]
    datlist_validate[[i]] <- datt[(i*num-num+1):(i*num),]
  }
  
  # identify the column for AA_k
  ind <- which(colnames(datt)==paste0("W",wi))
  
  # for each fold, get empirical risk for validation group
  # train = datlist_train[[1]]; validate = datlist_validate[[1]]
  empirical_risk_apply <- function(train,validate){
    
    # fit RF for OR 
    fit1.RF <- randomForest::randomForest(factor(Y)~.,data = train, maxnodes=mnds, ntree=ntrees) # 
    # get the importance score for each covariate
    imp1 = round(randomForest::importance(fit1.RF), 2)
    # get the names of ranked covariates (based on importance)
    imp_name1 = names(sort(imp1[,1],decreasing = T))
    
    # train
    newdat=list()
    for(i in seq(levels)){
      temp = train
      temp[,ind]=rep(as.character(i),nrow(train))#levels(temp[,ind])=as.factor(seq(levels))
      newdat[[i]] = temp
    }
    rm(temp,i)
    Qn10 = list()
    # generate Qn list
    for (j in seq(levels)) {
      Qn10[[j]] = predict(fit1.RF, newdata = newdat[[j]], type = "prob")[,2]
    }
    rm(newdat)
    
    #validate
    newdat=list()
    for(i in seq(levels)){
      temp = validate
      temp[,ind]=rep(as.character(i),nrow(validate))#levels(temp[,ind])=as.factor(seq(levels))
      newdat[[i]] = temp
    }
    rm(temp,i)
    Qn10_validate = list()
    # generate Qn list
    for (j in seq(levels)) {
      Qn10_validate[[j]] = predict(fit1.RF, newdata = newdat[[j]], type = "prob")[,2]
    }
    rm(newdat)
    
    # Qn_obs for validation dataset
    # # predict directly from train model
    # Qn_validate_obs = predict(fit1.RF, newdata = validate, type = "prob")[,2]
    # extract from validation prediction Qn10 table
    Qn_validate_obs = Reduce(cbind,Qn10_validate)[cbind(1:num,as.numeric(validate[,ind]))]
    # Qn_validate_obs == Qn_validate_obs1
    
    AA_k = train[,ind]
    
    datt_k = train[,-c(1,ind)]
    datt_kv = validate[,-c(1,ind)]
    
    # the triplets contains three list in it: (1)gn; (2)Qn; (3)QnStar.
    triplets = list()
    # set the initial neg log-lik to Infinity
    loglik <- Inf 
    Qnstar_validate_obs <- list()
    ic_cov_validate_TMLE <- list()
    ic_cov_validate_CTMLE <- list()
    QnStar_v_TMLEs <- list()
    QnStar_v_CTMLEs <- list()
    
    for(i in seq_along(selections)){
      # select top-ranked AA's for PS regression in training group
      datt_k1 = datt_k[,colnames(datt_k)%in%imp_name1[1:selections[i]]]
      # select top-ranked AA's for PS regression in validation group
      datt_k_v = datt_kv[,colnames(datt_kv)%in%imp_name1[1:selections[i]]]
      
      # fit a RF for propensity score (PS)
      fit_PS <- randomForest::randomForest(factor(AA_k) ~ ., data=datt_k1, maxnodes=mnds, ntree=ntrees)
      
      # get P(A = xxx | W) by calling predict
      gn_W <- predict(fit_PS, newdata = datt_k1, type = "prob")
      # truncation of gn at tolg (all gn estimates less than tolg(e.g. 0.01) will be replaced by tolg)
      gn_W[gn_W < tolg] <- tolg
      # generate gn list
      gn10 <- split(gn_W, col(gn_W))
      
      # get P(A = xxx | W) by calling predict for validation group (here we keep the matrix format for matching in latter step)
      gn_validate <- predict(fit_PS, newdata = datt_k_v, type = "prob")
      gn_validate[gn_validate < tolg] <- tolg
      
      gn10_validate <- split(gn_validate, col(gn_validate))
      
      # set-up the initial triplet (e.x. with selections=5 covariates in it)
      if(i==1){
        # fluctuate the initial Qn based on gn
        QnStar <- fluctuate_Q(Y = train[,1], A = AA_k, W = datt_k, Qn = Qn10, gn = gn10, a_0 = seq(levels))
        # TMLE use original Qn10 and gn10 to update Qnstar
        QnStar_v_TMLE <- pred_fluctuate_Q(Qn = Qn10_validate, gn = gn10_validate, a_0 = seq(levels),
                                          alpha = do.call("c", lapply(QnStar, function(x){return(x$eps)})))
      }else{
        # fluctuate the initial Qn based on gn
        QnStar <- fluctuate_Q(Y = train[,1], A = AA_k, W = datt_k, Qn = triplets[[i-1]][[2]], gn = gn10, a_0 = seq(levels))
        QnStar_v_TMLE <- pred_fluctuate_Q(Qn = Qn10_validate, gn = gn10_validate, a_0 = seq(levels),
                                          alpha = do.call("c", lapply(QnStar, function(x){return(x$eps)})))
      }
      
      # obtain the observed Qn to calculate empirical loss
      Qn_obs <- do.call("cbind",lapply(QnStar,FUN=function(x){return(x[[1]])}))
      Qn_obs <- Qn_obs[cbind(1:nrow(train),as.numeric(AA_k))]
      
      # negative log-likelihood ---empirical loss
      loglik1 <- -1/nrow(train)*sum(train[,1]*log(Qn_obs)+(1-train[,1])*log(1-Qn_obs))
      # if the neg-loglik is smaller than that of Qnstar_bar_(k-1), accept the candidate triplet
      if(loglik1<loglik){
        loglik = loglik1
        triplets[[i]] <- list(gn10, Qn10, do.call("list",lapply(QnStar,function(x){x[[1]]})),eps = do.call("c",lapply(QnStar,function(x){x[[2]]})))
      }
      
      # if the neg-loglik is larger than that of Qnstar_bar_(k-1), then set Qn_bar = Qnstar_bar_(k-1) and go back to steps above
      if(loglik1>loglik){
        tmp <- triplets[[i-1]] # reset the Qn_bar to the updated one and re-estimate Qnstar
        # re-estimate the Qnstar with the updated Qn
        QnStar1 <- fluctuate_Q(Y = train[,1], A = AA_k, W = datt_k, Qn = tmp[[3]], gn = gn10, a_0 = seq(levels))
        # obtain the observed Qn to calculate empirical loss
        Qn_obs <- do.call("cbind",lapply(QnStar1,FUN=function(x){return(x[[1]])}))
        Qn_obs <- Qn_obs[cbind(1:nrow(train),as.numeric(AA_k))]
        # negative log-likelihood ---empirical loss
        loglik1 <- -1/nrow(train)*sum(train[,1]*log(Qn_obs)+(1-train[,1])*log(1-Qn_obs))
        loglik = loglik1
        
        triplets[[i]] <- list(gn10, Qn10, do.call("list",lapply(QnStar1,function(x){x[[1]]})),eps = do.call("c",lapply(QnStar1,function(x){x[[2]]})))
      }
      # cat(loglik1)
      # rename the triplet elements
      names(triplets[[i]]) = c("gn","Qn","Qnstar","eps")
      
      # get the observed Qnstar for validation group
      # get the coefficient alpha_i for each 1/g term
      eps = triplets[[i]]$eps
      
      # calculate the Qn_star for validate data using different values: logit^{-1}(logit(Qn_validate_obs) + alpha_i*1/g(A=i|w))
      # gn_validate[cbind(seq(num),as.numeric(validate[,ind]))]: based on observed A, we choose the corresponding gn for the calculation
      # eps[as.numeric(validate[,ind])]):                        for the alpha_i, we choose the one corresponding to the treatment(validate[,ind])
      # from CTMLE (eps is estimated through CTMLE procedure)
      Qnstar_validate_obs[[i]] <- plogis(qlogis(Qn_validate_obs) + 1/gn_validate[cbind(seq(num),as.numeric(validate[,ind]))]*eps[as.numeric(validate[,ind])])
      
      #Qnstar_validation (CTMLE) update from Qnstar_validate
      Qnstar_v_CTMLE <- pred_fluctuate_Q(Qn = Qn10_validate, gn = gn10_validate, a_0 = seq(levels),alpha = eps)
      
      # calculate covariance matrix for each top
      ic_cov_validate_TMLE[[i]]  <- ic_cov_func(gn_v = gn10_validate, Qn_v = QnStar_v_TMLE, dat_v = validate)
      QnStar_v_TMLEs[[i]] <- QnStar_v_TMLE
      
      # calculate covariance matrix for each top
      ic_cov_validate_CTMLE[[i]] <- ic_cov_func(gn_v = gn10_validate, Qn_v = Qnstar_v_CTMLE, dat_v = validate)
      QnStar_v_CTMLEs[[i]] <- Qnstar_v_CTMLE
      
    }# end loop
    names(ic_cov_validate_TMLE) = selections
    names(ic_cov_validate_CTMLE) = selections
    names(QnStar_v_TMLEs) = selections
    names(QnStar_v_CTMLEs) = selections
    
    # output the empirical risk for each elements (with differnet k), they decreases as expected. (for one pair of train and validation group)
    logliks <- mapply(FUN = neg_loglik, Qn_obs = Qnstar_validate_obs,MoreArgs = list(Y=validate$Y))
    
    
    return(list(loglik = logliks,
                ic_cov_validate_TMLE  = ic_cov_validate_TMLE,
                ic_cov_validate_CTMLE = ic_cov_validate_CTMLE,
                QnStar_v_TMLEs = QnStar_v_TMLEs,
                QnStar_v_CTMLEs = QnStar_v_CTMLEs))
  }
  
  # we got empirical risks for all pairs of train and validation group
  empirical_results <- mapply(FUN = empirical_risk_apply, train = datlist_train, validate = datlist_validate,SIMPLIFY = F)
  empirical_risks <- do.call("cbind",lapply(empirical_results,function(x){x[[1]]})) # in each "empirical_results" list, the first element is empirical risk
  # calcuate the mean empirical risk for each k: we want the least one to be our choice
  meanRisk = rowMeans(empirical_risks)
  names(meanRisk) = selections
  
  # transform back to updated Qnstar for tmle
  tmp <- do.call("list",lapply(empirical_results,function(x){x[[4]]}))
  tmp1 <- lapply(tmp,function(x){lapply(x, function(t){unlist(lapply(t,mean))})})
  # point estimate for each top-number
  Q_TMLE <- list()
  for (i in seq(selections)) {
    Q_TMLE[[i]] <- colMeans(do.call("rbind",lapply(tmp1, function(x){x[[i]]})))
  }
  names(Q_TMLE) <- selections
  rm(tmp,tmp1)
  
  # transform back to updated Qnstar for both ctmle
  tmp <- do.call("list",lapply(empirical_results,function(x){x[[5]]}))
  tmp1 <- lapply(tmp,function(x){lapply(x, function(t){unlist(lapply(t,mean))})})
  # point estimate for each top-number
  Q_CTMLE <- list()
  for (i in seq(selections)) {
    Q_CTMLE[[i]] <- colMeans(do.call("rbind",lapply(tmp1, function(x){x[[i]]})))
  }
  names(Q_CTMLE) <- selections
  rm(tmp,tmp1)
  
  
  # the sorted list of empirical risk of all k's: the one with the least empirical risk will be chosen as the targeted k
  return(list(sortrisk = sort(meanRisk),
              cv_cov_TMLE  = do.call("list",lapply(empirical_results,function(x){x[[2]]})),
              cv_cov_CTMLE = do.call("list",lapply(empirical_results,function(x){x[[3]]})),
              cv_est_TMLE = Q_TMLE,
              cv_est_CTMLE = Q_CTMLE))
}
