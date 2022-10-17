# 

# ------------------------------------------
#   add method comparison with LASSO
# ------------------------------------------
library(glmnet, quietly = T)
library(ranger)
library(fastDummies)

# data sets
antibody <- c("vrc01", "vrc26", "1074", "PGT121", "PGT145")
s_max <- c(328, 229, 303, 313, 294)

dir.dat <- "C:/Users/yjin85/OneDrive - Emory University/Research/RF_CTMLE/Codes/20220419_HIV_Visuals/data/"
# final plots and tables
dir_save = paste0("C:/Users/yjin85/OneDrive - Emory University/Research/RF_CTMLE/Results/20220419_HIV_Visuals/")

dir.dat <- "/Users/yutong/Library/CloudStorage/OneDrive-EmoryUniversity/Research/RF_CTMLE/Codes/20220419_HIV_Visuals/data/"
dir_save = paste0("/Users/yutong/Library/CloudStorage/OneDrive-EmoryUniversity/Research/RF_CTMLE/Results/20220419_HIV_Visuals/")

i=1
for (i in seq_along(antibody)){
  dat <- get(load(paste0(dir.dat,"dat_",antibody[i],".rda")))
  
  # LASSO
  dat_site <- dat[,-1]
  dat_site <- do.call(cbind, apply(dat_site, 2, function(x){fastDummies::dummy_cols(x,remove_selected_columns = TRUE)}))
  
  # one-hot encoding dataset
  dat0 <- cbind(Y=dat[,1], dat_site)
  dat0 <- as.matrix(dat0)
  
  split <- floor(nrow(dat0)/2)
  set.seed(12345)
  sample <- sample(1:nrow(dat0),split)
  #training set (one-hot encoding)
  dat1 <- dat0[sample,]
  # testing dat (original)
  dat2 <- dat[-sample,]
  # testing dat (one-hot encoding)
  dat3 <- dat0[-sample,]
  
  # CV LASSO
  set.seed(134553)
  fit1<- cv.glmnet(x=dat1[,-1],y=dat1[,1],
                   alpha=1, # 0 for ridge, 1 for LASSO
                   family='binomial')
  plot(fit1)
  
  # the selected lambda (lse)
  fit1$lambda.1se
  # 0.05790046
  
  # the coef for LASSO features
  LASSO.est1 <- coef(fit1, s = fit1$lambda.1se)
  
  # how many features based on 0.05790046 lambda
  (coef_count <- sum(LASSO.est1[-1]!=0))  # 19
  
  # The detailed sites
  LASSO_coef <- cbind(LASSO.est1@Dimnames[[1]][LASSO.est1@i+1], LASSO.est1@x)
  
  # extract selected hxb2 sites
  tmp <- LASSO_coef[grep(LASSO_coef[,1], pattern = "hxb2."),1]
  # unique selected sites
  site_all <- do.call(c, lapply(strsplit(tmp, "[.]"), function(x){x[2]}))
  site_unique <- unique(site_all)
  # 
  # # features identified from cv-LASSO
  # var_names <- c(LASSO_coef[grep(LASSO_coef[,1], pattern = "length|num"),1],
  #                paste0("hxb2.", site_unique))
  # 
  # # fit glm for selected features
  # fit1.glm <- glm(Y~., 
  #                 data = dat2[,colnames(dat2) %in% c("Y", var_names)],
  #                 family = "binomial")
  # 
  # sum1 <- summary(fit1.glm)

  fit1.glm2 <- glm(Y~., 
                   data = as.data.frame(dat3[,colnames(dat3) %in% c("Y", LASSO_coef[-1,1])]),
                   family = "binomial")
  sum2 <- summary(fit1.glm2)
  # Call:
  #   glm(formula = Y ~ ., family = "binomial", data = as.data.frame(dat3[, 
  #                                                                       colnames(dat3) %in% c("Y", LASSO_coef[-1, 1])]))
  # 
  # Deviance Residuals: 
  #   Min       1Q   Median       3Q      Max  
  # -2.5669  -0.4777  -0.2908  -0.2167   2.7435  
  # 
  # Coefficients: (2 not defined because of singularities)
  # Estimate Std. Error z value Pr(>|z|)    
  # (Intercept)          -1.83644    0.73644  -2.494  0.01264 *  
  # hxb2.169..data_Q      1.00193    0.60564   1.654  0.09806 .  
  # hxb2.174..data_A      0.59783    0.46458   1.287  0.19815    
  # hxb2.396..data_T      0.08052    0.67335   0.120  0.90482    
  # hxb2.399..data_N      1.02462    0.56622   1.810  0.07036 .  
  # hxb2.409..data_D     -0.25327    0.70157  -0.361  0.71810    
  # hxb2.456..data_Other  2.99291    0.66031   4.533 5.83e-06 ***
  # hxb2.456..data_R           NA         NA      NA       NA    
  # hxb2.459..data_G     -1.90342    0.63756  -2.985  0.00283 ** 
  # hxb2.459..data_Other       NA         NA      NA       NA    
  # hxb2.461..data_Other  0.47350    0.44743   1.058  0.28994    
  # hxb2.463..data_Other  0.72262    0.42641   1.695  0.09014 .  
  # hxb2.471..data_Other  0.84529    0.48065   1.759  0.07864 .  
  # hxb2.513..data_gap    0.34183    0.68649   0.498  0.61852    
  # hxb2.588..data_Q      0.63808    0.47136   1.354  0.17584    
  # hxb2.602..data_Other  1.33707    0.60225   2.220  0.02641 *  
  # hxb2.853..data_A      1.02882    0.37582   2.738  0.00619 ** 
  #   ---
  # 
  # (Dispersion parameter for binomial family taken to be 1)
  # 
  # Null deviance: 369.47  on 412  degrees of freedom
  # Residual deviance: 244.68  on 398  degrees of freedom
  # AIC: 274.68
  # 
  # Number of Fisher Scoring iterations: 6
  
  
  
  # ------------------------------------------
  #   add method comparison with RF
  # ------------------------------------------
  # # training set (original)
  # dat1b <- dat[sample,]
  # full data (priginal; all used for training in RF)
  dat1b <- dat
  # fit2 <- randomForest::randomForest(Y~.,
  #                                    data = dat1b, maxnodes=10, ntree=1000)
  # 
  # imp <- fit2$importance
  # 
  # # get feature importance for all covariates
  # imp = sort(round(randomForest::importance(fit2), 4)[,1],decreasing = T)
  # # get the names of ranked covariates (based on importance)
  # imp = cbind(names(imp), imp)
  # 
  # # features identified from RF
  # var_names <- imp[1:10,1]
  # 
  # # fit glm for selected features
  # fit2.glm <- glm(Y~., 
  #                 data = dat2[,colnames(dat) %in% c("Y", var_names)],
  #                 family = "binomial")
  # 
  # sum2 <- summary(fit2.glm)
  
  
  fit2 <- ranger::ranger(Y~.,data = dat1b, max.depth=10, num.trees=1000,
                         importance = "impurity_corrected")
  imp2 <- ranger::importance_pvalues(fit2)
  
  # site_imp <- imp2[which(imp2[,2]<0.05/328),]
  
}


#------------------------------------------
# method comparison plot: CTMLE, LASSO, RF
#------------------------------------------
plotcomp_dat <- vector("list",3)
# i = 1 # vrc01

# boundary line
hline= 0.05

# LASSO results
tmp_dat <- cbind(site = as.numeric(site_unique),
                 logp = -log10(sum2$coefficients[-1,4]),
                 sig = as.numeric(sum2$coefficients[-1,4] < hline/length(site_unique)))
plotcomp_dat[[1]] <- as.matrix(tmp_dat)


# RF results
# AA positions
f_split <- strsplit(rownames(imp2),split = "[.]")
# sites being tested (328 sites / total=856)
AA_sites <- as.numeric(do.call("c",lapply(f_split, function(x){return(x[2])})))

tmp_dat <- cbind(site = AA_sites,
                 pval = -log10(imp2[,2]),
                 sig = as.numeric(imp2[,2] < hline/s_max[i]))
tmp_dat <- tmp_dat[!is.na(tmp_dat[,1]),]
plotcomp_dat[[2]] <- tmp_dat


# CTMLE results
load(file = paste0(dir_save, "plots_data.rda"))
tmp_dat <- plots_data[[i]]
rm(plots_data)

tmp_dat <- cbind(site = tmp_dat[,1],
                 pval = tmp_dat[,2],
                 sig = as.numeric(tmp_dat[,2] > -log10(hline/s_max[i])))
plotcomp_dat[[3]] <- tmp_dat


manhattan_ind <- function(plotcomp_dat, 
                          total_num=856, 
                          hline=0.05,
                          log_p_max = 15,
                          ...){
  
  plot(y=-log10(rep(1,total_num)), x=seq(total_num),
       xlim = c(0,total_num+100),ylim=c(0,log_p_max+1),
       xlab = "Env Residue (HXB2-referenced)", 
       ylab = "",
       # main = "AA designations", 
       col="white", pch=19, cex.lab=1.5,
       bty="n", xaxt="n", yaxt="n")
  title(ylab=expression(paste("-log"[10],"(p-value)",sep='')), line=1.5, cex.lab=1.5)

  # mtext(side=1, at=450, line=5.5, adj=0.5, cex=1.3,
  #       paste0("(",LETTERS[num],") ",Antibody[num]))
  # add x-axis
  # tck: length of the vertical bar above each value
  # mgp: distance between vertical bar and value
  axis(side = 1, at = seq(0,856,100), las=1, cex.axis=1.5, tck=-0.01, mgp=c(1,0.5,0), gap.axis=0.25)
  # add y-axis
  axis(side = 2, at = seq(0,log_p_max, 1), las=1, cex.axis=1.5, tck=-0.01, mgp=c(1,0.5,0), gap.axis=0.25)
  
  
  # color shade for each sub-region of gp120
  rect(131,-0.5,157,log_p_max, col=rgb(red=1, green=0.5, blue=0, alpha=0.2), border=NA) # v1
  rect(158,-0.5,196,log_p_max, col=rgb(red=1, green=0.5, blue=0, alpha=0.2), border=NA) # v2
  rect(275,-0.5,283,log_p_max, col=rgb(red=0.1, green=0, blue=1, alpha=0.2), border=NA) # loop d
  rect(296,-0.5,331,log_p_max, col=rgb(red=1, green=0.5, blue=0, alpha=0.2), border=NA) # v3
  rect(353,-0.5,357,log_p_max, col=rgb(red=0.1, green=0, blue=1, alpha=0.2), border=NA) # loop e
  rect(385,-0.5,418,log_p_max, col=rgb(red=0.1, green=0, blue=1, alpha=0.2), border=NA) # v4
  rect(460,-0.5,469,log_p_max, col=rgb(red=0.1, green=0, blue=1, alpha=0.2), border=NA) # v5
  
  
  # black text labels for subtypes
  text(144,log_p_max+0.5, pos=3, labels = "V1", cex=1, col = "chocolate1")
  text(177,log_p_max+0.5, pos=3, labels = "V2", cex=1, col = "chocolate1")
  text(279,log_p_max, pos=3, labels = "Loop D", cex=1, col = "slateblue1")
  text(313,log_p_max+0.5, pos=3, labels = "V3", cex=1, col = "chocolate1")
  text(355,log_p_max, pos=3, labels = "Loop E", cex=1, col = "slateblue1")
  text(401,log_p_max, pos=3, labels = "V4", cex=1, col = "slateblue1")
  text(464,log_p_max, pos=3, labels = "V5", cex=1, col = "slateblue1")
  
  
  
  # LASSO
  dat_pch <- plotcomp_dat[[1]]
  points(y = dat_pch[,2], x = dat_pch[,1], pch =15, col = "gray80", cex=0.7)
  dat_pch <- dat_pch[which(dat_pch[,3] == 1),]
  points(y = dat_pch[,2], x = dat_pch[,1], pch =15, col = "darkcyan", cex=0.7)
  
  # RF
  dat_pch <- plotcomp_dat[[2]]
  dat_pch[which(dat_pch==Inf, arr.ind = T)] = log_p_max
  points(y = dat_pch[,2], x = dat_pch[,1], pch =21, col = "gray80", cex=0.7)
  dat_pch <- dat_pch[which(dat_pch[,3] == 1),]
  points(y = dat_pch[,2], x = dat_pch[,1], pch =21, col = "forestgreen", cex=0.7)
  
  # CTMLE
  dat_pch <- plotcomp_dat[[3]]
  dat_pch[which(dat_pch==Inf, arr.ind = T)] = log_p_max
  points(y = dat_pch[,2], x = dat_pch[,1], pch =17, col = "gray80", cex=0.7)
  dat_pch <- dat_pch[which(dat_pch[,3] == 1),]
  points(y = dat_pch[,2], x = dat_pch[,1], pch =17, col = "brown", cex=0.7)
  
  
  # boundary line
  h1 <- -log10(hline/s_max[i])
  h2 <- -log10(hline/length(site_unique))
  segments(0, h1,856,h1, lty=2, lwd=2)
  segments(0, h2,856,h2, lty=3, lwd=2)
  text(850,h1, pos=3, labels = "CTMLE/RF-corrected", cex=0.7)
  text(850,h2, pos=1, labels = "LASSO-corrected", cex=0.7)
  
  legend(850, par('usr')[4]*3/4, bty="o",cex=0.7,y.intersp = 1,
         legend = c("LASSO", "RF","CTMLE"), 
         col=c("darkcyan", "forestgreen", "brown"), pch=c(15,21,17))
  
}

pdf(file = paste0(dir_save, "AA_designations_comparison_vrc01_",format(Sys.Date(), "%y%m%d"),".pdf"),width = 9,height = 6)
manhattan_ind(plotcomp_dat = plotcomp_dat)
dev.off()








#-------------------------------------
# examine categories on each dataset
#-------------------------------------
# vrc01
i=1
dat <- get(paste0("dat_",antibody[i]))
test <- apply(dat[,-c(1:18)], 2, function(x){length(table(x))})
table(test)
# 2  3  4  5  6  7  8  9 10 11 12 
# 68 82 55 37 31 23 17  7  4  3  1 
(s_max[i]-31-23-17-15)/s_max[i] 
# 0.9542683


# vrc26
i=2
dat <- get(paste0("dat_",antibody[i]))
test <- apply(dat[,-c(1:18)], 2, function(x){length(table(x))})
table(test)
# 2  3  4  5  6  7  8 
# 63 66 48 30 15  6  1 


# 1074
i=3
dat <- get(paste0("dat_",antibody[i]))
test <- apply(dat[,-c(1:18)], 2, function(x){length(table(x))})
table(test)
# 2  3  4  5  6  7  8  9 10 11 
# 73 85 51 43 24 16  6  3  1  1 

# PGT121
i=4
dat <- get(paste0("dat_",antibody[i]))
test <- apply(dat[,-c(1:18)], 2, function(x){length(table(x))})
table(test)
# 2  3  4  5  6  7  8  9 10 11 
# 66 89 55 42 23 22  9  5  1  1

# PGT145
i=5
dat <- get(paste0("dat_",antibody[i]))
test <- apply(dat[,-c(1:18)], 2, function(x){length(table(x))})
table(test)
# 2  3  4  5  6  7  8  9 
# 86 81 56 41 13 13  3  1 








