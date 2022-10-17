# options(echo=TRUE) # if you want to see commands in output file
# args=(commandArgs(TRUE))
# print(args)
# arguments = matrix(unlist(strsplit(args,"=")),ncol=2,byrow = T)
# 
# for (args_i in 1:length(args)) {
#   if(nchar(arguments[args_i,2])>5){
#     assign(arguments[args_i,1],arguments[args_i,2])
#   }
#   else{
#     assign(arguments[args_i,1],as.numeric(arguments[args_i,2]))
#   }
# }

rhos=0.75
sig=5
splt = 1
s=500
g=0.01

#############################
#####  parameter sets  ######
#############################

# AA levels for one particular position
# lel <- 2
lel <- 4
# lel <- 6


# probabilites for each group
p <- c(0.25,0.25,0.25,0.25)
# p <- c(0.85, 0.05, 0.05, 0.05)
# p <- c(0.315, 0.315, 0.32, 0.05)


# group_num <- 10
# simulations
sim <- (splt*group_num-(group_num-1)):(splt*group_num)

# # number of Cores used in simulation
# numCores = 16

# sample size (rows)
samples.size <- s

# selected AA
wis = c(10,85,87)

# indep. r.v. size (sub-cols)
indeps.size <- 0

# r.v. size (cols)
# rvs.size <- 30 + indeps.size
rvs.size <- 200



#############################
#####  directory paths  #####
#############################

# # local directory
# #results for TMLE different tops
# dir  = paste0("H:/Research/RF_TMLE/Results/20210222_New_Sims_3AA/Sample",samples.size,"/AA",wis,"/")
# # dir2 = paste0("H:/Research/RF_TMLE/Results/20210222_New_Sims_3AA/dat/")

# # path to read in data    e.g. w_ind and betas
fn  = "/Users/yutong/Library/CloudStorage/OneDrive-EmoryUniversity/Research/RF_CTMLE/Codes/20220822_Revision2/"
fn_func =  "/Users/yutong/Library/CloudStorage/OneDrive-EmoryUniversity/Research/RF_CTMLE/Codes/20220822_Revision2/"


# remote directory
# save data
dir  = paste0(path_save,"/Sample",samples.size,"/AA",wis,"/")
# dir2 = paste0(path_save,"/dat/")
# dir  = paste0("/home/yjin85/TMLE_RF/Results/20210222_New_Sims_3AA/Sample",samples.size,"/AA",wis,"/")
# dir2 = paste0("/home/yjin85/TMLE_RF/Results/20210222_New_Sims_3AA/dat/")

# read data
fn  = path_read
fn_func =  paste0(path_read,"/R_func/")

start.time = proc.time()
start.time



##############################
#####  required package  #####
##############################
library(MASS)
library(randomForest)
# library(drtmle)
library(parallel)
library(dummies)



#######################
#####  load data  #####
#######################
# make sure the betas are reasonable
load(paste0(fn,"beta_", sig,"_",lel,"lel.rda"))
load(paste0(fn,"w_ind_",sig,"_",lel,"lel.rda"))


#######################
#####  load func  #####
#######################

# simulation set-up function
source(paste0(fn_func,"Simulation_Setup.R"))



#####################
###  simulation 
#####################
# sample size: 500
a = proc.time()
# a vector for storing RF p-values
rf_imp_vec_500 <- matrix(NA, nrow=rvs.size, ncol = 1000)
for (sim in 1:1000){
  set.seed(sim)
  dat_W <- AA_generate(sample.size=500,rv.size=rvs.size,cov=var_mtx,prob = p)  # dat_W (uppercase W)
  dat = suppressWarnings(whole_dat(dat_W,level = lel, samples.size=500))
  # save(dat,file = paste0(dir2,"sim",sim,"_dat.rda"))
  
  fit_rf <- ranger::ranger(Y~.,data = dat, max.depth=10, num.trees=1000,
                           importance = "impurity_corrected")
  imp2 <- suppressWarnings(ranger::importance_pvalues(fit_rf))
  rf_imp_vec_500[,sim] <- imp2[,2]
}
proc.time()-a 

rownames(rf_imp_vec_500) <- paste0("W",1:rvs.size)

# sample size: 1000
a = proc.time()
# a vector for storing RF p-values
rf_imp_vec_1000 <- matrix(NA, nrow=rvs.size, ncol = 1000)
for (sim in 1:1000){
  set.seed(sim)
  dat_W <- AA_generate(sample.size=1000,rv.size=rvs.size,cov=var_mtx,prob = p)  # dat_W (uppercase W)
  dat = suppressWarnings(whole_dat(dat_W,level = lel, samples.size=1000))
  # save(dat,file = paste0(dir2,"sim",sim,"_dat.rda"))
  
  fit_rf <- ranger::ranger(Y~.,data = dat, max.depth=10, num.trees=1000,
                           importance = "impurity_corrected")
  imp2 <- suppressWarnings(ranger::importance_pvalues(fit_rf))
  rf_imp_vec_1000[,sim] <- imp2[,2]
}
proc.time()-a 

rownames(rf_imp_vec_1000) <- paste0("W",1:rvs.size)

# prepare simulation results table for target residues: 10, 85, 87
# rf_imp_vec_wis <- rbind(rf_imp_vec_500[wis,])
pct500 <- apply(rf_imp_vec_500, 1, function(x)(sum(x<0.05)/1000))
pct1000 <- apply(rf_imp_vec_1000, 1, function(x)(sum(x<0.05)/1000))


# load CTMLE results
dir21 = paste0("/Users/yutong/Library/CloudStorage/OneDrive-EmoryUniversity/Research/RF_CTMLE/Results/20210222_New_Sims_3AA/Sample500/")
dir22 = paste0("/Users/yutong/Library/CloudStorage/OneDrive-EmoryUniversity/Research/RF_CTMLE/Results/20210222_New_Sims_3AA/Sample1000/")

fn1 <- paste0(dir21,"sam",500,"_lel4_maxp0.25_AA",wis,"_pval_cv.rda")
fn2 <- paste0(dir22,"sam",1000,"_lel4_maxp0.25_AA",wis,"_pval_cv.rda")
# load simulation results
for (i in seq_along(fn1)) {
  load(fn1[i])
  load(fn2[i])
}

list_cv_500 = paste0("sam",500,"_lel4_maxp0.25_AA",wis,"_pval_cv")
list_cv_1000 = paste0("sam",1000,"_lel4_maxp0.25_AA",wis,"_pval_cv")

power_cal = function(p_table){
  pval = pval1 = p_table
  # Power
  pval[which(pval1 <0.05,arr.ind = T)]=1
  pval[which(pval1>=0.05,arr.ind = T)]=0
  Y = colMeans(pval)
  return(Y)
}

# calculate Power
Y_cv_500 = NULL
Y_cv_1000 = NULL
for (i in seq_along(list_cv_500)) {
  Y_cv_500 <- rbind(Y_cv_500, power_cal(p_table = get0(list_cv_500[i])))
  Y_cv_1000 <- rbind(Y_cv_1000, power_cal(p_table = get0(list_cv_1000[i])))
}


pcts <- rbind(c(n=500, Methods="raw RF", pct500[wis]),
              c(n=500, Methods="cv-CTMLE", Y_cv_500[,6]),
              c(n=1000, Methods="raw RF", pct1000[wis]),
              c(n=1000, Methods="cv-CTMLE", Y_cv_1000[,6]))



knitr::kable(pcts, format = "latex")



count <- apply(rf_imp_vec, 2, function(x)(sum(x<0.05)))
