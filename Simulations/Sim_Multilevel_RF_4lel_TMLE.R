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

# rhos=0.75
# sig=5
# splt = 1
# s=500
# g=0.01

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
# 
# # path to read in data    e.g. w_ind and betas
# fn  = "H:/Research/RF_TMLE/Codes/20210222_New_Sims_3AA/"
# fn_func =  "H:/Research/RF_TMLE/Codes/20210222_New_Sims_3AA/R_func/"


# remote directory
# save data
dir  = paste0(path_save,"/Sample",samples.size,"/AA",wis,"/")
# dir2 = paste0(path_save,"/dat/")
# dir  = paste0("/home/yjin85/TMLE_RF/Results/20210222_New_Sims_3AA/Sample",samples.size,"/AA",wis,"/")
# dir2 = paste0("/home/yjin85/TMLE_RF/Results/20210222_New_Sims_3AA/dat/")

# read data
# fn  = "/home/yjin85/TMLE_RF/Codes/20210222_New_Sims_3AA/"
# fn_func =  "/home/yjin85/TMLE_RF/Codes/20210222_New_Sims_3AA/R_func/"
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

# general wald test
source(paste0(fn_func,"wald_test.R"))

# fluctuation for Qbar
source(paste0(fn_func,"fluctuation_Q.R"))

# prediction of fluctuation with known alpha
source(paste0(fn_func,"fluctuation_pred.R"))

# C/TMLE using full dataset
source(paste0(fn_func,"AAk_full.R"))
# C/TMLE using cross validation
source(paste0(fn_func,"AAk_cv.R"))

# estimation of C/TMLE for full & cv
source(paste0(fn_func,"AAs_ctmle.R"))




#####################
###  simulation 
#####################


# simulation for four levels
one_sim = function(sim){
  #use correlated structure for AA
  set.seed(sim)
  dat_W <- AA_generate(sample.size=samples.size,rv.size=rvs.size,cov=var_mtx,prob = p)  # dat_W (uppercase W)
  dat = suppressWarnings(whole_dat(dat_W,level = lel))
  # save(dat,file = paste0(dir2,"sim",sim,"_dat.rda"))
  
  for (j in seq_along(wis)){
    # assign the variable 
    assign(paste0("sim",sim,"_sam",samples.size,"_lel",lel,"_maxp",round(max(p),2),"_AA",wis[j]), 
           OR_AAk_tmle_ctmle(datt1 = dat, s = c(5,10,50,100,200), w_i = wis[j], l = lel, n_fold=4))
    
    # save the data
    save(list=paste0("sim",sim,"_sam",samples.size,"_lel",lel,"_maxp",round(max(p),2),"_AA",wis[j]),
         file = paste0(dir[j],"sim",sim,"_sam",samples.size,"_lel",lel,"_maxp",round(max(p),2),"_AA",wis[j],".rda"))
  }
}




# lapply
a = proc.time()
infos = lapply(X=sim,FUN=one_sim)
proc.time()-a





