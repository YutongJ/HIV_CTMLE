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

# splt = 17
# splt_max=17
# anti = 2
# group_num = 20

# --------------
#   Parameters
# --------------
# Antibodies
antibody <- c("vrc26", "1074", "PGT121", "PGT145")
s_max <- c(229, 303, 313, 294)

# # numbers for each group
# group_num <- 20





# --------------
#   Library
# --------------
library(MASS)
library(randomForest)
# library(drtmle)
library(parallel)
library(dummies)
library(plyr)



# rm(list=ls())

# dir <- "H:/Research/RF_TMLE/Codes/20210511_HIV_Visuals/"
# dir_fun <- "H:/Research/RF_TMLE/Codes/20210511_HIV_Visuals/R_func/"
# dir_save <- paste0("H:/Research/RF_TMLE/Results/20210511_HIV_Visuals/", antibody[anti],"/AA/")

dir <- path_read
dir_fun <- paste0(path_read, "R_func/")
dir_save <- paste0(path_save, antibody[anti],"/AA/")
# dir <- "/home/yjin85/TMLE_RF/Codes/20210511_HIV_Visuals/"
# dir_fun <- "/home/yjin85/TMLE_RF/Codes/20210511_HIV_Visuals/R_func/"
# dir_save <- paste0("/home/yjin85/TMLE_RF/Results/20210511_HIV_Visuals/",antibody[anti],"/AA/")

# # --------------
# # Data cleaning
# # --------------
# # load cleaned dataset
# # odat_vrc01 <- read.csv(file = "H:/Research/RF_TMLE/Codes/20200514_HIV_QC/multiab_catnap_VRC01_24Mar2020_Revised.csv",stringsAsFactors = F)
# odat_vrc26 <- read.csv(file = paste0(dir,"data/slapnap_VRC26.08_11May2021.csv"),stringsAsFactors = F)
# odat_1074 <- read.csv(file = paste0(dir,"data/slapnap_10-1074_21May2021.csv"),header = TRUE, stringsAsFactors = F)
# odat_PGT121 <- read.csv(file = paste0(dir,"data/slapnap_PGT121_21May2021.csv"),header = TRUE, stringsAsFactors = F)
# odat_PGT145 <- read.csv(file = paste0(dir,"data/slapnap_PGT145_21May2021.csv"),header = TRUE, stringsAsFactors = F)
# 
# 
# # prepare dataset originally from catnap
# # input
# # odat: dataset originally from catnap
# # thred1: merging sub-level for each position with less than "thred1" AA. Default:30.
# # thred2: keep AA if min-level >= thred2. Default:30.
# # output:
# # Revised dataset with two thresholds:
# #   Y: outcome
# #   categorical variable: origin and subtype
# #   numerical variable
# #   AA information (min # of AA at each position >= thred2)
# #
# # odat=odat_vrc26;thred1=30; thred2 = 30
# prep.odat <- function(odat, thred1=30, thred2 = 30){
#   # remove "sequon_actual" columns
#   odat <- odat[,-grep(pattern = "sequon_actual",colnames(odat))]
# 
# 
#   # origin
#   geo <- odat[,grep(pattern = "geographic.region.of.origin.is.",colnames(odat),fixed = T)]
#   origin <- as.matrix(geo) %*% c(1:5)
#   origin <- plyr::mapvalues(origin,c(1:5),c("Americas", "Asia", "Europe", "NAfrica", "SAfrica"))
# 
#   # subtype
#   type <- odat[,grep(pattern = "subtype.is.",colnames(odat),fixed = T)]
#   subtype <- as.matrix(odat[,12:19]) %*% c(1:8)
#   subtype <- plyr::mapvalues(subtype,c(1:8),c("01_AE", "02_AG", "07_BC", "A1", "B", "C", "D", "Other"))
# 
#   odat_colname <- colnames(odat)[20:(ncol(odat)-15)]
#   # extract all AA residue subsets
#   AA_info <- do.call(rbind,strsplit(odat_colname, split="[.]"))
# 
#   # unique AA residues
#   AA_residue <- unique(AA_info[,2])
#   AA_pos_num <- length(AA_residue)
# 
# 
#   # Whether if there are some positions containing only one
#   for (i in seq(length(AA_residue))) {
#     ind <- grep(pattern = paste0("hxb2.",AA_residue[i],"."),colnames(odat),fixed = T)
#     temp <- odat[,ind]
# 
#     colsum <- colSums(temp)
#     if(sum(colsum==0)==(length(colsum)-1) | length(ind)==1){
#       print(AA_residue[i])
#     }
#   }
#   # No such positions
# 
#   # change the format of AA
#   dat_aa <- matrix(NA, nrow = nrow(odat), ncol = AA_pos_num)
#   for (i in seq_along(AA_residue)){
#     temp_AA_pos <- AA_info[AA_info[,2]==AA_residue[i],]
#     temp_colname <- apply(temp_AA_pos,1, function(x){paste(x,collapse = ".")})
# 
#     temp_odat <- odat[,match(temp_colname, colnames(odat))]
#     temp_AA <- as.matrix(temp_odat) %*% seq_along(temp_colname)
#     temp_AA <- plyr::mapvalues(temp_AA,seq_along(temp_colname),temp_AA_pos[,3])
#     dat_aa[,i] <- temp_AA
#   }
#   # change colnames of AA residue datasets
#   colnames(dat_aa) <- paste0("hxb2.",AA_residue)
# 
# 
# 
#   # outcome: "VRC01.ic50.censored"
#   Y <- odat$sens
# 
# 
# 
#   # -------------
#   # thresholding
#   # -------------
#   # newdat_aa <- matrix(NA,nrow = dim(dat_aa)[1], ncol=dim(dat_aa)[2])
#   newdat_aa <- NULL
#   ncols <- NULL
#   for (i in seq(ncol(dat_aa))){
#     temp <- dat_aa[,i]
#     nams <- names(which(table(temp) <= thred1))
#     temp[which(temp %in% nams)] <- "Other"
#     if (min(table(temp)) >= thred2){
#       newdat_aa <- cbind(newdat_aa, temp)
#       ncols <- c(ncols, i)
#     }
#   }
#   colnames(newdat_aa) <- colnames(dat_aa)[ncols]
#   rm(i,temp,nams,ncols)
# 
# 
#   # final Revised dataset
#   newdat <- cbind(Y=Y, origin, subtype, odat[,((ncol(odat)-14):ncol(odat))], newdat_aa)
# 
#   return(newdat)
# }
# 
# #
# dat_vrc26 <- prep.odat(odat_vrc26)
# dat_vrc26 <- dat_vrc26[complete.cases(dat_vrc26),]
# save(dat_vrc26, file = paste0(dir,"data/dat_vrc26.rda"))
# dat_1074  <- prep.odat(odat_1074)
# dat_1074 <- dat_1074[complete.cases(dat_1074),]
# save(dat_1074, file = paste0(dir,"data/dat_1074.rda"))
# dat_PGT121  <- prep.odat(odat_PGT121)
# dat_PGT121 <- dat_PGT121[complete.cases(dat_PGT121),]
# save(dat_PGT121, file = paste0(dir,"data/dat_PGT121.rda"))
# dat_PGT145  <- prep.odat(odat_PGT145)
# dat_PGT145 <- dat_PGT145[complete.cases(dat_PGT145),]
# save(dat_PGT145, file = paste0(dir,"data/dat_PGT145.rda"))




# load cleaned-dataset
dat_AA <- get(load(paste0(dir,"data/dat_", antibody[anti], ".rda")))

# wald test
source(paste0(dir_fun,"waldtest.R"))

# fluctuating the initial Qn_bar
source(paste0(dir_fun,"fluctuate_Q.R"))

# Given alpha, getting observed Qn_bar
source(paste0(dir_fun,"pred_fluctuate_Q.R"))

# calculation based on full data
source(paste0(dir_fun,"AAk_step2.R"))

# calculation based on cross-validation
source(paste0(dir_fun,"cv_triplets.R"))

# OR_AAk_tmle_ctmle
source(paste0(dir_fun,"OR_AAk_tmle_ctmle.R"))




#--------------------
#    data analysis 
#--------------------

# simulation for four levels
one_AA = function(AA_ind){
  # analysis one particular AA
  assign(paste0("AA_",AA_ind), OR_AAk_tmle_ctmle(datt1 = dat_AA, ind1=AA_ind, s = c(5,10,50,100,s_max[anti]), n_fold=5, g=0.01, q=0.0001))
    
  # save the data
  save(list=paste0("AA_",AA_ind), file = paste0(dir_save,"AA_",AA_ind,".rda"))
}



if (splt==1){
  sim <- 19:20
}


if (splt>1 & splt!= splt_max){
  # simulations
  sim <- (splt*group_num-(group_num-1)):(splt*group_num)
}
if(splt==splt_max){
  sim <- (splt*group_num-(group_num-1)):ncol(dat_AA)
}

# lapply
a = proc.time()
infos = lapply(X=sim,FUN=one_AA)
proc.time()-a

# # mclapply
# a = proc.time()
# infos = mclapply(19:346,one_AA,mc.cores = 5)
# proc.time()-a












