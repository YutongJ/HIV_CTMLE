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
antibody <- c("vrc01", "vrc26", "1074", "PGT121", "PGT145")
s_max <- c(328, 229, 303, 313, 294)


# directory to save result
# local path
dir  = paste0("H:/Research/RF_TMLE/Results/20220419_HIV_Visuals/",antibody[anti],"/AA/")
dir4 = paste0("H:/Research/RF_TMLE/Results/20220419_HIV_Visuals/",antibody[anti],"/")

# cluster path
dir  = paste0(path_save, antibody[anti],"/AA/")
dir4 = paste0(path_save, antibody[anti],"/")

# AA positions
pos <- 3:20

# cluster
pos <- 3:(2+s_max[anti])


# input:
# x   : the simulation results (e.g. "sim1_sam500_lel4_maxp0.25_AA87")
# ind : the indicator in x$estimates[[ii]] list: (ii was passed in from outsides, which indicates which results to use, e.g. TMLE5, TMLE10, ..., CTMLE)
# ind = 1, estimates (4 level for AA's)
# ind = 2, covariance matrix (4*4)
# ind = 3, p values for pairwise test

comlst <- function(x,ind){
  temp <- get0(x)
  res <- as.vector(temp[[1]][[ii]][[ind]])
  return(res)
}

pval_full <- NULL
pval_cv <- NULL
for (j in pos){
  load(paste0(dir,"AA_",j,".rda"))
  
  tmp <- paste0("AA_",j)
  pval_full <- rbind(pval_full,unlist(get0(tmp)[[2]]))
  pval_cv <- rbind(pval_cv,unlist(get0(tmp)[[3]]))
}

# remove the specific AA10 simulation results
rm(list = ls(pattern = "AA_"))

save(pval_full, file = paste0(dir4,"pval_full.rda"))
save(pval_cv, file = paste0(dir4,"pval_cv.rda"))



