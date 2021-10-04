
#' Manhattan Plot
#'
#' Fluctuating original Qn to get Qnstar.
#'
#' @param Y   The outcome vector
#' @param A   The trail
#' @param W   The covariates
#' @param Qn  A list of outcome regression estimates evaluated on observed data
#' @param gn  A list of propensity regression estimates evaluated on observed data
# @param a_0 A list of fixed treatment values
#' @param p_t0 Pr(T=0)
#' @param trt_r A ratio of P(T=0|W)/P(T=1|W) --- used in H(t,w)
#' @param trimQ value to truncate scaled Qn at. Default is tolQ.
#' @param step three steps in total for fluctuating S.
#' 
#' @return Qnstar: a fluctuated Qn.
#'
#' @examples
#' dat <- data_generate(123) 
#' fluctuate_Q(Y = dat$pS1, Delta=dat$Delta, A = dat$Trt, W = dat[,5:6], Qn = Qn1, gn = gn1, a_0 = c(0,1),step=1)
#'
#' @importFrom SuperLearner trimLogit
#' @importFrom stats predict glm
#'
#' @export
#' 

cat("\f")
rm(list=ls())

# pval = pval_cv[,6];hline=0.05;total_num=856;dir=dir_save;antibody = antibody[anti]
site_summary <- function(pval,hline=0.05, dir=dir_save, total_num=856,antibody){
  #----------------
  # P-value table
  #----------------
  
  # if p-value < 0.05/328  (Bonferroni corrected) --> 1, otherwise --> 0
  pval_sig <- ifelse(pval<(hline/length(pos)),1,0)
  
  # AA positions
  f_split <- strsplit(names(pval),split = ".",fixed = T)
  # sites being tested (328 sites / total=856)
  AA_sites <- as.numeric(do.call("c",lapply(f_split, function(x){return(x[2])})))
  # sites not being tested (NoT included in 328 sites / total=856)
  AA_sites_none <- seq(total_num)[-AA_sites]
  
  
  # add-on est(CI) information
  est_table <- matrix(NA, dimnames = NULL,
                      nrow = length(pval), 
                      ncol = max(unlist(lapply(est_list, function(x){length(x)}))))
  for(i in seq(pval)){
    temp <- est_list[[i]]
    est_table[i,1:length(temp)] <- temp
  }
  
  # first three columns:  (1) p-values for each AA site
  #                       (2) -log10(p_value)
  #                       (3) significant: Yes=1; No=0
  pval_cv_summary <- cbind(p_value = pval,
                           "-log10_pval" = -log10(pval),
                           significant = pval_sig)
  pval_cv_summary <- cbind(pval_cv_summary, est_table)
  
  # ordered list based on p-values
  pval_cv_tab <- pval_cv_summary[order(as.numeric(pval_cv_summary[,1]), decreasing = F),]
  # output ordered table as csv file
  # pval_cv_tab <- cbind(pval_cv_tab,rank = 1:)
  write.csv(pval_cv_tab, file = paste0(dir,"Ordered_Summary_Sites_Details_",antibody,".csv"), na="")
  # write.csv(pval_cv_tab, file = paste0(dir,"Ordered_Summary_Sites_Details_Transform.csv"), na="")

  
  #----------------
  # Manhattan Plot
  #----------------
  log_p = -log10(pval)
  log_p[which(log_p==Inf, arr.ind = T)] = max(log_p[log_p!=Inf])+1
  log_p_max <- max(log_p)
  
  pdf(file = paste0(dir, "AA_designations_",antibody,".pdf"),width = 15)
  # par(oma=c(0, 0, 0, 0))
  
  par(mar=c(5,6,4,1)+.1)
  
  plot(y=-log10(rep(1,total_num)), x=seq(total_num),
       xlim = c(0,total_num+100),ylim=c(0,max(log_p)+1),
       xlab = "Env Residue (HXB2-referenced)", 
       ylab = expression(paste("-log"[10],"(p-value)",sep='')),
       # main = "AA designations", 
       col="white", pch=19, cex.lab=1.2,
       bty="n", xaxt="n", yaxt="n")
  
  # add x-axis
  # tck: length of the vertical bar above each value
  # mgp: distance between vertical bar and value
  axis(side = 1, at = seq(0,856,50), las=1, cex.axis=1, tck=-0.01, mgp=c(1,0.5,0), gap.axis=0.25)
  # add y-axis
  axis(side = 2, at = seq(0,log_p_max, 1), las=1, cex.axis=1, tck=-0.01, mgp=c(1,0.5,0), gap.axis=0.25)
  
  
  # color shade for each sub-region of gp120
  rect(131,-0.5,157,log_p_max, col=rgb(red=0.5, green=0, blue=1, alpha=0.2), border=NA)
  rect(158,-0.5,196,log_p_max, col=rgb(red=0.5, green=0, blue=1, alpha=0.2), border=NA)
  rect(275,-0.5,283,log_p_max, col=rgb(red=0.5, green=0, blue=1, alpha=0.2), border=NA)
  rect(296,-0.5,331,log_p_max, col=rgb(red=0.5, green=0, blue=1, alpha=0.2), border=NA)
  rect(353,-0.5,357,log_p_max, col=rgb(red=0.5, green=0, blue=1, alpha=0.2), border=NA)
  rect(385,-0.5,418,log_p_max, col=rgb(red=0.5, green=0, blue=1, alpha=0.2), border=NA)
  rect(460,-0.5,469,log_p_max, col=rgb(red=0.5, green=0, blue=1, alpha=0.2), border=NA)

  # rect(131,-0.5,157,17, col=rgb(red=1, green=0, blue=0, alpha=0.3), border=NA)
  # rect(158,-0.5,196,17, col=rgb(red=1, green=0.5, blue=0, alpha=0.3), border=NA)
  # rect(275,-0.5,283,17, col=rgb(red=1, green=1, blue=0, alpha=0.3), border=NA)
  # rect(296,-0.5,331,17, col=rgb(red=0.5, green=1, blue=0, alpha=0.3), border=NA)
  # rect(353,-0.5,357,17, col=rgb(red=0, green=1, blue=1, alpha=0.3), border=NA)
  # rect(385,-0.5,418,17, col=rgb(red=0, green=0.5, blue=1, alpha=0.3), border=NA)
  # rect(460,-0.5,469,17, col=rgb(red=0.5, green=0, blue=1, alpha=0.3), border=NA)
  
  # black text labels for subtypes
  text(144,log_p_max, pos=3, labels = "V1", cex=1, col = "purple")
  text(177,log_p_max, pos=3, labels = "V2", cex=1, col = "purple")
  text(279,log_p_max, pos=3, labels = "Loop D", cex=1, col = "purple")
  text(313,log_p_max, pos=3, labels = "V3", cex=1, col = "purple")
  text(355,log_p_max, pos=3, labels = "Loop E", cex=1, col = "purple")
  text(401,log_p_max, pos=3, labels = "V4", cex=1, col = "purple")
  text(464,log_p_max, pos=3, labels = "V5", cex=1, col = "purple")
  
  # colorful text labels for subtypes
  # text(131,17, labels = "V1", col=rgb(red=1, green=0, blue=0), cex=0.6)
  # text(158,17, labels = "V2", col=rgb(red=1, green=0.5, blue=0), cex=0.6)
  # text(275,17, labels = "Loop D", col=rgb(red=0, green=0, blue=0), cex=0.6)
  # text(296,17, labels = "V3", col=rgb(red=0.5, green=1, blue=0), cex=0.6)
  # text(353,17, labels = "Loop E", col=rgb(red=0, green=1, blue=1), cex=0.6)
  # text(385,17, labels = "V4", col=rgb(red=0, green=0.5, blue=1), cex=0.6)
  # text(460,17, labels = "V5", col=rgb(red=0.5, green=0, blue=1), cex=0.6)
  
  
  
  # boundary line
  h = -log10(hline/length(pos))
  segments(0, h,856,h, lty=2, lwd=2)
  
  # sites not being tested
  points(y = rep(0,length(AA_sites_none)), x = AA_sites_none, pch = 21, col = "gray80")
  
  # sites being tested
  points(y = log_p, x = AA_sites, pch=20, col = "black")
  #------------------------------------------
  # Red-ish colors = signal peptide
  ind <- which(AA_sites>=1 & AA_sites<=30)
  points(y = log_p[ind], x = AA_sites[ind], pch=20, col = "red")
  # subset text label
  ind1 <- which(AA_sites>=1 & AA_sites<=30 & log_p > h)
  if (length(ind1) != 0) {
    text(AA_sites[ind1]-3, log_p[ind1],labels = AA_sites[ind1], pos = 4, cex = 0.5)
  }
  
  #------------------------------------------
  # Blue-ish/purple-ish colors = gp120
  ind <- which(AA_sites>=31 & AA_sites<=511)
  points(y = log_p[ind], x = AA_sites[ind], pch=20, col = "blue")
  # subset text label
  ind1 <- which(AA_sites>=31 & AA_sites<=511 & log_p > h)
  if (length(ind1) != 0) {
    text(AA_sites[ind1]-3, log_p[ind1],labels = AA_sites[ind1], pos = 4, cex = 0.5)
  }
  
  #------------------------------------------
  # Green-ish colors = gp41
  ind <- which(AA_sites>=512 & AA_sites<=856)
  points(y = log_p[ind], x = AA_sites[ind], pch=20, col = "chartreuse3")
  # subset text label
  ind1 <- which(AA_sites>=512 & AA_sites<=856 & log_p > h)
  if (length(ind1) != 0) {
    text(AA_sites[ind1]-3, log_p[ind1],labels = AA_sites[ind1], pos = 4, cex = 0.5)
  }
  
  legend(850, par('usr')[4]/2, bty="n",cex=1,
         legend = c("signal peptide", "gp120","gp41"), 
         col=c("red", "blue", "chartreuse3"), pch=c(19,19,19))

  dev.off()

}



# directory to save result
# local path
# load dat_AA.rda
dir0 = paste0("H:/Research/RF_TMLE/Codes/20210511_HIV_Visuals/data/")
# load pval_cv.rda
dir1 = paste0("H:/Research/RF_TMLE/Results/20210511_HIV_Visuals/")
# final plots and tables
dir_save = paste0("H:/Research/RF_TMLE/Results/20210511_HIV_Visuals/")

# --------------
#   Parameters
# --------------
# Antibodies
antibody <- c("vrc01", "vrc26", "1074", "PGT121", "PGT145")
Antibody <- c("VRC01", "VRC26.08", "10-1074", "PGT121", "PGT145")
s_max <- c(328, 229, 303, 313, 294)



plots_data <- vector(mode = "list", length = length(antibody))
for (anti in seq_along(antibody)){
  # AA positions
  pos <- 19:(18+s_max[anti])
  # load p-value
  # load(paste0(dir1,"pval_full.rda"))
  load(paste0(dir1,antibody[anti],"/pval_cv.rda"))
  
  # load AA position information
  dat_AA <- get(load(paste0(dir0,"dat_", antibody[anti], ".rda")))
  AA_position <- colnames(dat_AA)[pos]
  # rownames(pval_full) <- AA_position
  rownames(pval_cv) <- AA_position
  rm(AA_position)
  
  # prepare for one large plots
  plots_data[[anti]] <- manhattan_prep(pval = pval_cv[,6], log_p_max = 15)
  
  # load est(CI) list
  load(paste0(dir1,antibody[anti],"/est_list_transform.rda"))
  
  # # max column number for datialed info
  # max(unlist(lapply(est_list, function(x){length(x)})))
  
  site_summary(pval=pval_cv[,6], 
               hline=0.05, 
               dir=dir_save, 
               total_num=856,
               antibody = antibody[anti])
}






###################### June 10
# Big plot
######################
#
manhattan_prep <- function(pval,log_p_max=15){
  log_p = -log10(pval)
  # log_p[which(log_p==Inf, arr.ind = T)] = log_p_max
  log_p[which(log_p > log_p_max, arr.ind = T)] = log_p_max
  
  # AA positions
  f_split <- strsplit(names(pval),split = ".",fixed = T)
  # sites being tested (328 sites / total=856)
  AA_sites <- as.numeric(do.call("c",lapply(f_split, function(x){return(x[2])})))

  return(cbind(AA_sites, log_p))
  }
  

# plots_data <- vector(mode = "list", length = length(antibody))
# for (anti in seq_along(antibody)){
#   # AA positions
#   pos <- 19:(18+s_max[anti])
#   # load p-value
#   # load(paste0(dir1,"pval_full.rda"))
#   load(paste0(dir1,antibody[anti],"/pval_cv.rda"))
#   
#   # load AA position information
#   dat_AA <- get(load(paste0(dir0,"dat_", antibody[anti], ".rda")))
#   AA_position <- colnames(dat_AA)[pos]
#   # rownames(pval_full) <- AA_position
#   rownames(pval_cv) <- AA_position
#   
#   plots_data[[anti]] <- manhattan_prep(pval = pval_cv[,6], log_p_max = 15)
# }
# rm(AA_position,pos,anti,dat_AA,pval_cv)

# plot_dat=plots_data[[1]];total_num=856;hline=0.05;log_p_max=15
manhattan_ind <- function(plot_dat, 
                          total_num=856, 
                          hline=0.05,
                          log_p_max = 15,
                          num){
  
  log_p <- plot_dat[,2]
  plot(y=-log10(rep(1,total_num)), x=seq(total_num),
       xlim = c(0,total_num+100),ylim=c(0,log_p_max+1),
       xlab = "Env Residue (HXB2-referenced)", 
       ylab = expression(paste("-log"[10],"(p-value)",sep='')),
       # main = "AA designations", 
       col="white", pch=19, cex.lab=1.5,
       bty="n", xaxt="n", yaxt="n")
  mtext(side=1, at=450, line=5.5, adj=0.5, cex=1.3,
        paste0("(",LETTERS[num],") ",Antibody[num]))
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
  
  # boundary line
  h = -log10(hline/nrow(plot_dat))
  segments(0, h,856,h, lty=2, lwd=2)
  
  
  AA_sites <- as.numeric(plot_dat[,1])
  AA_sites_none <- seq(total_num)[-AA_sites]
  # sites not being tested
  points(y = rep(0,length(AA_sites_none)), x = AA_sites_none, pch = 20, col = "gray80")
  
  # sites being tested
  points(y = log_p, x = AA_sites, pch=46, col = "black")
  #------------------------------------------
  # Red-ish colors = signal peptide
  ind <- which(AA_sites>=1 & AA_sites<=30)
  points(y = log_p[ind], x = AA_sites[ind], pch=17, col = "firebrick2")
  # subset text label
  ind1 <- which(AA_sites>=1 & AA_sites<=30 & log_p > h)
  if (length(ind1) != 0) {
    text(AA_sites[ind1]-3, log_p[ind1],labels = AA_sites[ind1], pos = 4, cex = 1)
  }
  
  #------------------------------------------
  # Blue-ish/purple-ish colors = gp120
  ind <- which(AA_sites>=31 & AA_sites<=511)
  points(y = log_p[ind], x = AA_sites[ind], pch=18, col = "dodgerblue2")
  # subset text label
  ind1 <- which(AA_sites>=31 & AA_sites<=511 & log_p > h)
  if (length(ind1) != 0) {
    text(AA_sites[ind1]-3, log_p[ind1],labels = AA_sites[ind1], pos = 4, cex = 1)
  }
  
  #------------------------------------------
  # Green-ish colors = gp41
  ind <- which(AA_sites>=512 & AA_sites<=856)
  points(y = log_p[ind], x = AA_sites[ind], pch=20, col = "olivedrab3")
  # subset text label
  ind1 <- which(AA_sites>=512 & AA_sites<=856 & log_p > h)
  if (length(ind1) != 0) {
    text(AA_sites[ind1]-3, log_p[ind1],labels = AA_sites[ind1], pos = 4, cex = 1)
  }
  
  # legend(850, par('usr')[4]/2, bty="n",cex=1,
  #        legend = c("signal peptide", "gp120","gp41"), 
  #        col=c("red", "blue", "chartreuse3"), pch=c(19,19,19))
  # 
  # 
}




pdf(file = paste0(dir1, "AA_designations_all.pdf"),width = 12,height = 12)
# mar: c(bottom, left, top, right)
par(oma = c(2,1,1,0.1), mfrow = c(3,2), mar = c(6, 5, 1, 0.1))
for (i in seq_along(antibody)){
  manhattan_ind(plot_dat=plots_data[[i]],num=i)
}
# par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n',xlab="",ylab="")
legend("center", xpd = TRUE, horiz = FALSE, cex = 1.7, bty = 'n',
       x.intersp=1, y.intersp=1.2, #line space between each notation
       legend = c("signal peptide", "gp120","gp41"), 
       col=c("firebrick2", "dodgerblue2", "olivedrab3"), 
       pch=c(17,18,20))
# xpd = TRUE makes the legend plot to the figure
dev.off()







##### sig_site
# Antibodies
antibody <- c("vrc01", "vrc26", "PGT145", "PGT121", "1074")
Antibody <- c("VRC01", "VRC26.08", "PGT145", "PGT121", "10-1074")
s_max <- c(328, 229, 294, 313, 303)


list_antibody_ordered <- vector(mode = "list", length = length(antibody))
for (i in seq(length(antibody))) {
  list_antibody_ordered[[i]] <- read.csv(paste0(dir1,"Ordered_Summary_Sites_Details_",antibody[i],".csv"))
}
names(list_antibody_ordered) <- Antibody

lapply(list_antibody_ordered, function(x){return(c(sum(x$significant),
                                                   sum(x$significant)/nrow(x)))})

# $VRC01
# [1] 23.00000000  0.07012195
# 
# $VRC26.08
# [1] 3.00000000 0.01310044
# 
# $PGT145
# [1] 3.00000000 0.01020408
# 
# $PGT121
# [1] 8.00000000 0.02555911
# 
# $`10-1074`
# [1] 8.00000000 0.02640264



