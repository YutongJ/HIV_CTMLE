##############################
#####  required package  #####
##############################
# n <- 5; sigma <- 1.5; rho <- 0.5; t<-2
AR1 <- function(n,sigma,rho,t){
  delta <- 0:(n-1)
  out <- outer(delta,delta,"-")
  if (t!=0){
    out[1:t,] <- Inf
    out[,1:t] <- Inf
    diag(out) <- 0
    Sigma <- sigma^2*rho^abs(out)
  }
  else{
    diag(out) <- 0
    Sigma <- sigma^2*rho^abs(out)
  }
  return(Sigma)
}

var_mtx <- AR1(n = rvs.size, sigma = 0.9, rho = rhos, t=indeps.size)




# sample.size: the number of observations (rows)
# rv.size    : the number of AA (columns)
# cov        : AR-1 correlation matrix used for multivariate normal data generation
# prob       : probability for each level. Example: prob = c(0.1,0.4,0.3,0.2) which means: 0.1A, 0.4B,0.3C, 0.2D
AA_generate <- function(sample.size=samples.size,rv.size=rvs.size,cov=var_mtx,prob=c(0.1,0.4,0.3,0.2)){
  mu <- rep(0,rv.size)
  
  dat <- MASS::mvrnorm(sample.size,mu,Sigma=cov)  #mu <- rep(0,rvs.size);dat <- MASS::mvrnorm(200,mu,Sigma=var_mtx)
  dat = findInterval(dat,qnorm(cumsum(c(0,prob))),left.open=TRUE)
  
  # use qnorm(cumprob) as bars to separate them into 4 groups
  dat = matrix(as.character(dat), ncol=rv.size,byrow = FALSE)
  
  # dat <- apply(dat, 2, function(x){ as.numeric(x > 0)})
  colnames(dat) <- paste0("W",1:rv.size)
  dat <- as.data.frame(dat,stringsAsFactors = F)
  return(dat)
}
# dat_W = AA_generate(prob=c(0.5,0.5))



# datt   : datasets generated from function AA_generate()
# level  : the levels of AA
whole_dat = function(datt,level){
  # True outcome regression (OR) model
  Qbar <- function(datt,betas=beta,levels,wind=w_ind){
    return(plogis(as.matrix(datt[,match(paste0("W",rep(wind,each=levels),1:levels),colnames(datt))])%*%betas))
  }
  
  Dummy.dat = dummies::dummy.data.frame(datt)
  Qbar_obs <- Qbar(Dummy.dat,levels=level)
  # set.seed(6543)
  Y <- rbinom(samples.size, 1, Qbar_obs)
  dat <- data.frame(Y,datt,stringsAsFactors = F)
  return(dat)
}
# dat = whole_dat(dat_W,level = 2)
