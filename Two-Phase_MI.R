library(tidyverse)
library(mvtnorm)
library(MASS)

######################################
set.seed(210818)

Nsim <- 1000
n  <- 2000
n2 <- 500
######################################

beta <- c(1, 1, 1)
e_U <- c(sqrt(3),sqrt(3))
mx <- 0; sx <- 1; zrange <- 1; zprob <- .5

#######################################

print_step <- 100
res <- results_est <- list(0)
#######################################

#pre-allocate different values of sigma 
file_est  <- paste0('MI_eU',e_U[1],'eU',e_U[2],'_Zprob',zprob,"_beta",beta[2],"_betaZ",beta[3],"_N",n,"_n",n2,".txt")

## MI
#######################################

f_imp <- function(data, Nimpute) {
  #An imputation model
  modx   <- lm(X ~ X_tilde + Y + Z, subset = (R==1), data=data)
  
  #retrieves the residual standard error of the model
  #square it to a variance
  #retrieves the degrees of freedom of the residual standard error
  #random generates Nimpute number of chisquare based with degrees of freedom of the residual standard error
  
  #it is an unbiased estimation of the error variance:
  #rchisq() simulating the real ?
  sigmas <- sqrt(summary(modx)$sigma^2*rchisq(Nimpute, summary(modx)$df[2], ncp=0)/summary(modx)$df[2])
  
  #generates random multivariate normal based on the initial coefficients of the model 
  #scale the variance of the simulation by the residual variance and the variance-covariance matrix
  desmatrix <- rmvnorm(Nimpute, mean=modx$coeff, sigma=summary(modx)$sigma^2*summary(modx)$cov)
  #simulates the coefficients again based on the mle and the variance
  x.coeff.mi1 <- x.coeff.var.mi1 <- NULL
  for (k in 1:Nimpute) {
    # intercept, and other variables.
    # add a error term with mean zero and the true unbiased estimation of the error variance we simulated previously
    # inside the as.vector() is the initial estimation on X based on the coefficients simulated.
    pred.x   <- as.vector(desmatrix[k,] %*% t(cbind(1, data$X_tilde, data$Y , data$Z)) + rnorm(nrow(data), 0, sigmas[k]))
    # add predicted X to the dataframe
    x.impute <- ifelse(data$R==0, pred.x, data$X)
    dat <- data.frame(data, Ximp = x.impute)
    # refit the model using the imputed data
    mod.imputed <- lm(Y ~ Ximp + Z, data=dat)
    # save the coef and variance.
    x.coeff.mi1[k]     <- mod.imputed$coeff[2]
    x.coeff.var.mi1[k] <- summary(mod.imputed)$coeff[2,2]^2
  }
  # return it
  coefs <- mean(x.coeff.mi1)
  # why 
  vars  <- mean(x.coeff.var.mi1)+(Nimpute+1)*var(x.coeff.mi1)/Nimpute
  return(list(coefs=coefs, vars=vars))
}


#### Neyman allocation
NeymanAllocation <- function(data,qs,beta=c(1,1,1),type,binary=FALSE,n=1e6,missclas_prop=0.05,n2=500) {
  # how qs is specificed?
  #data <- gen_data(n, beta=beta)
  
  if (!binary) {
    
    if (type=='Yds') {
      # cut them by the quantiles, labelling it by which category it belongs to.
      Y1 <- cut(data$Y_tilde, breaks=c(-Inf, quantile(data$Y_tilde, probs=qs), Inf), labels=paste(1:(length(qs)+1), sep=','))
    }
    if (type=='Xds') {
      Y1 <- cut(data$X_tilde, breaks=c(-Inf, quantile(data$X_tilde, probs=qs), Inf), labels=paste(1:(length(qs)+1), sep=','))
    }
    if (type=='RS') {
      if (is.null(data$Z)) {
        #use the residuals of the simple model
        rs <- resid(lm(Y_tilde ~ X_tilde, data=data))
      } else {
        rs <- resid(lm(Y_tilde ~ X_tilde + Z, data=data))
      }
      # cut the residuals by the probability
      Y1 <- cut(rs, breaks=c(-Inf, quantile(rs, probs=qs), Inf), labels=paste(1:(length(qs)+1), sep=','))
    }
    if (type=='WRS') {
      # variance of the residual
      varX1 <- var(residuals(lm(X_tilde ~ Z, data=data))[data$Z==1])
      varX2 <- var(residuals(lm(X_tilde ~ Z, data=data))[data$Z==2])
      varX <- c(varX1, varX2)[data$Z]
      if (is.null(data$Z)) {
        wrs <- resid(lm(Y_tilde ~ X_tilde, data=data))
      } else {
        #times the residual by the variance of X.
        wrs <- resid(lm(Y_tilde ~ X_tilde + Z, data=data))*varX
      }
      Y1 <- cut(wrs, breaks=c(-Inf, quantile(wrs, probs=qs), Inf), labels=paste(1:(length(qs)+1), sep=','))
    }
    if (type=='IF') {
      if (is.null(data$Z)) {
        ifd <- data$X_tilde*resid(lm(Y_tilde ~ X_tilde, data=data))
      } else {
        ifd <- data$X_tilde*resid(lm(Y_tilde ~ X_tilde + Z, data=data))
      }
      Y1 <- cut(ifd, breaks=c(-Inf, quantile(ifd, probs=qs), Inf), labels=paste(1:(length(qs)+1), sep=','))
    }
    if (type=='SR') {
      S1   <- data$X_tilde*resid(lm(Y_tilde ~ X_tilde, data=data))
      S2   <- resid(lm(Y_tilde ~ X_tilde, data=data))
      S3   <- resid(lm(Y_tilde ~ X_tilde, data=data))^2
      if (is.null(data$Z)) {
        S1S2 <- resid(lm(S1 ~ S2+S3))
      } else {
        S4   <- data$Z*resid(lm(Y_tilde ~ X_tilde, data=data))
        S1S2 <- resid(lm(S1 ~ S2+S3+S4))
      }
      Y1 <- cut(S1S2, breaks=c(-Inf, quantile(S1S2, probs=qs), Inf), labels=paste(1:(length(qs)+1), sep=','))
    }
    
    mu1 <- residuals(lm(Y_tilde ~ X_tilde, data=data))
    if (!is.null(data$Z))
      mu1 <- residuals(lm(Y_tilde ~ X_tilde + Z, data=data))
    
    s1 <- sapply(1:(length(qs)+1), FUN=function(x) sd((data$X_tilde*mu1)[Y1==x]))
    
    # Get/retun proportions
    res <- table(Y1)*s1/sum(table(Y1)*s1)
    
    res2 <- cbind(round(res*n2), table(Y1))
    res2 <- res2[order(res2[,1], decreasing=TRUE),]
    while(any(res2[,1] > res2[,2]))
    {
      pos <- which(res2[,1] > res2[,2]); pos <- pos[1]
      res2. <- matrix(c(res2[res2[,1] < res2[,2],]), ncol=2)
      res2[1,1] <- res2.[1,1] + (res2[pos,1] - res2[pos,2])
      res2[pos,1] <- res2[pos,2]
      res2[res2[,1] < res2[,2],] <- res2.
    }
  }
  if (binary){
    
    # Generate data and calculate Influence function
    data <- gen_data(N, beta=beta, binary=binary, missclas_prop=missclas_prop)
    mu <- sum(data$Y_tilde==1)/length(data$Y_tilde)
    U <- data$X_tilde*(data$Y_tilde - mu)
    
    # Get SD for each stratum
    Y1 <- data$Y_tilde==1
    Y0 <- data$Y_tilde==0
    
    # Get SD for each stratum
    s1 <- sd(U[Y1])
    s2 <- sd(U[Y0])
    
    # Get influence function
    mu1 <- sum(data$Y_tilde==1)/length(data$Y_tilde)
    mu0 <- sum(data$Y_tilde==1)/length(data$Y_tilde)
    U1 <- data$X_tilde[Y1]*(data$Y_tilde[Y1] - mu1)
    U0 <- data$X_tilde[Y0]*(data$Y_tilde[Y0] - mu0)
    
    # Get SD for each stratum
    s1 <- sd(U[Y1])
    s0 <- sd(U[Y0])
    
    # Get/return proportions
    p1 <- mu / (mu + (1-mu)*s2/s1) 
    res <- c(p1, (1 - p1))
  }
  return(list(res=res, res2=res2))
}


func_cut <- function(x, qs) 
  cut(x, breaks=c(-Inf, quantile(x, probs=qs), Inf), labels=paste(1:(length(qs)+1), sep=','))
func_samp <- function(x,y,z=nrow(data)) {
  samps <- c()
  for (jj in 1:length(table(x))){
    samps <- c(samps, sample((1:z)[x==jj], y[jj]))
  }
  samps
}





############################################################
###################  Analysis  #############################
############################################################

for (nsim in 1:Nsim) {
  #
  simZ   <- rbinom(n, zrange, zprob)
  simX   <- (1-simZ)*rnorm(n, 0, 1) + simZ*rnorm(n, 0.5, 1)
  #error
  epsilon <- rnorm(n, 0, 1)
  #Y is added with epsilon
  simY    <- beta[1] + beta[2]*simX + beta[3]*simZ + epsilon
  #make X_star depends on it
  simX_tilde <- simX + rnorm(n, 0, e_U[1]*(simZ==0) + e_U[2]*(simZ==1))
  data <- data.frame(Y_tilde=simY, X_tilde=simX_tilde, Y=simY, X=simX, Z=simZ)
  
  if (nsim == 1){
    n_t <- 1e7
    simZ_t   <- rbinom(n_t, zrange, zprob)
    simX_t   <- (1-simZ_t)*rnorm(n_t, 0, 1) + simZ_t*rnorm(n_t, 0.5, 1)
    epsilon_t <- rnorm(n_t, 0, 1)
    simY_t    <- beta[1] + beta[2]*simX_t + beta[3]*simZ_t + epsilon_t
    simX_tilde_t <- simX_t + rnorm(n_t, 0, e_U[1]*(simZ_t==0) + e_U[2]*(simZ_t==1))
    data_t <- data.frame(Y_tilde=simY_t, X_tilde=simX_tilde_t, Y=simY_t, X=simX_t, Z=simZ_t)
    
    v1 <- sd(resid(lm(X ~ X_tilde + Z, data=data_t))[data_t$Z==0])
    v2 <- sd(resid(lm(X ~ X_tilde + Z, data=data_t))[data_t$Z==1])
    vt <- c(v1, v2)
  }
  
  ##### Designs
  ## SRS
  id_phase2 <- c(sample(n, n2))
  dat_srs   <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), X = ifelse(R==0, NA, X))
  res[[1]]  <- f_imp(data=dat_srs, Nimpute=20)$coefs
  
  ## SSRS
  id_phase2 <- c(sample((1:n)[data$Z==0], n2/2), sample((1:n)[data$Z==1], n2/2))
  dat_ssrs  <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), X = ifelse(R==0, NA, X))
  res[[2]]  <- f_imp(data=dat_ssrs, Nimpute=20)$coefs
  
  ## ODS
  order_Y   <- order(data$Y_tilde)
  id_phase2 <- c(order_Y[1:(n2/2)], order_Y[(n-n2/2+1):n])
  dat_ods   <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), X = ifelse(R==0, NA, X))
  res[[3]]  <- f_imp(data=dat_ods, Nimpute=20)$coefs
  
  ## RS
  order_Y   <- order(residuals(lm(Y_tilde ~ Z, data=data)))
  id_phase2 <- c(order_Y[1:(n2/2)], order_Y[(n-n2/2+1):n])
  dat_rs    <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), X = ifelse(R==0, NA, X))
  res[[4]]  <- f_imp(data=dat_rs, Nimpute=20)$coefs
  
  ## WRS
  order_Y   <- order(residuals(lm(Y_tilde ~ Z, data=data))*vt[data$Z+1])
  id_phase2 <- c(order_Y[1:(n2/2)], order_Y[(n-n2/2+1):n])
  dat_wrs   <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), X = ifelse(R==0, NA, X))
  res[[5]]  <- f_imp(data=dat_wrs, Nimpute=20)$coefs
  
  ## SRS
  order_Y   <- order(residuals(lm(Y_tilde ~ Z, data=data))*data$X_tilde)
  id_phase2 <- c(order_Y[1:(n2/2)], order_Y[(n-n2/2+1):n])
  dat_if   <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), X = ifelse(R==0, NA, X))
  res[[6]] <- f_imp(data=dat_if, Nimpute=20)$coefs
  
  #  ## SRS_2
  #  order_Y   <- order(residuals(lm(Y_tilde ~ X_tilde + Z, data=data))*data$X_tilde)
  #  id_phase2 <- c(order_Y[1:(n2/2)], order_Y[(n-n2/2+1):n])
  #  dat_if   <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), X = ifelse(R==0, NA, X))
  #  res[[7]] <- smle_MEXY(Y="Y", X="X", Z="Z", Y_tilde="Y_tilde", X_tilde="X_tilde",
  #                        Bspline=colnames(Bspline), data=dat_if, noSE=TRUE)
  
  ## SRS2
  qs3 <- c(.19,.81)
  IFnods  <- NeymanAllocation(data=data, q=qs3, beta=beta, type='IF', n2=n2)$res2[,1]
  IFnods3 <- IFnods[order(names(IFnods))]
  IFRS    <- data$X_tilde*resid(lm(Y_tilde ~ Z, data=data))
  IFRS_strata <- func_cut(IFRS, qs3)
  samps    <- func_samp(IFRS_strata, IFnods3)
  dat_srs2 <- data %>% mutate(R = ifelse(c(1:n) %in% samps, 1, 0), X = ifelse(R==0, NA, X))
  res[[7]] <- f_imp(data=dat_srs2, Nimpute=20)$coefs
  
  if (nsim == 1){
    results_conv <- matrix(NA, nrow=Nsim, ncol=length(res))
    for (i in 1:length(res)) {
      results_est[[i]]  <- matrix(NA, nrow=Nsim, ncol=3)
      colnames(results_est[[i]]) <- c("intercept", "X", "Z")
    }
  }
  
  ## Save and print
  for (i in 1:length(res)){
    results_est[[i]][nsim,] <- t(res[[i]])
  }
}

bias <- sapply(results_est, FUN = function(x) colMeans(x[1:Nsim,], na.rm=TRUE)) -
  matrix(rep(beta), length(results_est), nrow=length(beta))
se   <- sapply(results_est, FUN = function(x) apply(x[1:Nsim,], 2, var, na.rm=TRUE))
mse <- t(bias)^2 + t(se)
mse_ratio <- mse[1,2]/mse[,2]
lm_imp <- (list(round(bias*1000, 3), round(se*1000, 3), round(mse_ratio*1000,3)))
