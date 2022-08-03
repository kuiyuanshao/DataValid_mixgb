library(mixgb)
library(mice)
library(missRanger)
library(mitools)
library(tidyverse)
library(MASS)
resultFun <- function(Nsim, n, n2, beta, e_U, mx, sx, 
                      zrange, zprob, print_step, res, results_est, tool, method){
  file_est  <- paste0(tool, "_", method, 'MI_eU',e_U[1],'eU',e_U[2],'_Zprob',zprob,"_beta",
                      beta[2],"_betaZ",beta[3],"_N",n,"_n",n2,".txt")
  for (nsim in 1:Nsim) {
    set.seed(210818)
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
    dat_srs   <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), 
                                 X = ifelse(R==0, NA, X))
    if (tool == "mixgb"){
      res[[1]]  <- mixgbmi(data=dat_srs, N=20, method)$coefs
    }else if (tool == "mice"){
      res[[1]] <- micemi(data=data_srs, N=20, method = method)$coefs
    }else{
      res[[1]] <- mrangermi(data=data_srs, N=20)$coefs
    }
    
    ## SSRS
    id_phase2 <- c(sample((1:n)[data$Z==0], n2/2), sample((1:n)[data$Z==1], n2/2))
    dat_ssrs  <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), 
                                 X = ifelse(R==0, NA, X))
    if (tool == "mixgb"){
      res[[2]]  <- mixgbmi(data=dat_ssrs, N=20, method)$coefs
    }else if (tool == "mice"){
      res[[2]] <- micemi(data=data_ssrs, N=20, method = method)$coefs
    }else{
      res[[2]] <- mrangermi(data=data_ssrs, N=20)$coefs
    }
    
    ## ODS
    order_Y   <- order(data$Y_tilde)
    id_phase2 <- c(order_Y[1:(n2/2)], order_Y[(n-n2/2+1):n])
    dat_ods   <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), 
                                 X = ifelse(R==0, NA, X))
    if (tool == "mixgb"){
      res[[3]]  <- mixgbmi(data=dat_ods, N=20, method)$coefs
    }else if (tool == "mice"){
      res[[3]] <- micemi(data=data_ods, N=20, method = method)$coefs
    }else{
      res[[3]] <- mrangermi(data=data_ods, N=20)$coefs
    }
    
    ## RS
    order_Y   <- order(residuals(lm(Y_tilde ~ Z, data=data)))
    id_phase2 <- c(order_Y[1:(n2/2)], order_Y[(n-n2/2+1):n])
    dat_rs    <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), 
                                 X = ifelse(R==0, NA, X))
    if (tool == "mixgb"){
      res[[4]]  <- mixgbmi(data=dat_rs, N=20, method)$coefs
    }else if (tool == "mice"){
      res[[4]] <- micemi(data=data_rs, N=20, method = method)$coefs
    }else{
      res[[4]] <- mrangermi(data=data_rs, N=20)$coefs
    }
    
    ## WRS
    order_Y   <- order(residuals(lm(Y_tilde ~ Z, data=data))*vt[data$Z+1])
    id_phase2 <- c(order_Y[1:(n2/2)], order_Y[(n-n2/2+1):n])
    dat_wrs   <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), 
                                 X = ifelse(R==0, NA, X))
    if (tool == "mixgb"){
      res[[5]]  <- mixgbmi(data=dat_wrs, N=20, method)$coefs
    }else if (tool == "mice"){
      res[[5]] <- micemi(data=data_wrs, N=20, method = method)$coefs
    }else{
      res[[5]] <- mrangermi(data=data_wrs, N=20)$coefs
    }
    
    ## SRS
    order_Y   <- order(residuals(lm(Y_tilde ~ Z, data=data))*data$X_tilde)
    id_phase2 <- c(order_Y[1:(n2/2)], order_Y[(n-n2/2+1):n])
    dat_if   <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), 
                                X = ifelse(R==0, NA, X))
    if (tool == "mixgb"){
      res[[6]]  <- mixgbmi(data=dat_if, N=20, method)$coefs
    }else if (tool == "mice"){
      res[[6]] <- micemi(data=data_if, N=20, method = method)$coefs
    }else{
      res[[6]] <- mrangermi(data=data_if, N=20)$coefs
    }
    
    
    ## SRS2
    qs3 <- c(.19,.81)
    IFnods  <- NeymanAllocation(data=data, q=qs3, beta=beta, type='IF', n2=n2)$res2[,1]
    IFnods3 <- IFnods[order(names(IFnods))]
    IFRS    <- data$X_tilde*resid(lm(Y_tilde ~ Z, data=data))
    IFRS_strata <- func_cut(IFRS, qs3)
    samps    <- func_samp(IFRS_strata, IFnods3)
    dat_srs2 <- data %>% mutate(R = ifelse(c(1:n) %in% samps, 1, 0), 
                                X = ifelse(R==0, NA, X))
    if (tool == "mixgb"){
      res[[7]]  <- mixgbmi(data=dat_srs2, N=20, method)$coefs
    }else if (tool == "mice"){
      res[[7]] <- micemi(data=data_srs2, N=20, method = method)$coefs
    }else{
      res[[7]] <- mrangermi(data=data_srs2, N=20)$coefs
    }
    
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
  return (results_est)
}