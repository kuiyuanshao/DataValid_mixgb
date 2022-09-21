library(mixgb)
library(mice)
library(missRanger)
library(mitools)
library(tidyverse)
library(MASS)
resultFun2 <- function(Nsim, n, n2, beta, e_U, mx, sx, 
                      zrange, zprob, print_step, res, results_est, tool, method,
                      weights = T){
  res <- list (0)
  simZ   <- rbinom(n, zrange, zprob)
  simX   <- (1-simZ)*rnorm(n, 0, 1) + simZ*rnorm(n, 0.5, 1)
  #error
  epsilon <- rnorm(n, 0, 1)
  #Y is added with epsilon
  simY    <- beta[1] + beta[2]*simX + beta[3]*simZ + epsilon
  #make X_star depends on it
  simX_tilde <- simX + rnorm(n, 0, e_U[1]*(simZ==0) + e_U[2]*(simZ==1))
  data <- data.frame(Y_tilde=simY, X_tilde=simX_tilde, Y=simY, X=simX, Z=simZ)
  ##### Designs
  ## SRS
  id_phase2 <- c(sample(n, n2))
  dat_srs   <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), 
                               X = ifelse(R==0, NA, X))
  dat_srs$prob <- 1 / dim(dat_srs)[1]
  if (tool == "mixgb"){
    res[[1]]  <- mixgbmi(data=dat_srs, N=20, method, strategy = "srs")$coefs
  }else if (tool == "mice"){
    res[[1]] <- micemi(data=dat_srs, N=20, method = method)$coefs
  }else{
    res[[1]] <- mrangermi(data=dat_srs, N=20)$coefs
  }
  return (res)
}
resultFun <- function(Nsim, n, n2, beta, e_U, mx, sx, 
                      zrange, zprob, print_step, res, results_est, tool, method,
                      weights = T){
    res <- list (0)
    simZ   <- rbinom(n, zrange, zprob)
    simX   <- (1-simZ)*rnorm(n, 0, 1) + simZ*rnorm(n, 0.5, 1)
    #error
    epsilon <- rnorm(n, 0, 1)
    #Y is added with epsilon
    simY    <- beta[1] + beta[2]*simX + beta[3]*simZ + epsilon
    #make X_star depends on it
    simX_tilde <- simX + rnorm(n, 0, e_U[1]*(simZ==0) + e_U[2]*(simZ==1))
    data <- data.frame(Y_tilde=simY, X_tilde=simX_tilde, Y=simY, X=simX, Z=simZ)
    ##### Designs
    ## SRS
    id_phase2 <- c(sample(n, n2))
    dat_srs   <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), 
                                 X = ifelse(R==0, NA, X))
    dat_srs$prob <- 1 / dim(dat_srs)[1]
    if (tool == "mixgb"){
      res[[1]]  <- mixgbmi(data=dat_srs, N=20, method, strategy = "srs")$coefs
    }else if (tool == "mice"){
      res[[1]] <- micemi(data=dat_srs, N=20, method = method)$coefs
    }else{
      res[[1]] <- mrangermi(data=dat_srs, N=20)$coefs
    }
    
    ## SSRS
    id_phase2 <- c(sample((1:n)[data$Z==0], n2/2), sample((1:n)[data$Z==1], n2/2))
    dat_ssrs  <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), 
                                 X = ifelse(R==0, NA, X))
    dat_ssrs$prob <- NA
    dat_ssrs$prob[which(dat_ssrs$Z == 1)] <- 1 / length(which(dat_ssrs$Z == 1))
    dat_ssrs$prob[which(dat_ssrs$Z == 0)] <- 1 / length(which(dat_ssrs$Z == 0))

    if (tool == "mixgb"){
      res[[2]]  <- mixgbmi(data=dat_ssrs, N=20, method, strategy = "ssrs")$coefs
    }else if (tool == "mice"){
      res[[2]] <- micemi(data=dat_ssrs, N=20, method = method)$coefs
    }else{
      res[[2]] <- mrangermi(data=dat_ssrs, N=20)$coefs
    }
    
    ## ODS
    order_Y   <- order(data$Y_tilde)
    id_phase2 <- c(order_Y[1:(n2/2)], order_Y[(n-n2/2+1):n])
    dat_ods   <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), 
                                 X = ifelse(R==0, NA, X))
    
    if (tool == "mixgb"){
      res[[3]]  <- mixgbmi(data=dat_ods, N=20, method)$coefs
    }else if (tool == "mice"){
      res[[3]] <- micemi(data=dat_ods, N=20, method = method)$coefs
    }else{
      res[[3]] <- mrangermi(data=dat_ods, N=20)$coefs
    }
    
    
    ## RS
    order_Y   <- order(residuals(lm(Y_tilde ~ Z, data=data)))
    id_phase2 <- c(order_Y[1:(n2/2)], order_Y[(n-n2/2+1):n])
    dat_rs    <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), 
                                 X = ifelse(R==0, NA, X))
    dat_rs$resid_z <- residuals(lm(Y_tilde ~ Z, data=data))
    if (tool == "mixgb"){
      res[[4]]  <- mixgbmi(data=dat_rs, N=20, method)$coefs
    }else if (tool == "mice"){
      res[[4]] <- micemi(data=dat_rs, N=20, method = method)$coefs
    }else{
      res[[4]] <- mrangermi(data=dat_rs, N=20)$coefs
    }
    
    ## WRS
    order_Y   <- order(residuals(lm(Y_tilde ~ Z, data=data))*vt[data$Z+1])
    id_phase2 <- c(order_Y[1:(n2/2)], order_Y[(n-n2/2+1):n])
    dat_wrs   <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), 
                                 X = ifelse(R==0, NA, X))
    dat_wrs$weighted_resid <- residuals(lm(Y_tilde ~ Z, data=data))*vt[data$Z+1]

    if (tool == "mixgb"){
      res[[5]]  <- mixgbmi(data=dat_wrs, N=20, method)$coefs
    }else if (tool == "mice"){
      res[[5]] <- micemi(data=dat_wrs, N=20, method = method)$coefs
    }else{
      res[[5]] <- mrangermi(data=dat_wrs, N=20)$coefs
    }
    
    ## IF
    order_Y   <- order(residuals(lm(Y_tilde ~ Z, data=data))*data$X_tilde)
    id_phase2 <- c(order_Y[1:(n2/2)], order_Y[(n-n2/2+1):n])
    dat_if   <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), 
                                X = ifelse(R==0, NA, X))
    dat_if$score <- residuals(lm(Y_tilde ~ Z, data=data))*data$X_tilde
    
    if (tool == "mixgb"){
      res[[6]]  <- mixgbmi(data=dat_if, N=20, method)$coefs
    }else if (tool == "mice"){
      res[[6]] <- micemi(data=dat_if, N=20, method = method)$coefs
    }else{
      res[[6]] <- mrangermi(data=dat_if, N=20)$coefs
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
    dat_srs2$strata <- IFRS_strata
    dat_srs2$prob <- NA
    dat_srs2$prob[which(IFRS_strata == 1)] = 1 / sum(IFRS_strata == 1)
    dat_srs2$prob[which(IFRS_strata == 2)] = 1 / sum(IFRS_strata == 2)
    dat_srs2$prob[which(IFRS_strata == 3)] = 1 / sum(IFRS_strata == 3)
    dat_srs2$score <- IFRS
    
    if (tool == "mixgb"){
      res[[7]]  <- mixgbmi(data=dat_srs2, N=20, method, strategy = "srs2")$coefs
    }else if (tool == "mice"){
      res[[7]] <- micemi(data=dat_srs2, N=20, method = method)$coefs
    }else{
      res[[7]] <- mrangermi(data=dat_srs2, N=20)$coefs
    }
    
    ## Save 
    return (res)
}