library(mixgb)
library(mice)
library(missRanger)
library(mitools)
library(tidyverse)
library(MASS)
library(survey)
#pmm.type = 0, 1, 2
mixgbmi <- function(data, N, method, strategy = "else"){
  if (strategy == "srs"){
    design <- svydesign(id = ~ 1, prob = ~ prob, data = data) 
    lm <- svyglm(X ~ X_tilde + Y + Z, design = design, family = gaussian())
  }else if (strategy == "ssrs"){
    design <- svydesign(id = ~ 1, strata = ~ Z, prob = ~ prob, data = data)
    lm <- svyglm(X ~ X_tilde + Y + Z, design = design, family = gaussian())
  }else if (strategy == "srs2"){
    design <- svydesign(id = ~ 1, strata = ~ strata, prob = ~ prob, data = data)
    lm <- svyglm(X ~ X_tilde + Y + Z, design = design, family = gaussian())
  }else{
    lm <- lm(X ~ X_tilde + Y + Z, data = data)
  }
  data$resid[data$R == 1] <- residuals(lm)
  data$fitted[data$R == 1] <- lm$fitted.values
  data$fitted[data$R == 0] <- predict(lm, data[data$R == 0, ], type = "response")
  resid_data <- data[, -4]
  nround <- NULL
  for (i in 1:10){
    cv.results <- mixgb_cv(resid_data, verbose = F, response = "resid")
    nround[i] <- cv.results$best.nrounds
  }
  nround <- as.integer(mean(nround))
  midata <- mixgb(resid_data, m = N, nrounds = nround)
  modcoef <- modvar <- NULL
  for (i in 1:N){
    midata[[i]]$X <- NA
    yhatobs <- data$fitted[data$R == 1] + midata[[i]]$resid[data$R == 1]
    yhatmis <- data$fitted[data$R == 0] + midata[[i]]$resid[data$R == 0]
    yobs <- data$X[data$R == 1]
    idx <- mice::matchindex(yhatobs, yhatmis, k = 5)
    yhatmis <- yobs[idx]
    
    midata[[i]]$X[data$R == 1] <- data$X[data$R == 1]
    midata[[i]]$X[data$R == 0] <- yhatmis
    
    modi <- lm(Y ~ X + Z, data=midata[[i]])
    modcoef[i] <- modi$coeff[2]
    modvar[i] <- summary(modi)$coeff[2,2]^2
  }
  
  coefs <- mean(modcoef)
  vars <- mean(modvar) + (N + 1) * var(modcoef) / N
  return (list(coefs=coefs, vars=vars))
}
#logreg, logreg.boot, cart, rf
micemi <- function(data, N, method){
  midata <- mice(data, m = N, print=FALSE, method = method)
  fulldata <- complete(midata, 'all')
  modcoef <- modvar <- NULL
  for (i in 1:N){
    modi <- lm(Y ~ X + Z, data=fulldata[[i]])
    modcoef[i] <- coef(modi)[2]
    modvar[i] <- vcov(modi)[2, 2]
  }
  coefs <- mean(modcoef)
  vars <- mean(modvar) + (N + 1) * var(modcoef) / N
  return (list(coefs=coefs, vars=vars))
}

mrangermi <- function(data, N){
  midata <- replicate(n = N,missRanger(data = data, num.trees=100, num.threads=4, pmm.k=5), 
                      simplify=FALSE)
  fulldata <- imputationList(midata)
  modcoef <- modvar <- NULL
  for (i in 1:N){
    modi <- lm(Y ~ X + Z, data=fulldata[[i]])
    modcoef[i] <- coef(modi)[2]
    modvar[i] <- vcov(modi)[2, 2]
  }
  coefs <- mean(modcoef)
  vars <- mean(modvar) + (N + 1) * var(modcoef) / N
  return (list(coefs=coefs, vars=vars))
}


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
func_samp <- function(x,y,z=2000) {
  samps <- c()
  for (jj in 1:length(table(x))){
    samps <- c(samps, sample((1:z)[x==jj], y[jj]))
  }
  samps
}
cal <- function(results_est){
  temp <- NULL
  for (i in 1:Nsim){
    temp <- cbind(temp, unlist(results_est[, i]))
  }
  bias <- rowMeans(temp) - 1
  se   <- apply(temp, 1, var)
  
  bias_m <- rbind(bias, bias, bias)
  se_m <- rbind(se, se, se)
  
  mse <- t(bias_m)^2 + t(se_m)
  return (list(round(bias*1000, 3), round(se*1000, 3), round(mse[, 1]*1000,3)))
}
