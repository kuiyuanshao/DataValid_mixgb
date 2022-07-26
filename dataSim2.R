##############################################################################################################################
######################################################## Simulation ##########################################################
##############################################################################################################################
set.seed(2022)
library(purrr)
library(optimall)

simFun <- function(N, n){
  Z <- rbernoulli(N, p = 0.5)
  X <- (1 - mean(Z))/2 * rnorm(N, sd = 1)
  error <- rnorm(N)
  beta0 <- rnorm(N, mean = 1, sd = 1)
  beta1 <- 0.3
  beta2 <- 0.2
  Y <- beta0 + beta1*X + beta2*Z
  X_star <- X + rnorm(N, mean = 0, sd = 0.5)
  full_data <- data.frame(Y, X, Z, X_star)
  full_data$id <- 1:nrow(full_data)
  ############### Unvalidated Data #############
  unvalidated <- full_data
  unvalidated$X <- NA
  ############### Phase-2 Sampling #############
  ### Optimal Allocation
  sampling_design <- optimum_allocation(data = unvalidated, strata = "Z", 
                                      y = "X_star", nsample = n, method = "WrightII")
  entire <- sample_strata(data = unvalidated, strata = "Z", id = "id", 
                        design_data = sampling_design, design_strata = "strata",
                        n_allocated = "stratum_size")
  ############### Validated Data ###############
  entire$X[entire$sample_indicator == 1] <- full_data$X[entire$sample_indicator == 1]
  subset <- entire[entire$sample_indicator == 1, ]
  return (list(full_data, entire, subset))
}

