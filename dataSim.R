##############################################################################################################################
######################################################## Simulation ##########################################################
##############################################################################################################################
library(mvtnorm)
set.seed(2022)

n <- 50000
y1 <- as.factor(rbinom(n, size = 1, prob = c(0.4, 0.6)))

#Construct correlation
ind <- y1 == 0
y2 <- rep(0, n)
y2[ind] <- rpois(sum(ind), lambda = 50)
y2[!ind] <- rpois(sum(!ind), lambda = 40)

#Consturct a variance-covariance matrix for the 10 X covariates.
cov_m <- matrix(seq(5, 1, length.out = 100), 10)
cov_m[lower.tri(cov_m)] = t(cov_m)[lower.tri(cov_m)]
x_s <- rmvnorm(n, mean = rep(10, 10) + seq(10, 1), sigma = cov_m)
x_s[, 1] <- 0.2*y2
full_data <- data.frame(y1, y2, x1 = x_s[, 1], x2 = x_s[, 2], x3 = x_s[, 3],
                        x4 = x_s[, 4], x5 = x_s[, 5], x6 = x_s[, 6], x7 = x_s[, 7],
                        x8 = x_s[, 8], x9 = x_s[, 9], x10 = x_s[, 10])

rm(y1, y2, x_s)

############### Validated Data #############
R_ind <- sample(1:n, n/10)
validated <- full_data[R_ind, ]
############### Unvalidated Data ###########
unvalidated <- full_data[-R_ind, ]
#Adding measurement errors
n2 <- dim(unvalidated)[1]
unvalidated$y1[sample(1:n2, n2/10)] <- 
  as.factor(rbinom(n2/10, size = 1, prob = c(0.5, 0.5)))
unvalidated$y2 <- unvalidated$y2 + rnorm(n2, 0, sd = 10)
#Measurement errors in all 10 x variables are the same
for (i in seq(3:dim(unvalidated)[2])){
  unvalidated[, i] <- unvalidated[, i] + rnorm(n2, 0, sd = 5) 
}


