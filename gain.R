


gain <- function(data_x, batch_size, hint_rate, alpha, iterations){
  data_m <- 1 - is.na(data_x)
  
  nRow <- dim(data_x)[1]
  nCol <- dim(data_x)[2]
  
  h_dim <- integer(dim)
  
  norm_result <- normalize(data_x)
  norm_data <- norm_result$norm_data
  norm_parameters <- norm_result$norm_parameters
  
  norm_data_x <- data_x
  norm_data_x[is.na(norm_data_x)] <- 0
  
  X <- tf$placeholder(tf$float32, shape = c(NULL, nCol))
  M <- tf$placeholder(tf$float32, shape = c(NULL, nCol))
  H <- tf$placeholder(tf$float32, shape = c(NULL, nCol))
  
  D_W1 <- tf$Variable(xavier_init(c(nCol*2, h_dim)))
  D_b1 <- tf$Variable(tf$zeros(shape = c(h_dim)))
  
  D_W2 <- tf$Variable(xavier_init(c(h_dim, h_dim)))
  D_b2 = tf$Variable(tf$zeros(shape = c(h_dim)))
  
  D_W3 <- tf$Variable(xavier_init(c(h_dim, nCol)))
  D_b3 <- tf$Variable(tf$zeros(shape = c(nCol)))
  
  
  theta_D <- c(D_W1, D_W2, D_W3, D_b1, D_b2, D_b3)
  
  G_W1 <- tf$Variable(xavier_init(c(nCol*2, h_dim)))  
  G_b1 <- tf$Variable(tf$zeros(shape = c(h_dim)))
  
  G_W2 <- tf$Variable(xavier_init(c(h_dim, h_dim)))
  G_b2 <- tf$Variable(tf$zeros(shape = c(h_dim)))
  
  G_W3 <- tf$Variable(xavier_init(c(h_dim, nCol)))
  G_b3 <- tf$Variable(tf$zeros(shape = c(nCOl)))
  
  theta_G <- c(G_W1, G_W2, G_W3, G_b1, G_b2, G_b3)
  
  G_sample <- generator(X, M)
  
  # Combine with observed data
  Hat_X <- X * M + G_sample * (1 - M)
  
  # Discriminator
  D_prob <- discriminator(Hat_X, H)
  
  ## GAIN loss
  D_loss_temp <- -tf$reduce_mean(M * tf$log(D_prob + 1e-8) + 
                                   (1-M) * tf$log(1. - D_prob + 1e-8)) 
  
  G_loss_temp <- -tf$reduce_mean((1-M) * tf$log(D_prob + 1e-8))
  
  MSE_loss <- tf$reduce_mean((M * X - M * G_sample)**2) / tf$reduce_mean(M)
  
  D_loss <- D_loss_temp
  G_loss <- G_loss_temp + alpha * MSE_loss 
  
  ## GAIN solver
  D_solver <- tf$train$AdamOptimizer()$minimize(D_loss, var_list = theta_D)
  G_solver <- tf$train$AdamOptimizer()$minimize(G_loss, var_list = theta_G)
  
  ## Iterations
  sess <- tf$Session()
  sess$run(tf$global_variables_initializer())
  
  
  for (it in 1:iterations){ 
    # Sample batch
    batch_idx <- sample_batch_index(nRow, batch_size)
    X_mb <- norm_data_x[batch_idx, ]  
    M_mb <- data_m[batch_idx, ]  
    # Sample random vectors  
    Z_mb <- uniform_sampler(0, 0.01, batch_size, nCol) 
    # Sample hint vectors
    H_mb_temp <- binary_sampler(hint_rate, batch_size, nCol)
    H_mb <- M_mb * H_mb_temp
  
    # Combine random vectors with observed vectors
    X_mb <- M_mb * X_mb + (1-M_mb) * Z_mb 
  
    D_loss_curr <- sess$run(c(D_solver, D_loss_temp), 
                               feed_dict = dict(M = M_mb, X = X_mb, H = H_mb))[2]
    run_result <- sess$run(c(G_solver, G_loss_temp, MSE_loss),
                                              feed_dict = dict(X = X_mb, M = M_mb, H = H_mb))
    G_loss_curr <- run_result[2]
    
    MSE_loss_curr <- run_result[3]
  }
  
    ## Return imputed data      
  Z_mb <- uniform_sampler(0, 0.01, nRow, nCol) 
  M_mb <- data_m
  X_mb <- norm_data_x          
  X_mb <- M_mb * X_mb + (1-M_mb) * Z_mb 
  
  imputed_data <- sess$run(c(G_sample), feed_dict = dict(X = X_mb, M = M_mb))[1]
  
  imputed_data <- data_m * norm_data_x + (1 - data_m) * imputed_data
  
  # Renormalization
  imputed_data <- renormalize(imputed_data, norm_parameters)  
  
  # Rounding
  imputed_data <- rounding(imputed_data, data_x)  
  
  return (imputed_data)
}



generator <-  function(x, m){
  # Concatenate Mask and Data
  inputs <- tf$concat(values = c(x, m), axis = 1) 
  G_h1 <- tf$nn$relu(tf$matmul(inputs, G_W1) + G_b1)
  G_h2 <- tf$nn$relu(tf$matmul(G_h1, G_W2) + G_b2)   
  # MinMax normalized output
  G_prob <- tf$nn$sigmoid(tf$matmul(G_h2, G_W3) + G_b3) 
  return (G_prob)
}
# Discriminator
discriminator <- function(x, h){
  # Concatenate Data and Hint
  inputs <- tf$concat(values = c(x, h), axis = 1) 
  D_h1 <- tf$nn$relu(tf$matmul(inputs, D_W1) + D_b1)  
  D_h2 <- tf$nn$relu(tf$matmul(D_h1, D_W2) + D_b2)
  D_logit <- tf$matmul(D_h2, D_W3) + D_b3
  D_prob <- tf$nn$sigmoid(D_logit)
  return (D_prob)
  
}

