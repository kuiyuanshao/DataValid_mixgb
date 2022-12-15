



normalize <- function(data, parameters = NULL){
  '
  Normalize data in [0, 1] range.
  
  Args:
    - data: original data
  
  Returns:
    - norm_data: normalized data
    - norm_parameters: min_val, max_val for each feature for renormalization
  '
  
  nVar <- dim(data)[2]
  norm_data <- data
  
  if (is.null(parameters)){
    min_val <- max_val <- rep(0, nVar)
    for (i in 1:nVar){
      min_val[i] <- min(norm_data[, i])
      norm_data[, i] <- norm_data[, i] - min_val[i]
      max_val[i] <- max(norm_data[, i])
      norm_data[, i] <- norm_data[, i] / (max_val[i] + 1e-6)
    }
    norm_parameters <- list(min_val = min_val,
                            max_val = max_val)
  }else{
    min_val <- parameters$min_val
    max_val <- parameters$max_val
    
    for (i in 1:nVar){
      norm_data[, i] <- norm_data[, i] - min_val[i]
      norm_data[, i] <- norm_data[, i] / (max_val[i] + 1e-6)
    }
    norm_parameters <- parameters
  }
  return (list(norm_data = norm_data, norm_parameters = norm_parameters))
}

renormalize <- function(norm_data, norm_parameters){
  '
  Renormalize data from [0, 1] range to the original range.
  
  Args:
    - norm_data: normalized data
    - norm_parameters: min_val, max_val for each feature for renormalization
  
  Returns:
    - renorm_data: renormalized original data
  '
  min_val <- norm_parameters$min_val
  max_val <- norm_parameters$max_val
  
  nVar <- dim(norm_data)[2]
  renorm_data <- norm_data
  for (i in 1:nVar){
    renorm_data[, i] <- renorm_data[, i] * (max_val[i] + 1e-6)
    renorm_data[, i] <- renorm_data[, i] + min_val[i]
  }
  
  return (renorm_data)
}

rounding <- function(imputed_data, data_x){
  '
  Round imputed data for categorical variables.
  
  Args:
    - imputed_data: imputed data
    - data_x: original data with missing values
    
  Returns:
    - rounded_data: rounded imputed data
  '
  nVar <- dim(norm_data)[2]
  rounded_data <- imputed_data
  
  for (i in 1:nVar){
    temp <- data_x[!is.na(data_x[, i]), i]
    
    if (length(unique(temp)) < 20){
      rounded_data[, i] <- round(rounded_data[, i])
    }
  }
  
  return (rounded_data)
}

xavier_init <- function(size){
  '
  Xavier initialization.
  
  Args:
    - size: vector size
    
  Returns:
    - initialized random vector.
  '
  in_dim <- size[1]
  xavier_stddev <- 1 / tf$sqrt(in_dim / 2)
  
  return (tf$random_normal(shape = size, stddev = xavier_stddev))
}

uniform_sampler <- function(low, high, rows, cols){
  '
  Sample uniform random variables.
  
  Args:
    - low: low limit
    - high: high limit
    - rows: the number of rows
    - cols: the number of columns
    
  Returns:
    - uniform_random_matrix: generated uniform random matrix.
  '
  random_unif <- runif(min = low, max = high, n = rows * cols)
  unif_matrix <- matrix(random_unif, nrow = rows, ncol = cols, byrow = T)
  return (unif_matrix)
}



sample_batch_index <- function(total, batch_size){
  'Sample index of the mini-batch.
  
  Args:
    - total: total number of samples
    - batch_size: batch size
    
  Returns:
    - batch_idx: batch index
  '
  total_idx <- sample(total)
  batch_idx <- total_idx[1:batch_size]
  return (batch_idx)

}













