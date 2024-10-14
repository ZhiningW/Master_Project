rm(list = ls())
################################################################################

display_result <- function(est_lambda, est_psi, real_lambda, real_psi, n, upper_triangle = TRUE, time_to_run){
  p <- nrow(real_lambda)
  m <- ncol(real_lambda)
  if (p < m){
    print('Bad model')
    return(0)
  }
  
  
  # Permutation optimization
  after_perm <- find_optimal_permutation(real_lambda, est_lambda)
  est_lambda <- after_perm$optimal_N
  est_lambda[which(abs(est_lambda) <= 0.1)] <- 0
  
  # MSE
  MSE <- after_perm$min_mse
  
  # Compute the sparsity of the est_lambda
  if (upper_triangle){
    sparsity <- (sum(est_lambda == 0) - (1/2) * (m - 1)^2) / (p * m)
  }
  else {
    sparsity <- sum(est_lambda == 0) / (p * m)
  }
  
  # Compute TPR
  real_zero_positions <- real_lambda == 0
  est_zero_positions <- est_lambda == 0
  true_positives <- real_zero_positions & est_zero_positions
  total_zeros_in_real <- sum(real_zero_positions)
  correctly_predicted_zeros <- sum(true_positives)
  if (upper_triangle){
    # If we enforce the est_lambda to be upper triangle
    TPR <- (correctly_predicted_zeros - (1/2) * (m - 1)^2) / (total_zeros_in_real - - (1/2) * (m - 1)^2)
  }
  else{
    TPR <- correctly_predicted_zeros / total_zeros_in_real
  }
  
  # Compute FPR
  false_positives <- (!real_zero_positions) & est_zero_positions
  total_non_zeros_in_real <- sum(!real_zero_positions)  
  false_positive_count <- sum(false_positives)
  FPR <- false_positive_count / total_non_zeros_in_real
  
  # summarize the result
  result <- c(n, p, m, sparsity, MSE, TPR, FPR, time_to_run)
  return(result)
}

find_optimal_permutation <- function(M, N) {
  
  # M: The real matrix to compare to
  # N: The matrix need to be permutated
  
  
  p <- nrow(M)
  q <- ncol(M)
  
  # Generate all permutations of column indices
  perm <- gtools::permutations(q, q)
  
  min_mse <- Inf
  optimal_N <- N
  
  # Iterate over all permutations of columns of N
  for (i in 1:nrow(perm)) {
    permuted_N <- N[, perm[i, ]] # Permute columns of N
    
    # Calculate Mean Squared Error
    mse <- sum((M - permuted_N)^2) / (p * q)
    
    # Update minimum MSE and optimal N if necessary
    if (mse < min_mse) {
      min_mse <- mse
      optimal_N <- permuted_N
    }
  }
  
  return(list(min_mse = min_mse, optimal_N = optimal_N))
}

generate_sample <- function(n, p, real_lambda, real_psi){
  Y <- MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = real_lambda %*% t(real_lambda) + real_psi) 
  return(Y)
}

########################################################################

N <- c(50,100,200,400,1000)

# loading 1

real_lambda <- matrix(c(0.95,0,0.9,0,0.85,0,0,0.8,0,0.75,0,0.7), nrow = 6, ncol = 2, byrow = TRUE)
real_psi <- diag(diag(diag(rep(1,6)) - real_lambda %*% t(real_lambda)))

# loading 2
#p <- 12
#m <- 4
#real_lambda <- matrix(c(
#  0.8, 0, 0, 0,
#  0.8, 0, 0, 0,
#  0.8, 0, 0, 0,
#  0, 0.7, 0.7, 0,
#  0, 0.7, 0.7, 0,
#  0, 0.7, 0.7, 0,
#  0.6, 0.6, 0, 0,
#  0.6, 0.6, 0, 0,
#  0.6, 0.6, 0, 0,
#  0, 0, 0.5, 0.5,
#  0, 0, 0.5, 0.5,
#  0, 0, 0.5, 0.5
#), nrow = 12, ncol = 4, byrow = TRUE)
#real_psi <- diag(c(0.4, 0.4, 0.4, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.2, 0.3))

p <- nrow(real_lambda)
m <- col(real_lambda)

samples <- generate_sample(max(N),p, real_lambda, real_psi)

# Create an empty data frame with the specified column names
result.dataframe <- data.frame(
  n = numeric(),
  p = numeric(),
  m = numeric(),
  sparsity = numeric(),
  MSE = numeric(),
  TPR = numeric(),
  FPR = numeric(),
  timetorun = numeric()
)

# Display the empty data frame
print(df)


for (n in N){
  simulation_times <- 50
  M <- matrix(0, nrow = simulation_times, ncol = 8) # a matrix to store the simulation result and to calculate the mean
  for (i in 1:50){
    rows_to_extract <- sample(1:nrow(samples), n)
    Y <- samples[rows_to_extract, ,drop = FALSE]
    tictoc::tic()
    two_step <- factanal(factors = 2, covmat = cor(Y), rotation = 'varimax')
    est_lambda <- two_step$loadings
    est_psi <- two_step$uniquenesses
    time_to_run <- tictoc::toc()
    M[i,] <- display_result(est_lambda, est_psi, real_lambda, real_psi, n, upper_triangle = FALSE, time_to_run[[2]] - time_to_run[[1]])
  }
  aver_result <- colMeans(M)
  result.dataframe <- rbind(result.dataframe, data.frame(
    n = aver_result[1],
    p = aver_result[2],
    m = aver_result[3],
    sparsity = aver_result[4],
    MSE = aver_result[5],
    TPR = aver_result[6],
    FPR = aver_result[7],
    timetorun = aver_result[8]
  ))
}
#saveRDS(result.dataframe, "C://Users//zhini//desktop//study material//A. Research Project//Master_Project//Outputs//thesis//Simulation code//result_twostep_loading1.rds")
