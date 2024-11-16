rm(list = ls())
set.seed(123)
###############################################################################
loading3 <- function(p, m, sparsity){
  repeat {
    # Initialize an empty matrix
    M <- matrix(0, nrow = p, ncol = m)
    
    # Populate the non-upper triangle part of the matrix with values from N(0, 1)
    for (i in 1:p) {
      for (j in 1:min(i, m)) {
        M[i, j] <- rnorm(1)
      }
    }
    # Randomly set elements to zero to achieve desired sparsity
    non_upper_indices <- which(lower.tri(M, diag = TRUE), arr.ind = TRUE)
    num_non_upper_triangle <- nrow(non_upper_indices)
    num_to_zero <- round(num_non_upper_triangle * sparsity) 
    zero_indices <- non_upper_indices[sample(num_non_upper_triangle, num_to_zero), ]
    
    for (k in 1:nrow(zero_indices)) {
      M[zero_indices[k, 1], zero_indices[k, 2]] <- 0
    }
    
    # Check if there are any zero rows or zero columns
    if (!any(rowSums(M) == 0) && !any(colSums(M) == 0)) {
      break
    }
    print("finish generate loading matrix")
  }
  psi <- diag(0.2* rep(1,p))
  return(list(M,psi))
}

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
  
  # MSE
  MSE <- after_perm$min_mse * p*m / (p + sum(real_lambda != 0)) 
  + norm(diag(real_psi) - est_psi, type = "2")/(p + sum(real_lambda != 0))
  
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
    TPR <- (correctly_predicted_zeros - (1/2) * m * (m-1)) / (total_zeros_in_real - (1/2) * m * (m-1))
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
  
  
  r <- nrow(M)
  q <- ncol(M)
  
  # Generate all permutations of column indices
  perm <- gtools::permutations(q, q)
  
  min_mse <- Inf
  optimal_N <- N
  
  # Iterate over all permutations of columns of N
  for (i in 1:nrow(perm)) {
    permuted_N <- N[, perm[i, ]] # Permute columns of N
    
    # Calculate Mean Squared Error
    mse <- sum((M - permuted_N)^2) / (r*q)
    
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
set.seed(123)
N <- c(100, 200, 400, 1000, 2000, 5000)
P <- c(20,40,60,80)
m_collec <- c(4,6,8)
result.dataframe <- data.frame(
  n = numeric(),
  p = numeric(),
  m = numeric(),
  sparsity = numeric(),
  MSE = numeric(),
  TPR = numeric(),
  FPR = numeric(),
  timetorun = numeric(),
  true_sparsity = numeric()
)
simulation_times <- 20
sp <- 0.2
for (p in P){
  for (m in m_collec){
    for (n in N){
      if (n < p){
        next
      }
      matrix_generate <- loading3(p,m,sp)
      real_lambda <- matrix_generate[[1]]
      real_psi <- matrix_generate[[2]]
      samples <- generate_sample(10000,p, real_lambda, real_psi)
      
      Info <- matrix(0, nrow = simulation_times, ncol = 8) # a matrix to store the simulation result and to calculate the mean
      for (i in 1:simulation_times){
        rows_to_extract <- sample(1:nrow(samples), n)
        Y <- samples[rows_to_extract, ,drop = FALSE]
        #if (det(cor(Y)) < 1e-13){
        #  next
        #}
        tictoc::tic()
        tryCatch({
          two_step <- factanal(Y, factors = m, rotation = 'varimax')
        }, error = function(e){
          two_step <- factanal(Y, factors = m, rotation = 'varimax', start = rep(1,p))
        })
        est_lambda <- two_step$loadings
        est_lambda[which(abs(est_lambda) <= 0.05)] = 0
        est_psi <- two_step$uniquenesses
        time_to_run <- tictoc::toc()
        Info[i,] <- display_result(est_lambda, est_psi, real_lambda, real_psi, n, upper_triangle = FALSE, time_to_run[[2]] - time_to_run[[1]])
      }
      Info_nonzerorow <- Info[apply(Info, 1, function(row) any(row != 0)), , drop = FALSE]
      aver_result <- colMeans(Info_nonzerorow)
      result.dataframe <- rbind(result.dataframe, data.frame(
        n = n,
        p = p,
        m = m,
        sparsity = aver_result[4],
        MSE = aver_result[5],
        TPR = aver_result[6],
        FPR = aver_result[7],
        timetorun = aver_result[8],
        true_sparsity = sp
      ))
    }
  }
}

saveRDS(result.dataframe, "C://Users//zhini//desktop//study material//A. Research Project//Master_Project//Outputs//thesis//Simulation code//result_twostep_loading3_0.2sp.rds")
#readRDS("C://Users//zhini//desktop//study material//A. Research Project//Master_Project//Outputs//thesis//Simulation code//result_twostep_loading3.rds")
