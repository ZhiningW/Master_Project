rm(list = ls())
########## Functions to be used in the simulation ################

A <- function(lambda, psi) {
  # Input: lambda: loading matrix of size p*k we have now
  #        psi: the variance matrix of size p*p of common factor we have now
  # Output: matrix A of size k*p (See definition in section 3.1)
  inv_mat <- solve(lambda %*% t(lambda) + psi)
  result <- t(lambda) %*% inv_mat
  return(result)
}

B <- function(lambda, psi) {
  # Input: lambda: loading matrix of size p*k we have now
  #        psi: the variance matrix of size p*p of common factor we have now
  # Output: B (See definition in section 3.1)
  k <- ncol(lambda)
  inv_mat <- solve(lambda %*% t(lambda) + psi)
  result <- diag(rep(1, k)) - t(lambda) %*% inv_mat %*% lambda
  return(result)
}

tr <- function(M){
  # Calculate the trace of matrix M
  return(sum(diag(M)))
}

E_step <- function(Y,lambda,psi){
  # Input: Y: the sample matrix of size n*p
  #        lambda: loading matrix of size p*k in current iteration
  #        psi: the common variance matrix of size p*p in current iteration
  # Output: The expectation in currect iteration (See equation (3.1) in section 3.1)
  n <- nrow(Y)
  mat_A <- A(lambda, psi)
  mat_B <- B(lambda, psi)
  result <- - n/2 * log(det(psi)) - 1/2 * tr(Y %*% solve(psi) %*% t(Y)) 
  + tr(lambda %*% mat_A %*% t(Y) %*% Y %*% solve(psi)) 
  - 1/2 * tr(lambda %*% ( n * mat_B + mat_A %*% t(Y) 
                          %*% Y %*% t(mat_A)) %*% t(lambda) %*% solve(psi))
  return(result)
}

psi_update <- function(lambda,psi,Y,j){
  # Input: lambda: loading matrix of size p*k we have now
  #        psi: the variance matrix of size p*p of common factor we have now
  #        Y: the response matrix of size n*p
  #        j: the element in (j,j) entry of matrix psi to update
  # Output: Updated psi[j,j]
  
  n <- nrow(Y)
  mat_A <- A(lambda,psi)
  mat_B <- B(lambda, psi)
  lambda_j <- lambda[j, ,drop = FALSE]
  Y.j <-Y[ ,j,drop = FALSE] # The j-th column of Y
  
  # Start calculating
  first_term <- (1/n) * t(Y.j) %*% Y.j
  second_term <- (-2/n) * lambda_j %*% mat_A %*% t(Y) %*% Y.j
  third_term <- (1/n) * lambda_j %*% (n * mat_B + mat_A %*% t(Y) %*% Y %*% t(mat_A)) %*% t(lambda_j)
  return(first_term + second_term + third_term)
}

psi_valid <- function(lambda,psi,Y){
  # Input: lambda: loading matrix of size p*k we have now
  #        psi: the variance matrix of size p*p of common factor we have now
  #        Y: the response matrix of size n*p
  # Output: True if psi is good to update (i.e in the concave region), False otherwise
  
  n <- nrow(Y)
  p <- ncol(Y)
  mat_A <- A(lambda,psi)
  mat_B <- B(lambda, psi)
  constraint <- numeric(length=p)
  for (j in 1:p){
    lambda_j <- lambda[j,,drop=FALSE]
    
    # Start calculating
    first_term <- (2/n) * sum(Y[,j]*Y[,j])
    second_term <- -(4/n) * lambda_j %*% as.matrix(rowSums(tcrossprod(mat_A,Y) 
                                                           %*% Y[,j, drop=FALSE]),ncol=1)
    third_term <- (2/n) * lambda_j %*% (n * mat_B + tcrossprod(mat_A,Y) 
                                        %*% tcrossprod(Y , mat_A)) %*% t(lambda_j)
    constraint[j]<- first_term + second_term + third_term
  }
  return(diag(psi)<=constraint)
}



lambda_update <- function(Y, lambda, psi, q){
  # Input: lambda: loading matrix of size p*k we have now
  #        psi: the variance matrix of size p*p of common factor we have now
  #        Y: the response matrix of size n*p
  #        q: the objective row of lambda to update
  # Output: Updated lambda[q, ]
  n <- nrow(Y)
  mat_A <- A(lambda,psi)
  mat_B <- B(lambda,psi)
  Y.q <- Y[ , q, drop = FALSE]
  result <- t(solve(n * mat_B + mat_A %*% t(Y) %*% Y %*% t(mat_A)) %*% mat_A %*% t(Y) %*% Y.q)
  return(result)
}

Standard_EM_update <- function(Y, m, initial_lambda, initial_psi){
  #### Input:
  ###    Y: nxp sample matrix
  ###    m: the order of the factor model
  ###    initial_loading: the initialization of the loading matrix used for the algorithm
  ###    initial_psi: the initialization of the psi-matrix used for the algorithm
  #### Output:
  ###    expectation: the expectation during updating
  ###    loading_result: the loading matrix after updating
  ###    psi_result: the psi-matrix after updating
  
  p <- ncol(Y)
  
  epsilon <- 0.05 # error tolerance
  n_max_iter <- 1000 # maximum iteration number
  
  
  step <- 0 # record the number of iterations
  loading_diff <- 1000 # Set a big difference in case smaller than tolerance at very beginning
  
  # Variables to store result
  expectation <- numeric(length = n_max_iter) # to store the expectation during updating
  
  ## Initialization
  lambda <- initial_lambda
  psi <- initial_psi
  
  while(loading_diff >= epsilon){
    step <- step + 1
    expectation[step] <- E_step(Y,lambda,psi)
    
    ## Update psi elementwisely
    psi_update_position <- which(psi_valid(lambda,psi,Y))
    for (j in psi_update_position){
      psi[j,j] <- psi_update(lambda,psi,Y,j)
    }
    #print(diag(psi))
    
    
    ## Update loading matrix using current psi rowwisely
    lambda_new <- matrix(0, nrow = p, ncol = m)
    for (j in 1:p){
      lambda_new[j,] <- lambda_update(Y, lambda, psi,j)
    }
    lambda_new[upper.tri(lambda_new)] <- 0
    #print(lambda_new)
    loading_diff <- norm(lambda_new - lambda, type = 'F')
    
    ## Update parameters
    lambda <- lambda_new
    
    if (step > n_max_iter){
      cli::cli_alert_warning('Failed to converge!')
      break
    }
  }
  plot(expectation[which(expectation != 0)])
  iter_to_converge <- length(expectation[which(expectation != 0)])
  result <- list(expectation[which(expectation != 0)], lambda, psi, iter_to_converge)
  return(result)
}

generate_sample <- function(n, p, real_lambda, real_psi){
  Y <- MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = real_lambda %*% t(real_lambda) + real_psi) 
  return(Y)
}

display_result <- function(est_lambda, est_psi, real_lambda, real_psi, n, upper_triangle = TRUE, time_to_run){
  p <- nrow(real_lambda)
  m <- ncol(real_lambda)
  if (p < m){
    print('Bad model')
    return(0)
  }
  # Compute the sparsity of the est_lambda
  if (upper_triangle){
    sparsity <- (sum(est_lambda == 0))  / (p * m)
  }
  
  # Compute MSE
  MSE <- (norm(est_lambda - real_lambda, type = "F") + norm(est_psi - real_psi, type = "F"))/(p + sum(real_lambda != 0))
  
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

################################################################################
par(mfrow = c(2,2))
set.seed(123)
N <- c(100, 200, 400, 1000, 2000, 5000)
P <- c(20,40,60,80)
m_collec <- c(4,6,8)


# Create an empty data frame with the specified column names
result.dataframe <- data.frame(
  n = numeric(),
  p = numeric(),
  m = numeric(),
  sparsity = numeric(),
  MSE = numeric(),
  TPR = numeric(),
  FPR = numeric(),
  timetorun = numeric(),
  True_sparsity = numeric()
)

sp <- 0.5
simulation_times <- 30


for (p in P){
  for (m in m_collec){
    initial_lambda <- matrix(rep(1,p*m),nrow = p, ncol = m)
    initial_psi <- diag(rep(1,p))
    matrix_generate <- loading3(p,m,sp)
    real_lambda <- matrix_generate[[1]]
    real_psi <- matrix_generate[[2]]
    samples <- generate_sample(max(N),p, real_lambda, real_psi)
    
    for (n in N){
        M <- matrix(0, nrow = simulation_times, ncol = 8) # a matrix to store the simulation result and to calculate the mean
      for (i in 1:simulation_times){
        rows_to_extract <- sample(1:nrow(samples), n)
        Y <- samples[rows_to_extract, ,drop = FALSE]
        SEM_result <- Standard_EM_update(Y, m, initial_lambda, initial_psi)
        est_lambda <- SEM_result[[2]]
        est_lambda[which(abs(est_lambda)<=0.05)] <- 0
        est_psi <- SEM_result[[3]]
        M[i,] <- display_result(est_lambda, est_psi, real_lambda, real_psi, n, upper_triangle = TRUE, time_to_run[[2]] - time_to_run[[1]])
      }
      aver_result <- colMeans(M)
      result.dataframe <- rbind(result.dataframe, data.frame(
        n = n,
        p = p,
        m = m,
        sparsity = aver_result[4],
        MSE = aver_result[5],
        TPR = aver_result[6],
        FPR = aver_result[7],
        timetorun = aver_result[8],
        True_sparsity = sp
      ))
    }
  }

}

saveRDS(result.dataframe, "C://Users//zhini//desktop//study material//A. Research Project//Master_Project//Outputs//thesis//Simulation code//result_standardEM_loading3_0.5sp.rds")
#readRDS("C://Users//zhini//desktop//study material//A. Research Project//Master_Project//Outputs//thesis//Simulation code//result_standardEM_loading2.rds")
