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

E_step <- function(Y,lambda,psi,rho){
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
            - rho * (sum(lambda))
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

gradient <- function(Y, lambda, lambda.q, psi, q){
  # lambda.q: q-th row of the lambda
  # To calculate the gradient of q-th row of lambda
  n <- nrow(Y)
  mat_A <- A(lambda,psi)
  mat_B <- B(lambda, psi)
  coe <- 2/psi[q,q]
  result <- coe * ((n * mat_B + mat_A %*% t(Y) %*% Y %*% t(mat_A)) %*% t(lambda.q) 
                   - mat_A %*% t(Y) %*% Y[ , q, drop = FALSE])
  return(result)
}

proximal_method <- function(Y, lambda, psi, q, rho = 1){
  s <- 0.01 # the step size
  epsilon <- 0.01 # tolerance
  updated_lambda_q <- lambda[q, , drop = FALSE]
  error <- 10^5
  N_max <- 10^3
  i <- 1
  while (error > epsilon && i < N_max){
    grad <- gradient(Y, lambda, updated_lambda_q, psi, q)
    xi <- updated_lambda_q - s * t(grad)
    new_lambda.q <- sign(xi) * pmax(abs(xi) - rho * s, 0)
    error <- norm(new_lambda.q - updated_lambda_q, type = "2")
    updated_lambda_q <- new_lambda.q
    i <- i + 1
  }
  return <- updated_lambda_q
}

penalized_EM_algorithm <- function(N, real_lambda, real_psi, initial_lambda, initial_psi, rho){
  
  p <- nrow(real_lambda)
  m <- ncol(real_lambda)
  Y <- MASS::mvrnorm(n = N, mu = rep(0,p), Sigma = real_lambda %*% t(real_lambda) + real_psi) 
  
  epsilon <- 0.00001 # error tolerance
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
    expectation[step] <- E_step(Y, lambda, psi, rho)
    
    ## Update psi elementwisely
    psi_update_position <- which(psi_valid(lambda,psi,Y))
    for (j in psi_update_position){
      psi[j,j] <- psi_update(lambda,psi,Y,j)
    }
    #print(diag(psi))
    
    ## Update lambda rowwisely
    lambda_new <- matrix(0, nrow = p, ncol = m)
    for (q in 1:p){
      lambda_new[q,] <- proximal_method(Y, lambda, psi, q, rho)
    }
    loading_diff <- norm(lambda_new - lambda, type = 'F')
  }
  result <- list(expectation[which(expectation != 0)], lambda, psi)
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
#################################################################################
n <- 200
p <- 6
m <- 2
rho <- 1
# true_lambda <- matrix(rep(1,12), nrow = p, ncol = k, byrow = TRUE)
real_lambda <- matrix(c(0.95,0,0.9,0,0.85,0,0,0.8,0,0.75,0,0.7), nrow = p, ncol = m, byrow = TRUE)
# true_lambda <- matrix(c(0.9,0,0.9,0,0.8,0,0,0.8,0,0.7,0,0.7), nrow= p, ncol =m, byrow = TRUE)

real_psi <- diag(diag(diag(rep(1,p)) - real_lambda %*% t(real_lambda)))

initial_lambda <- matrix(rep(1,p*m),nrow = p, ncol = m)
initial_psi <- diag(rep(1,p))
penalized_EM_algorithm(n, real_lambda, real_psi, initial_lambda, initial_psi, rho)