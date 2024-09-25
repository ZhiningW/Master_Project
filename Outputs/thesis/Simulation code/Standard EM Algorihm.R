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

Standard_EM_update <- function(real_lambda, real_psi, N, initial_lambda, initial_psi){
  #### Input:
  ###    real_lambda: the real loading matrix
  ###    N: the sample size
  ###    initial_loading: the initialization of the loading matrix used for the algorithm
  ###    initial_psi: the initialization of the psi-matrix used for the algorithm
  #### Output:
  ###    expectation: the expectation during updating
  ###    loading_result: the loading matrix after updating
  ###    psi_result: the psi-matrix after updating
  
  p <- nrow(real_lambda)
  k <- ncol(real_lambda)
  Y <- MASS::mvrnorm(n = N, mu = rep(0,p), Sigma = real_lambda %*% t(real_lambda) + real_psi) 
  
  epsilon <- 0.0001 # error tolerance
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
    lambda_new <- matrix(0, nrow = p, ncol = k)
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
  result <- list(expectation[which(expectation != 0)], lambda, psi)
  return(result)
}
################################################################################

n <- 200
p <- 6
k <- 2

# true_lambda <- matrix(rep(1,12), nrow = p, ncol = k, byrow = TRUE)
true_lambda <- matrix(c(0.95,0,0.9,0,0.85,0,0,0.8,0,0.75,0,0.7), nrow = p, ncol = k, byrow = TRUE)
# true_lambda <- matrix(c(0.9,0,0.9,0,0.8,0,0,0.8,0,0.7,0,0.7), nrow= p, ncol = k, byrow = TRUE)

true_psi <- diag(diag(diag(rep(1,p)) - true_lambda %*% t(true_lambda)))

initial_lambda <- matrix(rep(1,p*k),nrow = p, ncol = k)
initial_psi <- diag(rep(1,p))


Standard_EM_update(true_lambda, true_psi, n, initial_lambda, initial_psi)

true_psi
true_lambda
