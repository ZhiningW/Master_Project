##### This code is written to simulation the EM algorithm in (penalized) FA models.
##### We will assume k is known in advanced.

###### try profvis

library(MASS)


###############################################################################

### Functions may be used in the simulation
rm(list = ls())
################### Functions may be used in the simulation ####################
## To calculate A in article
A <- function(lambda, psi) {
  # Input: lambda: loading matrix of size p*k we have now
  #        psi: the variance matrix of size p*p of common factor we have now
  # Output: A 
  inv_mat <- ginv(lambda %*% t(lambda) + psi)
  result <- t(lambda) %*% inv_mat
  return(result)
}

B <- function(lambda, psi) {
  # Input: lambda: loading matrix of size p*k we have now
  #        psi: the variance matrix of size p*p of common factor we have now
  # Output: B
  k <- ncol(lambda)
  inv_mat <- ginv(lambda %*% t(lambda) + psi)
  result <- diag(rep(1, k)) - t(lambda) %*% inv_mat %*% lambda
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
  lambda_j <- lambda[j,,drop=FALSE]
  
  # Start calculating
  first_term <- (1/n) * sum(Y[,j]*Y[,j])
  second_term <- -(2/n) * lambda_j %*% as.matrix(rowSums(tcrossprod(mat_A,Y) 
                                                         %*% Y[,j, drop=FALSE]),ncol=1)
  third_term <- (1/n) * lambda_j %*% (n * mat_B + tcrossprod(mat_A,Y) 
                                      %*% tcrossprod(Y , mat_A)) %*% t(lambda_j)
  return(first_term + second_term + third_term)
}


pem_E <- function(lambda,psi,Y,rho){
  Out <- npem_E(lambda,psi,Y)-1/2*rho*sum(abs(lambda))
  return(Out)
}



npem_E <- function(lambda,psi,Y){
  # Input: lambda: loading matrix of size p*k we have now
  #        psi: the variance matrix of size p*p of common factor we have now
  #        Y: the response matrix of size n*p
  # Output: non-penalized expectation (ignore the constant)
  
  n <- nrow(Y)
  p <- nrow(lambda)
  
  # Precompute A and B
  mat_A <- A(lambda, psi)
  mat_B <- B(lambda, psi)
  
  # first term of expectation
  first_term <- -(n/2) * sum(log(diag(psi)))
  
  # second term of expectation
  second_term <- 0
  for (i in 1:n){
    for (j in 1:p){
      yij <- Y[i,j]
      lambda_j <- t(as.matrix(lambda[j,]))
      Y_i <- t(as.matrix(Y[i,]))
      second_term <- second_term 
      + (yij^2 - 2 * yij * lambda_j %*% mat_A %*% t(Y_i))/psi[j,j]
    }
    
  }
  second_term <- -1/2 * second_term
  
  # third term of expectation
  third_term <- 0
  for (i in 1:n){
    for (j in 1:p){
      Y_i <- t(as.matrix(Y[i,]))
      lambda_j <- t(as.matrix(lambda[j,]))
      third_term <- third_term 
      + (lambda_j %*% (mat_B + mat_A %*% t(Y_i) %*% Y_i %*% t(mat_A)) %*% t(lambda_j))/psi[j,j]
    }
  }
  third_term <- -1/2 * third_term
  
  return(first_term + second_term + third_term)
}


loading_update <- function(lambda,psi,Y,j){
  # Input: lambda: loading matrix of size p*k we have now
  #        psi: the variance matrix of size p*p of common factor we have now
  #        Y: the response matrix of size n*p
  #        j: the j-th row of lambda to be updated
  # Output: Updated lambda[j,]
  
  n <- nrow(Y)
  mat_A <- A(lambda,psi)
  mat_B <- B(lambda, psi)
  k <- ncol(lambda)
  
  
  # first term of updating formula
  factor1 <- matrix(0, nrow = k, ncol = 1)
  for (i in 1:n){
    Y_i <- as.matrix(Y[i,])
    factor1 <- factor1 - 2 * Y[i,j] * mat_A %*% Y_i 
  }
  
  
  # second term of updating formula
  factor2 <- matrix(0, nrow = k, ncol = k)
  for (i in 1:n){
    Y_i <- as.matrix(Y[i,])
    factor2 <- factor2 + 2 * mat_B + 2 * mat_A %*% Y_i %*% t(Y_i) %*% t(mat_A)
  }
  
  
  return(- t(factor1) %*% t(solve(factor2)))
}
  
###############################################################################

set.seed(124)

n <- 5000  # Size of observations
p <- 6     # Dimension of each observation
k <- 2    # Order of FA model (dimension of latent variables)


real_loading <- matrix(c(0.95,0,0.9,0,0.85,0,0,0.8,0,0.75,0,0.7), nrow=p, ncol=k, byrow = TRUE)
psi_diag <- diag(diag(p)-tcrossprod(real_loading))
real_psi <- diag(psi_diag)

Y <- mvrnorm(n = n, mu = rep(0,p), Sigma = real_loading %*% t(real_loading) + real_psi)
# Real parameters of the FA model
# Generate the loading matrix
#sparsity_level <- 0.3  # The proportion of zero elements in the loading matrix
#non_zeros <- round(k * p * (1 - sparsity_level))
#poss_values_1 <- c(0.5, 0.6, 0.7, 0.8, 0.9, -0.5, -0.6, -0.7, -0.8, -0.9)
#loading_position <- sample(k * p, non_zeros)
#real_loading <- matrix(0, nrow = p, ncol = k)
#real_loading[loading_position] <- sample(poss_values_1, non_zeros, replace = TRUE)

# Ensure no all-zero rows or columns
#while (any(rowSums(real_loading) == 0) || any(colSums(real_loading) == 0)) {
#  loading_position <- sample(k * p, non_zeros)
#  real_loading <- matrix(0, nrow = p, ncol = k)
#  real_loading[loading_position] <- sample(poss_values_1, non_zeros, replace = TRUE)
#}

print(real_loading)

# Generate the real variance of the unique factors epsilon
poss_values_2 <- c(0.1, 0.2, 0.3, 0.4)
diagonal_values <- sample(poss_values_2, p, replace = TRUE)
real_psi <- diag(diagonal_values)

print(real_psi)

# Generate the factor matrix F, response matrix Y, and unique factor epsilon
F <- matrix(rnorm(n * k, mean = 0, sd = 1), nrow = n, ncol = k)
epsilon <- mvrnorm(n = n, mu = rep(0, p), Sigma = real_psi)
Y <- F %*% t(real_loading) + epsilon


###############################################################################

### Perform a FA using factanal() function in R to the response matrix Y (traditional
### two steps method).
fa_2steps_result <- factanal(x = Y, factors = k, scores = 'regression', rotation = 'varimax')
loading_2steps <- as.matrix(fa_2steps_result$loadings)
error_2steps <- norm(loading_2steps - real_loading,type = 'F')
print(loading_2steps)

###############################################################################
### Perform a non-penalized EM algorithm 

## Set the initial guess
npem_initial_loading <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1), 
                               nrow=p, ncol=k, byrow = TRUE)   # Take distinct initial values to see what happens
npem_initial_psi <- diag(rep(0.1,p))     # Take distinct initial values to see what happens

## End the iteration if the error between two steps is less than np_tolerance
## and regard as convergent
npem_tolerance <- 0.01 
n_max_iter <- 30 # Set a maximum iteration in case an infinity loop.

## Set the parameters for iteration
npem_step <- 0 # record the number of iterations
npem_loading_diff <- 1000 # Set a big difference in case smaller than tolerance at very beginning
npem_psi_diff <- 1000
npem_expectation <- numeric(length = n_max_iter)# Record expectations during iteration

# update iteratively
npem_loading_old <- npem_initial_loading
npem_psi_old <- npem_initial_psi

while(npem_loading_diff >= npem_tolerance && npem_psi_diff >= npem_tolerance){
  npem_step <- npem_step + 1
  npem_expectation[npem_step] <- npem_E(npem_loading_old,npem_psi_old,Y)
  ## Update psi elementwisely
  npem_psi_new <- matrix(0, nrow = p, ncol = p)
  for (j in 1:p){
    npem_psi_new[j,j] <- psi_update(npem_loading_old,npem_psi_old,Y,j)
  }
  npem_psi_diff <- norm(npem_psi_new - npem_psi_old, type = 'F') # calculate the difference between psi matrix
  
  ## Update loading matrix using current psi rowwisely
  
  npem_loading_new <- matrix(0, nrow = p, ncol = k)
  for (j in 1:p){
    npem_loading_new[j,] <- loading_update(npem_loading_old, npem_psi_new, Y, j)
  }
  npem_loading_diff <- norm(npem_loading_new - npem_loading_old, type = 'F')
  
  ## Update parameters
  npem_psi_old <- npem_psi_new
  npem_loading_old <- npem_loading_new
  
  if (npem_step > n_max_iter){
    cat('Failed to converge!')
    break
  }
  
}
## final result of non-penalized EM algorithm
npem_loading <- npem_loading_new
npem_psi <- npem_psi_new
plot(npem_expectation)

  


  





