rm(list = ls())
### Functions may be used in the simulation
library(MASS)
set.seed(123)
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

psi_valid <- function(lambda,psi,Y){
  # Input: lambda: loading matrix of size p*k we have now
  #        psi: the variance matrix of size p*p of common factor we have now
  #        Y: the response matrix of size n*p
  # Output: True if psi is good to iterate (i.e in the concave region), False otherwise
  
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


subgradient <- function(lambda,psi,Y,epsilon,j,rho){
  ### to find the subgradient
  n <- nrow(Y)
  mat_A <- A(lambda,psi)
  mat_B <- B(lambda,psi)
  lambda_j <- lambda[j, , drop = FALSE]
  first_term <- - as.matrix(rowSums(tcrossprod(mat_A,Y)%*%Y[,j, drop = FALSE]*2/psi[j,j]),ncol=1) 
  
  second_term <- 2 * ( n * mat_B + tcrossprod(mat_A,Y) %*% tcrossprod(Y,mat_A)) %*% t(lambda[j, ,drop= FALSE ]) /psi[j,j]
  third_term <- rho * sign(t(lambda_j)) 
  return(first_term+second_term+third_term)
}

 


subg_method <- function(lambda,psi,Y,epsilon,j,rho){
  ### Input:  lambda: loading matrix
  ###         psi: the variance matrix of size p*p of common factor we have now
  ###         Y: response 
  ###         epsilon: the torlance for the subgradient method regarded as convergent
  ###         j: the j-th row we want to update
  ###         rho: the penalization parameter
  ### Output: the j-th row after updating
  n <- nrow(Y)
  mat_A <- A(lambda,psi)
  mat_B <- B(lambda,psi)
  lambda_j_old <- lambda[j, ,drop= FALSE]
  error <- 10000
  iteration <- 0
  
  while(error>epsilon && iteration < 30){
    subg <- subgradient(lambda,psi,Y,epsilon,j,rho)
    t <- 1/ ((iteration + 1)*norm(subg,type='2'))# step size
    lambda_j_new <- lambda_j_old - t * t(subg)
    iteration <- iteration + 1
    error <- norm(lambda_j_new-lambda_j_old,type='2')
    lambda_j_old <- lambda_j_new
    iteration <- iteration + 1
  }
  return(lambda_j_old)
}


#####################################################
### Background settings
n <- 5000  # Size of observations
p <- 8     # Dimension of each observation
k <- 4     # Order of FA model (dimension of latent variables)

# Real parameters of the FA model
# Generate the loading matrix
sparsity_level <- 0.3  # The proportion of zero elements in the loading matrix
non_zeros <- round(k * p * (1 - sparsity_level))
poss_values_1 <- c(0.5, 0.6, 0.7, 0.8, 0.9, -0.5, -0.6, -0.7, -0.8, -0.9)
loading_position <- sample(k * p, non_zeros)
real_loading <- matrix(0, nrow = p, ncol = k)
real_loading[loading_position] <- sample(poss_values_1, non_zeros, replace = TRUE)

# Ensure no all-zero rows or columns
while (any(rowSums(real_loading) == 0) || any(colSums(real_loading) == 0)) {
  loading_position <- sample(k * p, non_zeros)
  real_loading <- matrix(0, nrow = p, ncol = k)
  real_loading[loading_position] <- sample(poss_values_1, non_zeros, replace = TRUE)
}


# Generate the real variance of the unique factors epsilon
poss_values_2 <- c(0.1, 0.2, 0.3, 0.4)
diagonal_values <- sample(poss_values_2, p, replace = TRUE)
real_psi <- diag(diagonal_values)


# Generate the factor matrix F, response matrix Y, and unique factor epsilon
F <- matrix(rnorm(n * k, mean = 0, sd = 1), nrow = n, ncol = k)
epsilon <- mvrnorm(n = n, mu = rep(0, p), Sigma = real_psi)
Y <- F %*% t(real_loading) + epsilon # generating Y directly??
#####################################################

rho <- 2 # penality parameter   #AIC / BIC 
## Set the initial guess
pem_initial_loading <- matrix(1, nrow=p, ncol=k)   # Take distinct initial values to see what happens
pem_initial_psi <- diag(rep(0.1,p))     # Take distinct initial values to see what happens

## End the iteration if the error between two steps is less than np_tolerance
## and regard as convergent
pem_tolerance <- 0.01 
n_max_iter <- 30
## Set the parameters for iteration
pem_step <- 0 # record the number of iterations
pem_loading_diff <- 1000 # Set a big difference in case smaller than tolerance at very beginning
pem_psi_diff <- 1000
pem_expectation <- numeric(length = n_max_iter) # Record expectations during iteration

# update iteratively
pem_loading_old <- pem_initial_loading
pem_psi_old <- pem_initial_psi

while(pem_loading_diff >= pem_tolerance){
  pem_step <- pem_step + 1
  pem_expectation[pem_step] <- pem_E(pem_loading_old,pem_psi_old,Y,rho)
  ## Update psi elementwisely
  # psi_goodtoupdate <- psi_valid(pem_loading_old,pem_psi_old,Y)
  # psi_update_position <- which(psi_goodtoupdate)
  pem_psi_new <- matrix(0, nrow = p, ncol = p)
  for (j in 1:p){
    pem_psi_new[j,j] <- psi_update(pem_loading_old,pem_psi_old,Y,j)
  }
  print(diag(pem_psi_new))
  ## Update loading matrix using current psi rowwisely
  
  pem_loading_new <- matrix(0, nrow = p, ncol = k)
  for (j in 1:p){
    pem_loading_new[j,] <- subg_method(pem_loading_old,pem_psi_new,Y,pem_tolerance,j,rho)
  }
  print(pem_loading_new)
  pem_loading_diff <- norm(pem_loading_new - pem_loading_old, type = 'F')
  
  ## Update parameters
  pem_loading_old <- pem_loading_new
  pem_psi_old <- pem_psi_new
  
  if (pem_step > n_max_iter){
    cli::cli_alert_warning('Failed to converge!')
    break
  }
  
}
## final result of penalized EM algorithm
pem_loading <- pem_loading_new
pem_psi <- pem_psi_new
plot(pem_expectation[-1])  
