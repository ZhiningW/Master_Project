#' The A matrix
#'
#' @param lambda loading matrix of size p*k we have now
#' @param psi the variance matrix of size p*p of common factor we have now
#' @name matrices
A <- function(lambda, psi) {
  inv_mat <- ginv(lambda %*% t(lambda) + psi)
  result <- t(lambda) %*% inv_mat
  return(result)
}

#' @rdname matrices
B <- function(lambda, psi) {
  k <- ncol(lambda)
  inv_mat <- ginv(lambda %*% t(lambda) + psi)
  result <- diag(rep(1, k)) - t(lambda) %*% inv_mat %*% lambda
  return(result)
}

#' E step
#'
#' @inheritParams matrices
#' @param Y the response matrix of size n*p
#' @return non-penalized expectation (ignore the constant)
npem_E <- function(lambda,psi,Y){

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
      + (yij^2 + 2 * yij * lambda_j %*% mat_A %*% t(Y_i))/psi[j,j]
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

#' @inheritParams npem_E
#' @param j  the element in (j,j) entry of matrix psi to update
#' @return Updated psi[j,j]
psi_update <- function(lambda,psi,Y,j){

  n <- nrow(Y)
  mat_A <- A(lambda,psi)
  mat_B <- B(lambda, psi)

  # Start calculating
  first_term <- (1/n) * sum(Y[,j]*Y[,j])
  second_term <- 0
  lambda_j <- t(as.matrix(lambda[j,]))
  for (i in 1:n){
    Y_i <- t(as.matrix(Y[i,]))
    second_term <- second_term + 2 * Y[i,j] * lambda_j %*% mat_A %*% t(Y_i)
  }
  second_term <- 1/n * second_term
  third_term <- 0
  for(i in 1:n){
    third_term <- third_term + lambda_j %*% (mat_B + mat_A %*% t(Y_i) %*% Y_i %*% t(mat_A)) %*% t(lambda_j)
  }
  third_term <- 1/n * third_term
  return(first_term + second_term + third_term)
}

#' @inheritParams npem_E
#' @param j the j-th row of lambda to be updated
#' @return Updated lambda[j,]
loading_update <- function(lambda,psi,Y,j){
  n <- nrow(Y)
  mat_A <- A(lambda,psi)
  mat_B <- B(lambda, psi)
  k <- ncol(lambda)


  # first term of updating formula
  factor1 <- matrix(0, nrow = k, ncol = 1)
  for (i in 1:n){
    Y_i <- as.matrix(Y[i,])
    factor1 <- factor1 + 2 * Y[i,j] * mat_A %*% Y_i
  }


  # second term of updating formula
  factor2 <- matrix(0, nrow = k, ncol = k)
  for (i in 1:n){
    Y_i <- as.matrix(Y[i,])
    factor2 <- factor2 + 2 * mat_B + 2 * mat_A %*% Y_i %*% t(Y_i) %*% t(mat_A)
  }


  return(- t(factor1) %*% t(solve(factor2)))
}


pem_E <- function(lambda,psi,Y,rho){
  Out <- npem_E(lambda,psi,Y)-1/2*rho*sum(abs(lambda))
  return(Out)
}

subgradient <- function(lambda,psi,Y,epsilon,j,rho){
  ### to find the subgradient
  mat_A <- A(lambda,psi)
  mat_B <- B(lambda,psi)
  lambda_j <- lambda[j,]
  first_term <- 0
  for (i in 1:n){
    first_term <- first_term + 2 * Y[i,j] * tcrossprod(mat_A ,Y[i, ,drop = FALSE])
  }
  first_term <- first_term/psi[j,j]
  second_term <- 0
  for (i in 1:n){
    # use tcrossprod
    second_term <- second_term+ (2 * mat_B + 2 * mat_A %*% t(Y[i, , drop=FALSE]) %*% Y[i, , drop=FALSE] %*% t(mat_A)) %*% t(lambda_j)
  }
  second_term <- second_term/psi[j,j]
  third_term <- rho * sign(t(lambda_j))
  return(first_term+second_term+third_term)
}



#' @inheritParams npem_E
#' @param j the j-th row we want to update
#' @param rho the penalization parameter
#' @return the j-th row after updating
subg_method <- function(lambda, psi, Y, epsilon, j, rho){
  n <- nrow(Y)
  mat_A <- A(lambda,psi)
  mat_B <- B(lambda,psi)
  lambda_j_old <- lambda[j, ,drop= FALSE]
  error <- 10000
  iteration <- 0
  while(error>epsilon){
    subg <- subgradient(lambda,psi,Y,epsilon,j,rho)
    t <- 1/ ((iteration + 1)*norm(subg,type='2'))# step size
    lambda_j_new <- lambda_j_old - t * t(subg)
    iteration <- iteration + 1
    error <- norm(lambda_j_new-lambda_j_old,type='2')
    lambda_j_old <- lambda_j_new
  }
  return(lambda_j_old)
}
