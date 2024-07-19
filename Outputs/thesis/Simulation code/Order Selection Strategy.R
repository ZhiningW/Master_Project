hypo_test_order_selection <- function(Y,lambda,psi,test_level){
  ### Input: Y: Sample matrix
  ###       lambda: the estimated loading matrix
  ###       psi: the estimated psi matrix
  ###       test_level: the level to make first type error, 0.05 is widely used
  ### Output: 
  ###       TRUE if we do not reject the null hypothesis, i.e. k is the best order 
  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(lambda)
  s <- 1/2 * (p-k)^2 - (p+k) # degree of freedom the chi-square distribution has
  cov_Y <- cov(Y)
  Sigma <- lambda %*% t(lambda) + psi
  TS <- sum(diag(solve(Sigma) %*% cov_Y)) - log(det(solve(Sigma) %*% cov_Y)) - p
  chi <- qchisq(1-test_level, s)
  
  return(chi > TS)
}


set.seed(123)
Y <- matrix(rnorm(100), nrow = 10, ncol = 10)
lambda <- matrix(rnorm(30), nrow = 10, ncol = 3)
psi <- diag(nrow(lambda))-lambda %*% t(lambda)

result <- hypo_test_order_selection(Y, lambda, psi, 0.05)
print(result)