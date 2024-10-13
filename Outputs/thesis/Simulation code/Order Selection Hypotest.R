generate_sample <- function(n, p, real_lambda, real_psi){
  Y <- MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = real_lambda %*% t(real_lambda) + real_psi) 
  return(Y)
}

hypo_TS <- function(Y,lambda,psi){
  ### Input: Y: Sample matrix
  ###       lambda: the estimated loading matrix
  ###       psi: the estimated psi matrix
  ###       test_level: the level to make first type error, 0.05 is widely used
  ### Output: 
  ###       Test statistics
  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(lambda)
  S <- cov(Y)
  print(dim(S))
  Sigma <- lambda %*% t(lambda) + psi
  print(Sigma)
  TS <- sum(diag(solve(Sigma) %*% S)) - log(det(solve(Sigma) %*% S)) - p
  return(TS)
}

################################################################################
set.seed(123)

n <- 200
p <- 12
m <- 4
real_lambda <- matrix(c(
  0.8, 0, 0, 0,
  0.8, 0, 0, 0,
  0.8, 0, 0, 0,
  0, 0.7, 0.7, 0,
  0, 0.7, 0.7, 0,
  0, 0.7, 0.7, 0,
  0.6, 0.6, 0, 0,
  0.6, 0.6, 0, 0,
  0.6, 0.6, 0, 0,
  0, 0, 0.5, 0.5,
  0, 0, 0.5, 0.5,
  0, 0, 0.5, 0.5
), nrow = 12, ncol = 4, byrow = TRUE)
s <- 1/2 * (p - m)^2 - 1/2 * (p + m)
real_psi <- diag(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2))
Y <- generate_sample(n, p, real_lambda, real_psi)
TS_collection <- numeric(p/2 - 1)
for (k in 2:p/2){
  mle_result <- factanal(factors = k, covmat = cor(Y), rotation = 'varimax')
  est_lambda <- mle_result$loadings
  print(est_lambda)
  est_psi <- mle_result$uniquenesses
  print(est_psi)
  TS <- hypo_TS(Y, est_lambda, est_psi)
  TS_collection[k-1] <- TS
}
plot(c(2,3,4,5,6), TS_collection)
abline(h = qchisq(0.95, s), col = "red", lty = 2)
