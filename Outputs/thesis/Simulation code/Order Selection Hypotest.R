generate_sample <- function(n, p, real_lambda, real_psi){
  Y <- MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = real_lambda %*% t(real_lambda) + real_psi) 
  return(Y)
}

tr <- function(M){
  # Calculate the trace of matrix M
  return(sum(diag(M)))
}

hypo_TS <- function(Y,lambda,psi){
  ### Input: Y: Sample matrix
  ###       lambda: the estimated loading matrix
  ###       psi: the estimated psi matrix
  ### Output: 
  ###       Test statistics
  n <- nrow(Y)
  p <- ncol(Y)
  m <- ncol(lambda)
  S <- cov(Y)
  print(psi)
  Sigma <- lambda %*% t(lambda) + diag(psi)
  print(det(Sigma))
  #print(Sigma)
  TS <- n * (log(det(Sigma)) + tr(solve(Sigma) %*% S) - log(det(S)) - p )
  s <- 1/2 * (p - m)^2 - 1/2 * (p + m)
  p_value <- 1 - pchisq(TS,s)
  return(p_value)
}

################################################################################
set.seed(123)

n <- 1000
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
real_psi <- diag(rep(0.15, 12))

p <- nrow(real_lambda)
m <- ncol(real_lambda)

Y <- generate_sample(n, p, real_lambda, real_psi)
TS_collection <- numeric(p/2 - 1)
for (k in 2:p/2){
  mle_result <- factanal(factors = k, covmat = cor(Y), rotation = 'varimax')
  est_lambda <- mle_result$loadings
  print(est_lambda)
  est_psi <- mle_result$uniquenesses
  #print(est_psi)
  TS <- hypo_TS(Y, est_lambda, est_psi)
  TS_collection[k-1] <- TS
}
plot(c(2,3,4,5,6), TS_collection)
s <- 1/2 * (p - m)^2 - 1/2 * (p + m)
abline(h = qchisq(0.95, s))
abline(h = 0.05, col = "red", lty = 2)
