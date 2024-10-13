generate_sample <- function(n, p, real_lambda, real_psi){
  Y <- MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = real_lambda %*% t(real_lambda) + real_psi) 
  return(Y)
}

tr <- function(M){
  # Calculate the trace of matrix M
  return(sum(diag(M)))
}

heuris <- function(lambda, psi){
  result <- tr(lambda %*% t(lambda)) / tr(lambda %*% t(lambda) + psi)
}
################################################################################
set.seed(123)
N <- c(50,100,200,400,600,1000)

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
real_psi <- diag(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2))
sample <- generate_sample(max(N), p, real_lambda, real_psi)
par(mfrow = c(3,2))

for (n in N){
  Y <- sample[sample(1:max(N), n, replace = TRUE), , drop = FALSE]
  TS_collection <- numeric(p/2 - 1)
  for (k in 2:p/2){
    mle_result <- factanal(factors = k, covmat = cor(Y), rotation = 'varimax')
    est_lambda <- mle_result$loadings
    print(est_lambda)
    est_psi <- mle_result$uniquenesses
    print(est_psi)
    TS <- heuris(est_lambda, est_psi)
    TS_collection[k-1] <- TS
  }
  plot(c(2,3,4,5,6), TS_collection, main = paste("Order Selection Procedure When n = ",n, "m = ", m), 
       xlab = "Potential order of the model", ylab = "Heuristic statistics" )
  abline(h = 0.8, col = "red", lty = 2)
}