rm(list = ls())
set.seed(123)
###############################################################################
loading1 <- function(){
  real_lambda <- matrix(c(0.95,0,0.9,0,0.85,0,0,0.8,0,0.75,0,0.7), nrow = 6, ncol = 2, byrow = TRUE)
  real_psi <- diag(diag(diag(rep(1,6)) - real_lambda %*% t(real_lambda)))
  result <- list(real_lambda,real_psi)
  return(result)
}

loading2 <- function(){
  real_lambda <- matrix(c(
    0.8, 0, 0, 0,
    0.8, 0, 0, 0,
    0.8, 0, 0, 0,
    0, 0.7, 0.7, 0,
    0, 0.7, 0.7, 0,
    0, 0.7, 0.7, 0,
    0.6, 0.6, 0, 0.7,
    0.6, 0.6, 0, 0.7,
    0.6, 0.6, 0, 0.7,
    0, 0, 0.5, 0.7,
    0, 0, 0.5, 0.7,
    0, 0, 0.5, 0.7
  ), nrow = 12, ncol = 4, byrow = TRUE)
  real_psi <- diag(c(0.4, 0.4, 0.4, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.2, 0.3))
  result <- list(real_lambda,real_psi)
  return(result)
}

generate_sample <- function(n, p, real_lambda, real_psi){
  Y <- MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = real_lambda %*% t(real_lambda) + real_psi) 
  return(Y)
}

#################################################################################

real_lambda <- loading2()[[1]]
real_psi <- loading2()[[2]]


p <- nrow(real_lambda)
m <- ncol(real_lambda)
N <- c(5000)

for (n in N) {
  Y <- generate_sample(n, p, real_lambda, real_psi)
  rho <- seq(0,0.1,length.out = 200)
  AIC_tuning <- numeric(length(rho))
  BIC_tuning <- numeric(length(rho))
  index <- 1
  for(i in rho){
    PFA <- fanc::fanc(Y, m, rho = i, gamma = Inf)
    out <- fanc::out(PFA, rho = i, gamma = Inf)
    AIC_tuning[index] <- out$criteria[1]
    BIC_tuning[index] <- out$criteria[2]
    index <- index + 1
    print(index)
  }
  print(AIC_tuning)
  print(BIC_tuning)
  par(mfrow = c(2,1))
  plot(rho,AIC_tuning, ylab = "AIC", main = paste(" Parameter Tuning Procedure by AIC (Model1, n = ", n,")"))
  plot(rho,BIC_tuning, ylab = "BIC", main = paste(" Parameter Tuning Procedure by BIC (Model1, n = ", n,")"))
  cat("the minimum AIC suggests rho = ", rho[which(AIC_tuning == min(AIC_tuning))], "\n")
  cat("the minimum BIC suggests rho = ", rho[which(BIC_tuning == min(BIC_tuning))], "\n")
}

