set.seed(123)
n <- 10000 # Number of observations
true_loading <- matrix(c(0.8,0,0.6,0,0,0,0.8,0,0.7,0.6),nrow=5,ncol=2)
true_psi <- diag(c(0.1,0.2,0.2,0.1,0.1))

samples <- MASS::mvrnorm(n = n, mu = rep(0, 5), 
                         Sigma = true_loading %*% t(true_loading) + true_psi)

# the function factanal has already used varimax method to get the optimal loading
fa_result <- factanal(factors = 2, covmat = cor(samples))
fa_result
