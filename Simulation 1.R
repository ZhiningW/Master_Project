set.seed(123)
n <- 1000000 # Number of observations
p <- 5   # Number of variables

# Simulate data: Factors with explicitly sparse structure
factor1 <- rnorm(n, 0, 1)
factor2 <- rnorm(n, 0, 1)

# Generate observed variables with noise
# Ensure sparse structure: Only a few variables load strongly on each factor
X <- matrix(NA, nrow=n, ncol=p)
X[,1] <- 0.8*factor1 + rnorm(n, 0, 1)  # Strong loading on factor 1
X[,2] <- 0.7*factor2 + rnorm(n, 0, 1)  # Strong loading on factor 2
X[,3] <- 0.8*factor1 + rnorm(n, 0, 1)  # Noise, no strong loading on any factor
X[,4] <- 0.5*factor1 - 0.6*factor2 + rnorm(n, 0, 1)  
X[,5] <- 0.6*factor2 + rnorm(n, 0, 1) 

fa_result <- factanal(factors = 2, covmat = cor(X))

rotated_fa <- varimax(fa_result$loadings)

print(rotated_fa$loadings)

