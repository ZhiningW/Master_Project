set.seed(123)
library(GPArotation)
library(tidyverse)
# True loading matrix
lambda_true <- matrix(c(0.8,0,0.6,0,-0.5,0,0.8,0,0.7,0.6),
                      nrow=5,ncol=2) 
# True psi matrix
psi_true <- diag(c(0.1,0.2,0.2,0.1,0.1))
i <- 1
n_vec <- c(50,100,200,400,800,1600,3200,6400,12800,20000,40000,100000, 500000)
samples_all <- MASS::mvrnorm(n = max(n_vec), mu = rep(0, 5), 
                         Sigma = lambda_true %*% t(lambda_true) + psi_true)
psi_estimated <- matrix(0, nrow = length(n_vec), ncol = 5)
for (n in n_vec){
#n <- 100 # Number of observations
# Generate samples from the multivariate normal distribution
samples <- samples_all[1:n, ]

# Traditional two-step procedure
# MLE: 
mle_result <- factanal(factors = 2, 
                       covmat = cor(samples), rotation = 'none')
print(mle_result)
psi_estimated[i,] <- mle_result$uniquenesses

i <- i + 1
#loading_unrotated <- mle_result$loadings # MLE of loading matrix
#print(loading_unrotated)

# Rotation technique: Varimax
#rotated_result <- Varimax(loading_unrotated)
#loading_rotated <- rotated_result$loadings # rotation result
#loading_rotated[abs(loading_rotated) < 0.05] = 0 # omit too small loadings
#print(loading_rotated)
}
psi_estimated
as.data.frame(psi_estimated)
as.data.frame(psi_estimated) |> mutate(index = 1:n()) |> pivot_longer(-index)
psi_dataframe <- as.data.frame(psi_estimated) |> mutate(index = 1:n(), n = n_vec) |> pivot_longer( c(-index, -n))
ggplot(psi_dataframe, aes(n, value, color = name)) + geom_line() + geom_hline(aes(yintercept = value), 
                                                                                  data = data.frame(value = diag(psi_true), name = paste0('V', 1:5)))+ 
                                                                                    facet_wrap(~name)

                                                                                  