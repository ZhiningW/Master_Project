test_that("model works", {
  et.seed(124)

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

  print(real_loading)

  # Generate the real variance of the unique factors epsilon
  poss_values_2 <- c(0.1, 0.2, 0.3, 0.4)
  diagonal_values <- sample(poss_values_2, p, replace = TRUE)
  real_psi <- diag(diagonal_values)

  print(real_psi)

  # Generate the factor matrix F, response matrix Y, and unique factor epsilon
  F <- matrix(rnorm(n * k, mean = 0, sd = 1), nrow = n, ncol = k)
  epsilon <- mvrnorm(n = n, mu = rep(0, p), Sigma = real_psi)
  Y <- F %*% t(real_loading) + epsilon


  ###############################################################################

  ### Perform a FA using factanal() function in R to the response matrix Y (traditional
  ### two steps method).
  fa_2steps_result <- factanal(x = Y, factors = k, scores = 'regression', rotation = 'varimax')
  loading_2steps <- as.matrix(fa_2steps_result$loadings)
  error_2steps <- norm(loading_2steps - real_loading,type = 'F')
  print(loading_2steps)

  ###############################################################################
  ### Perform a non-penalized EM algorithm

  ## Set the initial guess
  npem_initial_loading <- real_loading   # Take distinct initial values to see what happens
  npem_initial_psi <- diag(rep(0.1,p))     # Take distinct initial values to see what happens

  ## End the iteration if the error between two steps is less than np_tolerance
  ## and regard as convergent
  npem_tolerance <- 0.01
  n_max_iter <- 30 # Set a maximum iteration in case an infinity loop.

  ## Set the parameters for iteration
  npem_step <- 0 # record the number of iterations
  npem_loading_diff <- 1000 # Set a big difference in case smaller than tolerance at very beginning
  npem_psi_diff <- 1000
  npem_expectation <- numeric(length = n_max_iter)# Record expectations during iteration

  # update iteratively
  npem_loading_old <- npem_initial_loading
  npem_psi_old <- npem_initial_psi

  while(npem_loading_diff >= npem_tolerance && npem_psi_diff >= npem_tolerance){
    npem_step <- npem_step + 1
    npem_expectation[npem_step] <- npem_E(npem_loading_old,npem_psi_old,Y)
    ## Update psi elementwisely
    npem_psi_new <- matrix(0, nrow = p, ncol = p)
    for (j in 1:p){
      npem_psi_new[j,j] <- psi_update(npem_loading_old,npem_psi_old,Y,j)
    }
    npem_psi_diff <- norm(npem_psi_new - npem_psi_old, type = 'F') # calculate the difference between psi matrix

    ## Update loading matrix using current psi rowwisely

    npem_loading_new <- matrix(0, nrow = p, ncol = k)
    for (j in 1:p){
      npem_loading_new[j,] <- loading_update(npem_loading_old, npem_psi_new, Y, j)
    }
    npem_loading_diff <- norm(npem_loading_new - npem_loading_old, type = 'F')

    ## Update parameters
    npem_psi_old <- npem_psi_new
    npem_loading_old <- npem_loading_new

    if (npem_step > n_max_iter){
      cat('Failed to converge!')
      break
    }

  }
  ## final result of non-penalized EM algorithm
  npem_loading <- npem_loading_new
  npem_psi <- npem_psi_new
  plot(npem_expectation)




})
