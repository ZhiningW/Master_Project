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

loading1 <- function(){
  real_lambda <- matrix(c(0.95,0,0.9,0,0.85,0,0,0.8,0,0.75,0,0.7), nrow = 6, ncol = 2, byrow = TRUE)
  real_psi <- diag(diag(diag(rep(1,6)) - real_lambda %*% t(real_lambda)))
  result <- list(real_lambda,real_psi)
  return(result)
}
################################################################################
set.seed(123)
library("reshape2")
library("tidyverse")
library("patchwork")
N <- c(50,100,200,400,600,1000)
result.dataframe <- data.frame(
  n = numeric(),
  "2" = numeric(),
  "3" = numeric(),
  "4" = numeric(),
  "5" = numeric(),
  "6" = numeric()
)
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
real_psi <- diag(rep(0.3, 12))

p <- nrow(real_lambda)
m <- ncol(real_lambda)
for (n in N){
  
  Y <- generate_sample(n, p, real_lambda, real_psi)
  p_value <- numeric(p/2 - 1)
  for (k in 2:p/2){
    mle_result <- factanal(Y, factors = k,  rotation = 'varimax')
    est_lambda <- mle_result$loadings
    print(est_lambda)
    est_psi <- mle_result$uniquenesses
    #print(est_psi)
    TS <- mle_result$PVAL
    p_value[k-1] <- TS
  }
  result.dataframe <- rbind(result.dataframe, data.frame(
    n = n,
    "2" = p_value[1],
    "3" = p_value[2],
    "4" = p_value[3],
    "5" = p_value[4],
    "6" = p_value[5]
  ))
}
df_long <- melt(result.dataframe, id.vars = "n", variable.name = "order", value.name = "statistics")
plot2 <- ggplot(df_long, aes(x = order, y = statistics, color = as.factor(n), group = n)) +
  geom_line() +
  geom_point() +
  labs(title = "Setting 1b",
       x = "Order", 
       y = "p-value", 
       color = "n") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.05,
             color = "red",,
             linetype = "dashed") +
  colorspace::scale_color_discrete_qualitative() +
  scale_x_discrete(labels = c(2, 3, 4, 5, 6))
plot1
plot <- plot1 + plot2 + plot_layout(guides = "collect")
plot
ggsave("result_plot/OSresult_hypotest.pdf", height = 4) 
