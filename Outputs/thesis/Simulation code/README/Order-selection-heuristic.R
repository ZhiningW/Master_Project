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



###########################################################################
library(ggplot2)
library(reshape2)
library(patchwork)
N <- c(50,100,200,400,600,1000)
p <- nrow(real_lambda)
m <- ncol(real_lambda)
result.dataframe <- data.frame(
  n = numeric(),
  "2" = numeric(),
  "3" = numeric(),
  "4" = numeric(),
  "5" = numeric(),
  "6" = numeric()
)
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
  result.dataframe <- rbind(result.dataframe, data.frame(
    n = n,
    "2" = TS_collection[1],
    "3" = TS_collection[2],
    "4" = TS_collection[3],
    "5" = TS_collection[4],
    "6" = TS_collection[5]
  ))
}
df_long <- melt(result.dataframe, id.vars = "n", variable.name = "order", value.name = "statistics")
plot2 <- ggplot(df_long, aes(x = order, y = statistics, color = as.factor(n), group = n)) +
  geom_line() +
  geom_point() +
  labs(title = "Setting 1b", 
       x = "Order", 
       y = "Heuristic Statistics", 
       color = "n") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_hline(yintercept = 0.8,
             color = "red",,
             linetype = "dashed") +
  colorspace::scale_color_discrete_qualitative() +
  scale_x_discrete(labels = c(2, 3, 4, 5, 6))
plot1
plot <- plot1 + plot2 + plot_layout(guides = "collect")
plot
ggsave("result_plot/OSresult_heu.pdf", height = 4)  
