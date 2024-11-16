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
N <- c(50)

for (n in N) {
  Y <- generate_sample(n, p, real_lambda, real_psi)
  rho <- seq(0,0.1,length.out = 200)
  AIC_tuning <- numeric(length(rho))
  BIC_tuning <- numeric(length(rho))
  Sparsity <- numeric(length(rho))
  index <- 1
  for(i in rho){
    PFA <- fanc::fanc(Y, m, rho = i, gamma = Inf)
    out <- fanc::out(PFA, rho = i, gamma = Inf)
    AIC_tuning[index] <- out$criteria[1]
    BIC_tuning[index] <- out$criteria[2]
    est_lambda <- out$loadings
    est_lambda[abs(est_lambda) <= 0.05] <- 0
    Sparsity[index] <- (m * p - sum(est_lambda != 0))/(m*p)
    index <- index + 1
    print(index)
  }
  print(AIC_tuning)
  print(BIC_tuning)
  par(mfrow = c(2,1))
  plot(rho,AIC_tuning, ylab = "AIC")
  plot(rho,BIC_tuning, ylab = "BIC")
  
  cat("the minimum AIC suggests rho = ", rho[which(AIC_tuning == min(AIC_tuning))], "\n")
  cat("the minimum BIC suggests rho = ", rho[which(BIC_tuning == min(BIC_tuning))], "\n")
}
data <- data.frame(
  rho = rep(rho, 2),  # Repeat rho twice, once for AIC and once for BIC
  value = c(AIC_tuning, BIC_tuning),  # Combine AIC and BIC into a single vector
  metric = rep(c("AIC", "BIC"), each = length(rho))  # A new column to label AIC and BIC
)

# Plot using ggplot
plot2 <- ggplot(data, aes(x = rho, y = value, color = metric)) +
  geom_line() +
  geom_point() +
  labs(title = "Setting 2b",
       x = "rho", 
       y = "AIC or BIC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  colorspace::scale_color_discrete_qualitative()
plot2
ggsave("result_plot/PTresult_50.pdf", height = 4) 
