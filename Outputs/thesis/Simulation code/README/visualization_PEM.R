library("tidyverse")
DF1 <- readRDS("Outputs//thesis//Simulation code//result_PEM_loading3_0.2sp.rds")
DF2 <- readRDS("Outputs//thesis//Simulation code//result_PEM_loading3_0.3sp.rds")
DF3 <- readRDS("Outputs//thesis//Simulation code//result_PEM_loading3_0.5sp.rds")
DF <- bind_rows(DF1,DF2,DF3) |> 
  as_tibble()


DF |> ggplot(aes(n, MSE, color = as.factor(true_sparsity))) +
  geom_line(aes(group = interaction(p, m, true_sparsity)))  +
  facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  labs(color = "True-sparsity") +
  colorspace::scale_color_discrete_qualitative()

DF |> ggplot(aes(n, sparsity, color = as.factor(true_sparsity))) +
  geom_line(aes(group = interaction(p, m, true_sparsity)))  +
  facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  labs(color = "True-sparsity") +
  colorspace::scale_color_discrete_qualitative()

DF |> ggplot(aes(true_sparsity, sparsity, color = as.factor(n ))) +
  geom_line(aes(group = interaction(p, m, n)))  +
  facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(color = "n") +
  colorspace::scale_color_discrete_qualitative()+
  labs( y = "Sparsity recovered", x = "True sparsity")

DF |> ggplot(aes(n, 2 * (TPR/(TPR+FPR) * TPR) /  (TPR/(TPR+FPR) + TPR), color = as.factor(true_sparsity))) +
  geom_line(aes(group = interaction(p, m, true_sparsity)))  +
  facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  labs( y = "F1 score") +
  labs(color = "True-sparsity") +
  colorspace::scale_color_discrete_qualitative()

DF |> ggplot(aes(n, TPR, color = as.factor(true_sparsity))) +
  geom_line(aes(group = interaction(p, m, true_sparsity)))  +
  facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  labs(color = "True-sparsity")+
  colorspace::scale_color_discrete_qualitative()

DF |> ggplot(aes(n, FPR, color = as.factor(true_sparsity))) +
  geom_line(aes(group = interaction(p, m, true_sparsity)))  +
  facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  labs(color = "True-sparsity")+
  colorspace::scale_color_discrete_qualitative()
