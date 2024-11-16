library("tidyverse")
library("patchwork")
theme_set(theme_bw(base_size = 14))
DF1 <- readRDS("Outputs//thesis//Simulation code//result_twostep_loading2.rds")
DF1$method <- "Two-step"
DF2 <- readRDS("Outputs//thesis//Simulation code//result_standardEM_loading2.rds")
DF2$method <- "Standard EM"
DF3 <- readRDS("Outputs//thesis//Simulation code//result_PEM_loading2.rds")
DF3$method <- "Penalized EM"
DF <- bind_rows(DF1,DF2,DF3) |> 
  as_tibble()

DF |> ggplot(aes(n, MSE, color = as.factor(method))) +
  geom_line(aes(group = interaction(p, m, method)))  +
  #facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  ggtitle("MSE v.s. sample size (Model 2)") +
  labs(color = "Method") +
  colorspace::scale_color_discrete_qualitative()

ggsave("result_plot/MSE-compar.pdf")


DF |> ggplot(aes(n, sparsity , color = as.factor(method))) +
  geom_line(aes(group = interaction(p, m, method)))  +
  #facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  ggtitle("Sparsity Recovered v.s. sample size (Model 2)") +
  labs(color = "Method") +
  colorspace::scale_color_discrete_qualitative()


DF |> ggplot(aes(n, TPR , color = as.factor(method))) +
  geom_line(aes(group = interaction(p, m, method)))  +
  #facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  ggtitle("TPR v.s. n by Different Methods (Model 2)")

DF |> ggplot(aes(n, FPR , color = as.factor(method))) +
  geom_line(aes(group = interaction(p, m, method)))  +
  #facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  ggtitle("FPR v.s. n by Different Methods (Model2)")
