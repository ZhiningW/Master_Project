library("tidyverse")
library("patchwork")
DF1 <- readRDS("Outputs//thesis//Simulation code//result_twostep_loading1.rds")
DF1$method <- "Two-step"
DF2 <- readRDS("Outputs//thesis//Simulation code//result_standardEM_loading1.rds")
DF2$method <- "Standard EM"
DF3 <- readRDS("Outputs//thesis//Simulation code//result_PEM_loading1.rds")
DF3$method <- "Penalized EM"
DF_model1 <- bind_rows(DF1,DF2,DF3) |> 
  as_tibble()
DF4 <- readRDS("Outputs//thesis//Simulation code//result_twostep_loading2.rds")
DF4$method <- "Two-step"
DF5 <- readRDS("Outputs//thesis//Simulation code//result_standardEM_loading2.rds")
DF5$method <- "Standard EM"
DF6 <- readRDS("Outputs//thesis//Simulation code//result_PEM_loading2.rds")
DF6$method <- "Penalized EM"
DF_model2 <- bind_rows(DF4,DF5,DF6) |> 
  as_tibble()


plotMSE1 <- DF_model1 |> ggplot(aes(n, MSE, color = as.factor(method))) +
  geom_line(aes(group = interaction(p, m, method)))  +
  #facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  ggtitle("Setting 2a") +
  labs(color = "Method") +
  colorspace::scale_color_discrete_qualitative() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)) 

plotMSE2 <- DF_model2 |> ggplot(aes(n, MSE, color = as.factor(method))) +
  geom_line(aes(group = interaction(p, m, method)))  +
  #facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  ggtitle("Setting 2b") +
  labs(color = "Method") +
  colorspace::scale_color_discrete_qualitative() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

MSE_plot <- plotMSE1 + plotMSE2 + plot_layout(guides = "collect")
MSE_plot
ggsave("result_plot/MSE-compar.pdf", height = 4)

#####################

plotspar1 <- DF_model1 |> ggplot(aes(n, sparsity , color = as.factor(method))) +
  geom_line(aes(group = interaction(p, m, method)))  +
  #facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  ggtitle("Setting 2a") +
  labs(color = "Method") +
  colorspace::scale_color_discrete_qualitative() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

plotspar2 <- DF_model2 |> ggplot(aes(n, sparsity , color = as.factor(method))) +
  geom_line(aes(group = interaction(p, m, method)))  +
  #facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  ggtitle("Setting 2b")+
  labs(color = "Method") +
  colorspace::scale_color_discrete_qualitative() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

spar_plot <- plotspar1 + plotspar2 + plot_layout(guides = "collect")
spar_plot
ggsave("result_plot/sparsity-compar.pdf", height = 4)

######################

plotTPR1 <- DF_model1 |> ggplot(aes(n, TPR , color = as.factor(method))) +
  geom_line(aes(group = interaction(p, m, method)))  +
  #facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  ggtitle("Setting 2a") + 
  labs(color = "Method") +
  colorspace::scale_color_discrete_qualitative() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

plotTPR2 <- DF_model2 |> ggplot(aes(n, TPR , color = as.factor(method))) +
  geom_line(aes(group = interaction(p, m, method)))  +
  #facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  ggtitle("Setting 2b") +
  labs(color = "Method") +
  colorspace::scale_color_discrete_qualitative() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
TPR_plot <- plotTPR1 + plotTPR2 + plot_layout(guides = "collect")
TPR_plot
ggsave("result_plot/TPR-compar.pdf", height = 4)

############################
plotFPR1 <- DF_model1 |> ggplot(aes(n, FPR , color = as.factor(method))) +
  geom_line(aes(group = interaction(p, m, method)))  +
  #facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  ggtitle("Setting 2a") +
  ylim(c(0,0.3)) +
  labs(color = "Method") +
  colorspace::scale_color_discrete_qualitative() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

plotFPR2 <- DF_model2 |> ggplot(aes(n, FPR , color = as.factor(method))) +
  geom_line(aes(group = interaction(p, m, method)))  +
  #facet_grid(p~m, labeller = label_both)  +
  geom_point() +
  ggtitle("Setting 2b") +
  ylim(c(0,0.3)) +
  labs(color = "Method") +
  colorspace::scale_color_discrete_qualitative() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

FPR_plot <- plotFPR1 + plotFPR2 + plot_layout(guides = "collect")
FPR_plot
ggsave("result_plot/FPR-compar.pdf", height = 4)
