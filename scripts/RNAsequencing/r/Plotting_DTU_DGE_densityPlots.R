library(ggplot2)
library(cowplot)


upVdown_subset

upVdown_dtu_subset




dge_density <- ggplot(data = upVdown_subset) + 
  geom_density(aes(x = estimate)) +
  theme_classic() +
  xlab("Estimate") +
  xlab(expression(log[2](CPM))) +
  ylab("Density") +
  geom_vline(xintercept = c(-1,1), linetype="dotted", 
             color = "red")
  
dge_density


dtu_density <- ggplot(data = upVdown_dtu_subset) + 
  geom_density(aes(x = estimate)) +
  theme_classic() +
  xlab("Estimate") +
  xlab(expression(log[2](Counts))) +
  ylab("Density") + 
  geom_vline(xintercept = c(-1,1), linetype="dotted", 
                               color = "red")



plot_grid(dge_density, dtu_density,
          nrow = 1, ncol = 2,
          labels = c("A)", "B)"),
          label_size = 11)



