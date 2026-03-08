
sample_data(physeq_obj)$Day <- as.factor(sample_data(physeq_obj)$Day)

alpha_plot_time <- plot_richness(physeq_obj, x="Day", measures=c("Shannon", "Simpson"), color="Day") +
  geom_boxplot(alpha=0.6, outlier.shape = NA) + # Adds the boxplot layer
  geom_jitter(width=0.2, size=2, alpha=0.8) +   # Overlays individual mouse data points
  theme_bw() + 
  theme(text = element_text(size=14), legend.position="none") + 
  labs(title = "Longitudinal Stability of Gut Alpha Diversity",
       x = "Clinical Timepoint (Day)",
       y = "Alpha Diversity Measure")

ggsave(filename = "figures/alpha_diversity_longitudinal.png", 
       plot = alpha_plot_time, 
       width = 8, 
       height = 6, 
       dpi = 300)

alpha_plot_time