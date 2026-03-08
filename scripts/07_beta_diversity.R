library(phyloseq); library(ggplot2)

physeq_rel <- transform_sample_counts(physeq_obj, function(x) x / sum(x) )

ord <- ordinate(physeq_rel, method="PCoA", distance="bray")

beta_plot <- plot_ordination(physeq_rel, ord, color="Day") +
  geom_point(size=4, alpha=0.8) + # Creates large, slightly transparent dots
  theme_bw() +
  theme(text = element_text(size=14)) +
  ggtitle("Beta Diversity: Gut Microbiome Compositional Shift (PCoA)")

ggsave(filename = "figures/beta_diversity_pcoa.png", 
       plot = beta_plot, 
       width = 8, 
       height = 6, 
       dpi = 300)

beta_plot