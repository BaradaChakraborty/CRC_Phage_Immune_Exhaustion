library(phyloseq); library(ggplot2)

physeq_rel <- transform_sample_counts(physeq_obj, function(x) x / sum(x))

physeq_phylum <- tax_glom(physeq_rel, taxrank = "Phylum")

abundance_plot <- plot_bar(physeq_phylum, x="Day", fill="Phylum") +
  geom_bar(stat="identity", color="black", size=0.1) + # Adds clean black outlines to the segments
  theme_bw() +
  theme(text = element_text(size=14),
        legend.position = "right") +
  labs(title = "Microbiome Composition by Phylum Over Time",
       x = "Clinical Timepoint (Day)",
       y = "Relative Abundance (100%)")

ggsave(filename = "figures/relative_abundance_phylum.png", 
       plot = abundance_plot, 
       width = 10, 
       height = 6, 
       dpi = 300)

abundance_plot