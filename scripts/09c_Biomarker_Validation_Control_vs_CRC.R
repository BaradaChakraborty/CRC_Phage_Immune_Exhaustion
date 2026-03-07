# STEP 09c: BIOMARKER VALIDATION (CONTROL VS CRC)
library(curatedMetagenomicData)
library(ggplot2)
library(ggpubr)
zeller_tse <- curatedMetagenomicData("ZellerG_2014.relative_abundance", dryrun = FALSE)[[1]]
# Extract F. nucleatum abundance and Study Condition
fuso_idx <- grep("s__Fusobacterium_nucleatum", rownames(zeller_tse))
fuso_abundance <- as.numeric(assay(zeller_tse)[fuso_idx, ])
comp_data <- data.frame(
  Condition = factor(zeller_tse$study_condition, levels = c("control", "CRC")),
  F_nucleatum = fuso_abundance
)
biomarker_plot <- ggplot(comp_data, aes(x = Condition, y = F_nucleatum, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_bw() +
  scale_fill_manual(values = c("control" = "#4575b4", "CRC" = "#d73027")) +
  # Add the Wilcoxon test to prove statistical significance
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5) +
  labs(title = "Fusobacterium nucleatum Enrichment in CRC",
       subtitle = "Validation across Zeller et al. (2014) clinical cohort",
       x = "Clinical Status",
       y = "Relative Abundance (%)")

# Save the evidence
ggsave("figures/fuso_biomarker_validation.png", biomarker_plot, width=6, height=6)
print("Biomarker validation plot complete.")