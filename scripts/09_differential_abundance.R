# ---------------------------------------------------------
# STEP 13: DIFFERENTIAL ABUNDANCE (DESeq2)
# ---------------------------------------------------------
# Load required libraries
library(phyloseq)
library(DESeq2)
library(ggplot2)

# 1. Clean the data: Remove any samples where the 'Day' metadata is missing (NA)
# This strictly prevents the DESeq2 "cannot contain NA" error
physeq_clean <- subset_samples(physeq_obj, !is.na(Day))

# 2. Create a clear binary condition on the CLEANED object
# We group Days 0-9 as "Early" and everything else as "Late"
sample_data(physeq_clean)$Time_Phase <- ifelse(as.numeric(as.character(sample_data(physeq_clean)$Day)) < 10, "Early", "Late")
sample_data(physeq_clean)$Time_Phase <- factor(sample_data(physeq_clean)$Time_Phase, levels = c("Early", "Late"))

# 3. Convert clean phyloseq object to DESeq2 format and run the core algorithm
diagdds <- phyloseq_to_deseq2(physeq_clean, ~ Time_Phase)
diagdds <- DESeq(diagdds, test="Wald", fitType="parametric")

# 4. Extract the statistical results and filter for strict significance
res <- results(diagdds, cooksCutoff = FALSE)
sigtab <- res[which(res$padj < 0.01), ] # Keep only highly significant shifts

# 5. Bind the biological names back to the math results so we know WHO shifted
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(physeq_clean)[rownames(sigtab), ], "matrix"))

# 6. Plot the Log2 Fold Change of significant bacteria
deseq_plot <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  geom_point(size=4, alpha=0.8) + 
  geom_hline(yintercept=0, linetype="dashed", color="black") + # Adds a baseline at zero
  theme_bw() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust=1, size=12),
        text = element_text(size=14)) +
  ggtitle("Statistically Significant Shifts: Late vs. Early Timepoints")

# 7. Save and display
ggsave("figures/deseq2_log2fc.png", deseq_plot, width=12, height=8, dpi=300)
deseq_plot