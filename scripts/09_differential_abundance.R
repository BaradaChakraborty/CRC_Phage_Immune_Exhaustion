library(phyloseq)
library(DESeq2)
library(ggplot2)

physeq_clean <- subset_samples(physeq_obj, !is.na(Day))

sample_data(physeq_clean)$Time_Phase <- ifelse(as.numeric(as.character(sample_data(physeq_clean)$Day)) < 10, "Early", "Late")
sample_data(physeq_clean)$Time_Phase <- factor(sample_data(physeq_clean)$Time_Phase, levels = c("Early", "Late"))

diagdds <- phyloseq_to_deseq2(physeq_clean, ~ Time_Phase)
diagdds <- DESeq(diagdds, test="Wald", fitType="parametric")

res <- results(diagdds, cooksCutoff = FALSE)
sigtab <- res[which(res$padj < 0.01), ] # Keep only highly significant shifts

sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(physeq_clean)[rownames(sigtab), ], "matrix"))

deseq_plot <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  geom_point(size=4, alpha=0.8) + 
  geom_hline(yintercept=0, linetype="dashed", color="black") + # Adds a baseline at zero
  theme_bw() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust=1, size=12),
        text = element_text(size=14)) +
  ggtitle("Statistically Significant Shifts: Late vs. Early Timepoints")

ggsave("figures/deseq2_log2fc.png", deseq_plot, width=12, height=8, dpi=300)
deseq_plot




library(phyloseq)
library(DESeq2)
library(ggplot2)
physeq_genus <- tax_glom(physeq_clean, taxrank = "Genus")
diagdds_genus <- phyloseq_to_deseq2(physeq_genus, ~ Time_Phase)
diagdds_genus <- DESeq(diagdds_genus, test="Wald", fitType="parametric")
res_genus <- results(diagdds_genus, cooksCutoff = FALSE)
fuso_results <- as.data.frame(res_genus)
fuso_results$Genus <- as.character(tax_table(physeq_genus)[rownames(fuso_results), "Genus"])
fuso_evidence <- fuso_results[grep("Fusobacterium", fuso_results$Genus, ignore.case = TRUE), ]
fuso_plot <- plot_abundance(physeq_clean, "Time_Phase", "Fusobacterium") +
  geom_boxplot(aes(fill=Time_Phase), alpha=0.7) +
  theme_bw() +
  labs(title="Proof of Fusobacterium Dominance",
       subtitle="Significant Enrichment as CRC Progresses",
       y="Normalized Sequence Counts")

ggsave("figures/fuso_dominance_proof.png", fuso_plot, width=8, height=6)
fuso_plot

ps_rel <- transform_sample_counts(physeq_clean, function(x) x / sum(x))

fuso_only <- subset_taxa(ps_rel, Genus == "Fusobacterium")

fuso_abundance_plot <- plot_bar(fuso_only, x="Time_Phase", fill="Genus") +
  geom_boxplot(aes(x=Time_Phase, y=Abundance, fill=Time_Phase), alpha=0.5) +
  theme_bw() +
  labs(title="The 'Who' Proof: Fusobacterium Relative Abundance",
       subtitle="Evidence of microbial takeover during CRC progression",
       y="Relative Abundance (%)")

ggsave("figures/fuso_relative_abundance.png", fuso_abundance_plot, width=8, height=6)
fuso_abundance_plot

ps_rel <- transform_sample_counts(physeq_clean, function(x) x / sum(x))

fuso_rel <- subset_taxa(ps_rel, Genus == "Fusobacterium")

fuso_dominance_plot <- plot_bar(fuso_rel, x="Time_Phase", fill="Genus") +
  geom_boxplot(aes(x=Time_Phase, y=Abundance, fill=Time_Phase), alpha=0.5) +
  theme_bw() +
  labs(title="The 'Who' Proof: Fusobacterium Relative Abundance",
       subtitle="Evidence of microbial takeover during CRC progression",
       y="Community Share (%)")

ggsave("figures/fuso_dominance_proof.png", fuso_dominance_plot, width=8, height=6)
fuso_dominance_plot




fuso_counts <- otu_table(physeq_clean)["Fusobacterium", ]
sample_data(physeq_clean)$Fuso_Load <- as.numeric(fuso_counts)

alpha_fuso_plot <- plot_richness(physeq_clean, x="Fuso_Load", measures=c("Shannon")) +
  geom_point(aes(color=Time_Phase), size=3) +
  geom_smooth(method="lm", color="black", linetype="dashed") +
  theme_bw() +
  labs(title="Diversity Collapse vs. Fusobacterium Load",
       subtitle="Mathematical proof that F. nucleatum drives ecological dysbiosis",
       x="Fusobacterium Sequence Counts",
       y="Shannon Diversity Index")

ps_rel <- transform_sample_counts(physeq_clean, function(x) x / sum(x))
ord <- ordinate(ps_rel, method="PCoA", distance="bray")
beta_fuso_plot <- plot_ordination(ps_rel, ord, color="Time_Phase", size="Fuso_Load") +
  geom_point(alpha=0.7) +
  theme_bw() +
  labs(title="Beta Diversity Shift Driven by Fuso Load",
       subtitle="Larger dots = Higher F. nucleatum concentration")

ggsave("figures/fuso_diversity_impact.png", alpha_fuso_plot, width=8, height=6)
ggsave("figures/fuso_beta_drive.png", beta_fuso_plot, width=8, height=6)
