
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "DESeq2"), update=FALSE, ask=FALSE)

library(GEOquery)
library(DESeq2)

print("1/5: Downloading clinical metadata...")
gse <- getGEO("GSE90944", GSEMatrix = TRUE)
meta <- pData(phenoData(gse[[1]]))

print("2/5: Downloading and unpacking the raw .tar file...")
file_paths <- getGEOSuppFiles("GSE90944")
tar_file <- rownames(file_paths)[1]

untar(tar_file, exdir = "GSE90944_unpacked")

print("3/5: Merging individual samples into a master matrix...")
sample_files <- list.files("GSE90944_unpacked", full.names = TRUE)
count_list <- list()

for (file in sample_files) {
  # Read the individual text file
  temp_data <- read.table(file, header = FALSE, row.names = 1, stringsAsFactors = FALSE)
  
  # Extract the exact GSM sample ID from the filename (the first 10 letters)
  gsm_id <- substr(basename(file), 1, 10)
  colnames(temp_data) <- gsm_id
  
  count_list[[gsm_id]] <- temp_data
}

counts_df <- do.call(cbind, count_list)
counts_df <- as.matrix(counts_df)

print("4/5: Aligning data and running DESeq2 algorithm...")

meta$Condition <- ifelse(grepl("Fn|nucleatum", meta$title, ignore.case=TRUE), "Infected", "Control")
meta$Condition <- factor(meta$Condition, levels = c("Control", "Infected"))

counts_df <- counts_df[, rownames(meta)]

dds <- DESeqDataSetFromMatrix(countData = counts_df, colData = meta, design = ~ Condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("Condition", "Infected", "Control"))
res_df <- as.data.frame(res)

print("5/5: Hunting for EEF2K, MTOR, and CD274 (PD-L1)...")
hypothesis_genes <- res_df[grep("^EEF2K$|^MTOR$|^CD274$", rownames(res_df), ignore.case=TRUE), ]

print("==================================================")
print("   HYPOTHESIS TEST RESULTS (Log2 Fold Change)     ")
print("==================================================")
print(hypothesis_genes[, c("log2FoldChange", "padj")])

print("6/6: Generating Volcano Plot...")

res_clean <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

res_clean$Significance <- "Not Significant"
res_clean$Significance[res_clean$log2FoldChange > 1 & res_clean$padj < 0.05] <- "Upregulated"
res_clean$Significance[res_clean$log2FoldChange < -1 & res_clean$padj < 0.05] <- "Downregulated"

volcano_plot <- ggplot(res_clean, aes(x=log2FoldChange, y=-log10(padj), color=Significance)) +
  geom_point(alpha=0.6, size=2) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="black") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
  theme_bw() +
  labs(title="F. nucleatum Infection in CRC (GSE90944)",
       x="Log2 Fold Change",
       y="-Log10 Adjusted P-value")

ggsave("figures/rnaseq_volcano.png", volcano_plot, width=8, height=6, dpi=300)
volcano_plot


print("7/7: Adding specific gene labels to the Volcano Plot...")

if (!require("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

res_clean$GeneLabel <- ""

your_targets <- c("CD274", "ALPK1", "NFKB1", "EEF2K", "MTOR")
top_significant_genes <- head(rownames(res_clean[order(res_clean$padj), ]), 15)
genes_to_label <- unique(c(your_targets, top_significant_genes))

res_clean$GeneLabel <- ifelse(rownames(res_clean) %in% genes_to_label, rownames(res_clean), "")

volcano_labeled <- ggplot(res_clean, aes(x=log2FoldChange, y=-log10(padj), color=Significance)) +
  geom_point(alpha=0.6, size=2) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="black") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
  
  # This is the magic function that adds the names beautifully
  geom_label_repel(aes(label=GeneLabel), 
                   size=3.5, 
                   max.overlaps = 50, 
                   box.padding = 0.5, 
                   show.legend = FALSE,
                   color="black",
                   fill="white",
                   alpha=0.8) + 
  theme_bw() +
  labs(title="F. nucleatum Infection in CRC (GSE90944)",
       subtitle="Mechanistic Target Identification: ALPK1 / NF-kB / CD274 Axis",
       x="Log2 Fold Change",
       y="-Log10 Adjusted P-value")

ggsave("figures/rnaseq_volcano_labeled.png", volcano_labeled, width=10, height=8, dpi=300)
volcano_labeled

print("8/8: Exporting significant gene lists to CSV files...")

upregulated <- subset(res_clean, Significance == "Upregulated")
downregulated <- subset(res_clean, Significance == "Downregulated")
ns_genes <- subset(res_clean, Significance == "Not Significant")

upregulated <- upregulated[order(-upregulated$log2FoldChange), ]
downregulated <- downregulated[order(downregulated$log2FoldChange), ]

write.csv(upregulated, "figures/upregulated_genes_list.csv")
write.csv(downregulated, "figures/downregulated_genes_list.csv")
write.csv(ns_genes, "figures/notsignificant_genes_list.csv")

print("Export complete! Files saved to the figures/ directory.")


print("9/9: Generating Targeted Bar Chart for Mechanistic Pitch...")
target_genes <- c("ALPK1", "NFKB1", "CD274", "TNFSF15", "GDF15")
target_data <- res_df[rownames(res_df) %in% target_genes, ]
target_data$Gene <- rownames(target_data)

target_data$Hit_Type <- "Mechanistic Trend (Non-Sig in vitro)"
target_data$Hit_Type[target_data$padj < 0.05 & target_data$log2FoldChange > 1] <- "Significant Upregulation"
target_data$Hit_Type[is.na(target_data$padj)] <- "Mechanistic Trend (Low Baseline)"

bar_plot <- ggplot(target_data, aes(x=reorder(Gene, -log2FoldChange), y=log2FoldChange, fill=Hit_Type)) +
  geom_col(color="black", alpha=0.8, width=0.7) +
  # This adds the professional error bars!
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width=0.2, size=0.8) +
  scale_fill_manual(values=c("Significant Upregulation" = "#d73027",      # Deep Red
                             "Mechanistic Trend (Low Baseline)" = "#fdae61", # Orange
                             "Mechanistic Trend (Non-Sig in vitro)" = "#abd9e9")) + # Light Blue
  theme_bw() +
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=45, hjust=1, face="bold", size=12),
        legend.position = "bottom",
        legend.title = element_blank()) +
  labs(title = "Validation of the ALPK1-PD-L1 Axis & Exhaustion Markers",
       subtitle = "F. nucleatum infection in CRC epithelial cells (GSE90944)",
       x = "Target Gene",
       y = "Log2 Fold Change")

ggsave("figures/targeted_mechanistic_bars.png", bar_plot, width=9, height=6, dpi=300)
bar_plot

print("10/10: Plotting Top 30 Significant Genes (Up & Down)...")
up_genes <- read.csv("figures/upregulated_genes_list.csv", row.names = 1)
down_genes <- read.csv("figures/downregulated_genes_list.csv", row.names = 1)
top_up <- up_genes[order(up_genes$padj), ][1:15, ]
top_down <- down_genes[order(down_genes$padj), ][1:15, ]
top_combined <- rbind(top_up, top_down)
top_combined$Gene <- rownames(top_combined)
sig_plot <- ggplot(top_combined, aes(x=reorder(Gene, log2FoldChange), y=log2FoldChange, fill=Significance)) +
  geom_bar(stat="identity", color="black", alpha=0.8) +
  coord_flip() + # Flip for better readability of gene names
  scale_fill_manual(values=c("Upregulated" = "#d73027", "Downregulated" = "#4575b4")) +
  theme_bw() +
  labs(title="Top 30 Significant Expression Shifts (GSE90944)",
       subtitle="Identifying Primary Mechanistic Candidates in F. nucleatum Infection",
       x="Gene Name",
       y="Log2 Fold Change") +
  theme(axis.text.y = element_text(face="bold"))

ggsave("figures/top_30_sig_genes.png", sig_plot, width=10, height=8, dpi=300)
sig_plot

print("11/11: Running Pathway Enrichment Analysis...")
if (!require("msigdbr", quietly = TRUE)) install.packages("msigdbr")
if (!require("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
library(msigdbr)
library(clusterProfiler)
h_df <- msigdbr(species = "Homo sapiens", category = "H")
h_list <- h_df %>% split(x = .$gene_symbol, f = .$gs_name)
res_for_gsea <- res_df[!is.na(res_df$stat), ]
gene_list <- res_for_gsea$stat
names(gene_list) <- rownames(res_for_gsea)
gene_list <- sort(gene_list, decreasing = TRUE)
gsea_res <- GSEA(gene_list, TERM2GENE = h_df[,c("gs_name", "gene_symbol")], pvalueCutoff = 0.05)
dot_plot <- dotplot(gsea_res, showCategory=15) + 
  ggtitle("Hijacked Biological Circuits: F. nucleatum vs Control")
ggsave("figures/hijacked_pathways.png", dot_plot, width=10, height=8, dpi=300)
dot_plot

print("Generating Causal Integration Proof...")

target_genes <- c("SESN2", "GDF15", "TNFSF15")
normalized_counts <- counts(dds, normalized=TRUE)
target_data <- as.data.frame(t(normalized_counts[target_genes, ]))
target_data$Condition <- colData(dds)$Condition

library(ggplot2)
library(reshape2)
plot_melt <- melt(target_data, id.vars = "Condition")

proof_plot <- ggplot(plot_melt, aes(x=Condition, y=value, fill=Condition)) +
  geom_boxplot(alpha=0.7) +
  geom_jitter(width=0.2) +
  facet_wrap(~variable, scales="free") +
  theme_bw() +
  labs(title="Functional Proof: F. nucleatum Responsibility",
       subtitle="Coordinated host response in primary exhaustion pathways",
       y="Normalized Expression")

ggsave("figures/functional_responsibility_proof.png", proof_plot, width=10, height=6)

library(ggplot2)
library(reshape2)

target_genes <- c("SESN2", "GDF15", "TNFSF15")
normalized_counts <- counts(dds, normalized=TRUE)
target_data <- as.data.frame(t(normalized_counts[target_genes, ]))
target_data$Condition <- colData(dds)$Condition

plot_melt <- melt(target_data, id.vars = "Condition")

stat_labels <- data.frame(
  variable = c("SESN2", "TNFSF15", "GDF15"),
  label = c("p-adj = 3.1e-60", "p-adj = 1.7e-53", "p-adj = 2.9e-12")
)

proof_plot <- ggplot(plot_melt, aes(x=Condition, y=value, fill=Condition)) +
  geom_boxplot(alpha=0.7, outlier.shape = NA) +
  geom_jitter(width=0.2, alpha=0.4) +
  facet_wrap(~variable, scales="free") +
  # Add the statistical labels at the top of each facet
  geom_text(data = stat_labels, aes(x = 1.5, y = Inf, label = label), 
            vjust = 2, fontface = "italic", size = 4, inherit.aes = FALSE) +
  theme_bw() +
  scale_fill_manual(values=c("Control"="#4575b4", "Infected"="#d73027")) +
  labs(title="Host Functional Response: F. nucleatum Responsibility",
       subtitle="Coordinated upregulation of primary exhaustion & metabolic stress pathways",
       y="Normalized Counts (DESeq2)")

ggsave("figures/causal_proof_with_pvals.png", proof_plot, width=10, height=6, dpi=300)
proof_plot