# ---------------------------------------------------------
# HYPOTHESIS TEST: F. nucleatum -> mTOR -> Immune Exhaustion
# ---------------------------------------------------------

# STEP 1: LOAD THE LIBRARIES
# (We install the GEO downloader if you don't have it)
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "DESeq2"), update=FALSE, ask=FALSE)

library(GEOquery)
library(DESeq2)

# STEP 1: DOWNLOAD METADATA
print("1/5: Downloading clinical metadata...")
gse <- getGEO("GSE90944", GSEMatrix = TRUE)
meta <- pData(phenoData(gse[[1]]))

# STEP 2: DOWNLOAD AND UNPACK THE RAW DATA (.tar archive)
print("2/5: Downloading and unpacking the raw .tar file...")
file_paths <- getGEOSuppFiles("GSE90944")
tar_file <- rownames(file_paths)[1]

# This is the fix! We unpack the archive into a new folder
untar(tar_file, exdir = "GSE90944_unpacked")

# STEP 3: GLUE THE 6 SAMPLES INTO ONE MASTER TABLE
print("3/5: Merging individual samples into a master matrix...")
sample_files <- list.files("GSE90944_unpacked", full.names = TRUE)
count_list <- list()

# Loop through all 6 unpacked files and read them
for (file in sample_files) {
  # Read the individual text file
  temp_data <- read.table(file, header = FALSE, row.names = 1, stringsAsFactors = FALSE)
  
  # Extract the exact GSM sample ID from the filename (the first 10 letters)
  gsm_id <- substr(basename(file), 1, 10)
  colnames(temp_data) <- gsm_id
  
  count_list[[gsm_id]] <- temp_data
}

# Bind them side-by-side and force them into a clean numeric matrix
counts_df <- do.call(cbind, count_list)
counts_df <- as.matrix(counts_df)

# STEP 4: ALIGN THE METADATA AND RUN DESEQ2
print("4/5: Aligning data and running DESeq2 algorithm...")

# Create our "Infected" vs "Control" groups
meta$Condition <- ifelse(grepl("Fn|nucleatum", meta$title, ignore.case=TRUE), "Infected", "Control")
meta$Condition <- factor(meta$Condition, levels = c("Control", "Infected"))

# Ensure the columns in our math table perfectly match the rows in our metadata table
counts_df <- counts_df[, rownames(meta)]

# Run the core statistics
dds <- DESeqDataSetFromMatrix(countData = counts_df, colData = meta, design = ~ Condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("Condition", "Infected", "Control"))
res_df <- as.data.frame(res)

# STEP 5: PROVE THE HYPOTHESIS
print("5/5: Hunting for EEF2K, MTOR, and CD274 (PD-L1)...")
hypothesis_genes <- res_df[grep("^EEF2K$|^MTOR$|^CD274$", rownames(res_df), ignore.case=TRUE), ]

print("==================================================")
print("   HYPOTHESIS TEST RESULTS (Log2 Fold Change)     ")
print("==================================================")
print(hypothesis_genes[, c("log2FoldChange", "padj")])

# STEP 6: VISUALIZATION (VOLCANO PLOT)
print("6/6: Generating Volcano Plot...")

# Remove NA values so ggplot doesn't throw an error
res_clean <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

# Create a new column to label significant genes (padj < 0.05 AND absolute log2FC > 1)
res_clean$Significance <- "Not Significant"
res_clean$Significance[res_clean$log2FoldChange > 1 & res_clean$padj < 0.05] <- "Upregulated"
res_clean$Significance[res_clean$log2FoldChange < -1 & res_clean$padj < 0.05] <- "Downregulated"

# Generate the plot
volcano_plot <- ggplot(res_clean, aes(x=log2FoldChange, y=-log10(padj), color=Significance)) +
  geom_point(alpha=0.6, size=2) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="black") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
  theme_bw() +
  labs(title="F. nucleatum Infection in CRC (GSE90944)",
       x="Log2 Fold Change",
       y="-Log10 Adjusted P-value")

# Save and display
ggsave("figures/rnaseq_volcano.png", volcano_plot, width=8, height=6, dpi=300)
volcano_plot


# STEP 7: ADDING SMART GENE LABELS (ggrepel)
print("7/7: Adding specific gene labels to the Volcano Plot...")

# 1. Install and load the ggrepel package
if (!require("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

# 2. Create an empty column for our labels
res_clean$GeneLabel <- ""

# 3. Define exactly who we want to label!
# We want your specific mechanism targets, plus we will grab the top 15 most significant genes
your_targets <- c("CD274", "ALPK1", "NFKB1", "EEF2K", "MTOR")
top_significant_genes <- head(rownames(res_clean[order(res_clean$padj), ]), 15)
genes_to_label <- unique(c(your_targets, top_significant_genes))

# 4. Fill the label column ONLY if the gene is in our target list
res_clean$GeneLabel <- ifelse(rownames(res_clean) %in% genes_to_label, rownames(res_clean), "")

# 5. Generate the upgraded, labeled plot
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

# 6. Save and display the final masterpiece
ggsave("figures/rnaseq_volcano_labeled.png", volcano_labeled, width=10, height=8, dpi=300)
volcano_labeled

# ---------------------------------------------------------
# STEP 8: EXPORTING GENE LISTS TO CSV
# ---------------------------------------------------------
print("8/8: Exporting significant gene lists to CSV files...")

# Create separate tables for each category
upregulated <- subset(res_clean, Significance == "Upregulated")
downregulated <- subset(res_clean, Significance == "Downregulated")
ns_genes <- subset(res_clean, Significance == "Not Significant")

# Order the Upregulated and Downregulated genes by fold change so the most extreme are at the top
upregulated <- upregulated[order(-upregulated$log2FoldChange), ]
downregulated <- downregulated[order(downregulated$log2FoldChange), ]

# Save them as permanent spreadsheet files in your data or figures folder
write.csv(upregulated, "figures/upregulated_genes_list.csv")
write.csv(downregulated, "figures/downregulated_genes_list.csv")
# Note: We usually don't export the Not Significant genes to save space, but here it is:
write.csv(ns_genes, "figures/notsignificant_genes_list.csv")

print("Export complete! Files saved to the figures/ directory.")

# ---------------------------------------------------------
# STEP 9: TARGETED MECHANISTIC BAR CHART

print("9/9: Generating Targeted Bar Chart for Mechanistic Pitch...")
target_genes <- c("ALPK1", "NFKB1", "CD274", "TNFSF15", "GDF15")
# (We use the original res_df so we don't lose CD274 due to its NA p-value)
target_data <- res_df[rownames(res_df) %in% target_genes, ]
target_data$Gene <- rownames(target_data)

# Create a custom category to color-code the bars professionally
target_data$Hit_Type <- "Mechanistic Trend (Non-Sig in vitro)"
target_data$Hit_Type[target_data$padj < 0.05 & target_data$log2FoldChange > 1] <- "Significant Upregulation"
target_data$Hit_Type[is.na(target_data$padj)] <- "Mechanistic Trend (Low Baseline)"

# Generate the bar chart with Standard Error (lfcSE) whiskers
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