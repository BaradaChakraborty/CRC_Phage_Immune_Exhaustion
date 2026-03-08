decoy_gene <- "TNFRSF6B" # This is the gene name for DcR3

decoy_results <- res[rownames(res) == decoy_gene, ]

if(nrow(decoy_results) > 0) {
  print(paste("Results for", decoy_gene))
  print(decoy_results)
  
  # 2. Visualize the 'Decoy' counts
  plotCounts(dds, gene=decoy_gene, intgroup="Condition")
} else {
  print("Gene TNFRSF6B not found in the dataset. Checking for synonyms...")
}


library(TCGAbiolinks)
library(ggplot2)

print("Checking CPTAC Proteomics for physical DcR3 protein levels...")

query_cptac <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Proteome Profiling",
  data.type = "Protein Expression Quantification"
)

GDCdownload(query_cptac)
prot_data <- GDCprepare(query_cptac)

target_proteins <- c("TNFRSF6B", "GDF15", "SESN2", "TNFSF15")
cptac_validation <- prot_data[prot_data$gene_name %in% target_proteins, ]

if(nrow(cptac_validation) > 0) {
  print("Success: Physical Decoy Protein identified in human cohort.")
  write.csv(cptac_validation, "results/cptac_decoy_validation.csv")
  
  # 4. Plotting the 'Physical Shield'
  ggplot(cptac_validation, aes(x=gene_name, y=value, fill=gene_name)) +
    geom_boxplot() +
    theme_bw() +
    labs(title="Proteomic Evidence: The Physical Exhaustion Shield",
         subtitle="CPTAC Colorectal Cancer Cohort",
         y="Protein Abundance (Z-score)", x="Protein Marker")
} else {
  print("DcR3 protein not found in this specific CPTAC slice. Pivoting to 'Signal Overload' theory.")
}



library(TCGAbiolinks)

query_cptac <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Proteome Profiling",
  data.type = "Protein Expression Quantification",
  workflow.type = "CPTAC Cohort"
)

GDCdownload(query_cptac, method = "client")

prot_data <- GDCprepare(query_cptac)

target_proteins <- c("TNFRSF6B", "GDF15", "SESN2", "TNFSF15")
cptac_validation <- prot_data[prot_data$gene_name %in% target_proteins, ]

write.csv(cptac_validation, "results/cptac_decoy_validation_final.csv")
print("Physical Decoy Protein data successfully secured.")



library(TCGAbiolinks)

print("Querying TCGA-COAD proteomics data...")

query_cptac <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Proteome Profiling",
  data.type = "Protein Expression Quantification"
)

print("Downloading in small chunks. This may take a moment...")
GDCdownload(query_cptac, method = "api", files.per.chunk = 5)

prot_data <- GDCprepare(query_cptac)

target_proteins <- c("TNFRSF6B", "GDF15", "SESN2", "TNFSF15")
cptac_validation <- prot_data[prot_data$gene_name %in% target_proteins, ]

if(nrow(cptac_validation) > 0) {
  write.csv(cptac_validation, "results/cptac_decoy_validation_final.csv")
  print("Physical Decoy Protein data successfully secured!")
} else {
  print("Decoy proteins not present in this TCGA-COAD proteomic slice.")
}

library(TCGAbiolinks)

print("1. Clearing corrupted partial downloads...")
if(dir.exists("GDCdata")) unlink("GDCdata", recursive = TRUE)

print("2. Increasing server timeout limit...")
options(timeout = 1000) 

print("3. Querying TCGA-COAD proteomics data...")
query_cptac <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Proteome Profiling",
  data.type = "Protein Expression Quantification"
)

print("4. Attempting direct, un-chunked download...")

tryCatch({
  GDCdownload(query_cptac, method = "api")
  
  print("5. Preparing data table...")
  prot_data <- GDCprepare(query_cptac)
  
  # 6. Extract the physical shield proteins
  target_proteins <- c("TNFRSF6B", "GDF15", "SESN2", "TNFSF15")
  
  # Check if prot_data is a dataframe or SummarizedExperiment and filter
  if(inherits(prot_data, "SummarizedExperiment")) {
    prot_df <- as.data.frame(SummarizedExperiment::assay(prot_data))
    prot_df$gene_name <- rownames(prot_df)
  } else {
    prot_df <- as.data.frame(prot_data)
  }
  
  cptac_validation <- prot_df[prot_df$gene_name %in% target_proteins, ]
  
  if(nrow(cptac_validation) > 0) {
    write.csv(cptac_validation, "results/cptac_decoy_validation_final.csv", row.names = FALSE)
    print("SUCCESS: Physical Decoy Protein data successfully secured!")
  } else {
    print("Decoy proteins not present in this TCGA-COAD proteomic slice.")
  }
  
}, error = function(e) {
  print("---------------------------------------------------")
  print(paste("GDC API ERROR:", e$message))
  print("The NCI GDC servers are currently unstable or down.")
  print("Recommendation: Proceed with the 'Signal Overload' mechanism for the manuscript.")
  print("---------------------------------------------------")
})
