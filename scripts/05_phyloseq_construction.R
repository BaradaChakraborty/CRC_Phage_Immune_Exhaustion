# STEP 9: PHYLOSEQ OBJECT CONSTRUCTION
library(dada2); library(phyloseq); library(tidyverse)
# We need clinical metadata for these mouse samples.
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject, 1, 1)
subject <- substr(subject, 2, 999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
rownames(samdf) <- samples.out
# Construct the master phyloseq object!
# This fuses the ASV counts, the clinical metadata, and the biological taxonomy
physeq_obj <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                       sample_data(samdf), 
                       tax_table(taxa))

# Print the object summary to verify it was built successfully
print(physeq_obj)