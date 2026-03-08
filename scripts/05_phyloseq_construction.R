library(dada2); library(phyloseq); library(tidyverse)
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject, 1, 1)
subject <- substr(subject, 2, 999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
rownames(samdf) <- samples.out
physeq_obj <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                       sample_data(samdf), 
                       tax_table(taxa))

print(physeq_obj)