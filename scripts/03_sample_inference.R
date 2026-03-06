# STEP 5: SAMPLE INFERENCE
library(dada2); library(tidyverse)
dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
dadaRs <- dada(filtRs, err = errR, multithread = FALSE)
print("Forward Read Inference for Sample 1:")
print(dadaFs[[1]])
# STEP 6: MERGING PAIRED-END READS
# Merge the forward and reverse reads together
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merged data frame from the very first sample
head(mergers[[1]])
# STEP 7: CONSTRUCT SEQUENCE TABLE & REMOVE CHIMERAS
# This creates a matrix: Rows are Samples, Columns are unique DNA sequences
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Remove Chimeras (PCR artifacts)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)

# Calculate the percentage of our data that was biologically real (non-chimeric)
real_proportion <- sum(seqtab.nochim) / sum(seqtab)
print(paste("Proportion of clean sequences retained:", round(real_proportion * 100, 2), "%"))