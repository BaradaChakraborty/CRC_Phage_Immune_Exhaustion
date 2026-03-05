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