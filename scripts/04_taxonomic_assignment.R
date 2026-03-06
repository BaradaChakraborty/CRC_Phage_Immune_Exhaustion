# STEP 8: TAXONOMIC ASSIGNMENT WITH SILVA
library(dada2); library(tidyverse)
options(timeout = 600)
silva_url <- "https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
download.file(silva_url, destfile = "data/silva_nr99_v138.1_train_set.fa.gz", mode = "wb")
print("Assigning taxonomy... This may take a few minutes depending on your computer's RAM.")
taxa <- assignTaxonomy(seqtab.nochim, "data/silva_nr99_v138.1_train_set.fa.gz", multithread = FALSE)
# Inspect the taxonomic assignments
# We create a copy to print and remove the massive DNA sequences from the row names 
# so the matrix fits cleanly in your R console
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
# Search the Genus column of your taxonomy table for 'Fusobacterium'
# grep() searches for a specific text pattern and ignores NAs
fuso_search <- taxa[grep("Fusobacterium", taxa[, "Genus"], ignore.case = TRUE), ]

# Print the result
print(fuso_search)