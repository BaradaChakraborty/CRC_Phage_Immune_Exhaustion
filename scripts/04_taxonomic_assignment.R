library(dada2); library(tidyverse)
options(timeout = 600)
silva_url <- "https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
download.file(silva_url, destfile = "data/silva_nr99_v138.1_train_set.fa.gz", mode = "wb")
print("Assigning taxonomy... This may take a few minutes depending on your computer's RAM.")
taxa <- assignTaxonomy(seqtab.nochim, "data/silva_nr99_v138.1_train_set.fa.gz", multithread = FALSE)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
fuso_search <- taxa[grep("Fusobacterium", taxa[, "Genus"], ignore.case = TRUE), ]

print(fuso_search)