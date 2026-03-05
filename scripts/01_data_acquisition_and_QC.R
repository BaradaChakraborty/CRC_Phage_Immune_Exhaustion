renv::init()
# STEP 1: ENVIRONMENT SETUP & PACKAGE INSTALLATION
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("dada2", "phyloseq"))
install.packages("tidyverse")
library(dada2)
library(phyloseq)
library(tidyverse)
print("All core microbiome packages successfully installed and loaded!")
url <- "https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip"
download.file(url, destfile = "data/miseqsopdata.zip", mode = "wb")
options(timeout = 600)
download.file(url, destfile = "data/miseqsopdata.zip", mode = "wb")
unzip("data/miseqsopdata.zip", exdir = "data")
path <- "data/MiSeq_SOP"
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
plotQualityProfile(fnFs[1:2])
dir.create("figures", showWarnings = FALSE)
forward_quality_plot <- plotQualityProfile(fnFs[1:2])
ggsave(filename = "figures/forward_read_quality.png", 
       plot = forward_quality_plot, 
       width = 8, 
       height = 6, 
       dpi = 300)
print("Figure successfully saved to the 'figures' folder at 300 DPI!")
