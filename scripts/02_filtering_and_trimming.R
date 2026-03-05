# STEP 3: FILTERING AND TRIMMING
library(dada2); library(tidyverse)

path <- "data/MiSeq_SOP"
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Execute the Filtering and Trimming mathematically
out <- filterAndTrim(fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs, 
                     truncLen = c(240, 160), # Cut Forward at 240, Reverse at 160
                     maxN = 0,               # Discard any read with an unknown letter (N)
                     maxEE = c(2, 2),        # Maximum Expected Errors allowed
                     truncQ = 2,             # Truncate at the first instance of a crash in quality
                     rm.phix = TRUE,         # Remove PhiX (a control virus used in the sequencing machine)
                     compress = TRUE, 
                     multithread = FALSE)    # MUST be FALSE on Windows computers!

# View the results (How many reads survived the filter?)
head(out)

# 1. Train the machine-learning model on the forward and reverse reads
# multithread = FALSE is required for Windows
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)

# 2. Plot the estimated error rates for the forward reads
# This visualizes the probability of every possible nucleotide mutation
error_plot <- plotErrors(errF, nominalQ = TRUE)

# 3. Save the error plot to your figures directory
ggsave(filename = "figures/error_rates_forward.png", 
       plot = error_plot, 
       width = 10, 
       height = 8, 
       dpi = 300)

# Display the plot in RStudio
error_plot

