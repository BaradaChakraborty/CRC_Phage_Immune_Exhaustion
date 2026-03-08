library(dada2); library(tidyverse)

path <- "data/MiSeq_SOP"
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs, 
                     truncLen = c(240, 160), # Cut Forward at 240, Reverse at 160
                     maxN = 0,               # Discard any read with an unknown letter (N)
                     maxEE = c(2, 2),        # Maximum Expected Errors allowed
                     truncQ = 2,             # Truncate at the first instance of a crash in quality
                     rm.phix = TRUE,         # Remove PhiX (a control virus used in the sequencing machine)
                     compress = TRUE, 
                     multithread = FALSE)    # MUST be FALSE on Windows computers!

head(out)

errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)

error_plot <- plotErrors(errF, nominalQ = TRUE)

ggsave(filename = "figures/error_rates_forward.png", 
       plot = error_plot, 
       width = 10, 
       height = 8, 
       dpi = 300)

error_plot

