#### ---- Directory ---- ####
# In fastq both forward and reverse are located together

#### ---- Path and prepare demultiplexed seqs ---- ####
pathF <- "FWD" # Path to forward reads
pathR <- "REV" #Path to reverse reads

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(pathF, pattern="_R1-trimmed.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(pathR, pattern="_R2-trimmed.fq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#### ---- Quality profiles ---- #### 
qualityF = plotQualityProfile(fnFs[1:2])
saveRDS(qualityF, "qualityprofileF.rds")
qualityR = plotQualityProfile(fnRs[1:2])
saveRDS(qualityR, "qualityprofileR.rds")

#### ---- Filter and trim reads ---- ####

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(pathF, "filtered", paste0(sample.names, "_F_filt.fq.gz"))
filtRs <- file.path(pathR, "filtered", paste0(sample.names, "_R_filt.fq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


# set maxEE to 2 as it is a better measure than average quality score
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     maxN=0, maxEE=c(2,2), truncQ=2,
                     rm.phix=TRUE, compress=TRUE, 
                     multithread=5, verbose = TRUE) # 14 min 

#Run each sample at a time
for(i in seq_along(fnFs)) {
  cat("Sample", i, "\n")
  cat("fnF:", fnFs[[i]], "-- fnR:", fnRs[[i]], "\n")
  out <- filterAndTrim(fnFs[[i]], filtFs[[i]], fnRs[[i]], filtRs[[i]], maxEE=c(2,2),
                       maxN=0, rm.phix=TRUE, compress=TRUE, multithread=5)
}


head(out)

tt = assignTaxonomy(seqtab.nochim, "silva_132.18s.99_rep_set.dada2.fa.gz", 
                    tryRC = TRUE, verbose = TRUE, multithread = 10)