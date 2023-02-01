
#########################
### DADA2-denoise ADs ###
#########################

#### ---- Packages ---- ####

library("dplyr")
library("Rcpp")
library("tidyr")
library("dada2"); packageVersion("dada2")


#### ---- Directory ---- ####
# In fastq both forward and reverse are located together

#### ---- Path and prepare demultiplexed seqs ---- ####
pathF <- "/home/perbec/PChem0002/Jyllinghavn_metagenomes2021/AD_amplicon/AD-amplicons-demuxed/demult_AD_Forward" # Path to forward reads
#pathR <- "REV/" #Path to reverse reads

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(pathF, pattern="_1_trimmed.fq.gz", full.names = TRUE))
#fnRs <- sort(list.files(pathR, pattern="_R2-trimmed.fq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#### ---- Quality profiles ---- #### 
qualityF = plotQualityProfile(fnFs[1:2])
saveRDS(qualityF, "qualityprofileF.rds")
#qualityR = plotQualityProfile(fnRs[1:2])
#saveRDS(qualityR, "qualityprofileR.rds")

#### ---- Filter and trim reads ---- ####

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(pathF, "filtered", paste0(sample.names, "_F_filt.fq.gz"))
#filtRs <- file.path(pathR, "filtered", paste0(sample.names, "_R_filt.fq.gz"))
names(filtFs) <- sample.names
#names(filtRs) <- sample.names


# set maxEE to 2 as it is a better measure than average quality score
out <- filterAndTrim(fnFs, filtFs, 
                     maxN=0, maxEE=2, truncQ=2,
                     rm.phix=TRUE, compress=TRUE, 
                     multithread=8, verbose = TRUE) # 14 min 
head(out)

# calculate percentage passing filter and trim (out_perc)
out_perc = as.data.frame(out)
out_perc$perc_pass = out_perc$reads.out/out_perc$reads.in*100
saveRDS(out_perc, "out_perc.rds")

##### ---- Learn error rates ---- ####

#use a parametric error model (err) to learn error rates in dataset.
#Learns by alternating estimation of err rates and inference of sample composition until solution. Starts with initial guess, 
#for which max possible error rate in the data set are used 
#(error rate if only most abundant seq is correct and all the rest are errors)

# Make error models #10.37
errF <- learnErrors(filtFs, multithread=16)
save.image(file='errF.RData')
#errR <- learnErrors(filtRs, multithread=TRUE, verbose = TRUE, nbases = 1e9)
#dada2:::checkConvergence(errF)
#dada2:::checkConvergence(errR)


# visualize the estimated error rates
#The error rates for each possible transition (A→C, A→G, …) are shown.
# Points are the observed error rates for each consensus quality score.
#The black line shows the estimated error rates after convergence of the machine-learning algorithm.
#The red line shows the error rates expected under the nominal definition of the Q-score

errFplot = plotErrors(errF, nominalQ=TRUE)
#errRplot = plotErrors(errR, nominalQ=TRUE)
#  If the plotted error model does not look like a good fit, try increasing the nreads parameter to see if the fit improves.
#saveRDS(errRplot, "errRplot.rds")
saveRDS(errFplot, "errFplot.rds")

#### --- Dereplication ---- #####

# Combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”
#remove those that is not on the list:
#The filter removed all reads: /home/perbec/PChem0002/Jyllinghavn_metagenomes2021/AD_amplicon/AD-amplicons-demuxed/demult_AD_Forward/filtered/JH21-25_F_filt.fq.gz not written.
#The filter removed all reads: /home/perbec/PChem0002/Jyllinghavn_metagenomes2021/AD_amplicon/AD-amplicons-demuxed/demult_AD_Forward/filtered/Seawater-T1-2_F_filt.fq.gz not written.
#The filter removed all reads: /home/perbec/PChem0002/Jyllinghavn_metagenomes2021/AD_amplicon/AD-amplicons-demuxed/demult_AD_Forward/filtered/JH21-16_F_filt.fq.gz not written.
#The filter removed all reads: /home/perbec/PChem0002/Jyllinghavn_metagenomes2021/AD_amplicon/AD-amplicons-demuxed/demult_AD_Forward/filtered/Seawater-T2-3_F_filt.fq.gz not written.
#The filter removed all reads: /home/perbec/PChem0002/Jyllinghavn_metagenomes2021/AD_amplicon/AD-amplicons-demuxed/demult_AD_Forward/filtered/Succession-T2-7_F_filt.fq.gz not written.

filtFs_new <- filtFs[!names(filtFs) %in% c("JH21-25", "Seawater-T1-2","JH21-16","Seawater-T2-3","Succession-T2-7")]

derepF <- derepFastq(filtFs_new,  verbose=TRUE)
#derepR <- derepFastq(filtRs, verbose=TRUE)
## TOO COMPUTATIONAL HEAVY on com - use sif 

# Name the derep-class objects by the sample names
sample.names_new <- sample.names[!names(sample.names) %in% c("JH21-25", "Seawater-T1-2","JH21-16","Seawater-T2-3","Succession-T2-7")]

names(derepF) <- sample.names_new
#names(derepR) <- sample.names


#### ---- Sample inference ---- ####
dadaFs <- dada(derepF, err=errF, multithread=TRUE)
#dadaRs <- dada(derepR, err=errR, multithread=TRUE)


#### ---- Sequence table ---- ####
#mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(dadaFs)

#### ---- Get sequencetable with no chim ---- ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=8, verbose=TRUE)
saveRDS(seqtab.nochim, "seqtab.nochim.rds")

#### ---- Statistics ---- ####

sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
write.table(track, file = "track_dada2.tsv", sep="\t")
