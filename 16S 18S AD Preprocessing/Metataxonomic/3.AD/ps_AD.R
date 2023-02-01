setwd("/Users/pernillekjersgaardbech/Documents/JH21_Biofilm/AD_amplicons/")
library(rRDPData)
library(debar)
library(tidyverse)
library("devtools")
library(Rcpp)
library("dada2")
library(plyr)
library("Biostrings")
library("phyloseq")
library(ggplot2)
library( "DESeq2" )
library(vegan)
library(venn)

#Load in the ASV AD table
seqtab.nochim.raw <- readRDS('seqtab.nochim.rds_AD.rds')
Meta_data_chlaC_mean <- as.data.frame(readRDS('../Metagenomes/function_annotations/Meta_data_chlaC_mean'), )

str(seqtab.nochim.raw)

samples.out <- rownames(seqtab.nochim.raw)

subject <- sapply(strsplit(samples.out, "-"), `[`, 1)
timepoint <- sapply(strsplit(samples.out, "-"), `[`, 2)
day <- as.numeric(sapply(strsplit(timepoint, "T"), `[`, 2))
strain <- sapply(strsplit(samples.out, "JH21-"), `[`, 2)
position <- sapply(strsplit(samples.out, "Succession-"), `[`, 2)

chla <- left_join(metadata, Meta_data_chlaC_mean, "Timepoint")
chla_suc <- as.numeric(chla$ChlaC_mean_suc)
chla_ass <- as.numeric(chla$ChlaC_mean_ass)


sposition <- sapply(strsplit(position, "-"), `[`, 2)
sample <- samples.out


metadata <- data.frame(Subject=subject, Timepoint=timepoint, Day=day, Strain=strain, 
                       Position=position, Sample=sample, Chla_suc = chla_suc, Chla_ass = chla_ass)
rownames(metadata) <- samples.out
metadata$Day[1:9] <- 8
metadata$Day[23:25] <- 0
metadata$Day[10:17] <- 8

metadata$Rep <- paste0(metadata$Subject,"-",metadata$Timepoint)
str(metadata)



#Remove false positives with Decomtam ----
library(decontam)
packageVersion("decontam") 

dim(seqtab.nochim.raw) # blank are column 17:22
vector_for_decontam <- c(rep(FALSE, c(16)), rep(TRUE, c(6)),  rep(FALSE, c(3)),  rep(TRUE, c(2)), rep(FALSE, c(132)))

contam_df <- isContaminant(t(seqtab.nochim.raw), neg=vector_for_decontam)

table(contam_df$contaminant) # identified 16 as contaminants

# getting vector holding the identified contaminant IDs
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])


# making new count table
seqtab.nochim.raw_decontam <- seqtab.nochim.raw[!row.names(seqtab.nochim.raw) %in% contam_asvs, ]

#remove blanks

#Phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim.raw_decontam, taxa_are_rows=FALSE), 
               sample_data(metadata))


#remove blanks
ps <- subset_samples(ps, Timepoint != "Negative" & Subject != "NegativeBioelement")


dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("AD_ASV", seq(ntaxa(ps)))
ps

saveRDS(ps, "ps_AD")



