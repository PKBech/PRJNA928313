
########################
#Clustering AD ASV at 99% 
#######################

## Packages ##
library(tidyverse)
library("devtools")
library("Biostrings")
library("phyloseq")
library(vegan)
library(biomeUtils) # For clustering
library(DECIPHER)
library(speedyseq)

#### Load data ####

#AD
load('16S 18S AD Preprocessing/Metataxonomic/3.AD/ps.AD.wTree.RData') 
ps_AD_filtered <- filter_taxa(ps_AD_reduced, function (x) {sum(x > 10) > 3}, prune=TRUE)

#Cluster ADs at 99%
dna <- refseq(ps_AD_filtered)

nproc <- 1 # Increase to use multiple processorss
set.seed(41)
aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
d <- DECIPHER::DistanceMatrix(aln, processors = nproc)

#DECIPHER package cluster algorithm 
set.seed(41)
clusters <- DECIPHER::TreeLine(myDistMatrix=d,
                               method = "complete",
                               cutoff = 0.01, # use `cutoff = 0.01` for a 99% OTU
                               type = "clusters",
                               processors = nproc)

# Merge clusters with the AD PS (Using speedyseq)
ps_AD_filtered_OBU99 <- speedyseq::merge_taxa_vec(ps_AD_filtered, group = clusters$cluster)

save(ps_AD_filtered_OBU99, 
     file = "16S 18S AD Preprocessing/Metataxonomic/3.AD/ps.AD.wTree_OBU99.RData")