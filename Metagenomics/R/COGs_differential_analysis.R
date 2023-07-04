
library("dplyr")
library("phyloseq")
#library("ADImpute")
#library(rRDPData)
#library(debar)
library(tidyverse)
library("devtools")
library(Rcpp)
library(plyr)
library("phyloseq")
library(ggplot2)
library( "DESeq2")
library(vegan)
#library(venn)


##### FUNCTIONALITY ######
### Create merged function dataframe


Coverage <- read.table("Metagenomics/Data/JH21_Co-assembly-GENE-COVERAGES.txt", sep = "\t", header = TRUE)
COG_cat <- read.table("Metagenomics/Data/COG20_CATEGORY.txt", sep = "\t", header = TRUE)
COG_Func <- read.csv("Metagenomics/Data/COG20_FUNCTION.csv", sep = "\t", header = TRUE)
head(as.data.frame(COG_cat))#KEGG <- read.table("KEGG_Class_functions.txt", sep = "\t", header = TRUE)
#Pfam <- read.table("Pfam.txt", sep = "\t", header = TRUE)
#Kofam <- read.csv("KOfam_functions.csv", sep = "\t", header = TRUE)
#colnames(Kofam) <- c("taxon", "source.kofam_FUNC","accession.kofam","kofam_FUNC","e_value.kofam_CAT")


colnames(Coverage)
Coverage <- Coverage[order(Coverage$key),]
data <- merge(Coverage,COG_cat, by.x="key", by.y="gene_callers_id", all.x = TRUE)
data <- merge(data,COG_Func, by.x="key", by.y="gene_callers_id", all.x = TRUE)
colnames(data) <- c("key","T2","T3","T4","T5","T6","T7","T9",
                    "T10","T11","T12","source.COG_CATEGORY","accession.COG","COG_CATEGORY",
                    "e_value.COG_CAT","source.COG_FUNCTION","accession.COG_FUNCTION","COG_FUNCTION","e_value.COG_FUNCTION")

head(data)
#write.csv(data,file = "Merged_Functions.csv")
# colnames(KEGG) <- c("gene_callers_id","source.KEGG","accession.KEGG","KEGG","e_value.KEGG")
# data <- merge(data,KEGG, by.x="key", by.y="gene_callers_id", all.x = TRUE)
# Pfam <- Pfam[order(Pfam$gene_callers_id),]
# Pfam <- Pfam[!duplicated(Pfam$gene_callers_id),]
# colnames(Pfam) <- c("gene_callers_id","source.Pfam","accession.Pfam","Pfam","e_value.Pfam")
# Pfam <- na.omit(Pfam)
# data <- merge(data,Pfam, by.x="key", by.y="gene_callers_id", all.x = TRUE)
#write.csv(data,file = "Merged_Functions.csv")



####### FUNCRIONS ALL #####
#Remove rows containing NAs


#cov.data <- data %>% na.omit 


cov.data <- data[,2:11]
str(cov.data)
rownames(cov.data) <- data$key

# Subset functions of genecalls from COGs, KEGG, and Pfam

Functions <- data[,c("key","COG_CATEGORY","COG_FUNCTION")]

Functions$COG_CAT <- data$COG_CATEGORY
#Create "human readable" categories for COG from http://clovr.org/docs/clusters-of-orthologous-groups-cogs/ & https://img.jgi.doe.gov/docs/COG.pdf
Functions$COG_CAT <- str_replace_all(Functions$COG_CAT,
                                     c("D"="Cell cycle control, cell division, chromosome partitioning",
                                       "M"="Cell wall/membrane/envelope biogenesis",
                                       "N"="Cell motility",
                                       "O"="Post-translational modification, protein turnover, and chaperones",
                                       "T"="Signal transduction mechanisms",
                                       "U"="Intracellular trafficking, secretion, and vesicular transport",
                                       "V"="Defense mechanisms",
                                       "W"="Extracellular structures",
                                       "Y"="Nuclear structure",
                                       "Z"="Cytoskeleton",
                                       "A"="RNA processing and modification",
                                       "B"="Chromatin structure and dynamics",
                                       "J"="Translation, ribosomal structure and biogenesis",
                                       "K"="Transcription",
                                       "L"="Replication, recombination and repair",
                                       "C"="Energy production and conversion",
                                       "E"="Amino acid transport and metabolism",
                                       "F"="Nucleotide transport and metabolism",
                                       "G"="Carbohydrate transport and metabolism",
                                       "H"="Coenzyme transport and metabolism",
                                       "I"="Lipid transport and metabolism",
                                       "P"="Inorganic ion transport and metabolism",
                                       "Q"="Secondary metabolites biosynthesis, transport, and catabolism",
                                       "R"="General function prediction only",
                                       "S"="Function unknown",
                                       "X"="Mobilome: prophages, transposons"))


## Metadata ###
sample_names <- c("T2","T3","T4","T5","T6",
                  "T7","T9", "T10","T11", "T12")
Phases <- c("pre-early","early","early","early", "early",
            "late","late","late","late", "late")
metadata <- data.frame(t(rbind(sample_names, Phases) ) )

rownames(metadata) <- metadata$sample_names

str(cov.data)
rownames(cov.data)
str(metadata)
str(Functions)
rownames(Functions) <- Functions$key
# Create Phyloseq object
tax <- as.matrix(Functions)

d <- phyloseq(otu_table(cov.data,taxa_are_rows=TRUE),
              tax_table(tax),
              sample_data(metadata))

#Exclude T2 from phyloseq
d <- subset_samples(d, Phases != "pre-early")
sample_data(d)$Phases <- factor(sample_data(d)$Phases, levels = c("early", "late"))
#Subset to specific functions (taxa)
#

unique(data$COG_CATEGORY)

### "Q"="Secondary metabolites biosynthesis, transport, and catabolism" #######
d_Q <- subset_taxa(d, COG_CATEGORY == "Secondary metabolites biosynthesis, transport and catabolism")


# Create DEseq object# using deseq instead of ANCOM-BC because power < 10
#ds2 <- phyloseq_to_deseq2(d, ~ Phases)
ds_Q2 <- phyloseq_to_deseq2(d_Q, ~ Phases)
# Run DESeq2 analysis (all taxa at once!) (cannot run this) - divide into specific function before
#dds <- DESeq(ds2)
dds_Q <- DESeq(ds_Q2)

# Investigate results
deseq.results_Q <- as.data.frame(results(dds_Q))
deseq.results_Q$taxon <- rownames(results(dds_Q))
str(Functions)
Functions_Q <- Functions %>% filter(COG_CATEGORY == "Secondary metabolites biosynthesis, transport and catabolism")
Functions_Q <- Functions_Q[-1,]
resultsNames(dds_Q)
deseq.results_Q$Functional_Category <- Functions_Q$COG_CAT[match(deseq.results_Q$taxon,Functions_Q$key)]
deseq.results_Q$COG <- Functions_Q$COG_FUNCTION[match(deseq.results_Q$taxon,Functions_Q$key)]
#deseq.results$KEGG <- Functions$KEGG[match(deseq.results$taxon,Functions$key)]
#deseq.results$Pfam <- Functions$Pfam[match(deseq.results$taxon,Functions$key)]
# Sort (arrange) by pvalue and effect size
deseq.results_Q <- deseq.results_Q %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

#Sigs for early stage
sig_early <- deseq.results_Q %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_Q %>%
  filter(log2FoldChange > 0, padj < 0.01)

sig_early_unique <- sig_early[-which(sig_early$COG %in% sig_late$COG),]
unique(sig_early_unique$COG)
arrange(sig_early_unique, desc(baseMean) )


unique_sec_early <- deseq.results_K_10cov <- deseq.results_K %>% filter(baseMean > 5) %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

sig_early_unique %>% dplyr::count(COG, sort = TRUE)

# Print the result table
# Let us only show significant hits
# library("knitr")
# knitr::kable(deseq.results_Q %>%
#         filter(padj < 0.01),
#       digits = 2)
sig <- deseq.results_Q %>%
  filter(padj < 0.01)
str(sig)

### "N"="Cell motility", #######
d_N <- subset_taxa(d, COG_CATEGORY == "Cell motility")


# Create DEseq object# using deseq instead of ANCOM-BC because power < 10
#ds2 <- phyloseq_to_deseq2(d, ~ Phases)
ds_N2 <- phyloseq_to_deseq2(d_N, ~ Phases)
# Run DESeq2 analysis (all taxa at once!) (cannot run this) - divide into specific function before
#dds <- DESeq(ds2)
dds_N <- DESeq(ds_N2)

# Investigate results
deseq.results_N <- as.data.frame(results(dds_N))
deseq.results_N$taxon <- rownames(results(dds_N))
str(Functions)
Functions_N <- Functions %>% filter(COG_CATEGORY == "Cell motility")
Functions_N <- Functions_N[-1,]
resultsNames(dds_N)
deseq.results_N$Functional_Category <- Functions_N$COG_CAT[match(deseq.results_N$taxon,Functions_N$key)]
deseq.results_N$COG <- Functions_N$COG_FUNCTION[match(deseq.results_N$taxon,Functions_N$key)]
#deseq.results$KEGG <- Functions$KEGG[match(deseq.results$taxon,Functions$key)]
#deseq.results$Pfam <- Functions$Pfam[match(deseq.results$taxon,Functions$key)]
# Sort (arrange) by pvalue and effect size
deseq.results_N <- deseq.results_N %>%
  arrange(pvalue, log2FoldChange, Functional_Category)
# Print the result table
# Let us only show significant hits
# library("knitr")
# knitr::kable(deseq.results_N %>%
#         filter(padj < 0.01),
#       digits = 2)
sig <- deseq.results_N %>%
  filter(padj < 0.01)

### "K"="Transcription", #######
#Look for	K T	Accessory gene regulator protein AgrB:  quorum sensing system is responsible for the regulation of the expression of virulence factor genes. 

d_K <- subset_taxa(d, COG_CATEGORY == "Transcription")


# Create DEseq object# using deseq instead of ANCOM-BC because power < 10
#ds2 <- phyloseq_to_deseq2(d, ~ Phases)
ds_K2 <- phyloseq_to_deseq2(d_K, ~ Phases)
# Run DESeq2 analysis (all taxa at once!) (cannot run this) - divide into specific function before
#dds <- DESeq(ds2)
dds_K <- DESeq(ds_K2)

# Investigate results
deseq.results_K <- as.data.frame(results(dds_K))
deseq.results_K$taxon <- rownames(results(dds_K))
str(Functions)
Functions_K <- Functions %>% filter(COG_CATEGORY == "Transcription")
Functions_K <- Functions_K[-1,]
resultsNames(dds_K)
deseq.results_K$Functional_Category <- Functions_K$COG_CAT[match(deseq.results_K$taxon,Functions_K$key)]
deseq.results_K$COG <- Functions_K$COG_FUNCTION[match(deseq.results_K$taxon,Functions_K$key)]
#deseq.results$KEGG <- Functions$KEGG[match(deseq.results$taxon,Functions$key)]
#deseq.results$Pfam <- Functions$Pfam[match(deseq.results$taxon,Functions$key)]
# Sort (arrange) by pvalue and effect size
deseq.results_K <- deseq.results_K %>%
  arrange(pvalue, log2FoldChange, Functional_Category)


#Sigs for early stage
sig_early <- deseq.results_K %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_K %>%
  filter(log2FoldChange > 0, padj < 0.01)

sig_early_unique <- sig_early[-which(sig_early$COG %in% sig_late$COG),]
unique(sig_early_unique$COG)
arrange(sig_early_unique, desc(baseMean) )
# Print the result table
# Let us only show significant hits
# library("knitr")
# knitr::kable(deseq.results_K %>%
#         filter(padj < 0.01),
#       digits = 2)
sig <- deseq.results_K %>%
  filter(padj < 0.01)

### "D"="Cell cycle control, cell division, chromosome partitioning" #######
d_D <- subset_taxa(d, COG_CATEGORY == "Cell cycle control, cell division, chromosome partitioning")


# Create DEseq object# using deseq instead of ANCOM-BC because power < 10
#ds2 <- phyloseq_to_deseq2(d, ~ Phases)
ds_D2 <- phyloseq_to_deseq2(d_D, ~ Phases)
# Run DESeq2 analysis (all taxa at once!) (cannot run this) - divide into specific function before
#dds <- DESeq(ds2)
dds_D <- DESeq(ds_D2)

# Investigate results
deseq.results_D <- as.data.frame(results(dds_D))
deseq.results_D$taxon <- rownames(results(dds_D))
str(Functions)
Functions_D <- Functions %>% filter(COG_CATEGORY == "Cell cycle control, cell division, chromosome partitioning")
Functions_D <- Functions_D[-1,]
resultsNames(dds_D)
deseq.results_D$Functional_Category <- Functions_D$COG_CAT[match(deseq.results_D$taxon,Functions_D$key)]
deseq.results_D$COG <- Functions_D$COG_FUNCTION[match(deseq.results_D$taxon,Functions_D$key)]
#deseq.results$KEGG <- Functions$KEGG[match(deseq.results$taxon,Functions$key)]
#deseq.results$Pfam <- Functions$Pfam[match(deseq.results$taxon,Functions$key)]
# Sort (arrange) by pvalue and effect size
deseq.results_D <- deseq.results_D %>%
  arrange(pvalue, log2FoldChange, Functional_Category)
# Print the result table
# Let us only show significant hits
# library("knitr")
# knitr::kable(deseq.results_D %>%
#         filter(padj < 0.01),
#       digits = 2)
sig <- deseq.results_D %>%
  filter(padj < 0.01)

### "T"="Signal transduction mechanisms" #######


#Look for	K T	Accessory gene regulator protein AgrB:  quorum sensing system is responsible for the regulation of the expression of virulence factor genes. 

d_T <- subset_taxa(d, COG_CATEGORY == "Signal transduction mechanisms")


# Create DEseq object# using deseq instead of ANCOM-BC because power < 10
#ds2 <- phyloseq_to_deseq2(d, ~ Phases)
ds_T2 <- phyloseq_to_deseq2(d_T, ~ Phases)
# Run DESeq2 analysis (all taxa at once!) (cannot run this) - divide into specific function before
#dds <- DESeq(ds2)
dds_T <- DESeq(ds_T2)

# Investigate results
deseq.results_T <- as.data.frame(results(dds_T))
deseq.results_T$taxon <- rownames(results(dds_T))
str(Functions)
Functions_T <- Functions %>% filter(COG_CATEGORY == "Signal transduction mechanisms")
Functions_T <- Functions_T[-1,]
resultsNames(dds_T)
deseq.results_T$Functional_Category <- Functions_T$COG_CAT[match(deseq.results_T$taxon,Functions_T$key)]
deseq.results_T$COG <- Functions_T$COG_FUNCTION[match(deseq.results_T$taxon,Functions_T$key)]
#deseq.results$KEGG <- Functions$KEGG[match(deseq.results$taxon,Functions$key)]
#deseq.results$Pfam <- Functions$Pfam[match(deseq.results$taxon,Functions$key)]
# Sort (arrange) by pvalue and effect size
deseq.results_T <- deseq.results_T %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

#Find unique functions for each group:
str(deseq.results_T)

#Sigs for early stage
sig_early <- deseq.results_T %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_T %>%
  filter(log2FoldChange > 0, padj < 0.01)

sig_early_unique <- sig_early[-which(sig_early$COG %in% sig_late$COG),]

# Print the result table
# Let us only show significant hits
# library("knitr")
# knitr::kable(deseq.results_T %>%
#         filter(padj < 0.01),
#       digits = 2)
sig <- deseq.results_T %>%
  filter(padj < 0.01)


### "S"="Function unknown": QUORUM SENSING #######


#Look for	K T	Accessory gene regulator protein AgrB:  quorum sensing system is responsible for the regulation of the expression of virulence factor genes. 

d_S <- subset_taxa(d, COG_CATEGORY == "Function unknown")


# Create DEseq object# using deseq instead of ANCOM-BC because power < 10
#ds2 <- phyloseq_So_deseq2(d, ~ Phases)
ds_S2 <- phyloseq_to_deseq2(d_S, ~ Phases)
# Run DESeq2 analysis (all taxa at once!) (cannot run this) - divide into specific function before
#dds <- DESeq(ds2)
dds_S <- DESeq(ds_S2)

# Investigate results
deseq.results_S <- as.data.frame(results(dds_S))
deseq.results_S$taxon <- rownames(results(dds_S))
str(Functions)
Functions_S <- Functions %>% filter(COG_CATEGORY == "Function unknown")
Functions_S <- Functions_S[-1,]
resultsNames(dds_S)
deseq.results_S$Functional_Category <- Functions_S$COG_CAT[match(deseq.results_S$taxon,Functions_S$key)]
deseq.results_S$COG <- Functions_S$COG_FUNCTION[match(deseq.results_S$taxon,Functions_S$key)]
#deseq.results$KEGG <- Functions$KEGG[match(deseq.results$taxon,Functions$key)]
#deseq.results$Pfam <- Functions$Pfam[match(deseq.results$taxon,Functions$key)]
# Sort (arrange) by pvalue and effect size
deseq.results_S <- deseq.results_S %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

#Find unique functions for each group:
str(deseq.results_S)

#Sigs for early stage
# sig_early <- deseq.results_S %>%
#   filter(log2FoldChange < 0, padj < 0.01)
# 
# sig_late <- deseq.results_S %>%
#   filter(log2FoldChange > 0, padj < 0.01)
# 
# sig_early_unique <- sig_early[-which(sig_early$COG %in% sig_late$COG),]

# Print the result table
# Let us only show significant hits
# library("knitr")
# knitr::kable(deseq.results_S %>%
#         filter(padj < 0.01),
#       digits = 2)
sig <- deseq.results_S %>%
  filter(padj < 0.01)




### "P"="Inorganic ion transport and metabolism", #######

d_P <- subset_taxa(d, COG_CATEGORY == "Inorganic ion transport and metabolism")


# Create DEseq object# using deseq instead of ANCOM-BC because power < 10
#ds2 <- phyloseq_Po_deseq2(d, ~ Phases)
ds_P2 <- phyloseq_to_deseq2(d_P, ~ Phases)
# Run DESeq2 analysis (all taxa at once!) (cannot run this) - divide into specific function before
#dds <- DESeq(ds2)
dds_P <- DESeq(ds_P2)

# Investigate results
deseq.results_P <- as.data.frame(results(dds_P))
deseq.results_P$taxon <- rownames(results(dds_P))
str(Functions)
Functions_P <- Functions %>% filter(COG_CATEGORY == "Inorganic ion transport and metabolism")
Functions_P <- Functions_P[-1,]
resultsNames(dds_P)
deseq.results_P$Functional_Category <- Functions_P$COG_CAT[match(deseq.results_P$taxon,Functions_P$key)]
deseq.results_P$COG <- Functions_P$COG_FUNCTION[match(deseq.results_P$taxon,Functions_P$key)]
#deseq.results$KEGG <- Functions$KEGG[match(deseq.results$taxon,Functions$key)]
#deseq.results$Pfam <- Functions$Pfam[match(deseq.results$taxon,Functions$key)]
# Sort (arrange) by pvalue and effect size
deseq.results_P <- deseq.results_P %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

#Find unique functions for each group:
str(deseq.results_P)

#Sigs for early stage
# sig_early <- deseq.results_P %>%
#   filter(log2FoldChange < 0, padj < 0.01)
# 
# sig_late <- deseq.results_P %>%
#   filter(log2FoldChange > 0, padj < 0.01)
# 
# sig_early_unique <- sig_early[-which(sig_early$COG %in% sig_late$COG),]

# Print the result table
# Let us only show significant hits
# library("knitr")
# knitr::kable(deseq.results_P %>%
#         filter(padj < 0.01),
#       digits = 2)
sig <- deseq.results_P %>%
  filter(padj < 0.01)





### "X"="Mobilome: prophages, transposons" #######

d_X <- subset_taxa(d, COG_CATEGORY == "Mobilome: prophages, transposons")


# Create DEseq object# using deseq instead of ANCOM-BC because power < 10
#ds2 <- phyloseq_Xo_deseq2(d, ~ Phases)
ds_X2 <- phyloseq_to_deseq2(d_X, ~ Phases)
# Run DESeq2 analysis (all taxa at once!) (cannot run this) - divide into specific function before
#dds <- DESeq(ds2)
dds_X <- DESeq(ds_X2)

# Investigate results
deseq.results_X <- as.data.frame(results(dds_X))
deseq.results_X$taxon <- rownames(results(dds_X))
str(Functions)
Functions_X <- Functions %>% filter(COG_CATEGORY == "Mobilome: prophages, transposons")
Functions_X <- Functions_X[-1,]
resultsNames(dds_X)
deseq.results_X$Functional_Category <- Functions_X$COG_CAT[match(deseq.results_X$taxon,Functions_X$key)]
deseq.results_X$COG <- Functions_X$COG_FUNCTION[match(deseq.results_X$taxon,Functions_X$key)]
#deseq.results$KEGG <- Functions$KEGG[match(deseq.results$taxon,Functions$key)]
#deseq.results$Pfam <- Functions$Pfam[match(deseq.results$taxon,Functions$key)]
# Sort (arrange) by pvalue and effect size
deseq.results_X <- deseq.results_X %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

#Find unique functions for each group:
str(deseq.results_X)

#Sigs for early stage
sig_early <- deseq.results_X %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_X %>%
  filter(log2FoldChange > 0, padj < 0.01)

sig_early_unique <- sig_early[-which(sig_early$COG %in% sig_late$COG),]

# Print the result table
# Let us only show significant hits
# library("knitr")
# knitr::kable(deseq.results_X %>%
#         filter(padj < 0.01),
#       digits = 2)
sig <- deseq.results_X %>%
  filter(padj < 0.01)






### "D"="Cell cycle control, cell division, chromosome partitioning" #######
d_D <- subset_taxa(d, COG_CATEGORY == "Cell cycle control, cell division, chromosome partitioning")

# Create DEseq object# using deseq instead of ANCOM-BC because power < 10
#ds2 <- phyloseq_to_deseq2(d, ~ Phases)
ds_D2 <- phyloseq_to_deseq2(d_D, ~ Phases)
# Run DESeq2 analysis (all taxa at once!) (cannot run this) - divide into specific function before
#dds <- DESeq(ds2)
dds_D <- DESeq(ds_D2)

# Investigate results
deseq.results_D <- as.data.frame(results(dds_D))
deseq.results_D$taxon <- rownames(results(dds_D))
str(Functions)
Functions_D <- Functions %>% filter(COG_CATEGORY == "Cell cycle control, cell division, chromosome partitioning")
Functions_D <- Functions_D[-1,]
resultsNames(dds_D)
deseq.results_D$Functional_Category <- Functions_D$COG_CAT[match(deseq.results_D$taxon,Functions_D$key)]
deseq.results_D$COG <- Functions_D$COG_FUNCTION[match(deseq.results_D$taxon,Functions_D$key)]
#deseq.results$KEGG <- Functions$KEGG[match(deseq.results$taxon,Functions$key)]
#deseq.results$Pfam <- Functions$Pfam[match(deseq.results$taxon,Functions$key)]
# Sort (arrange) by pvalue and effect size
deseq.results_D <- deseq.results_D %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

#Find unique functions for each group:
str(deseq.results_D)

#Sigs for early stage
sig_early <- deseq.results_D %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_D %>%
  filter(log2FoldChange > 0, padj < 0.01)

sig_early_unique <- sig_early[-which(sig_early$COG %in% sig_late$COG),]

# Print the result table
# Let us only show significant hits
# library("knitr")
# knitr::kable(deseq.results_D %>%
#         filter(padj < 0.01),
#       digits = 2)
sig <- deseq.results_D %>%
  filter(padj < 0.01)





###  "C"="Energy production and conversion" #######
d_C <- subset_taxa(d, COG_CATEGORY == "Energy production and conversion")

# Create DEseq object# using deseq instead of ANCOM-BC because power < 10
#ds2 <- phyloseq_to_Ceseq2(d, ~ Phases)
ds_C2 <- phyloseq_to_deseq2(d_C, ~ Phases)
# Run DESeq2 analysis (all taxa at once!) (cannot run this) - divide into specific function before
#dds <- DESeq(ds2)
dds_C <- DESeq(ds_C2)

# Investigate results
deseq.results_C <- as.data.frame(results(dds_C))
deseq.results_C$taxon <- rownames(results(dds_C))
str(Functions)
Functions_C <- Functions %>% filter(COG_CATEGORY == "Energy production and conversion")
Functions_C <- Functions_C[-1,]
resultsNames(dds_C)
deseq.results_C$Functional_Category <- Functions_C$COG_CAT[match(deseq.results_C$taxon,Functions_C$key)]
deseq.results_C$COG <- Functions_C$COG_FUNCTION[match(deseq.results_C$taxon,Functions_C$key)]
#deseq.results$KEGG <- Functions$KEGG[match(deseq.results$taxon,Functions$key)]
#deseq.results$Pfam <- Functions$Pfam[match(deseq.results$taxon,Functions$key)]
# Sort (arrange) by pvalue and effect size
deseq.results_C <- deseq.results_C %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

#Find unique functions for each group:
str(deseq.results_C)

#Sigs for early stage
sig_early <- deseq.results_C %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_C %>%
  filter(log2FoldChange > 0, padj < 0.01)

sig_early_unique <- sig_early[-which(sig_early$COG %in% sig_late$COG),]
sig_early_unique
# Print the result table
# Let us only show significant hits
# library("knitr")
# knitr::kable(deseq.results_C %>%
#         filter(padj < 0.01),
#       digits = 2)
sig <- deseq.results_C %>%
  filter(padj < 0.01)









###  "G"="Carbohydrate transport and metabolism"#######
d_G <- subset_taxa(d, COG_CATEGORY == "Carbohydrate transport and metabolism")

# Create DEseq object# using deseq instead of ANCOM-BC because power < 10
#ds2 <- phyloseq_to_Geseq2(d, ~ Phases)
ds_G2 <- phyloseq_to_deseq2(d_G, ~ Phases)
# Run DESeq2 analysis (all taxa at once!) (cannot run this) - divide into specific function before
#dds <- DESeq(ds2)
dds_G <- DESeq(ds_G2)

# Investigate results
deseq.results_G <- as.data.frame(results(dds_G))
deseq.results_G$taxon <- rownames(results(dds_G))
str(Functions)
Functions_G <- Functions %>% filter(COG_CATEGORY == "Carbohydrate transport and metabolism")
Functions_G <- Functions_G[-1,]
resultsNames(dds_G)
deseq.results_G$Functional_Category <- Functions_G$COG_CAT[match(deseq.results_G$taxon,Functions_G$key)]
deseq.results_G$COG <- Functions_G$COG_FUNCTION[match(deseq.results_G$taxon,Functions_G$key)]
#deseq.results$KEGG <- Functions$KEGG[match(deseq.results$taxon,Functions$key)]
#deseq.results$Pfam <- Functions$Pfam[match(deseq.results$taxon,Functions$key)]
# Sort (arrange) by pvalue and effect size
deseq.results_C <- deseq.results_C %>%
  arrange(pvalue, log2FoldChange, Functional_Category)



#Find unique functions for each group:
str(deseq.results_G)

#Sigs for early stage
sig_early <- deseq.results_G %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_G %>%
  filter(log2FoldChange > 0, padj < 0.01)

sig_early_unique <- sig_early[-which(sig_early$COG %in% sig_late$COG),]
sig_early_unique
# Print the result table
# Let us only show significant hits
# library("knitr")
# knitr::kable(deseq.results_G %>%
#         filter(padj < 0.01),
#       digits = 2)
sig <- deseq.results_G %>%
  filter(padj < 0.01)












###  "O"="Post-translational modification, protein turnover, and chaperones" #######
d_O <- subset_taxa(d, COG_CATEGORY == "Posttranslational modification, protein turnover, chaperones")

# Create DEseq object# using deseq instead of ANCOM-BC because power < 10
#ds2 <- phyloseq_to_Oeseq2(d, ~ Phases)
ds_O2 <- phyloseq_to_deseq2(d_O, ~ Phases)
# Run DESeq2 analysis (all taxa at once!) (cannot run this) - divide into specific function before
#dds <- DESeq(ds2)
dds_O <- DESeq(ds_O2)

# Investigate results
deseq.results_O <- as.data.frame(results(dds_O))
deseq.results_O$taxon <- rownames(results(dds_O))
str(Functions)
Functions_O <- Functions %>% filter(COG_CATEGORY == "Posttranslational modification, protein turnover, chaperones")
Functions_O <- Functions_O[-1,]
resultsNames(dds_O)
deseq.results_O$Functional_Category <- Functions_O$COG_CAT[match(deseq.results_O$taxon,Functions_O$key)]
deseq.results_O$COG <- Functions_O$COG_FUNCTION[match(deseq.results_O$taxon,Functions_O$key)]
#deseq.results$KEGG <- Functions$KEGG[match(deseq.results$taxon,Functions$key)]
#deseq.results$Pfam <- Functions$Pfam[match(deseq.results$taxon,Functions$key)]
# Sort (arrange) by pvalue and effect size
deseq.results_O <- deseq.results_O %>%
  arrange(pvalue, log2FoldChange, Functional_Category)




### "L"="Replication, recombination and repair"#######
d_L <- subset_taxa(d, COG_CATEGORY == "Replication, recombination and repair")

# Create DEseq object# using deseq instead of ANCOM-BC because power < 10
#ds2 <- phyloseq_to_Leseq2(d, ~ Phases)
ds_L2 <- phyloseq_to_deseq2(d_L, ~ Phases)
# Run DESeq2 analysis (all taxa at once!) (cannot run this) - divide into specific function before
#dds <- DESeq(ds2)
dds_L <- DESeq(ds_L2)

# Investigate results
deseq.results_L <- as.data.frame(results(dds_L))
deseq.results_L$taxon <- rownames(results(dds_L))
str(Functions)
Functions_L <- Functions %>% filter(COG_CATEGORY == "Replication, recombination and repair")
Functions_L <- Functions_L[-1,]
resultsNames(dds_L)
deseq.results_L$Functional_Category <- Functions_L$COG_CAT[match(deseq.results_L$taxon,Functions_L$key)]
deseq.results_L$COG <- Functions_L$COG_FUNCTION[match(deseq.results_L$taxon,Functions_L$key)]
#deseq.results$KEGG <- Functions$KEGG[match(deseq.results$taxon,Functions$key)]
#deseq.results$Pfam <- Functions$Pfam[match(deseq.results$taxon,Functions$key)]
# Sort (arrange) by pvalue and effect size
deseq.results_L <- deseq.results_L %>%
  arrange(pvalue, log2FoldChange, Functional_Category)




### "W"="Extracellular structures"#######
d_W <- subset_taxa(d, COG_CATEGORY == "Extracellular structures")

# Create DEseq object# using deseq instead of ANCOM-BC because power < 10
#ds2 <- phyloseq_to_Weseq2(d, ~ Phases)
ds_W2 <- phyloseq_to_deseq2(d_W, ~ Phases)
# Run DESeq2 analysis (all taxa at once!) (cannot run this) - divide into specific function before
#dds <- DESeq(ds2)
dds_W <- DESeq(ds_W2)

# Investigate results
deseq.results_W <- as.data.frame(results(dds_W))
deseq.results_W$taxon <- rownames(results(dds_W))
str(Functions)
Functions_W <- Functions %>% filter(COG_CATEGORY == "Extracellular structures")
Functions_W <- Functions_W[-1,]
resultsNames(dds_W)
deseq.results_W$Functional_Category <- Functions_W$COG_CAT[match(deseq.results_W$taxon,Functions_W$key)]
deseq.results_W$COG <- Functions_W$COG_FUNCTION[match(deseq.results_W$taxon,Functions_W$key)]
#deseq.results$KEGG <- Functions$KEGG[match(deseq.results$taxon,Functions$key)]
#deseq.results$Pfam <- Functions$Pfam[match(deseq.results$taxon,Functions$key)]
# Sort (arrange) by pvalue and effect size
deseq.results_W <- deseq.results_W %>%
  arrange(pvalue, log2FoldChange, Functional_Category)






### "V"="Defense mechanisms"#######
d_V <- subset_taxa(d, COG_CATEGORY == "Defense mechanisms")

# Create DEseq object# using deseq instead of ANCOM-BC because power < 10
#ds2 <- phyloseq_to_Veseq2(d, ~ Phases)
ds_V2 <- phyloseq_to_deseq2(d_V, ~ Phases)
# Run DESeq2 analysis (all taxa at once!) (cannot run this) - divide into specific function before
#dds <- DESeq(ds2)
dds_V <- DESeq(ds_V2)

# Investigate results
deseq.results_V <- as.data.frame(results(dds_V))
deseq.results_V$taxon <- rownames(results(dds_V))
str(Functions)
Functions_V <- Functions %>% filter(COG_CATEGORY == "Defense mechanisms")
Functions_V <- Functions_V[-1,]
resultsNames(dds_V)
deseq.results_V$Functional_Category <- Functions_V$COG_CAT[match(deseq.results_V$taxon,Functions_V$key)]
deseq.results_V$COG <- Functions_V$COG_FUNCTION[match(deseq.results_V$taxon,Functions_V$key)]
#deseq.results$KEGG <- Functions$KEGG[match(deseq.results$taxon,Functions$key)]
#deseq.results$Pfam <- Functions$Pfam[match(deseq.results$taxon,Functions$key)]
# Sort (arrange) by pvalue and effect size
deseq.results_V <- deseq.results_V %>%
  arrange(pvalue, log2FoldChange, Functional_Category)





####r plot differential test ####

library(ggrepel)
library(EnhancedVolcano)
str(deseq.results_Q)
deseq.results_Q_10cov <- deseq.results_Q %>% filter(baseMean > 5) %>% #Sekventeret mere end fem gange for at komme med 
  arrange(pvalue, log2FoldChange, Functional_Category)

EnhancedVolcano(deseq.results_Q_10cov,
                lab = deseq.results_Q_10cov$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
                pCutoff = 10e-3,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))


sig_early <- deseq.results_Q_10cov %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_Q_10cov %>%
  filter(log2FoldChange > 0, padj < 0.01)

unique_sec_early <- unique(sig_early)[-which(unique(sig_early$Functional_Category) %in% unique(sig_late$Functional_Category)),]
arrange(unique_transc_early, desc(baseMean))

unique_sec_early %>% dplyr::count(COG, sort = TRUE)

deseq.results_Q_10cov_Imidazolonepropionase <- filter(deseq.results_Q_10cov, grepl('Imidazolonepropionase|mycothiol', COG))
deseq.results_Q_10cov_Qselected <- filter(deseq.results_Q_10cov, grepl('Imidazolonepropionase|mycothiol|DinB|EntF', COG)) # selected for final plot

EnhancedVolcano(deseq.results_Q_10cov_Imidazolonepropionase,
                lab = deseq.results_Q_10cov_Imidazolonepropionase$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))

##Imidazolonepropionase or related amidohydrolase (HutI) (PDB:3OOQ) 19
#Bacillithiol/mycothiol S-transferase BstA/DinB, DinB/YfiT family (unrelated to E. coli DinB) (DinB) (PDB:2QE9) (PUBMED:22059487;24821014) 15


deseq.results_Q_10cov_mycothiol <- filter(deseq.results_Q_10cov, grepl('mycothiol', COG))
p = EnhancedVolcano(deseq.results_Q_10cov_mycothiol,
                    lab = deseq.results_Q_10cov_mycothiol$COG,
                    x = 'log2FoldChange',
                    y = 'padj',
                    #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
                    pCutoff = 10e-2,
                    labSize = 3,
                    title = '',
                    legendPosition = 'bottom',
                    xlim = c(-10,10))

p + theme_bw()

####plots with secondary met functions #####
# 
# deseq.results_Q_Bacillithiol <- filter(deseq.results_Q, grepl('Bacillithiol', COG)) #to cope with endogenous or exogenous toxic metabolites and oxidants.
# deseq.results_Q_Polyketide <- filter(deseq.results_Q, grepl('polyketide', COG))
# deseq.results_Q_Virulence <- filter(deseq.results_Q, grepl('Virulence', COG))
# deseq.results_Q_epoxidase <- filter(deseq.results_Q, grepl('1,2-phenylacetyl-CoA epoxidase', COG))
# deseq.results_Q_siderophore <- filter(deseq.results_Q, grepl('siderophore', COG))
# 
# p1 <- EnhancedVolcano(deseq.results_Q_Bacillithiol,
#                 lab = deseq.results_Q_Bacillithiol$COG,
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
#                 pCutoff = 10e-2,
#                 labSize = 3,
#                 title = '',
#                 legendPosition = 'bottom',
#                 xlim = c(-10,10))
# p2 <- EnhancedVolcano(deseq.results_Q_Polyketide,
#                       lab = deseq.results_Q_Polyketide$COG,
#                       x = 'log2FoldChange',
#                       y = 'padj',
#                       #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
#                       pCutoff = 10e-2,
#                       labSize = 3,
#                       title = '',
#                       legendPosition = 'bottom',
#                       xlim = c(-10,10))
# p3 <- EnhancedVolcano(deseq.results_Q_Virulence,
#                       lab = deseq.results_Q_Virulence$COG,
#                       x = 'log2FoldChange',
#                       y = 'padj',
#                       #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
#                       pCutoff = 10e-2,
#                       labSize = 3,
#                       title = '',
#                       legendPosition = 'bottom',
#                       xlim = c(-10,10))
# p4 <- EnhancedVolcano(deseq.results_Q_epoxidase,
#                       lab = deseq.results_Q_epoxidase$COG,
#                       x = 'log2FoldChange',
#                       y = 'padj',
#                       #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
#                       pCutoff = 10e-2,
#                       labSize = 3,
#                       title = '',
#                       legendPosition = 'bottom',
#                       xlim = c(-10,10))

# p5 <- EnhancedVolcano(deseq.results_Q_siderophore,
#                       lab = deseq.results_Q_siderophore$COG,
#                       x = 'log2FoldChange',
#                       y = 'padj',
#                       #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
#                       pCutoff = 10e-2,
#                       labSize = 3,
#                       title = '',
#                       legendPosition = 'bottom',
#                       xlim = c(-10,10))

#plot_grid(p1,p2,p3,p4)

deseq.results_N_10cov <- deseq.results_N %>% filter(baseMean > 5) %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

sig_early <- deseq.results_N_10cov %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_N_10cov %>%
  filter(log2FoldChange > 0, padj < 0.01)

unique_N_early <- unique(sig_early)[-which(unique(sig_early$Functional_Category) %in% unique(sig_late$Functional_Category)),]

unique_N_early %>% dplyr::count(COG, sort = TRUE) #Nothing with coverage > 5x

unique_N_late <- unique(sig_late)[-which(unique(sig_late$Functional_Category) %in% unique(sig_early$Functional_Category)),]
unique_N_late %>% dplyr::count(COG, sort = TRUE)
deseq.results_N_EPS_Flag <- filter(deseq.results_N_10cov, grepl('Flagellar|Flagellin', COG))

EnhancedVolcano(deseq.results_N_EPS_Flag,
                lab = deseq.results_N_EPS_Flag$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))

# deseq.results_N_EPS_Flag <- filter(deseq.results_N, grepl('hook', COG))
# 
# EnhancedVolcano(deseq.results_N_EPS_Flag,
#                      lab = deseq.results_N_EPS_Flag$COG,
#                      x = 'log2FoldChange',
#                      y = 'padj',
#                      #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
#                      pCutoff = 10e-2,
#                      labSize = 3,
#                      title = '',
#                      legendPosition = 'bottom',
#                      xlim = c(-10,10))

#deseq.results_K_AgrB <- filter(deseq.results_K, grepl('accessory gene regulator', COG))

###plots with transcriptional factors #####

deseq.results_K_10cov <- deseq.results_K %>% filter(baseMean > 5) %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

sig_early <- deseq.results_K_10cov %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_K_10cov %>%
  filter(log2FoldChange > 0, padj < 0.01)

unique_transc_early <- unique(sig_early)[-which(unique(sig_early$Functional_Category) %in% unique(sig_late$Functional_Category)),]
arrange(unique_transc_early, desc(baseMean))

unique_transc_early %>% dplyr::count(COG, sort = TRUE)

unique_K_late <- unique(sig_late)[-which(unique(sig_late$Functional_Category) %in% unique(sig_early$Functional_Category)),]
unique_K_late %>% dplyr::count(COG, sort = TRUE) #RpoB/RpoC


EnhancedVolcano(deseq.results_K_10cov,
                lab = deseq.results_K_10cov$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('AraC-type DNA-binding domain and AraC-containing proteins (AraC) (PDB:1BL0)'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))

deseq.results_K_sigma <- filter(deseq.results_K_10cov, grepl('2Q1Z', COG))
#In some bacterial species, such as Clostridium botulinum, this sigma factor may be necessary for sporulation.[5] 
#Stress response --> sporolation
EnhancedVolcano(deseq.results_K_sigma,
                lab = deseq.results_K_sigma$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))

left_join(deseq.results_K_sigma, deseq.results_Pfam, "taxon")

##Cell division

deseq.results_D_10cov <- deseq.results_D %>% filter(baseMean > 5) %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

sig_early <- deseq.results_D_10cov %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_D_10cov %>%
  filter(log2FoldChange > 0, padj < 0.01)

unique_cell_early <- unique(sig_early)[-which(unique(sig_early$Functional_Category) %in% unique(sig_late$Functional_Category)),]
arrange(unique_cell_early, desc(baseMean)) 
unique_cell_early %>% dplyr::count(COG, sort = TRUE) #Nothing


unique_cell_late <- unique(sig_late)[-which(unique(sig_late$Functional_Category) %in% unique(sig_early$Functional_Category)),]
unique_cell_late %>% dplyr::count(COG, sort = TRUE) #Chromosome segregation ATPase Smc (Smc) (PDB:5XG3) 13


EnhancedVolcano(deseq.results_D_10cov,
                lab = deseq.results_D_10cov$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('AraC-type DNA-binding domain and AraC-containing proteins (AraC) (PDB:1BL0)'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))



#Signal trans#
deseq.results_T_10cov <- deseq.results_T %>% filter(baseMean > 5) %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

sig_early <- deseq.results_T_10cov %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_T_10cov %>%
  filter(log2FoldChange > 0, padj < 0.01)

unique_signalT_early <- unique(sig_early)[-which(unique(sig_early$Functional_Category) %in% unique(sig_late$Functional_Category)),]
arrange(unique_signalT_early, desc(baseMean))

unique_signalT_early %>% dplyr::count(COG, sort = TRUE)

unique_signalT_late <- unique(sig_late)[-which(unique(sig_late$Functional_Category) %in% unique(sig_early$Functional_Category)),]
unique_signalT_late %>% dplyr::count(COG, sort = TRUE) ## Serine/threonine protein kinase (SPS1) (PDB:6G4J)
#Serine/Threonine Kinase receptors play a role in the regulation of cell proliferation, programmed cell death (apoptosis), cell differentiation, and embryonic development.

#Nucleotide-binding universal stress protein,  UspA family (UspA) (PDB:1JMV):
#The primary function of this superfamily is to protect the organism from environmental stress 
#such as exposure to UV light, which may induce genes containing the USP domain in order to 
#protect the DNA and more generally the cell from further damage.[2] During bacterial starvation 
#the USP genes upregulated will often arrest cell growth and promote its metabolism to adapt to sparse nutrients.[2]

deseq.results_T_10cov_stress <- filter(deseq.results_T_10cov, grepl('stress', COG))
EnhancedVolcano(deseq.results_T_10cov_stress,
                lab = deseq.results_T_10cov_stress$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))

left_join(deseq.results_T_10cov_stress, Kofam) #eukaryotic-like serine/threonine-protein kinase

deseq.results_T_10cov_SPS1 <- filter(deseq.results_T_10cov, grepl('SPS1', COG))
EnhancedVolcano(deseq.results_T_10cov_SPS1,
                lab = deseq.results_T_10cov_SPS1$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))

left_join(deseq.results_T_10cov_SPS1, Kofam) #

EnhancedVolcano(deseq.results_T,
                lab = deseq.results_T$COG,
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = c('Bacteriophytochrome (light-regulated signal transduction histidine kinase) (PDB:2VEA)'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))


# deseq.results_T_light <- filter(deseq.results_T, grepl('light', COG))
# EnhancedVolcano(deseq.results_T_light,
#                 lab = deseq.results_T_light$COG,
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
#                 pCutoff = 10e-2,
#                 labSize = 3,
#                 title = '',
#                 legendPosition = 'bottom',
#                 xlim = c(-10,10))



#Signal trans#
deseq.results_S_10cov <- deseq.results_S %>% filter(baseMean > 5) %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

sig_early <- deseq.results_S_10cov %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_S_10cov %>%
  filter(log2FoldChange > 0, padj < 0.01)

unique_Uknown_early <- unique(sig_early)[-which(unique(sig_early$Functional_Category) %in% unique(sig_late$Functional_Category)),]
arrange(unique_Uknown_early, desc(baseMean))
unique_Uknown_early %>% dplyr::count(COG, sort = TRUE)

unique_Uknown_late <- unique(sig_late)[-which(unique(sig_late$Functional_Category) %in% unique(sig_early$Functional_Category)),]
unique_Uknown_late %>% dplyr::count(COG, sort = TRUE) ## nothing


EnhancedVolcano(deseq.results_S,
                lab = deseq.results_S$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('Bacteriophytochrome (light-regulated signal transduction histidine kinase) (PDB:2VEA)'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))

# Uncharacterized conserved protein, PKD repeat domain (PDB:2KZW) 68
deseq.results_S_10cov_stress <- filter(deseq.results_S_10cov, grepl('2KZW', COG))
EnhancedVolcano(deseq.results_S_10cov_stress,
                lab = deseq.results_S_10cov_stress$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))

#str(Kofam)
#Kofam$taxon <- as.character(Kofam$taxon)
#left_join(deseq.results_S_10cov_stress,Kofam, "taxon") ### UNKNOWN PROTEIN! FORGET ABOUT IT!

#deseq.results_S <- filter(deseq.results_S, grepl(' sensing', COG)) No qourum sensing
deseq.results_P_10cov <- deseq.results_P %>% filter(baseMean > 5) %>%
  arrange(pvalue, log2FoldChange, Functional_Category)


sig_early <- deseq.results_P_10cov %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_P_10cov %>%
  filter(log2FoldChange > 0, padj < 0.01)

unique_P_early <- unique(sig_early)[-which(unique(sig_early$Functional_Category) %in% unique(sig_late$Functional_Category)),]
arrange(unique_transc_early, desc(baseMean))
unique_P_early %>% dplyr::count(COG, sort = TRUE)

unique_P_late <- unique(sig_late)[-which(unique(sig_late$Functional_Category) %in% unique(sig_early$Functional_Category)),]
unique_P_late %>% dplyr::count(COG, sort = TRUE) ## nothing

EnhancedVolcano(deseq.results_P_10cov ,
                lab = deseq.results_P_10cov$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('Bacteriophytochrome (light-regulated signal transduction histidine kinase) (PDB:2VEA)'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))

#Outer membrane receptor protein, Fe transport (CirA) (PDB:1BY3) 34
deseq.results_P_1BY3 <- filter(deseq.results_P_10cov, grepl('1BY3|Fe |siderophore', COG))
deseq.results_P_MgtA <- filter(deseq.results_P_10cov, grepl('MgtA', COG)) # Magnesium-transporting ATPase (P-type) (MgtA) (PDB:1IWO) 34

#left_join(deseq.results_P_1BY3, deseq.results_Pfam, "taxon")


EnhancedVolcano(deseq.results_P_1BY3,
                lab = deseq.results_P_1BY3$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))

#Iron uptake dominant function in the early stage - iron is scarch - competition for iron
deseq.results_P_iron <- filter(deseq.results_P_10cov, grepl('Fe |siderophore', COG))
EnhancedVolcano(deseq.results_P_iron,
                lab = deseq.results_P_iron$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))

deseq.results_P_Rhodanese <- filter(deseq.results_P_10cov, grepl('Rhodanese', COG))
EnhancedVolcano(deseq.results_P_Rhodanese,
                lab = deseq.results_P_Rhodanese$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('2-keto-4-pentenoate hydratase (MhpD) (PDB:1SV6)'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))


# #Look for iron siderophore
# deseq.results_P_siderophore <- filter(deseq.results_P, grepl('siderophore', COG)) #No qourum sensing
# 
# EnhancedVolcano(deseq.results_P_siderophore ,
#                 lab = deseq.results_P_siderophore$COG,
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 #selectLab = c('Bacteriophytochrome (light-regulated signal transduction histidine kinase) (PDB:2VEA)'),
#                 pCutoff = 10e-2,
#                 labSize = 3,
#                 title = '',
#                 legendPosition = 'bottom',
#                 xlim = c(-10,10))

deseq.results_X_10cov <- deseq.results_X %>% filter(baseMean > 5) %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

sig_early <- deseq.results_X_10cov %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_X_10cov %>%
  filter(log2FoldChange > 0, padj < 0.01)

unique_mobilome_early <- unique(sig_early)[-which(unique(sig_early$Functional_Category) %in% unique(sig_late$Functional_Category)),]
unique_mobilome_early %>% dplyr::count(COG, sort = TRUE)

unique_mobilome_late <- unique(sig_late)[-which(unique(sig_late$Functional_Category) %in% unique(sig_early$Functional_Category)),]
unique_mobilome_late %>% dplyr::count(COG, sort = TRUE)

EnhancedVolcano(deseq.results_X_10cov,
                lab = deseq.results_X_10cov$COG,
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = c('Transposase'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))

deseq.results_X_Retron <- filter(deseq.results_X_10cov, grepl('Retron', COG)) 
#left_join(deseq.results_X_Retron,Kofam, "taxon")


deseq.results_X_Phage <- filter(deseq.results_X_10cov, grepl('Phage', COG)) 


EnhancedVolcano(deseq.results_X_Retron, #function as anti-phage defense systems
                lab = deseq.results_X_Retron$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('Transposase'),
                pCutoff = 10e-3,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))

EnhancedVolcano(deseq.results_X_Phage, #No difference in phage abundance
                lab = deseq.results_X_Phage$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('Transposase'),
                pCutoff = 10e-3,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))






deseq.results_D_10cov <- deseq.results_D %>% filter(baseMean > 5) %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

sig_early <- deseq.results_D_10cov %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_D_10cov %>%
  filter(log2FoldChange > 0, padj < 0.01)

unique_D_early <- unique(sig_early)[-which(unique(sig_early$Functional_Category) %in% unique(sig_late$Functional_Category)),]
unique_D_early %>% dplyr::count(COG, sort = TRUE) #Nothing really

unique_D_late <- unique(sig_late)[-which(unique(sig_late$Functional_Category) %in% unique(sig_early$Functional_Category)),]
unique_D_late %>% dplyr::count(COG, sort = TRUE) #Nothing


EnhancedVolcano(deseq.results_D_10cov,
                lab = deseq.results_D_10cov$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('Transposase'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))

# deseq.results_D_sporu <- filter(deseq.results_D, grepl('sporulation', COG)) 
# 
# 
# EnhancedVolcano(deseq.results_D_sporu,
#                 lab = deseq.results_D_sporu$COG,
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 #selectLab = c('Transposase'),
#                 pCutoff = 10e-2,
#                 labSize = 3,
#                 title = '',
#                 legendPosition = 'bottom',
#                 xlim = c(-10,10))

deseq.results_C_10cov <- deseq.results_C %>% filter(baseMean > 5) %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

sig_early <- deseq.results_C_10cov %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_C_10cov %>%
  filter(log2FoldChange > 0, padj < 0.01)

unique_C_early <- unique(sig_early)[-which(unique(sig_early$Functional_Category) %in% unique(sig_late$Functional_Category)),]
unique_C_early %>% dplyr::count(COG, sort = TRUE) #Nothing


deseq.results_C_10cov_flavoprotein <- filter(deseq.results_C_10cov, grepl('flavoprotein', COG)) 


EnhancedVolcano(deseq.results_C_10cov_flavoprotein,
                lab = deseq.results_C_10cov_flavoprotein$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('Transposase'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))


# EnhancedVolcano(deseq.results_C,
#                 lab = deseq.results_C$COG,
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 #selectLab = c('Transposase'),
#                 pCutoff = 10e-2,
#                 labSize = 3,
#                 title = '',
#                 legendPosition = 'bottom',
#                 xlim = c(-10,10))
deseq.results_G_10cov <- deseq.results_G %>% filter(baseMean > 5) %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

sig_early <- deseq.results_G_10cov %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_G_10cov %>%
  filter(log2FoldChange > 0, padj < 0.01)

unique_G_early <- unique(sig_early)[-which(unique(sig_early$Functional_Category) %in% unique(sig_late$Functional_Category)),]
unique_G_early %>% dplyr::count(COG, sort = TRUE)

unique_G_late <- unique(sig_late)[-which(unique(sig_late$Functional_Category) %in% unique(sig_early$Functional_Category)),]
unique_G_late %>% dplyr::count(COG, sort = TRUE)

deseq.results_G_10cov_Chitodextrinase <- filter(deseq.results_G_10cov, grepl('Chitodextrinase', COG)) 
deseq.results_G_10cov_GH31 <- filter(deseq.results_G_10cov, grepl('GH31', COG)) 



EnhancedVolcano(deseq.results_G_10cov_Chitodextrinase,
                lab = deseq.results_G_10cov_Chitodextrinase$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('Transposase'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))


EnhancedVolcano(deseq.results_G_10cov_GH31,
                lab = deseq.results_G_10cov_GH31$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('Transposase'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))


deseq.results_O_10cov <- deseq.results_O %>% filter(baseMean > 5) %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

sig_early <- deseq.results_O_10cov %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_O_10cov %>%
  filter(log2FoldChange > 0, padj < 0.01)

unique_O_early <- unique(sig_early)[-which(unique(sig_early$Functional_Category) %in% unique(sig_late$Functional_Category)),]
unique_O_early %>% dplyr::count(COG, sort = TRUE)


deseq.results_O_10cov_Nad <- filter(deseq.results_O_10cov, grepl('NAD', COG)) #Nothing found

EnhancedVolcano(deseq.results_O_10cov,
                lab = deseq.results_O_10cov$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('Transposase'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))



deseq.results_L_10cov <- deseq.results_L %>% filter(baseMean > 5) %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

sig_early <- deseq.results_L_10cov %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_L_10cov %>%
  filter(log2FoldChange > 0, padj < 0.01)

unique_L_early <- unique(sig_early)[-which(unique(sig_early$Functional_Category) %in% unique(sig_late$Functional_Category)),]
unique_L_early %>% dplyr::count(COG, sort = TRUE) #Nothing


deseq.results_W_10cov <- deseq.results_W %>% filter(baseMean > 5) %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

sig_early <- deseq.results_W %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_W %>%
  filter(log2FoldChange > 0, padj < 0.01)

unique_W_early <- unique(sig_early)[-which(unique(sig_early$Functional_Category) %in% unique(sig_late$Functional_Category)),]
unique_W_early %>% dplyr::count(COG, sort = TRUE)

#Nothing found

EnhancedVolcano(deseq.results_W,
                lab = deseq.results_W$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('Transposase'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10))



deseq.results_V_10cov <- deseq.results_V %>% filter(baseMean > 5) %>%
  arrange(pvalue, log2FoldChange, Functional_Category)

sig_early <- deseq.results_V_10cov %>%
  filter(log2FoldChange < 0, padj < 0.01)

sig_late <- deseq.results_V_10cov %>%
  filter(log2FoldChange > 0, padj < 0.01)

unique_V_early <- unique(sig_early)[-which(unique(sig_early$Functional_Category) %in% unique(sig_late$Functional_Category)),]
unique_V_early %>% dplyr::count(COG, sort = TRUE)

#Multidrug resistance efflux pump EmrA (EmrA) (PDB:4TKO) maybe a tedency? 
deseq.results_O_10cov_Vselected <- filter(deseq.results_V_10cov, grepl('EmrA|SalY|AcrB', COG))  # The final v-cat genes 

#left_join(deseq.results_O_10cov_EmrA, Kofam, "taxon")
#left_join(deseq.results_O_10cov_AmpC, Kofam, "taxon")


EnhancedVolcano(deseq.results_O_10cov_Vselected,
                lab = deseq.results_O_10cov_Vselected$COG,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('Transposase'),
                pCutoff = 10e-2,
                labSize = 3,
                title = '',
                legendPosition = 'bottom',
                xlim = c(-10,10),
                ylim = c(0,3))

# EnhancedVolcano(deseq.results_O_10cov_AmpC,
#                 lab = deseq.results_O_10cov_AmpC$COG,
#                 x = 'log2FoldChange',
#                 y = 'padj',
#                 #selectLab = c('Transposase'),
#                 pCutoff = 10e-2,
#                 labSize = 3,
#                 title = '',
#                 legendPosition = 'bottom',
#                 xlim = c(-10,10))


### Summary ###
#Imidazolonepropionase|mycothiol dominant in early (Q)

#RpoE: #In some bacterial species, such as Clostridium botulinum, this sigma factor may be necessary for sporulation.[5] 
#Stress response --> sporolation (K)

## Serine/threonine protein kinase (SPS1) (PDB:6G4J) (T) dominant in LATE
#Serine/Threonine Kinase receptors play a role in the regulation of cell proliferation, programmed cell death (apoptosis), cell differentiation, and embryonic development.

#Nucleotide-binding universal stress protein,  UspA family (UspA) (PDB:1JMV): dominant in early
#The primary function of this superfamily is to protect the organism from environmental stress 
#such as exposure to UV light, which may induce genes containing the USP domain in order to 
#protect the DNA and more generally the cell from further damage.[2] During bacterial starvation 
#the USP genes upregulated will often arrest cell growth and promote its metabolism to adapt to sparse nutrients.[2]

#Outer membrane receptor protein, Fe transport (CirA) (PDB:1BY3) 34 in early (P)
#Magnesium-transporting ATPase (P-type) (MgtA) (PDB:1IWO) 34 in late (P)

#Retrons highly dominant in late function as anti-phage defense systems (X)

#Multidrug resistance efflux pump EmrA (EmrA) (PDB:4TKO) maybe a tendency in early, and SalY only present in early og A (V)

# And betalactamase maybe


#### Final plot ####

# Final df
deseq.results_O_10cov_Vselected$COG_Cat = "Defense mechanisms (V)"
deseq.results_Q_10cov_Qselected$COG_Cat = "Sec. metabolites biosynthesis, transport and catabolism (Q)"
deseq.results_N_EPS_Flag$COG_Cat = "Cell motility (N)"
deseq.results_K_sigma$COG_Cat = "Transcription (K)"
deseq.results_T_10cov_stress$COG_Cat = "Signal transduction mechanisms (T)"
deseq.results_P_iron$COG_Cat = "Inorganic ion transport and metabolism (P)"

Purplerain <- c("#e0b62b","#f25a13", "#2b55e0", "#240785", "#9111ab", "#e33963", "#e35b39", "#e3b039",
                "#580a82", "#060896", "#0d4582", "#c439e3", "#4d44ab")
scales::show_col(Purplerain)


final_sig = rbind(deseq.results_O_10cov_Vselected, deseq.results_Q_10cov_Qselected, 
                  deseq.results_N_EPS_Flag, deseq.results_K_sigma, 
                  deseq.results_T_10cov_stress, deseq.results_P_iron)
write_csv(final_sig, file = "Metagenomics/Data/final_sig.COG.siginearly.csv")

unique(final_sig$COG)

final_sig$color[final_sig$COG_Cat == "Defense mechanisms (V)"] <- "#e33963"
final_sig$color[final_sig$COG_Cat == "Sec. metabolites biosynthesis, transport and catabolism (Q)"] <- "#9111ab"
final_sig$color[final_sig$COG_Cat == "Cell motility (N)"] <- "#0d4582"
final_sig$color[final_sig$COG_Cat == "Transcription (K)"] <- "#e0b62b"
final_sig$color[final_sig$COG_Cat == "Signal transduction mechanisms (T)"] <-  "#2b55e0"
final_sig$color[final_sig$COG_Cat == "Inorganic ion transport and metabolism (P)"] <-  "#580a82"
final_sig$color[final_sig$padj > 0.1] <-  "grey"

key.legend = final_sig$color

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements


names(key.legend)[key.legend == "#e33963"] <- "Defense mechanisms (V)"
names(key.legend)[key.legend == "#9111ab"] <- "Sec. metabolites biosynthesis, transport and catabolism (Q)"
names(key.legend)[key.legend == "#0d4582"] <- "Cell motility (N)"
names(key.legend)[key.legend == "#e0b62b"] <- "Transcription (K)"
names(key.legend)[key.legend ==  "#2b55e0"] <-  "Signal transduction mechanisms (T)"
names(key.legend)[key.legend ==  "#580a82"] <-"Inorganic ion transport and metabolism (P)"
names(key.legend)[key.legend ==  "grey"] <-"Not Significant"




final_sig_plot <- EnhancedVolcano(final_sig,
                lab = NA,
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = c('Transposase'),
                pCutoff = 10e-2,
                colAlpha = 0.8,
                pointSize = 4,
                caption = "",
                subtitle = "",
                labSize = 3,
                colCustom = key.legend,
                title = '',
                axisLabSize = 5,
                legendPosition = 'bottom',
                borderColour = "black",
                #xlim = c(-8,8),
                ylim = c(0,3)) +
  theme_bw(base_size = 12)+
  scale_x_continuous(
    breaks=seq(-8,8, 1))+
  theme(panel.grid.minor = element_blank(), legend.position = "right",
        axis.text = element_text(color = "black", face = "bold"))
library(svglite)
ggsave(file="Metagenomics/Figures/COGs_diff.svg", plot=final_sig_plot, width=8.3, height=5)

 saveRDS(final_sig, "Metagenomics/Figures/final_sig.RData")
