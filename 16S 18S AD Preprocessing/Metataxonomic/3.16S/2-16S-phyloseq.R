## 2. Creation of Phyloseq object ###


#### ---- Import data ---- ####
library("dada2")
library("phyloseq")
library("tidyr")

#Cleaned metadata
metadata <- read.csv("AmpliconAnalysis/16S/metadata-16s.csv", row.names = 1)

# load ASV table
seqtab <- readRDS("Dada2/16S/seqtab.nochim.rds") # file from dada2 

#Load taxanomy
#taxa_silva138 <- readRDS("Data/tax_silva138.rds") # file from assign taxonomy in dada2 october2020
taxa_silva138.1 <- readRDS("Dada2/16S/tt.plus.rds") # file from assign taxonomy and add species in dada2 with 138.1 version from september 2021

#### ---- Remove short length sequences ---- ####

# Distribution length of asv table
table(nchar(getSequences(seqtab)))
hist(nchar(getSequences(seqtab)))
abline(v=380) # cut off for short sequences

# Remove short length sequences from asv table
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 381:444] # remove everything below 381 bp length, they are in columns
table(nchar(getSequences(seqtab2)))
hist(nchar(getSequences(seqtab2)))

# Distribution of length of taxonomy and remove short length seqs (remove same length sequences in the taxonomy file) change to tax file you want
table(nchar(getSequences(taxa_silva138.1)))
hist(nchar(getSequences(taxa_silva138.1)))
taxa_silva138.2 <- taxa_silva138.1[nchar(rownames(taxa_silva138.1)) %in% 381:444,] # here they are rows
table(nchar(getSequences(taxa_silva138.2)))

#### ---- Standard files + add ASV names ---- ####

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab2)
asv_headers <- vector(dim(seqtab2)[2], mode="character")

for (i in 1:dim(seqtab2)[2]) {
  asv_headers[i] <- paste(">16S_ASV", i, sep="_")
}

# Making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
#write(asv_fasta, "AmpliconAnalysis/16S/16S_ASVs_sequences.fa")

# count table:
asv_tab <- t(seqtab2)
row.names(asv_tab) <- sub(">", "", asv_headers)
#write.table(asv_tab, "AmpliconAnalysis/16S/ASVs_tab.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
colnames(taxa_silva138.2) <- ranks
rownames(taxa_silva138.2) <- gsub(pattern=">", replacement="", x=asv_headers)

#write.table(taxa_silva138.2, "Data/ASVs_taxa_silva138.2.tsv", sep = "\t", quote=F, col.names=NA)

#### ---- Phyloseq object ---- ####

ASV = otu_table(asv_tab, taxa_are_rows = TRUE)
TAX = tax_table(taxa_silva138.2)
MET = sample_data(metadata)
SEQS <- Biostrings::readDNAStringSet("AmpliconAnalysis/16S/16S_ASVs_sequences.fa", format = "fasta")


taxa_names(TAX)
taxa_names(ASV)

sample_names(ASV)
sample_names(MET)


phys <- phyloseq(ASV, TAX, MET, SEQS)
phys = phys %>%
  prune_taxa(taxa_sums(.) > 0, .) # 78767 ASVs


sample_names(phys)
rank_names(phys)
sample_variables(phys)

#Save phyloseqobject 
save(phys, file = "phys.RData")

