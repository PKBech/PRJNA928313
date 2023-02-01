## 2. Creation of Phyloseq object with 18S data ###


#### ---- Import data ---- ####
library("dada2")
library("phyloseq")
library("tidyr")

#Cleaned metadata
metadata <- read.csv("Metadata/metadata-16s.csv", row.names = 1)


# load ASV table
seqtab <- readRDS("Dada2/18S/seqtab.nochim.rds") # file from dada2 

#Load taxanomy
taxa_silva138.1 <- readRDS("Dada2/18S/tt.silva.18S.rds") # file from assign taxonomy and add species in dada2 with 138.1 version from september 2021
taxa_rdp <- readRDS("Dada2/18S/tt.rdp.18S.rds") # file from assign taxonomy and add species in dada2 with 138.1 version from september 2021

#### ---- Remove short length sequences ---- ####

# Distribution length of asv table
table(nchar(getSequences(seqtab)))
hist(nchar(getSequences(seqtab)))

#### ---- Standard files + add ASV names ---- ####

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")

for (i in 1:dim(seqtab)[2]) {
  asv_headers[i] <- paste(">18S_ASV", i, sep="_")
}

# Making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
#write(asv_fasta, "AmpliconAnalysis/18S/ASVs_18S_fasta.fa")

# count table:
asv_tab <- t(seqtab)
row.names(asv_tab) <- sub(">", "", asv_headers)
#write.table(asv_tab, "AmpliconAnalysis/18S/ASVs_18S_tab.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("Kingdom", "phylum", "class", "order", "family", "genus", "species", "strain")
colnames(taxa_silva138.1) <- ranks
rownames(taxa_silva138.1) <- gsub(pattern=">", replacement="", x=asv_headers)

#write.table(taxa_silva138.2, "Data/ASVs_taxa_silva138.2.tsv", sep = "\t", quote=F, col.names=NA)

#### ---- Phyloseq object ---- ####

ASV = otu_table(asv_tab, taxa_are_rows = TRUE)
TAX = tax_table(taxa_silva138.1)
MET = sample_data(metadata)
REF <- refseq(Biostrings::readDNAStringSet("AmpliconAnalysis/18S/ASVs_18S_fasta.fa", format = "fasta"))

taxa_names(TAX)
taxa_names(ASV)
taxa_names(REF)

sample_names(ASV)
sample_names(MET)


phys.18s <- phyloseq(ASV, TAX, MET, REF)
phys.18s = phys.18s %>%
  prune_taxa(taxa_sums(.) > 0, .) # 16726 ASVs


sample_names(phys.18s)
rank_names(phys.18s)
sample_variables(phys.18s)

#Save phyloseqobject 
saveRDS(phys.18s, "AmpliconAnalysis/18S/phys.18s.rds")

# Import again
#phys <- readRDS("Data/phys.rds")

