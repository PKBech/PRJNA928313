# Cleaning and subset PS
library("ggplot2")
library("phyloseq")
library(decontam); packageVersion("decontam")
library(dplyr)
library(tidyr)

#### ---- Library sizes - w. negs and true ---- ####
#Library sizes/number of reads in each sample as a function of whether that sample was a sample or neg-control

# put sample_data into a ggplot-friendly data.frame
sample.data = as.data.frame(sample_data(phys))

#Add read counts to sadata as a variable
sample.data$lib.size = sample_sums(phys)

# Order the sadata after the library sizes
sample.data = sample.data[order(sample.data$lib.size),]

#Make index of the order
sample.data$Index = seq(nrow(sample.data))

# Plot library sizes 
ggplot(sample.data, aes(x=Index, y=lib.size, color=as.factor(timepoint))) +
  geom_point() + 
  theme_bw() +
  facet_grid(~element.type)

#### ---- decontam ---- ####
sample_data(phys)$is.neg <- sample_data(phys)$sample.type == "negative"
contamdf.prev <- isContaminant(phys, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

# try with higher decontam level 
contamdf.prev.05 <- isContaminant(phys, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev.05$contaminant)
phys.decon <- prune_taxa(!contamdf.prev.05$contaminant, phys)

#### ---- Filtering and subsetting ---- ####

# Remove eukaryotes, archea, chloroplasts and mitochondria
phys.decon.noEuk <- phys.decon %>%
  subset_taxa( domain == "Bacteria" & family  != "mitochondria" & class   != "Chloroplast") %>%
  prune_taxa(taxa_sums(.) > 0, .)

#### ---- Sample depths (should low read sample out)--- ####

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(phys.decon.noEuk), sampleID = sample_names(phys.decon.noEuk))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

#Max and min
min(sample_sums(phys.decon.noEuk))
max(sample_sums(phys.decon.noEuk)) 

ggplot(sample_sum_df, aes(x= reorder(sampleID, -sum), y= sum))+
  geom_bar(color= "black", fill = "indianred", stat = "identity") +
  ggtitle("Distriution of sample sequencing depth") +
  xlab("sample.ID") +
  ylab("Read counts") +
  theme(axis.text.x = element_text(angle=90, vjust = 0.4))+ 
  geom_hline(yintercept = 1000, linetype = "dashed", size = 1)

sample_sum_df = sample_sum_df[order(sample_sum_df$sum),]

#Rare curve
library("vegan")
library("MicEco")

lol = rcurve(phys.decon.noEuk, subsamp = seq(from = 1000, to = 100000, by = 1000))

ggplot(lol, aes(Reads, Richness, group = Sample, color = as.factor(timepoint))) +
  theme_bw() +
  geom_line() +
  facet_grid(~element.type)

#Remove negatives
phys.decon.noEuk.noneg = phys.decon.noEuk %>%
  subset_samples(sample.type == "sample") %>%
  prune_taxa(taxa_sums(.) > 0, .)
  
#rcurve.nonegs = rcurve(phys.decon.noEuk.noneg, subsamp = seq(from = 10, to = 10000000, by = 100))

ggplot(rcurve.nonegs, aes(Reads, Richness, group = Sample, color = as.factor(timepoint))) +
  theme_bw() +
  geom_line() +
  facet_grid(~element.type)

saveRDS(phys.decon.noEuk.noneg, "AmpliconAnalysis/16S/PS-16S-decon.noEuk.noneg.rds")

# Cleaning low abundant ASVs out - remove based on elbox plot

#Find read cut-off for generas
df.ASV= phys.decon.noEuk.noneg %>%
  psmelt()

df.ASV.sum = df.ASV %>%
  group_by(OTU) %>%
  summarise(Total.abundance = sum(Abundance))

#Make relative abundance
df.ASV.sum$RA = (100*(df.ASV.sum$Total.abundance/sum(df.ASV.sum$Total.abundance)))

#Elbow plot
ggplot(df.ASV.sum, aes(x= reorder(OTU, -Total.abundance), y = (Total.abundance))) +
  geom_bar(color= "black", fill = "indianred", stat = "identity") +
  ggtitle("Distribution of reads on ASVs") +
  xlab("ASV") +
  ylab("Total reads") +
  theme(axis.text.x = element_text(angle=90, vjust = 0.4))

# Take a look in table
ASV.totalsum = df.ASV.sum[order(df.ASV.sum$Total.abundance),]

# Percentage ASVs with reads <100 account for of total reads = 1.355%
100*sum(df.ASV.sum$Total.abundance[df.ASV.sum$Total.abundance<100])/sum(df.ASV.sum$Total.abundance)

#Remove these ASVs from phyloseq object with less than 100 reads
ps.ASV.reduced <- filter_taxa(phys.decon.noEuk.noneg, 
                             function (x) {sum(x > 100) > 0}, prune=TRUE)
save(ps.ASV.reduced, file = "ps.ASV.reduced.RData")
# Reduced ps.genus has 3661 taxa (removed 49'640 generas)

#Check plot with more than 100 reads
df.ASV.reduced  = ps.ASV.reduced  %>%
  psmelt()
df.ASV.reduced.sum = df.ASV.reduced %>%
  group_by(OTU) %>%
  summarise(Total.abundance = sum(Abundance))
ggplot(df.ASV.reduced.sum, 
       aes(x= reorder(OTU, -Total.abundance), y = log10(Total.abundance))) +
  geom_bar(color= "black", fill = "indianred", stat = "identity") +
  ggtitle("Distribution of reads on ASV") +
  xlab("ASV") +
  ylab("Reads in total") +
  theme(axis.text.x = element_text(angle=90, vjust = 0.4)) +
  geom_hline(yintercept = 2, linetype = "dashed", size = 1, color = "red")
# Maybe a little rough?

# Try only remove the ASVs that only count one read in less than two samples
phys.decon.noEuk.noneg.reduced = 
  filter_taxa(phys.decon.noEuk.noneg, function (x) {sum(x > 1) > 1}, prune=TRUE)
# Removed 37'547 ASVs - now 15'754 ASVs

save(phys.decon.noEuk.noneg.reduced, 
     file = "phys.decon.noEuk.noneg.reduced.RData")
