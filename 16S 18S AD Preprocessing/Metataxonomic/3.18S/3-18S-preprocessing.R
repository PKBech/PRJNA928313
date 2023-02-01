# Cleaning and subset PS

library("ggplot2")
library("phyloseq")
library(decontam); packageVersion("decontam")
library("tidyr")

phys.18s = readRDS("AmpliconAnalysis/18S/phys.18s.rds")

#### ---- Library sizes - w. negs and true ---- ####
#Library sizes/number of reads in each sample as a function of whether that sample was a sample or neg-control

# put sample_data into a ggplot-friendly data.frame
sample.data = as.data.frame(sample_data(phys.18s))

#Add read counts to sadata as a variable
sample.data$lib.size = sample_sums(phys.18s)

# Order the sadata after the library sizes
sample.data = sample.data[order(sample.data$lib.size),]

#Make index of the order
sample.data$Index = seq(nrow(sample.data))

# Plot library sizes 
ggplot(sample.data, aes(x=Index, y=lib.size, color=as.factor(sample.type))) +
  geom_point() + 
  theme_bw() +
  facet_grid(~element.type)

#### ---- decontam ---- ####
sample_data(phys.18s)$is.neg <- sample_data(phys.18s)$sample.type == "negative"
contamdf.prev <- isContaminant(phys.18s, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

# try with higher decontam level 
contamdf.prev.05 <- isContaminant(phys.18s, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev.05$contaminant)
phys.18s.decon <- prune_taxa(!contamdf.prev.05$contaminant, phys.18s)

#### ---- Filtering and subsetting ---- ####

# Remove eukaryotes, archea, chloroplasts and mitochondria
# phys.decon.noEuk <- phys.18s.decon %>%
#   subset_taxa( domain == "Bacteria" & family  != "mitochondria" & class   != "Chloroplast") %>%
#   prune_taxa(taxa_sums(.) > 0, .)

#### ---- Sample depths (should low read sample out)--- ####

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(phys.18s.decon), sampleID = sample_names(phys.18s.decon))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

#Max and min
min(sample_sums(phys.18s.decon))
max(sample_sums(phys.18s.decon)) 

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

rcurve.18s = rcurve(phys.18s.decon, 
                    subsamp = seq(from = 1, to = 100000, by = 1000))

ggplot(rcurve.18s[rcurve.18s$element.type != "bryozoan",], 
       aes(Reads, Richness, group = Sample, color = as.factor(timepoint))) +
  theme_bw() +
  geom_line() +
  facet_grid(~element.type)

#Remove negatives
phys.18s.decon.noneg = phys.18s.decon %>%
  subset_samples(sample.type == "sample") %>%
  prune_taxa(taxa_sums(.) > 0, .)

# rcurve.18s.nonegs = rcurve(phys.18s.decon.noneg, 
#                            subsamp = seq(from = 0, to = 10000000, by = 10000))

# ggplot(rcurve.nonegs, aes(Reads, Richness, group = Sample, color = as.factor(timepoint))) +
#   theme_bw() +
#   geom_line() +
#   facet_grid(~element.type)

saveRDS(phys.18s.decon.noneg, "AmpliconAnalysis/18S/PS-18S-decon.noneg.rds")
