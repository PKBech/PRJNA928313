# Analysis of untargeted LC-MS metabolome data collected across the succession of marine biofilm 

## Pre-processing of LC-MS/MS data
Raw HRMS data was converted to mzML using MSConvert (ProteoWizard) and preprocessed using MZmine 3 [74]. Molecular networking (MN) was all completed within the GNPS platform [75], which includes: Feature Based Molecular Networking [76] and Ion Identity Molecular Networking [77]. Ion Identity MN workflow can be found here: (To be filled out). Formula predictions, compound class annotations and NPclass annotations were run through SIRIUS/v5.5.7 and include: SIRIUS, ZODIAC, CSI:FingerID and CANOPUS. All data run in this study can be found in GNPS MassIVE (To be filled out). 

Formula predictions (SIRIUS) and the feature LC-MS table was imported and analysed further in R. 

### Load dependencies

```
library("dplyr")
library("phyloseq")
library(tidyverse)
library("devtools")
library("Biostrings")
library(ggplot2)
library(vegan)
library(ggvenn)
library(zCompositions)
library(compositions)
library(ggbeeswarm)
library(FSA)
```
### Clean-up input tables and create phyloseq object 
```
##Load feature table
#'Row ID' in the quant table and 'ID in the SIRIUS table are identical and can be linked.
quat_table <- read.csv("Metabolome/Data/Mzmine_GNPS output_1000 filtered_quant.csv", sep = ",", header = TRUE)
quat_table <- as.data.frame(quat_table)
str(quat_table)
quat_table <- quat_table[, -100]

##Load feature table with predictions
SERIUS_table <- read.csv("Metabolome/Data/SIRIUS_predictions.csv", sep = ";", header = TRUE)
SERIUS_table <- as.data.frame(SERIUS_table)

str(SERIUS_table)
colnames(SERIUS_table)[1] <- colnames(quat_table)[1]
#Join tables based in "row.ID"
met_table <- as.data.frame(left_join(quat_table, SERIUS_table)) 
met_table$row.m.z <- round(met_table$row.m.z, digits = 5)
#Remove metabolite "128.1066" (empty row)
met_table <- met_table[!met_table$row.m.z == "128.10663",]


#Wrangle data to get a featuretable, annotation table and metadata for import to phyloseq
#rownames(met_table) <- make.unique(as.character(met_table$row.m.z))
met_table$row.m.z[362] = 482.324691
rownames(met_table) <- (as.character(met_table$row.m.z))

colnames(met_table)[14:99]
#Make feature dataframe
feature_abundance <- met_table[,14:99]
#remove blancks 
feature_abundance <- feature_abundance[rowSums(feature_abundance[,c(1:3, 6, 66)])==0,]
#remove blanck cols
feature_abundance <- feature_abundance[,-c(1:3, 6, 66)]


#Make metadata table 
metabolome_sample_data <- read.csv("Metabolome/Data/sample_data2.csv", sep = ",", header = TRUE)
#match feature table to metabolome_sample_data
metabolome_sample_data <- metabolome_sample_data[match(colnames(feature_abundance), metabolome_sample_data$sample_name),]
#add readable names to metadata table
rownames(metabolome_sample_data) <- metabolome_sample_data$sample_name_2

#Get sample names from metadata to feature table
colnames(feature_abundance) <- metabolome_sample_data$sample_name_2

#Make annotation ('metabolome_tax_table') dataframe ("phyloseq tax table")
str(met_table[,100:118])
metabolome_annotations_table <- met_table[,100:118]
#match rownames with feature table
metabolome_annotations_table <- metabolome_annotations_table[rownames(feature_abundance) %in% rownames(metabolome_annotations_table),]
#Change NA's in annotation table to "Unclassified"

colnames(metabolome_annotations_table)

metabolome_annotations_table$mass <- rownames(metabolome_annotations_table)

#Exclude proberbilities for now
metabolome_annotations_table_1 <- metabolome_annotations_table %>% dplyr::select(mass, molecularFormula,adduct,NPC.pathway,
                                        NPC.superclass, NPC.class, ClassyFire.most.specific.class,
                                        ClassyFire.level.5,ClassyFire.subclass, ClassyFire.class,
                                        ClassyFire.superclass, ClassyFire.all.classifications) %>% dplyr::mutate_all(~ifelse(is.na(.),"Unclassified",.))



##Create Phyloseq Object
Phylo_metabolome <- phyloseq(otu_table(feature_abundance, taxa_are_rows = TRUE), 
                       tax_table(as.matrix(metabolome_annotations_table_1)), sample_data(metabolome_sample_data))


summary(as.data.frame(tax_table(Phylo_metabolome))$NPC.pathway == "Unclassified")
(data.frame(sample_data(Phylo_metabolome))) %>% dplyr::group_by(Phase) %>% dplyr::summarise(., n=n())

# Remove outliers 
outliers <- names(which(colSums(feature_abundance) <= 50000))
Phylo_metabolome <-  subset_samples(Phylo_metabolome, !sample_name_2 %in% outliers) 

#Filter samples in ps so they match rownames of AD amplicon samples processed in Figure 1 scripts
AD_filter_norm_toMetabolome <- read.csv("Metabolome/Data/AD_filter_norm_toMetabolome.csv", sep = ",", header = TRUE)
rownames(AD_filter_norm_toMetabolome) = AD_filter_norm_toMetabolome$X
Phylo_metabolome <- subset_samples(Phylo_metabolome, sample_name_2 %in% rownames(AD_filter_norm_toMetabolome)) 

```
### Venn diagram 
```
Phylo_metabolome_dat <- data.frame(t(otu_table(Phylo_metabolome)))
Phylo_metabolome_dat$sample_name_2 <- rownames(Phylo_metabolome_dat)
Phylo_metabolome_metadat <- data.frame(sample_data(Phylo_metabolome))
Phylo_metabolome_metadat_2 <- left_join(Phylo_metabolome_dat, Phylo_metabolome_metadat)
colnames(Phylo_metabolome_metadat_2)[1:795] <- rownames(otu_table(Phylo_metabolome))

#Change to long format
Phylo_metabolome_metadat_2_long <- Phylo_metabolome_metadat_2 %>% 
  gather(mass, abundance, -c(sample_name_2, sample_name, Replicate, Day, Time_Cat, Phase, Phase_major))
str(Phylo_metabolome_metadat_2)

#Clean up 
Phylo_metabolome_metadat_2_long <- Phylo_metabolome_metadat_2_long %>% 
  dplyr::filter(abundance != 0)  %>% 
  dplyr::group_by(Phase, mass) %>% 
  dplyr::summarize(n = n(), mean_abundance = mean(abundance), total_abundance = sum(abundance))


#Find unique masses per group: "Phase" 
library(ggvenn)
d<-list(Early = filter(Phylo_metabolome_metadat_2_long, Phase == "Early")$mass, 
        Peak = filter(Phylo_metabolome_metadat_2_long, Phase == "Peak")$mass, 
        Late = filter(Phylo_metabolome_metadat_2_long, Phase == "Late")$mass)


str(d)
##venn plot ######
Figure_4A <- ggvenn(
  d, 
  fill_color = c("#240785", "#e0b62b", "#f21395"),
  stroke_size = 0.5, set_name_size = 7
)
Figure_4A

ggsave(file="Metabolome/Figures/Figure_4A.svg", plot=Figure_4A, width=6.5, height=6.5)

```

![Venn diagram with features distributed between the three successional phases](https://github.com/PKBech/PRJNA928313/blob/main/Metabolome/Figures/Figure_4A.png)

