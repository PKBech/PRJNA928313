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
Venn diagram with features distributed between the three successional phases
![Venn diagram with features distributed between the three successional phases](https://github.com/PKBech/PRJNA928313/blob/main/Metabolome/Figures/Figure_4A.png)


### Removing primary metabolites

Filter all features that is less than 300 m/z, and features annotated as fatty acids to look for presence of secondary metabolites dynamics
```
Phylo_metabolome <- subset_taxa(Phylo_metabolome, mass > 300 & NPC.pathway != "Fatty acids")
Phylo_metabolome <- subset_taxa(Phylo_metabolome, NPC.pathway != "Fatty acids" & NPC.superclass != "Lignans")
Phylo_metabolome_NRPs <- subset_taxa(Phylo_metabolome, NPC.pathway == "Amino acids and Peptides")
```

### Investigating Zero-inflation
```
zPatterns(as.data.frame(otu_table(Phylo_metabolome)), label = 0)
```
More than 80% zeros are observed across all samples in the metabolome

### Normalization before beta-diversity analysis
We use center-log-transformation of the data before counting the distance matrix (Sisk-Hackworth & Kelley 2020; https://doi.org/10.1093/nargab/lqaa079)
```
#Impute pseudo count for zero's
Phylo_metabolome_pseudo_counts <- cmultRepl(as.data.frame(otu_table(Phylo_metabolome)), z.warning = 0.9)
#Phylo_metabolome_NRPs_pseudo_counts <- cmultRepl(as.data.frame(otu_table(Phylo_metabolome_NRPs)), z.warning = 0.9)

#Centered log ratio transform
library(compositions)
Phylo_metabolome_pseudo_counts_clr_dat <- clr(Phylo_metabolome_pseudo_counts)


#create new phyloseq object
Phylo_metabolome_pseudo_counts_clr <- phyloseq(otu_table(Phylo_metabolome_pseudo_counts_clr_dat, taxa_are_rows = TRUE), 
                             tax_table(as.matrix(metabolome_annotations_table_1)), sample_data(Phylo_metabolome))


#Plot PCoA with physeq package
main_col <- c("#240785",  "#f21395", "#e0b62b")
Phylo_metabolome_pseudo_counts_clr.ord_PCoA <- ordinate(Phylo_metabolome_pseudo_counts_clr, "PCoA", "euclidean")
p2 = plot_ordination(Phylo_metabolome_pseudo_counts_clr, Phylo_metabolome_pseudo_counts_clr.ord_PCoA, type="samples", color="Phase", axes = c(1,2)) 

pcoa_pseudo_counts_clr_euclidean_p =p2 + geom_point(size=7.5) + 
  theme_minimal() + scale_fill_manual(values = main_col) +
  scale_color_manual(values = main_col)

pcoa_pseudo_counts_clr_euclidean_p



#Beta-disper ####
## Calculate multivariate dispersions 
bdisp_metabolome_pseudo_counts_clr <- betadisper(vegdist(t(Phylo_metabolome_pseudo_counts_clr_dat), "euclidean"), 
                                                 factor(paste(data.frame(sample_data(Phylo_metabolome_pseudo_counts_clr))$Phase,sep=" ")))

plot(bdisp_metabolome_pseudo_counts_clr)

bdisp_metabolome_pseudo_counts_clr_day <- betadisper(vegdist(t(Phylo_metabolome_pseudo_counts_clr_dat), "euclidean"), 
                                                 factor(paste(data.frame(sample_data(Phylo_metabolome_pseudo_counts_clr))$Day,sep=" ")))



## Permutation test for F to test if betadisp is significant
permutest(bdisp_metabolome_pseudo_counts_clr, pairwise = TRUE, permutations = 999)
permutest(bdisp_metabolome_pseudo_counts_clr_day, pairwise = TRUE, permutations = 999)


## Tukey's Honest Significant Differences
bdisp.HSD <- TukeyHSD(bdisp_metabolome_pseudo_counts_clr)
plot(bdisp.HSD)

##Beta-disper it significant between late vs. early and peak and late.
#Therefore not sure if the significant different between group running PERMANOVA is
#caused by the high variance in the late phase or differences explained by time


#PERMANOVA (Adonis) #####
Metabolome_pseudo_counts_clr_dat_dist_eucledian = as.matrix(vegdist(t(data.frame(otu_table(Phylo_metabolome_pseudo_counts_clr))), "euclidean"))
metadata_pseudo_counts_clr_dat_dist_eucledian <- data.frame(sample_data(Phylo_metabolome_pseudo_counts_clr))
metadata_pseudo_counts_clr_dat_dist_eucledian$Phase <- factor(metadata_pseudo_counts_clr_dat_dist_eucledian$Phase, levels=c("Early", "Peak", "Late"))


PERMANOVA <- adonis2(Metabolome_pseudo_counts_clr_dat_dist_eucledian ~ Phase * Replicate, metadata_pseudo_counts_clr_dat_dist_eucledian)
PERMANOVA


#PCOA plot with vegan/betadisper #####
centroids.metabolome = data.frame(bdisp_metabolome_pseudo_counts_clr$centroids[,1:2])
centroids.metabolome$Phase = rownames(centroids.metabolome)
points.metabolome = data.frame(bdisp_metabolome_pseudo_counts_clr$vectors[,1:2], Phase = metadata_pseudo_counts_clr_dat_dist_eucledian$Phase)

PCoA1 <- bdisp_metabolome_pseudo_counts_clr$eig[1]/sum(bdisp_metabolome_pseudo_counts_clr$eig)*100
PCoA1
PCoA2 <- bdisp_metabolome_pseudo_counts_clr$eig[2]/sum(bdisp_metabolome_pseudo_counts_clr$eig)*100
PCoA1
PCoA1 + PCoA2

main_col <- c("#240785", "#e0b62b", "#f21395")

PCoA.metabolome.plot <- ggplot(centroids.metabolome, aes(PCoA1, PCoA2, color = Phase)) +
  theme_bw(base_size = 12)+
  geom_point(data = points.metabolome, aes(PCoA1, PCoA2, color = Phase), size = 1.5, alpha = 0.5) +
  geom_segment(aes(x = centroids.metabolome[1,1] , y = centroids.metabolome[1,2] , xend = PCoA1, yend = PCoA2), 
               data = points.metabolome[points.metabolome$Phase == "Early",], size =1, alpha = 0.3, colour = "#240785") +
  geom_segment(aes(x =  centroids.metabolome[2,1] , y = centroids.metabolome[2,2], xend = PCoA1, yend = PCoA2), 
               data = points.metabolome[points.metabolome$Phase == "Late",], size =1, alpha = 0.3, colour = "#f21395" ) +
  geom_segment(aes(x =  centroids.metabolome[3,1] , y = centroids.metabolome[3,2], xend = PCoA1, yend = PCoA2), 
               data = points.metabolome[points.metabolome$Phase == "Peak",], size =1, alpha = 0.3, colour = "#e0b62b") +
  geom_point(size = 4)+
  stat_ellipse(data = points.metabolome, aes(group=Phase, color = Phase), size = 1, alpha = 0.3)+
  labs(x = "PCoA1 [23.86%]", y = "\nPCoA2 [10.84%]")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = main_col) + 
  scale_color_manual(values = main_col)

PCoA.metabolome.plot

ggsave(file="Metabolome/Figures/Figure_4B.svg", plot=PCoA.metabolome.plot, width=4.9, height=3.5)
```

Beta-diversity
![PCoA](https://github.com/PKBech/PRJNA928313/blob/main/Metabolome/Figures/Figure_4B.png)






