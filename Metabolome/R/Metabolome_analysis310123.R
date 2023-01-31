#Packages
library("dplyr")
library("phyloseq")
library("ADImpute")
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

#Create dataframe, clean-up and phyloseq object#####
##Load feature table
#'Row ID' in the quant table and 'ID in the SIRIUS table are identical and can be linked.
quat_table <- read.csv("Metabolome/Data/Mzmine_GNPS output_1000 filtered_quant.csv", sep = ",", header = TRUE)
quat_table <- as.data.frame(quat_table)
str(quat_table)
quat_table <- quat_table[, -100]


quat_table[round(quat_table$row.m.z, digits = 4) %in% c(244.1696), ]

SERIUS_table <- read.csv("Metabolome/Data/SIRIUS_predictions.csv", sep = ";", header = TRUE)
SERIUS_table <- as.data.frame(SERIUS_table)

str(SERIUS_table)
colnames(SERIUS_table)[1] <- colnames(quat_table)[1]
#Join tables based in "row.ID"
met_table <- as.data.frame(left_join(quat_table, SERIUS_table)) 
met_table$row.m.z <- round(met_table$row.m.z, digits = 5)
#Remove metabolite "128.1066" (empty row)
met_table <- met_table[!met_table$row.m.z == "128.10663",]
# met_table[round(met_table$row.m.z) == 244,]
# met_table[met_table$row.m.z = "430.2443",]


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

dim(as.data.frame(tax_table(Phylo_metabolome)))[1]

##Distribution of the different NPC pathways across the metabolome
as.data.frame(tax_table(Phylo_metabolome)) %>% 
  dplyr::group_by(NPC.pathway) %>% 
  dplyr::summarise(total = n()) %>%
  mutate(total_pct = total/sum(total)*100)

#Venn diagram ####

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

# 
# tmp <- venn(d, intersection=TRUE)
# isect <- attr(tmp, "intersection")
# Early_unique <- isect$Early
# Peak_unique <- isect$Peak
# Late_unique <- isect$Late
# Early_Peak_unique <- isect$`Early:Peak`
# 
# subset_taxa(Phylo_metabolome, mass %in% c(Early_unique, Peak_unique, Late_unique)) %>% 
#   #transform_sample_counts(., function(x) x/sum(x)) %>%  
#   plot_bar(., "sample_name_2", fill="NPC.superclass") + facet_grid(Phase~Day, scales="free")
# 
# #Early
# subset_taxa(Phylo_metabolome, mass %in% c(Early_unique, Peak_unique, Late_unique)) %>% 
#   subset_samples(., Phase == "Early") %>%
#   filter_taxa(., function (x) {sum(x > 0) > 0}, prune=TRUE) %>%
#   #transform_sample_counts(., function(x) x/sum(x)) %>%  
#   plot_bar(., "sample_name_2", fill="NPC.pathway") + facet_grid(Phase~Day, scales="free")
# 
# #Peak
# subset_taxa(Phylo_metabolome, mass %in% c(Early_unique, Peak_unique, Late_unique) ) %>% 
#   subset_samples(., Phase == c("Peak") ) %>%
#   #filter_taxa(., "mass" >= 300) %>%
#   filter_taxa(., function (x) {sum(x > 0) > 0}, prune=TRUE) %>%
#   #transform_sample_counts(., function(x) x/sum(x)) %>%  
#   plot_bar(., "sample_name_2", fill="NPC.pathway") + facet_grid(Phase~Day, scales="free")
# 
# 
# #Late
# subset_taxa(Phylo_metabolome, mass %in% c(Early_unique, Peak_unique, Late_unique) ) %>% 
#   subset_samples(., Phase == c("Late") ) %>%
#   #filter_taxa(., "mass" >= 300) %>%
#   filter_taxa(., function (x) {sum(x > 0) > 0}, prune=TRUE) %>%
#   #transform_sample_counts(., function(x) x/sum(x)) %>%  
#   plot_bar(., "sample_name_2", fill="NPC.pathway") + facet_grid(Phase~Day, scales="free")
# 



#FILTERING ######
#Filter all features that is less than 300 m/z, and features annotated as fatty acids to look for presence of secondary metabolites dynamics
Phylo_metabolome <- subset_taxa(Phylo_metabolome, mass > 300 & NPC.pathway != "Fatty acids")
Phylo_metabolome <- subset_taxa(Phylo_metabolome, NPC.pathway != "Fatty acids" & NPC.superclass != "Lignans")
Phylo_metabolome_NRPs <- subset_taxa(Phylo_metabolome, NPC.pathway == "Amino acids and Peptides")


#Zero's across dataset #####
library(zCompositions)

zPatterns(as.data.frame(otu_table(Phylo_metabolome)), label = 0)
#zPatterns(as.data.frame(otu_table(Phylo_metabolome_NRPs)), label = 0)



#Clr transformation ####
# Laura Sisk-Hackworth, Scott T Kelley, 
# An application of compositional data analysis to multiomic time-series data, 
# NAR Genomics and Bioinformatics, Volume 2, Issue 4, December 2020, lqaa079, 
# https://doi.org/10.1093/nargab/lqaa079


#Impute pseudo count for zero's
Phylo_metabolome_pseudo_counts <- cmultRepl(as.data.frame(otu_table(Phylo_metabolome)), z.warning = 0.9)
#Phylo_metabolome_NRPs_pseudo_counts <- cmultRepl(as.data.frame(otu_table(Phylo_metabolome_NRPs)), z.warning = 0.9)

#Centered log ratio transform
library(compositions)
Phylo_metabolome_pseudo_counts_clr_dat <- clr(Phylo_metabolome_pseudo_counts)

# #heatmap
# Phylo_metabolome_pseudo_counts_unfilt_clr_dat <- as.data.frame(clr(Phylo_metabolome_pseudo_counts_unfilt))
# 
# 
# Phylo_metabolome_pseudo_counts_unfilt_clr_dat_2 <- as.data.frame(lapply(as.data.frame(Phylo_metabolome_pseudo_counts_clr_dat), function(x) ifelse(x < 0, 0, x)))
# rownames(Phylo_metabolome_pseudo_counts_unfilt_clr_dat_2) <- rownames(Phylo_metabolome_pseudo_counts_clr_dat)
# colnames(Phylo_metabolome_pseudo_counts_unfilt_clr_dat_2) <- colnames(Phylo_metabolome_pseudo_counts_clr_dat)
# 
# heatmap(as.matrix(Phylo_metabolome_pseudo_counts_unfilt_clr_dat_2), scale = "none") 

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

# p2 = plot_ordination(Phylo_metabolome_pseudo_counts_clr, Phylo_metabolome_pseudo_counts_clr.ord_PCoA, type="samples", color="Replicate", axes = c(1,2))
# 
# pcoa_pseudo_counts_clr_euclidean_p = p2 + geom_point(size=7.5) +
#   theme_minimal() #+ scale_fill_manual(values = main_col) +
#   #scale_color_manual(values = main_col)
# 
# pcoa_pseudo_counts_clr_euclidean_p


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

#Alpha-diversity #####
Observed_total_metabolome_dat <- data.frame(t(estimateR((round(t(otu_table(Phylo_metabolome)) ))))) #estimate richness
colnames(Observed_total_metabolome_dat)[1] <- "Obs_features"
Observed_total_metabolome_dat <- Observed_total_metabolome_dat[-c(2:5)]
Observed_total_metabolome_dat$NPC.pathway <- rep("All NPC Pathway groups")
Observed_total_metabolome_dat$sample_name_2 <- rownames(Observed_total_metabolome_dat)
Observed_total_metabolome_dat <- left_join(Observed_total_metabolome_dat, data.frame(sample_data(Phylo_metabolome)))

#Alkaloids
Phylo_metabolome_Alkaloids <- subset_taxa(Phylo_metabolome, NPC.pathway == "Alkaloids")

Observed_metabolome_Alkaloids_dat <- data.frame(t(estimateR((round(t(otu_table(Phylo_metabolome_Alkaloids)) ))))) #estimate richness
colnames(Observed_metabolome_Alkaloids_dat)[1] <- "Obs_features"
Observed_metabolome_Alkaloids_dat <- Observed_metabolome_Alkaloids_dat[-c(2:5)]
Observed_metabolome_Alkaloids_dat$NPC.pathway <- rep("Alkaloids")
Observed_metabolome_Alkaloids_dat$sample_name_2 <- rownames(Observed_metabolome_Alkaloids_dat)
Observed_metabolome_Alkaloids_dat <- left_join(Observed_metabolome_Alkaloids_dat, data.frame(sample_data(Phylo_metabolome)))


#Amino_acids_and_Peptides
Phylo_metabolome_Amino_acids_and_Peptides <- subset_taxa(Phylo_metabolome, NPC.pathway == "Amino acids and Peptides")

Observed_metabolome_Amino_acids_and_Peptides_dat <- data.frame(t(estimateR((round(t(otu_table(Phylo_metabolome_Amino_acids_and_Peptides)) ))))) #estimate richness
colnames(Observed_metabolome_Amino_acids_and_Peptides_dat)[1] <- "Obs_features"
Observed_metabolome_Amino_acids_and_Peptides_dat <- Observed_metabolome_Amino_acids_and_Peptides_dat[-c(2:5)]
Observed_metabolome_Amino_acids_and_Peptides_dat$NPC.pathway <- rep("Amino acids and Peptides")
Observed_metabolome_Amino_acids_and_Peptides_dat$sample_name_2 <- rownames(Observed_metabolome_Amino_acids_and_Peptides_dat)
Observed_metabolome_Amino_acids_and_Peptides_dat <- left_join(Observed_metabolome_Amino_acids_and_Peptides_dat, data.frame(sample_data(Phylo_metabolome)))


#Fatty_acids
# Phylo_metabolome_Fatty_acids <- subset_taxa(Phylo_metabolome, NPC.pathway == "Fatty acids")
# 
# Observed_metabolome_Fatty_acids_dat <- data.frame(t(estimateR((round(t(otu_table(Phylo_metabolome_Fatty_acids)) ))))) #estimate richness
# colnames(Observed_metabolome_Fatty_acids_dat)[1] <- "Obs_features"
# Observed_metabolome_Fatty_acids_dat <- Observed_metabolome_Fatty_acids_dat[-c(2:5)]
# Observed_metabolome_Fatty_acids_dat$NPC.pathway <- rep("Fatty acids")
# Observed_metabolome_Fatty_acids_dat$sample_name_2 <- rownames(Observed_metabolome_Fatty_acids_dat)
# Observed_metabolome_Fatty_acids_dat <- left_join(Observed_metabolome_Fatty_acids_dat, data.frame(sample_data(Phylo_metabolome)))
# 

#Polyketides
Phylo_metabolome_Polyketides <- subset_taxa(Phylo_metabolome, NPC.pathway == "Polyketides")

Observed_metabolome_Polyketides_dat <- data.frame(t(estimateR((round(t(otu_table(Phylo_metabolome_Polyketides)) ))))) #estimate richness
colnames(Observed_metabolome_Polyketides_dat)[1] <- "Obs_features"
Observed_metabolome_Polyketides_dat <- Observed_metabolome_Polyketides_dat[-c(2:5)]
Observed_metabolome_Polyketides_dat$NPC.pathway <- rep("Polyketides")
Observed_metabolome_Polyketides_dat$sample_name_2 <- rownames(Observed_metabolome_Polyketides_dat)
Observed_metabolome_Polyketides_dat <- left_join(Observed_metabolome_Polyketides_dat, data.frame(sample_data(Phylo_metabolome)))


#Shikimates_and_Phenylpropanoids
Phylo_metabolome_Shikimates_and_Phenylpropanoids <- subset_taxa(Phylo_metabolome, NPC.pathway == "Shikimates and Phenylpropanoids")

Observed_metabolome_Shikimates_and_Phenylpropanoids_dat <- data.frame(t(estimateR((round(t(otu_table(Phylo_metabolome_Shikimates_and_Phenylpropanoids)) ))))) #estimate richness
colnames(Observed_metabolome_Shikimates_and_Phenylpropanoids_dat)[1] <- "Obs_features"
Observed_metabolome_Shikimates_and_Phenylpropanoids_dat <- Observed_metabolome_Shikimates_and_Phenylpropanoids_dat[-c(2:5)]
Observed_metabolome_Shikimates_and_Phenylpropanoids_dat$NPC.pathway <- rep("Shikimates and Phenylpropanoids")
Observed_metabolome_Shikimates_and_Phenylpropanoids_dat$sample_name_2 <- rownames(Observed_metabolome_Shikimates_and_Phenylpropanoids_dat)
Observed_metabolome_Shikimates_and_Phenylpropanoids_dat <- left_join(Observed_metabolome_Shikimates_and_Phenylpropanoids_dat, data.frame(sample_data(Phylo_metabolome)))





#Terpenoids
Phylo_metabolome_Terpenoids <- subset_taxa(Phylo_metabolome, NPC.pathway == "Terpenoids")

Observed_metabolome_Terpenoids_dat <- data.frame(t(estimateR((round(t(otu_table(Phylo_metabolome_Terpenoids)) ))))) #estimate richness
colnames(Observed_metabolome_Terpenoids_dat)[1] <- "Obs_features"
Observed_metabolome_Terpenoids_dat <- Observed_metabolome_Terpenoids_dat[-c(2:5)]
Observed_metabolome_Terpenoids_dat$NPC.pathway <- rep("Terpenoids")
Observed_metabolome_Terpenoids_dat$sample_name_2 <- rownames(Observed_metabolome_Terpenoids_dat)
Observed_metabolome_Terpenoids_dat <- left_join(Observed_metabolome_Terpenoids_dat, data.frame(sample_data(Phylo_metabolome)))



# Observed_total_metabolome_dat <- rbind(Observed_total_metabolome_dat, Observed_metabolome_Alkaloids_dat, 
#                                        Observed_metabolome_Amino_acids_and_Peptides_dat, Observed_metabolome_Fatty_acids_dat, 
#                                        Observed_metabolome_Polyketides_dat, Observed_metabolome_Terpenoids_dat)

Observed_total_metabolome_dat <- rbind(Observed_total_metabolome_dat, Observed_metabolome_Alkaloids_dat, Observed_metabolome_Shikimates_and_Phenylpropanoids_dat,
                                       Observed_metabolome_Amino_acids_and_Peptides_dat, 
                                       Observed_metabolome_Polyketides_dat, Observed_metabolome_Terpenoids_dat)



Observed_total_metabolome_dat$NPC.pathway <- factor(Observed_total_metabolome_dat$NPC.pathway, levels = c("All NPC Pathway groups","Alkaloids","Amino acids and Peptides","Shikimates and Phenylpropanoids","Polyketides","Terpenoids"))
str(Observed_total_metabolome_dat)


#Statistics
p_Obs.stats = Observed_total_metabolome_dat %>%
  dplyr::group_by(Day, NPC.pathway,  Phase) %>%
  dplyr::summarise(mean.obs = mean(Obs_features),
                   sd.obs = sd(Obs_features),
                   counts = n(), zscore = abs((Obs_features-mean(Obs_features))/sd(Obs_features)))

#Find outliers per NPC.pathway per Day based on z-scores

p_Obs.stats %>% filter(abs(zscore) >= 3)

#Make levels of phase
p_Obs.stats$Phase = factor(p_Obs.stats$Phase, levels= c("Early", "Peak", "Late"))

library(ggbeeswarm)

main_col <- c("#240785",  "#f21395", "#e0b62b")

Figure_4D <- p_Obs.stats %>% 
  filter(NPC.pathway =="Amino acids and Peptides") %>% 
  ggplot(aes(x = Day, y = mean.obs)) + 
  geom_quasirandom(data = Observed_total_metabolome_dat %>% 
                     filter(NPC.pathway =="Amino acids and Peptides"), 
                   aes(y = Obs_features, col = Phase),size = 3, dodge.width=0, alpha = 0.6, stroke = NA) + 
  geom_errorbar(aes(ymax = mean.obs+sd.obs, ymin = mean.obs-sd.obs, col = Phase), width = 0, linewidth = 0.7, alpha = 0.2) +
  geom_line(size = 1, alpha = 0.6, col = "black") +
  labs(x= "\nTime (days)", y = "\n Feature richness \n", 
       #subtitle = "Biofilm - Metabolome"
       ) +
  scale_x_continuous(breaks = c(10,15,23,29,44,57,71,85,99,113))+
  theme_bw(base_size = 12) +
  #facet_wrap(.~NPC.pathway, scales = "free_y") +
  scale_color_manual(values = main_col)+
  scale_fill_manual(values = main_col) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right",
        axis.text.x = element_text(face = "bold", colour ="black", angle = 45, vjust =1, hjust = 0.5),
        axis.title = element_text(face = "bold", colour ="black", vjust =1),
        axis.text = element_text(face = "bold", colour ="black"))

Figure_4D

ggsave(file="Metabolome/Figures/Figure_4D.svg", plot=Figure_4D,width=4.9, height=3.5)


#Test if late is significant between the days using dunntest
groups <- unique(Observed_total_metabolome_dat$NPC.pathway)

# Initialize list to store results
dunnTest_results_perDay <- list()

# Loop through each group
for (i in 1:length(groups)) {
  # Subset data for current group
  group_data <- subset(Observed_total_metabolome_dat, Observed_total_metabolome_dat$NPC.pathway == groups[i])
  # Perform Dunn test for current group
  test_result <- dunnTest(Obs_features ~ as.factor(Day), data = group_data, method = "bh")
  # Store dunnTest_results in list with group name as name if P.adj is less than 0.05
  dunnTest_results_perDay[[groups[i]]] <- test_result
  names(dunnTest_results_perDay)[i] <- as.character(groups[i])
}

data.frame(dunnTest_results_perDay$`Amino acids and Peptides`$res) %>% filter(P.adj < 0.05)


#Test multible comparisons differences between the days within each NPC.pathway
# Get unique groups
groups <- unique(Observed_total_metabolome_dat$NPC.pathway)

# Initialize list to store results
dunnTest_results_perDay <- list()

library(FSA)

# Loop through each group
for (i in 1:length(groups)) {
  # Subset data for current group
  group_data <- subset(Observed_total_metabolome_dat, Observed_total_metabolome_dat$NPC.pathway == groups[i])
  # Perform Dunn test for current group
  test_result <- dunnTest(Obs_features ~ as.factor(Day), data = group_data, method = "bh")
  # Store dunnTest_results in list with group name as name if P.adj is less than 0.05
  dunnTest_results_perDay[[groups[i]]] <- test_result
  names(dunnTest_results_perDay)[i] <- as.character(groups[i])
}

data.frame(dunnTest_results_perDay$`Amino acids and Peptides`$res) %>% filter(P.adj < 0.05)


#Differential analysis#####
#Kruskal-Wallis: Differences between the features between the early/peak and the late (grouping=Phase_major)

#Make long-format
Phylo_metabolome_Early.Peak <- subset_samples(Phylo_metabolome, Phase_major == "Early")
Phylo_metabolome_Late <- subset_samples(Phylo_metabolome, Phase_major == "Late")

Phylo_metabolome_Early.Peak_dat <- as.data.frame(otu_table(Phylo_metabolome_Early.Peak))
Phylo_metabolome_Early.Peak_dat$m.z <- rownames(Phylo_metabolome_Early.Peak_dat)
Phylo_metabolome_Late_dat <- as.data.frame(otu_table(Phylo_metabolome_Late))
Phylo_metabolome_Late_dat$m.z <- rownames(Phylo_metabolome_Late_dat)

#Make long format
Phylo_metabolome_Early.Peak_dat_long <- Phylo_metabolome_Early.Peak_dat %>% 
  gather(key = "variable", value = "value", -m.z)

Phylo_metabolome_Early.Peak_dat_long$Phase_major <- rep("Early")

Phylo_metabolome_Late_dat_long <- Phylo_metabolome_Late_dat %>% 
  gather(key = "variable", value = "value", -m.z)

Phylo_metabolome_Late_dat_long$Phase_major <- rep("Late")

Phylo_metabolome_dat_long <- rbind(Phylo_metabolome_Early.Peak_dat_long, 
                                                     Phylo_metabolome_Late_dat_long)


Phylo_metabolome_dat_long$Phase_major <- factor(Phylo_metabolome_dat_long$Phase_major , levels=c("Early", "Late"))
head(Phylo_metabolome_dat_long)

#Make Kruskal-wallis test for each feature


results_kruskal.test <- by(Phylo_metabolome_dat_long, Phylo_metabolome_dat_long$m.z, 
                           function(x) {
                             test <- kruskal.test(x$value ~ as.factor(x$Phase_major))
                             data.frame(statistic = test$statistic, p.value = test$p.value)
                           })
#make dataframe
results_kruskal.test <- do.call("rbind", results_kruskal.test)
results_kruskal.test$m.z <- rownames(results_kruskal.test)

#P.adjusted values
results_kruskal.test$p.adj <- p.adjust(results_kruskal.test$p.value,method = "BH")

results_kruskal.test_sig <- results_kruskal.test %>% filter(p.adj < 0.05)

#metadata table
Phylo_metabolome_metadata <- as.data.frame(tax_table(Phylo_metabolome))[4:6]
Phylo_metabolome_metadata$m.z <- rownames(Phylo_metabolome_metadata) 
#join with metadata table
results_kruskal.test_sig <- left_join(results_kruskal.test_sig, Phylo_metabolome_metadata) 

#stat summary 
#calculate mean
Phylo_metabolome_dat_long_mean <- Phylo_metabolome_dat_long %>% group_by(Phase_major, m.z) %>% 
  dplyr::summarise(mean=mean(value)) 

Phylo_metabolome_dat_long_mean_max <- Phylo_metabolome_dat_long_mean %>% group_by(m.z) %>% 
  dplyr::summarise(mean = max(mean))


summary_mean_met <- left_join(Phylo_metabolome_dat_long_mean_max, Phylo_metabolome_dat_long_mean)

#join with table telling if it is significant in the early/peak or late
results_kruskal.test_sig <- left_join(results_kruskal.test_sig, summary_mean_met) 

results_kruskal.test_sig

results_kruskal.test_sig %>% group_by(Phase_major) %>% 
  dplyr::summarise(n=n())

results_kruskal.test_sig %>% group_by(Phase_major, NPC.pathway,m.z) %>% 
  dplyr::summarise(n=n()) %>% head(.,20)


## Mirrored beeswarm plot 
#Filter significant features from the phyloseq object:

Phylo_metabolome_sig <- subset_taxa(Phylo_metabolome, mass %in% results_kruskal.test_sig$m.z)
Phylo_metabolome_sig_Early.Peak <- subset_taxa(Phylo_metabolome, mass %in% results_kruskal.test_sig$m.z) %>% subset_samples(., Phase_major =="Early")
Phylo_metabolome_sig_Late <- subset_taxa(Phylo_metabolome, mass %in% results_kruskal.test_sig$m.z) %>% subset_samples(., Phase_major =="Late")


Phylo_metabolome_sig_Early.Peak_otutab <- as.data.frame(otu_table(Phylo_metabolome_sig_Early.Peak))
Phylo_metabolome_sig_Early.Peak_otutab$m.z <- rownames(Phylo_metabolome_sig_Early.Peak_otutab)
Phylo_metabolome_sig_Early.Peak_otutab_tax <- as.data.frame(tax_table(Phylo_metabolome_sig_Early.Peak)) %>% 
  mutate(., m.z = rownames(.)) %>% left_join(Phylo_metabolome_sig_Early.Peak_otutab,.)



Phylo_metabolome_sig_Late_otutab <- as.data.frame(otu_table(Phylo_metabolome_sig_Late))
Phylo_metabolome_sig_Late_otutab$m.z <- rownames(Phylo_metabolome_sig_Late_otutab)
Phylo_metabolome_sig_Late_otutab_tax <- as.data.frame(tax_table(Phylo_metabolome_sig_Late)) %>% 
  mutate(., m.z = rownames(.)) %>% left_join(Phylo_metabolome_sig_Late_otutab,.)

#Change to long format
Phylo_metabolome_sig_Early.Peak_otutab_tax_long <- Phylo_metabolome_sig_Early.Peak_otutab_tax %>% 
  gather(sample, values, -c(m.z, mass, molecularFormula, adduct, NPC.pathway, 
                            NPC.superclass, NPC.class, ClassyFire.most.specific.class,
                            ClassyFire.level.5,ClassyFire.subclass,ClassyFire.class,
                            ClassyFire.superclass, ClassyFire.all.classifications))

Phylo_metabolome_sig_Early.Peak_otutab_tax_long$Phase_Major <- rep("Early/Peak")

Phylo_metabolome_sig_Late_otutab_tax_long <- Phylo_metabolome_sig_Late_otutab_tax %>% 
  gather(sample, values, -c(m.z, mass, molecularFormula, adduct, NPC.pathway, 
                            NPC.superclass, NPC.class, ClassyFire.most.specific.class,
                            ClassyFire.level.5,ClassyFire.subclass,ClassyFire.class,
                            ClassyFire.superclass, ClassyFire.all.classifications))


Phylo_metabolome_sig_Late_otutab_tax_long$Phase_Major <- rep("Late")

#Join both dataframes

Phylo_metabolome_sig_otutab_tax_long <- rbind(Phylo_metabolome_sig_Early.Peak_otutab_tax_long, 
                                              Phylo_metabolome_sig_Late_otutab_tax_long)

Phylo_metabolome_sig_otutab_tax_long$Phase_Major <- factor(Phylo_metabolome_sig_otutab_tax_long$Phase_Major, levels= c("Early/Peak", "Late"))

#Change Early/Peak values to negative values to mirror values across intercept:

Phylo_metabolome_sig_otutab_tax_long$values <- ifelse(Phylo_metabolome_sig_otutab_tax_long$Phase_Major =="Early/Peak", -1*Phylo_metabolome_sig_otutab_tax_long$values, Phylo_metabolome_sig_otutab_tax_long$values)
#Statistics
p_values.stats = Phylo_metabolome_sig_otutab_tax_long %>%
  dplyr::group_by(Phase_Major, m.z) %>%
  dplyr::summarise(mean.values = mean(values),
                   sd.values = sd(values),
                   counts = n())

colnames(Phylo_metabolome_sig_Late_otutab_tax)

Phylo_metabolome_sig_otutab_tax_long_new <- left_join(Phylo_metabolome_sig_otutab_tax_long, p_values.stats[2:3])

#Make levels of phase
p_values.stats$Phase_Major = factor(p_values.stats$Phase_Major, levels= c("Early/Peak", "Late"))
p_values.stats <- left_join(p_values.stats, Phylo_metabolome_sig_Late_otutab_tax[c(50,54)] , by = "m.z")


p_values.stats <- left_join(p_values.stats, results_kruskal.test_sig[c(1,3:4)])
p_values.stats$p.adj <- round(p_values.stats$p.adj, digits = 3)
p_values.stats$p.adj <- ifelse(p_values.stats$p.adj==0, "p < 0.001", p_values.stats$p.adj)
p_values.stats$p.adj <- ifelse(p_values.stats$p.adj!="p < 0.001", paste("p = ", p_values.stats$p.adj, sep = ""), p_values.stats$p.adj)
p_values.stats$statistic <- round(p_values.stats$statistic, digits = 1)
p_values.stats$statistic <- paste("chi-square = ", p_values.stats$statistic, sep = "")
p_values.stats$kruskal.test <- paste(p_values.stats$statistic,p_values.stats$p.adj, sep = ", ")

p_values.stats$kruskal.test <- ifelse(p_values.stats$Phase_Major=="Early/Peak", NA, p_values.stats$kruskal.test)
p_values.stats$p.adj <- ifelse(p_values.stats$Phase_Major=="Early/Peak", NA, p_values.stats$p.adj)


  library(ggbeeswarm)

NPC.pathway_col <- c("#c439e3",  "#e3b039", "#e35b39", "#a8195a", "#240785", "#969696")


Figure_4C <- p_values.stats %>%
  ggplot(aes(y =reorder(m.z, -mean.values), x = mean.values)) + 
  geom_vline(xintercept=0) +
  geom_quasirandom(data = Phylo_metabolome_sig_otutab_tax_long_new,
                 aes(x = values, shape = Phase_Major, col=NPC.pathway),
                 size = 1.2,alpha = 0.7,  dodge.width=.4, width = 0.8) +
  theme_bw(base_size = 10) +
  geom_point(aes(shape = Phase_Major), size = 2, position=position_dodge(width=0.7)) +
  geom_errorbar(aes(xmax = mean.values+sd.values, xmin = mean.values-sd.values), 
                position=position_dodge2(0.7), width = 0.7, linewidth = 0.7, alpha=0.5) +  
  xlab("LC-MS Peak Area") + ylab(expression(italic("m/z"))) + labs(color='NPC Pathway', shape = 'Phase') + 
  scale_color_manual(values = NPC.pathway_col) +
  geom_text(aes(label=p.adj), position=position_dodge(width=0.9), hjust=-0.5, vjust=-0.9, size = 3) +
  #scale_fill_manual(values = main_col_phase_major) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size = 15), 
        #axis.title = element_text(colour = "black", face = "bold"),
        
        #axis.text = element_text(colour = "black", face = "bold"),
        #legend.position = "none",
        #legend.background = element_rect(linetype = 1, size = 0.5, colour = "lightgrey"),
        #axis.text.x = element_text(angle = 45, vjust =1, hjust = 0.5)
  )
Figure_4C

ggsave(file="Metabolome/Figures/Figure_4C.svg", plot=Figure_4C, width=11, height=6.5)


