
###########################
## Metabolome ordination ##
###########################

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
library(IMIFA) # Pareto scaling

##Load feature table
#'Row ID' in the quant table and 'ID in the SIRIUS table are identical and can be linked.
quat_table <- read.csv("Figure5-Metabolome/Mzmine_GNPS output_1000 filtered_quant.csv", sep = ",", header = TRUE)
quat_table$row.m.z
quat_table <- as.data.frame(quat_table)
str(quat_table)

SERIUS_table <- read.csv("Figure5-Metabolome/SIRIUS_predictions.csv", sep = ";", header = TRUE)
SERIUS_table <- as.data.frame(SERIUS_table)

colnames(SERIUS_table)[1] <- colnames(quat_table)[1]
#Join tables based in "row.ID"
met_table <- as.data.frame(left_join(quat_table, SERIUS_table))
#round row.m.zto 4 decimals
met_table$row.m.z <- round(met_table$row.m.z, digits = 5)
#Remove metabolite "128.1066" (empty row)
met_table <- met_table[!met_table$row.m.z == "128.1066",]



#Wrangle data to get a featuretable, annotation table and metadata for import to phyloseq
rownames(met_table) <- make.unique(as.character(met_table$row.m.z))
colnames(met_table)[14:99]
#Make feature dataframe
feature_abundance <- met_table[,14:99]
#remove blancks 
feature_abundance <- feature_abundance[rowSums(feature_abundance[,c(1:3, 6, 66)])==0,]
#remove blanck cols
feature_abundance <- feature_abundance[,-c(1:3, 6, 66)]




# # Remove metabolites not present in 3 of 81 of the samples
# feature_abundance_curated <- feature_abundance[apply(feature_abundance, 1, FUN = function(x){sum(x == 0)}) < 78,]
# #remove empty columns
# which(colSums(feature_abundance_curated)==0)
# feature_abundance_curated <- feature_abundance_curated[,-13]
#   # post blank removal
# dim(feature_abundance_curated)
# # post inflation removal
# feature_abundance <- feature_abundance_curated





#Make metadata table 
metabolome_sample_data <- read.csv("Figure5-Metabolome/sample_data2.csv", sep = ",", header = TRUE)
#match feature table to metabolome_sample_data
metabolome_sample_data <- metabolome_sample_data[match(colnames(feature_abundance), metabolome_sample_data$sample_name),]
#add readable names to metadata table
rownames(metabolome_sample_data) <- metabolome_sample_data$sample_name_2

#Get sample names from metadata to feature table
colnames(feature_abundance) <- metabolome_sample_data$sample_name_2

#Make annotation ('metabolome_tax_table') dataframe ("phyloseq tax table")
str(met_table[,101:119])
metabolome_annotations_table <- met_table[,101:119]
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




#set levels of Phase, Cat days
sample_data(Phylo_metabolome)$Phase <- factor(sample_data(Phylo_metabolome)$Phase, level=c("Early", "Peak", "Late"))


# Remove outliers 
outliers <- names(which(colSums(feature_abundance) <= 50000))
Phylo_metabolome <-  subset_samples(Phylo_metabolome, !sample_name_2 %in% outliers) 

#Filter samples in ps so they match rownames of AD amplicon samples processed in Figure 1 scripts
AD_filter_norm_toMetabolome <- read.csv("Figure5-Metabolome/AD_filter_norm_toMetabolome.csv", sep = ",", header = TRUE)
rownames(AD_filter_norm_toMetabolome) = AD_filter_norm_toMetabolome$X
Phylo_metabolome <- subset_samples(Phylo_metabolome, sample_name_2 %in% rownames(AD_filter_norm_toMetabolome)) 

summary(as.data.frame(tax_table(Phylo_metabolome)) == "Unclassified")

#Quick look at barplot with classifications before filtering an normalization of the data
colnames(tax_table(Phylo_metabolome))
plot_bar(Phylo_metabolome, "sample_name_2", fill="NPC.pathway") + facet_wrap(~Phase, scales="free_x")
Phylo_metabolome_filter <- filter_taxa(Phylo_metabolome, function (x) {sum(x > 0) > 10}, prune=TRUE) ####  we loose some of the pseudanes at this points. Not present ?
plot_bar(Phylo_metabolome_filter, "sample_name_2", fill="NPC.pathway") + facet_wrap(~Phase, scales="free_x")



subset_taxa(Phylo_metabolome) %>% plot_bar("sample_name_2", fill="NPC.superclass") + facet_grid(NPC.pathway~factor(Day), scales="free")
subset_taxa(Phylo_metabolome) %>% plot_bar("sample_name_2", fill="NPC.class") + facet_grid(NPC.pathway~factor(Day), scales="free") +
  theme(legend.position = "none" )
subset_taxa(Phylo_metabolome_norm, NPC.pathway == "Polyketides") %>% plot_bar("sample_name_2", fill="mass") + 
  facet_grid(NPC.pathway~factor(Day), scales="free") 
as.data.frame(tax_table(subset_taxa(Phylo_metabolome_norm, NPC.pathway == "Polyketides") ))[,1:6]



subset_taxa(Phylo_metabolome_filter, NPC.pathway == "Terpenoids") %>% plot_bar("sample_name_2", fill="mass") + 
  facet_grid(NPC.pathway~factor(Day), scales="free") 

subset_taxa(Phylo_metabolome_filter, NPC.pathway == "Shikimates and Phenylpropanoids") %>% plot_bar("sample_name_2", fill="mass") + 
  facet_grid(NPC.pathway~factor(Day), scales="free") 
as.data.frame(tax_table(subset_taxa(Phylo_metabolome, NPC.pathway == "Shikimates and Phenylpropanoids") ))[,1:6]
#Phenylpropanoids at day 29 and 44: worth looking up

subset_taxa(Phylo_metabolome, NPC.pathway == "Amino acids and Peptides" & mass >= 300) %>% plot_bar("sample_name_2", fill="mass") + 
  facet_grid(NPC.class~factor(Day), scales="free") 
#macrolides mass: 505.2542 common in late phase

as.data.frame(tax_table(subset_taxa(Phylo_metabolome, NPC.pathway == "Amino acids and Peptides") ))[,1:6]

subset_taxa(Phylo_metabolome, NPC.superclass == "Anthranilic acid alkaloids") %>% plot_bar("sample_name_2", fill="mass") + facet_grid(~Day, scales="free")
subset_taxa(Phylo_metabolome, NPC.pathway == "Alkaloids") %>% plot_bar("sample_name_2", fill="NPC.superclass") + facet_wrap(NPC.superclass~Day, scales="free")


#Normalize to relative abundances
Phylo_metabolome_norm = transform_sample_counts(Phylo_metabolome_filter, function(x) x/sum(x)*100000)
plot_bar(Phylo_metabolome_norm, "sample_name_2", fill="NPC.class") + facet_grid(~Day, scales="free")

#Clr transformation ####
# Laura Sisk-Hackworth, Scott T Kelley, 
# An application of compositional data analysis to multiomic time-series data, 
# NAR Genomics and Bioinformatics, Volume 2, Issue 4, December 2020, lqaa079, 
# https://doi.org/10.1093/nargab/lqaa079

library(zCompositions)
zPatterns(as.data.frame(otu_table(Phylo_metabolome_filter)), label = 0)
#Impute pseudo count for zero's
Phylo_metabolome_pseudo_counts <- cmultRepl(as.data.frame(otu_table(Phylo_metabolome)), z.warning = 0.99)


#Centered log ratio transform
library(compositions)
Phylo_metabolome_pseudo_counts_clr_dat <- clr(Phylo_metabolome_pseudo_counts)

#create new phyloseq object
Phylo_metabolome_pseudo_counts_clr <- phyloseq(otu_table(Phylo_metabolome_pseudo_counts_clr_dat, taxa_are_rows = TRUE), 
                             tax_table(as.matrix(metabolome_annotations_table_1)), sample_data(metabolome_sample_data))


#Plot PCoA and NMDS
main_col <- c("#240785",  "#f21395", "#e0b62b")
Phylo_metabolome_pseudo_counts_clr.ord_PCoA <- ordinate(Phylo_metabolome_pseudo_counts_clr, "PCoA", "euclidean")
p2 = plot_ordination(Phylo_metabolome_pseudo_counts_clr, Phylo_metabolome_pseudo_counts_clr.ord_PCoA, type="samples", color="Phase", axes = c(1,2)) 

pcoa_pseudo_counts_clr_euclidean_p =p2 + geom_point(size=7.5) + 
  theme_minimal() + scale_fill_manual(values = main_col) +
  scale_color_manual(values = main_col)

pcoa_pseudo_counts_clr_euclidean_p


Phylo_metabolome_pseudo_counts_clr.ord_PCoA <- ordinate(Phylo_metabolome_pseudo_counts_clr, "PCoA", "euclidean")
p2 = plot_ordination(Phylo_metabolome_pseudo_counts_clr, Phylo_metabolome_pseudo_counts_clr.ord_PCoA, type="samples", color="Replicate", axes = c(1,2)) 

pcoa_pseudo_counts_clr_euclidean_p = p2 + geom_point(size=7.5) + 
  theme_minimal() #+ scale_fill_manual(values = main_col) +
  scale_color_manual(values = main_col)

pcoa_pseudo_counts_clr_euclidean_p

# Phylo_metabolome_pseudo_counts_clr.ord_NMDS <- ordinate(Phylo_metabolome_pseudo_counts_clr, "NMDS", "euclidean")
# p2 = plot_ordination(Phylo_metabolome_pseudo_counts_clr, Phylo_metabolome_pseudo_counts_clr.ord_NMDS, type="samples", color="Phase", axes = c(1,2)) 
# 
# NMDS_pseudo_counts_clr_euclidean_p =p2 + geom_point(size=7.5) + 
#   theme_minimal() + scale_fill_manual(values = main_col) +
#   scale_color_manual(values = main_col)
# 
# NMDS_pseudo_counts_clr_euclidean_p

#Look at the other NPC.pathway separately
unique(as.data.frame(tax_table(Phylo_metabolome_pseudo_counts_clr))$NPC.pathway)

#Amino_acids_and_Peptides
Phylo_metabolome_pseudo_counts_clr_Amino_acids_and_Peptides <- subset_taxa(Phylo_metabolome_pseudo_counts_clr, NPC.pathway == "Amino acids and Peptides")

Phylo_metabolome_pseudo_counts_clr_Amino_acids_and_Peptides.ord_PCoA <- ordinate(Phylo_metabolome_pseudo_counts_clr_Amino_acids_and_Peptides, "PCoA", "euclidean")
p2 = plot_ordination(Phylo_metabolome_pseudo_counts_clr_Amino_acids_and_Peptides, Phylo_metabolome_pseudo_counts_clr_Amino_acids_and_Peptides.ord_PCoA, type="samples", color="Phase", axes = c(1,2)) 

pcoa_pseudo_counts_clr_euclidean_Amino_acids_and_Peptides_p =p2 + geom_point(size=7.5) + 
  theme_minimal() + scale_fill_manual(values = main_col) +
  scale_color_manual(values = main_col)

pcoa_pseudo_counts_clr_euclidean_Amino_acids_and_Peptides_p


#Terpenoids
Phylo_metabolome_pseudo_counts_clr_Terpenoids <- subset_taxa(Phylo_metabolome_pseudo_counts_clr, NPC.pathway == "Terpenoids")

Phylo_metabolome_pseudo_counts_clr_Terpenoids.ord_PCoA <- ordinate(Phylo_metabolome_pseudo_counts_clr_Terpenoids, "PCoA", "euclidean")
p2 = plot_ordination(Phylo_metabolome_pseudo_counts_clr_Terpenoids, Phylo_metabolome_pseudo_counts_clr_Terpenoids.ord_PCoA, type="samples", color="Phase", axes = c(1,2)) 

pcoa_pseudo_counts_clr_euclidean_Terpenoids_p =p2 + geom_point(size=7.5) + 
  theme_minimal() + scale_fill_manual(values = main_col) +
  scale_color_manual(values = main_col)

pcoa_pseudo_counts_clr_euclidean_Terpenoids_p


#Alkaloids
Phylo_metabolome_pseudo_counts_clr_Alkaloids <- subset_taxa(Phylo_metabolome_pseudo_counts_clr, NPC.pathway == "Alkaloids")

Phylo_metabolome_pseudo_counts_clr_Alkaloids.ord_PCoA <- ordinate(Phylo_metabolome_pseudo_counts_clr_Alkaloids, "PCoA", "euclidean")
p2 = plot_ordination(Phylo_metabolome_pseudo_counts_clr_Alkaloids, Phylo_metabolome_pseudo_counts_clr_Alkaloids.ord_PCoA, type="samples", color="Phase", axes = c(1,2)) 

pcoa_pseudo_counts_clr_euclidean_Alkaloids_p =p2 + geom_point(size=7.5) + 
  theme_minimal() + scale_fill_manual(values = main_col) +
  scale_color_manual(values = main_col)

pcoa_pseudo_counts_clr_euclidean_Alkaloids_p

#Fatty_acids
Phylo_metabolome_pseudo_counts_clr_Fatty_acids <- subset_taxa(Phylo_metabolome_pseudo_counts_clr, NPC.pathway == "Fatty acids")

Phylo_metabolome_pseudo_counts_clr_Fatty_acids.ord_PCoA <- ordinate(Phylo_metabolome_pseudo_counts_clr_Fatty_acids, "PCoA", "euclidean")
p2 = plot_ordination(Phylo_metabolome_pseudo_counts_clr_Fatty_acids, Phylo_metabolome_pseudo_counts_clr_Fatty_acids.ord_PCoA, type="samples", color="Phase", axes = c(1,2)) 

pcoa_pseudo_counts_clr_euclidean_Fatty_acids_p =p2 + geom_point(size=7.5) + 
  theme_minimal() + scale_fill_manual(values = main_col) +
  scale_color_manual(values = main_col)

pcoa_pseudo_counts_clr_euclidean_Fatty_acids_p

#Shikimates_and_Phenylpropanoids
Phylo_metabolome_pseudo_counts_clr_Shikimates_and_Phenylpropanoids <- subset_taxa(Phylo_metabolome_pseudo_counts_clr, NPC.pathway == "Shikimates and Phenylpropanoids")

Phylo_metabolome_pseudo_counts_clr_Shikimates_and_Phenylpropanoids.ord_PCoA <- ordinate(Phylo_metabolome_pseudo_counts_clr_Shikimates_and_Phenylpropanoids, "PCoA", "euclidean")
p2 = plot_ordination(Phylo_metabolome_pseudo_counts_clr_Shikimates_and_Phenylpropanoids, Phylo_metabolome_pseudo_counts_clr_Shikimates_and_Phenylpropanoids.ord_PCoA, type="samples", color="Phase", axes = c(1,2)) 

pcoa_pseudo_counts_clr_euclidean_Shikimates_and_Phenylpropanoids_p =p2 + geom_point(size=7.5) + 
  theme_minimal() + scale_fill_manual(values = main_col) +
  scale_color_manual(values = main_col)

pcoa_pseudo_counts_clr_euclidean_Shikimates_and_Phenylpropanoids_p


#Fatty_acids
Phylo_metabolome_pseudo_counts_clr_Fatty_acids <- subset_taxa(Phylo_metabolome_pseudo_counts_clr, NPC.pathway == "Fatty acids")

Phylo_metabolome_pseudo_counts_clr_Fatty_acids.ord_PCoA <- ordinate(Phylo_metabolome_pseudo_counts_clr_Fatty_acids, "PCoA", "euclidean")
p2 = plot_ordination(Phylo_metabolome_pseudo_counts_clr_Fatty_acids, Phylo_metabolome_pseudo_counts_clr_Fatty_acids.ord_PCoA, type="samples", color="Phase", axes = c(1,2)) 

pcoa_pseudo_counts_clr_euclidean_Fatty_acids_p =p2 + geom_point(size=7.5) + 
  theme_minimal() + scale_fill_manual(values = main_col) +
  scale_color_manual(values = main_col)

pcoa_pseudo_counts_clr_euclidean_Fatty_acids_p


#Carbohydrates
Phylo_metabolome_pseudo_counts_clr_Carbohydrates <- subset_taxa(Phylo_metabolome_pseudo_counts_clr, NPC.pathway == "Carbohydrates")

Phylo_metabolome_pseudo_counts_clr_Carbohydrates.ord_PCoA <- ordinate(Phylo_metabolome_pseudo_counts_clr_Carbohydrates, "PCoA", "euclidean")
p2 = plot_ordination(Phylo_metabolome_pseudo_counts_clr_Carbohydrates, Phylo_metabolome_pseudo_counts_clr_Carbohydrates.ord_PCoA, type="samples", color="Phase", axes = c(1,2)) 

pcoa_pseudo_counts_clr_euclidean_Carbohydrates_p =p2 + geom_point(size=7.5) + 
  theme_minimal() + scale_fill_manual(values = main_col) +
  scale_color_manual(values = main_col)

pcoa_pseudo_counts_clr_euclidean_Carbohydrates_p


#Polyketides
Phylo_metabolome_pseudo_counts_clr_Polyketides <- subset_taxa(Phylo_metabolome_pseudo_counts_clr, NPC.pathway == "Polyketides")

Phylo_metabolome_pseudo_counts_clr_Polyketides.ord_PCoA <- ordinate(Phylo_metabolome_pseudo_counts_clr_Polyketides, "PCoA", "euclidean")
p2 = plot_ordination(Phylo_metabolome_pseudo_counts_clr_Polyketides, Phylo_metabolome_pseudo_counts_clr_Polyketides.ord_PCoA, type="samples", color="Phase", axes = c(1,2)) 

pcoa_pseudo_counts_clr_euclidean_Polyketides_p =p2 + geom_point(size=7.5) + 
  theme_minimal() + scale_fill_manual(values = main_col) +
  scale_color_manual(values = main_col)

pcoa_pseudo_counts_clr_euclidean_Polyketides_p


### Differential abundances using metaboDiff
library(MetaboDiff)

#Remove features with too many zero's
Phylo_metabolome_filter <- filter_taxa(Phylo_metabolome, function (x) {sum(x > 0) > 10}, prune=TRUE) ####  we loose some of the pseudanes at this points. Not present ?


#Remove samples with more than 80% zero values 
# Calculate the number of zero values in each column
zero_counts <- colSums(as.data.frame(otu_table(Phylo_metabolome_filter))==0)
total <- dim(as.data.frame(otu_table(Phylo_metabolome_filter)))[1]
zero_counts[order(zero_counts)]/total

# Keep only the columns with more than 80% zero values
features_more_than_80 <- colnames(as.data.frame(otu_table(Phylo_metabolome_filter))[, zero_counts/nrow(as.data.frame(otu_table(Phylo_metabolome_filter))) < 0.8])
sample_data(Phylo_metabolome_filter)
Phylo_metabolome_filter <- subset_samples(Phylo_metabolome_filter, sample_name_2 %in% features_more_than_80)

#Replace zero's to NAs
assay_filter <- as.data.frame(lapply(as.data.frame(otu_table(Phylo_metabolome_filter)), function(x) ifelse(x == 0 , NA, x)))
rownames(assay_filter) <- rownames(as.data.frame(otu_table(Phylo_metabolome_filter)))
colnames(assay_filter) <- colnames(as.data.frame(otu_table(Phylo_metabolome_filter)))
dim(as.data.frame(otu_table(Phylo_metabolome_filter)))

#assay - a matrix containing the relative metabolic measurements
#assay <- as.matrix(as.data.frame(otu_table(Phylo_metabolome)))
assay <- as.matrix(assay_filter)
assay[1:5,1:5]


#rowData - a dataframe containing the available metabolite annotation
rowData <- (data.frame(tax_table(Phylo_metabolome_filter)))
rownames(rowData) <- rownames(assay)
str(rowData)
#colData - a dataframe containing sample metadata
colData <- data.frame(sample_data(Phylo_metabolome_filter))
#Fix order
colData$Phase <- factor(colData$Phase, levels = c("Late", "Early", "Peak") )
colData$Phase_major <- factor(colData$Phase, levels = c("Late", "Early") )

color_phase = c( "#f21395", "#240785", "#e0b62b")

met <- create_mae(assay,rowData,colData)


na_heatmap(met,
           group_factor="Phase",
           label_colors=color_phase)



met = knn_impute(met,cutoff=0.4)

#Normalization
met <- normalize_met(met)

#Quality control of normalization
quality_plot(met,
             group_factor=c("Phase"),
             label_colors=color_phase)



met = diff_test(met,
                group_factor=c("Phase","Phase_major"))   

str(metadata(met), max.level=2)

par(mfrow=c(1,2))
volcano_plot(met, 
             group_factor="Phase_major",
             label_colors=color_phase,
             dm_cutoff=0.5,
             p_adjust = FALSE)


df <- met@metadata[["ttest_Phase_major_Early_vs_Late"]]
metabolome_annotations_table_1$row.ID <- rownames(metabolome_annotations_table_1)
colnames(df)[1] <- "row.ID"
df_tax <- left_join(df, metabolome_annotations_table_1)
str(df_tax)

library(EnhancedVolcano)
EnhancedVolcano(df_tax,
                lab = c(ifelse(df_tax$adj_pval<10e-3, df_tax$NPC.pathway, 
                                "")),
                boxedLabels = TRUE,
                col = c('grey0', 'grey0', 'red3', 'red3'),
                x = 'dm',
                y = 'adj_pval',
                xlim = c(-5,5),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                labSize = 3,
                pointSize = c(ifelse(df$adj_pval<10e-3, 5, 1)),
                #title = 'Differential Abundance of metabolites related to feed',
                #subtitle = 'Comparison between CTRL, PRO, and SYN',
                xlab = bquote(~Log[2]~ 'Fold Change'),
                #.legend = c('NS','P & Log2 FC'),
                legendLabels = c('NS',expression(p-value~and~log[2]~FC)),
                legendPosition = 'bottom')



# #remove outlier "Succession-T9-9"
# Phylo_metabolome_pseudo_counts_clr <- subset_samples(Phylo_metabolome_pseudo_counts_clr, sample_name_2 != "Succession-T9-9")
# Phylo_metabolome_pseudo_counts_clr.ord_NMDS <- ordinate(Phylo_metabolome_pseudo_counts_clr, "NMDS", "euclidean")
# p2 = plot_ordination(Phylo_metabolome_pseudo_counts_clr, Phylo_metabolome_pseudo_counts_clr.ord_NMDS, type="samples", color="Phase", axes = c(1,2)) 
# 
# NMDS_pseudo_counts_clr_euclidean_p =p2 + geom_point(size=7.5) + 
#   theme_minimal() + scale_fill_manual(values = main_col) +
#   scale_color_manual(values = main_col)
# 
# NMDS_pseudo_counts_clr_euclidean_p




#Beta-disper
## Calculate multivariate dispersions 
bdisp_metabolome_pseudo_counts_clr <- betadisper(vegdist(t(Phylo_metabolome_pseudo_counts_clr_dat), "euclidean"), 
                                                 factor(paste(data.frame(sample_data(Phylo_metabolome_pseudo_counts_clr))$Phase,sep=" ")))


## Permutation test for F to test if betadisp is significant
permutest(bdisp_metabolome_pseudo_counts_clr, pairwise = TRUE, permutations = 999)

## Tukey's Honest Significant Differences
bdisp.HSD <- TukeyHSD(bdisp_metabolome_pseudo_counts_clr)
plot(bdisp.HSD)

##Beta-disper it significant between late vs. early and peak and late.
#Therefore not sure if the significant different between group running PERMANOVA is
#caused by the high variance in the late phase or differences explained by time

#PERMANOVA (Adonis)
Metabolome_pseudo_counts_clr_dat_dist_eucledian = as.matrix(vegdist(t(data.frame(otu_table(Phylo_metabolome_pseudo_counts_clr))), "euclidean"))
metadata_pseudo_counts_clr_dat_dist_eucledian <- data.frame(sample_data(Phylo_metabolome_pseudo_counts_clr))
metadata_pseudo_counts_clr_dat_dist_eucledian$Phase <- factor(metadata_pseudo_counts_clr_dat_dist_eucledian$Phase, levels=c("Early", "Peak", "Late"))


PERMANOVA <- adonis2(Metabolome_pseudo_counts_clr_dat_dist_eucledian ~ Replicate * Phase * Time_Cat, metadata_pseudo_counts_clr_dat_dist_eucledian)
PERMANOVA
















#Quick look for pseudanes/quinolones####
#look for Pseudane annotations
colSums(otu_table(Phylo_metabolome))

colnames(as.data.frame(tax_table(Phylo_metabolome_filter_norm)))
unique(as.data.frame(tax_table(Phylo_metabolome_filter_norm))$ClassyFire.class)
df <- as.data.frame(tax_table(Phylo_metabolome_filter_norm))
sapply(colnames(df), function(x) grep("Pseu", df[,x]))
as.data.frame(tax_table(Phylo_metabolome)) %>% filter(str_detect(NPC.superclass, 'lactams'))
rownames(otu_table(Phylo_metabolome)) == 244
as.data.frame(otu_table(Phylo_metabolome)) %>% filter(str_detect(rownames(.), '244.16'))

subset_taxa(Phylo_metabolome, NPC.superclass == "Anthranilic acid alkaloids") %>% plot_bar("sample_name_2", fill="mass") + facet_wrap(~Phase, scales="free")
subset_taxa(Phylo_metabolome, NPC.superclass == "β-lactams") %>% plot_bar("sample_name_2", fill="molecularFormula") + facet_wrap(~Phase, scales="free")


###alpha-diversity #####
Observed_total_metabolome_dat <- data.frame(t(estimateR((round(t(otu_table(Phylo_metabolome)) ))))) #estimate richness
colnames(Observed_total_metabolome_dat)[1] <- "Obs_features"
Observed_total_metabolome_dat <- Observed_total_metabolome_dat[-c(2:5)]
Observed_total_metabolome_dat$NPC.pathway <- rep("All")
Observed_total_metabolome_dat$sample_name_2 <- rownames(Observed_total_metabolome_dat)
Observed_total_metabolome_dat <- left_join(Observed_total_metabolome_dat, data.frame(sample_data(Phylo_metabolome_norm)))

#Alkaloids
Phylo_metabolome_Alkaloids <- subset_taxa(Phylo_metabolome, NPC.pathway == "Alkaloids")

Observed_metabolome_Alkaloids_dat <- data.frame(t(estimateR((round(t(otu_table(Phylo_metabolome_Alkaloids)) ))))) #estimate richness
colnames(Observed_metabolome_Alkaloids_dat)[1] <- "Obs_features"
Observed_metabolome_Alkaloids_dat <- Observed_metabolome_Alkaloids_dat[-c(2:5)]
Observed_metabolome_Alkaloids_dat$NPC.pathway <- rep("Alkaloids")
Observed_metabolome_Alkaloids_dat$sample_name_2 <- rownames(Observed_metabolome_Alkaloids_dat)
Observed_metabolome_Alkaloids_dat <- left_join(Observed_metabolome_Alkaloids_dat, data.frame(sample_data(Phylo_metabolome_norm)))


#Amino_acids_and_Peptides
Phylo_metabolome_Amino_acids_and_Peptides <- subset_taxa(Phylo_metabolome, NPC.pathway == "Amino acids and Peptides")

Observed_metabolome_Amino_acids_and_Peptides_dat <- data.frame(t(estimateR((round(t(otu_table(Phylo_metabolome_Amino_acids_and_Peptides)) ))))) #estimate richness
colnames(Observed_metabolome_Amino_acids_and_Peptides_dat)[1] <- "Obs_features"
Observed_metabolome_Amino_acids_and_Peptides_dat <- Observed_metabolome_Amino_acids_and_Peptides_dat[-c(2:5)]
Observed_metabolome_Amino_acids_and_Peptides_dat$NPC.pathway <- rep("Amino_acids and Peptides")
Observed_metabolome_Amino_acids_and_Peptides_dat$sample_name_2 <- rownames(Observed_metabolome_Amino_acids_and_Peptides_dat)
Observed_metabolome_Amino_acids_and_Peptides_dat <- left_join(Observed_metabolome_Amino_acids_and_Peptides_dat, data.frame(sample_data(Phylo_metabolome_norm)))


#Fatty_acids
Phylo_metabolome_Fatty_acids <- subset_taxa(Phylo_metabolome, NPC.pathway == "Fatty acids")

Observed_metabolome_Fatty_acids_dat <- data.frame(t(estimateR((round(t(otu_table(Phylo_metabolome_Fatty_acids)) ))))) #estimate richness
colnames(Observed_metabolome_Fatty_acids_dat)[1] <- "Obs_features"
Observed_metabolome_Fatty_acids_dat <- Observed_metabolome_Fatty_acids_dat[-c(2:5)]
Observed_metabolome_Fatty_acids_dat$NPC.pathway <- rep("Fatty acids")
Observed_metabolome_Fatty_acids_dat$sample_name_2 <- rownames(Observed_metabolome_Fatty_acids_dat)
Observed_metabolome_Fatty_acids_dat <- left_join(Observed_metabolome_Fatty_acids_dat, data.frame(sample_data(Phylo_metabolome_norm)))


#Polyketides
Phylo_metabolome_Polyketides <- subset_taxa(Phylo_metabolome, NPC.pathway == "Polyketides")

Observed_metabolome_Polyketides_dat <- data.frame(t(estimateR((round(t(otu_table(Phylo_metabolome_Polyketides)) ))))) #estimate richness
colnames(Observed_metabolome_Polyketides_dat)[1] <- "Obs_features"
Observed_metabolome_Polyketides_dat <- Observed_metabolome_Polyketides_dat[-c(2:5)]
Observed_metabolome_Polyketides_dat$NPC.pathway <- rep("Polyketides")
Observed_metabolome_Polyketides_dat$sample_name_2 <- rownames(Observed_metabolome_Polyketides_dat)
Observed_metabolome_Polyketides_dat <- left_join(Observed_metabolome_Polyketides_dat, data.frame(sample_data(Phylo_metabolome_norm)))


# #Shikimates_and_Phenylpropanoids
# Phylo_metabolome_Shikimates_and_Phenylpropanoids <- subset_taxa(Phylo_metabolome, NPC.pathway == "Shikimates and Phenylpropanoids")
# 
# Observed_metabolome_Shikimates_and_Phenylpropanoids_dat <- data.frame(t(estimateR((round(t(otu_table(Phylo_metabolome_Shikimates_and_Phenylpropanoids)) ))))) #estimate richness
# colnames(Observed_metabolome_Shikimates_and_Phenylpropanoids_dat)[1] <- "Obs_features"
# Observed_metabolome_Shikimates_and_Phenylpropanoids_dat <- Observed_metabolome_Shikimates_and_Phenylpropanoids_dat[-c(2:5)]

#Terpenoids
Phylo_metabolome_Terpenoids <- subset_taxa(Phylo_metabolome, NPC.pathway == "Terpenoids")

Observed_metabolome_Terpenoids_dat <- data.frame(t(estimateR((round(t(otu_table(Phylo_metabolome_Terpenoids)) ))))) #estimate richness
colnames(Observed_metabolome_Terpenoids_dat)[1] <- "Obs_features"
Observed_metabolome_Terpenoids_dat <- Observed_metabolome_Terpenoids_dat[-c(2:5)]
Observed_metabolome_Terpenoids_dat$NPC.pathway <- rep("Terpenoids")
Observed_metabolome_Terpenoids_dat$sample_name_2 <- rownames(Observed_metabolome_Terpenoids_dat)
Observed_metabolome_Terpenoids_dat <- left_join(Observed_metabolome_Terpenoids_dat, data.frame(sample_data(Phylo_metabolome_norm)))



Observed_total_metabolome_dat <- rbind(Observed_total_metabolome_dat, Observed_metabolome_Alkaloids_dat, 
                                       Observed_metabolome_Amino_acids_and_Peptides_dat, Observed_metabolome_Fatty_acids_dat, 
                                       Observed_metabolome_Polyketides_dat, Observed_metabolome_Terpenoids_dat)


# ggplot(Observed_total_metabolome_dat, aes(x=Day, y=Obs_features, fill=factor(Day))) +
#   geom_boxplot() + facet_wrap(~NPC.pathway, scales="free_y")  + geom_jitter(shape=16, position=position_jitter(0.2))
# 

#Statistics
p_Obs.stats = Observed_total_metabolome_dat %>%
  dplyr::group_by(Day, NPC.pathway,  Phase) %>%
  dplyr::summarise(mean.obs = mean(Obs_features),
                   sd.obs = sd(Obs_features),
                   counts = n())
#Make levels of phase
p_Obs.stats$Phase = factor(p_Obs.stats$Phase, levels= c("Early", "Peak", "Late"))

library(ggbeeswarm)


p_Obs.stats %>% 
  ggplot(aes(x = Day, y = mean.obs)) + 
  geom_quasirandom(data = Observed_total_metabolome_dat, 
                   aes(y = Obs_features, col = Phase),size = 3, dodge.width=0, alpha = 0.6, stroke = NA) + 
  geom_errorbar(aes(ymax = mean.obs+sd.obs, ymin = mean.obs-sd.obs, col = Phase), width = 0, linewidth = 0.7, alpha = 0.2) +
  geom_line(size = 1, alpha = 0.6, col = "black") +
  labs(x= "\nTime (days)", y = "\n Feature richness \n", subtitle = "Biofilm - Metabolome" )+
  scale_x_continuous(breaks = c(10,15,23,29,44,57,71,85,99,113))+
  theme_bw(base_size = 10) +
  facet_wrap(.~NPC.pathway, scales = "free_y") +
  scale_color_manual(values = color_phase)+
  scale_fill_manual(values = color_phase) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(colour = "black", face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(colour = "black", face = "bold"),
        legend.position = "none",
        #legend.background = element_rect(linetype = 1, size = 0.5, colour = "lightgrey"),
        axis.text.x = element_text(angle = 45, vjust =1, hjust = 0.5))



### beta-diversity with jaccard distance ######
GP.ord <- ordinate(Phylo_metabolome, "PCoA", "jaccard")
p2 = plot_ordination(Phylo_metabolome, GP.ord, type="samples", color="Phase", axes = c(1,2)) 
pcoa_jaccard_p =p2 + geom_point(size=7.5) + 
  theme_minimal() + scale_fill_manual(values = color_phase) +
  scale_color_manual(values = color_phase)
pcoa_jaccard_p

GP.ord_Alkaloids <- ordinate(Phylo_metabolome_Alkaloids, "PCoA", "jaccard")
p2 = plot_ordination(Phylo_metabolome_Alkaloids, GP.ord_Alkaloids, type="samples", color="Phase", axes = c(1,2)) 
pcoa_jaccard_p_Alkaloids =p2 + geom_point(size=7.5) + 
  theme_minimal() + scale_fill_manual(values = color_phase) +
  scale_color_manual(values = color_phase)
pcoa_jaccard_p_Alkaloids

GP.ord_Terpenoids <- ordinate(Phylo_metabolome_Terpenoids, "PCoA", "jaccard")
p2 = plot_ordination(Phylo_metabolome_Terpenoids, GP.ord_Terpenoids, type="samples", color="Phase", axes = c(1,2)) 
pcoa_jaccard_p_Terpenoids =p2 + geom_point(size=7.5) + 
  theme_minimal() + scale_fill_manual(values = color_phase) +
  scale_color_manual(values = color_phase)
pcoa_jaccard_p_Terpenoids

GP.ord_Fatty_acids <- ordinate(Phylo_metabolome_Fatty_acids, "PCoA", "jaccard")
p2 = plot_ordination(Phylo_metabolome_Fatty_acids, GP.ord_Fatty_acids, type="samples", color="Phase", axes = c(1,2)) 
pcoa_jaccard_p_Fatty_acids =p2 + geom_point(size=7.5) + 
  theme_minimal() + scale_fill_manual(values = color_phase) +
  scale_color_manual(values = color_phase)
pcoa_jaccard_p_Fatty_acids






#Beta-disper ----
## Calculate multivariate dispersions from all 18 combinations

metabolome_dist_jaccard_bdisp_stages <- betadisper(vegdist(t(as.matrix(otu_table(Phylo_metabolome))), "jaccard"), 
                                              factor(paste(data.frame(sample_data(Phylo_metabolome))$Phase,sep=" ")))
#PCoA
plot(metabolome_dist_jaccard_bdisp_stages)
#Test betadisper significans
## Permutation test for F
permutest(metabolome_dist_jaccard_bdisp_stages, pairwise = TRUE, permutations = 999)
#beta-disper is significant

# metabolome (use data from betadisper)
centroids.metabolome = data.frame(metabolome_dist_jaccard_bdisp_stages$centroids[,1:2])
centroids.metabolome$Phase = rownames(centroids.metabolome)
centroids.metabolome$Phase <- factor(centroids.metabolome$Phase, levels = c("Early", "Peak", "Late"))
points.metabolome = data.frame(metabolome_dist_jaccard_bdisp_stages$vectors[,1:2], Phase = data.frame(sample_data(Phylo_metabolome_norm))$Phase)
points.metabolome$Phase <- factor(points.metabolome$Phase, levels = c("Early", "Peak", "Late"))


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
  labs(x = "PCoA1 [8.2%]", y = "\nPCoA2 [2.7%]", title = "Metabolome")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = main_col) + 
  scale_color_manual(values = main_col)

PCoA.metabolome.plot






Observed_Phylo_metabolome_dat$Sample <- rownames(Observed_Phylo_metabolome_dat)
#load Observed_ADand16Sand18S_filter_norm_tab from DynamicsRichness_AD99pct.R and merge metabolome richness to dataframe
Observed_ADand16Sand18SandFeatures_filter_norm_tab <- Observed_ADand16Sand18S_filter_norm_tab %>% left_join(Observed_metabolome_dat) %>% filter(Obs_features>=0) 



#transform to longformat
alphadiversity_tab = Observed_ADand16Sand18SandFeatures_filter_norm_tab %>% 
  gather(key="Measure", value="Obs", "Shannon_AD", "Shannon_16S", "Shannon_18S", "Obs_18S", "Obs_16S", "Obs_AD", "Obs_features") %>%
  select(-c(Day,Timepoint,Rep_tec,Sample)) %>%
  separate(col=Measure, into=c('Diversity', 'Measure'), sep='_')





### Filtering with elbow method#####
# #Filter very low abundant metabolites using elbow method
# #read sum per genus 
# feature_sum <- rowSums(data.frame(otu_table(Phylo_metabolome)))
# #data frame
# feature_sum_dat=data.frame(N=names(feature_sum), val=feature_sum)
# #Ordering in decreasing order
# feature_sum_dat=feature_sum_dat[order(feature_sum_dat$val,decreasing = T),]
# 
# #Plot log10 of all feature sums
# plot(log10(feature_sum_dat$val))



#Normalize to relative abundances

Phylo_metabolome_filter_norm = transform_sample_counts(Phylo_metabolome, function(x) x/sum(x))
#And remove empty rows

colSums(otu_table(Phylo_metabolome_filter_norm))
plot_bar(Phylo_metabolome_filter_norm, "Day", fill="ClassyFire.superclass", facet_grid=~Phase)
















##### MetaboDiff #####
library(MetaboDiff)
#For metaboDiff to run metabolites remove samples with more than 80% zero values 



# Calculate the number of zero values in each column
zero_counts <- colSums(as.data.frame(otu_table(Phylo_metabolome))==0)
total <- dim(as.data.frame(otu_table(Phylo_metabolome)))[1]
zero_counts[order(zero_counts)]/total*100

# Keep only the columns with more than 80% zero values
features_more_than_80 <- colnames(as.data.frame(otu_table(Phylo_metabolome_filter))[, zero_counts/nrow(as.data.frame(otu_table(Phylo_metabolome_filter))) < 0.8])
sample_data(Phylo_metabolome_filter)
Phylo_metabolome_filter <- subset_samples(Phylo_metabolome_filter, sample_name_2 %in% features_more_than_80)



#Replace zero's to NAs
assay_filter <- as.data.frame(lapply(as.data.frame(otu_table(Phylo_metabolome_filter)), function(x) ifelse(x ==0, NA, x)))
rownames(assay_filter) <- rownames(as.data.frame(otu_table(Phylo_metabolome_filter)))
colnames(assay_filter) <- colnames(as.data.frame(otu_table(Phylo_metabolome_filter)))
dim(as.data.frame(otu_table(Phylo_metabolome_filter)))

#assay - a matrix containing the relative metabolic measurements
#assay <- as.matrix(as.data.frame(otu_table(Phylo_metabolome)))
assay <- as.matrix(assay_filter)
assay[1:5,1:5]
dim(assay)

#rowData - a dataframe containing the available metabolite annotation
rowData <- (data.frame(tax_table(Phylo_metabolome_filter)))
rownames(rowData) <- rownames(assay)
str(rowData)
#colData - a dataframe containing sample metadata
colData <- data.frame(sample_data(Phylo_metabolome_filter))
#Fix order
colData$Phase <- factor(colData$Phase, levels = c("Late", "Early", "Peak") )
colData$Phase_major <- factor(colData$Phase, levels = c("Late", "Early") )

color_phase = c( "#f21395", "#240785", "#e0b62b")
color_Phase_major = c( "#f21395", "#240785")


met <- create_mae(assay,rowData,colData)


na_heatmap(met,
           group_factor="Phase_major",
           label_colors=color_Phase_major)


met = knn_impute(met,cutoff=0.4)


#Outlier heatmap
outlier_heatmap(met,
                group_factor="Phase",
                label_colors=color_phase,
                k=2)


#Normalization
met <- normalize_met(met)



#Quality control of normalization
quality_plot(met,
             group_factor="Phase",
             label_colors=color_phase)





















## Remove outlier Succession-T9-9 
#Phylo_metabolome <- subset_samples(Phylo_metabolome, sample_name_2!="Succession-T9-9")

#Filter samples in ps so they match rownames of AD amplicon samples processed in Figure 1 scripts
AD_filter_norm_toMetabolome <- read.csv("Figure5-Metabolome/AD_filter_norm_toMetabolome.csv", sep = ",", header = TRUE)
rownames(AD_filter_norm_toMetabolome) = AD_filter_norm_toMetabolome$X
Phylo_metabolome <- subset_samples(Phylo_metabolome, sample_name_2 %in% rownames(AD_filter_norm_toMetabolome)) 



#Remove features with too many missing values
Phylo_metabolome_filter <- filter_taxa(Phylo_metabolome, function (x) {sum(x > 0) > 0}, prune=TRUE) ####  we loose some of the pseudanes at this points. Not present ?


#Normalize to total sum
Phylo_metabolome_filter_norm = transform_sample_counts(Phylo_metabolome_filter, function(x) x/sum(x))

#look for Pseudane annotations
colSums(otu_table(Phylo_metabolome_filter_norm))

colnames(as.data.frame(tax_table(Phylo_metabolome_filter_norm)))
unique(as.data.frame(tax_table(Phylo_metabolome_filter_norm))$ClassyFire.class)
df <- as.data.frame(tax_table(Phylo_metabolome_filter_norm))
sapply(colnames(df), function(x) grep("Pseu", df[,x]))
as.data.frame(tax_table(Phylo_metabolome_filter_norm)) %>% filter(str_detect(NPC.superclass, 'Pseu'))

sample_data(Phylo_metabolome_filter_norm)
###Barplot
#NPC.pathway
subset_taxa(Phylo_metabolome_filter_norm, NPC.pathway == "Alkaloids") %>% plot_bar("Day", fill="ClassyFire.superclass", facet_grid=~Phase)

plot_bar(Phylo_metabolome_filter, "Day", fill="ClassyFire.superclass", facet_grid=~Phase)

#Alpha diversity (observed richness)#####
Observed_metabolome_dat <- data.frame(t(estimateR((round(t(otu_table(Phylo_metabolome_filter)) ))))) #estimate richness
colnames(Observed_metabolome_dat)[1] <- "Obs_features"
Observed_metabolome_dat <- Observed_metabolome_dat[-c(2:5)]
Observed_metabolome_dat$Sample <- rownames(Observed_metabolome_dat)
#load Observed_ADand16Sand18S_filter_norm_tab from DynamicsRichness_AD99pct.R and merge metabolome richness to dataframe
Observed_ADand16Sand18SandFeatures_filter_norm_tab <- Observed_ADand16Sand18S_filter_norm_tab %>% left_join(Observed_metabolome_dat) %>% filter(Obs_features>=0) 



#transform to longformat
alphadiversity_tab = Observed_ADand16Sand18SandFeatures_filter_norm_tab %>% 
  gather(key="Measure", value="Obs", "Shannon_AD", "Shannon_16S", "Shannon_18S", "Obs_18S", "Obs_16S", "Obs_AD", "Obs_features") %>%
  select(-c(Day,Timepoint,Rep_tec,Sample)) %>%
  separate(col=Measure, into=c('Diversity', 'Measure'), sep='_')

color_phase = c("#240785", "#e0b62b", "#f21395")



#Statistics
p_Obs.stats = alphadiversity_tab %>%
  dplyr::group_by(Time, Subject, Measure, Diversity, Phase) %>%
  dplyr::summarise(mean.obs = mean(Obs),
                   sd.obs = sd(Obs),
                   counts = n())
#Make levels of phase
p_Obs.stats$Phase = factor(p_Obs.stats$Phase, levels= c("Early", "Mid", "Late"))

# test = p_Obs.stats %>%
#   filter(Subject != "Seawater" & Diversity != "Shannon") 

#Plot: observed ricness only biofilm
p_Obs.stats %>% filter(Subject != "Seawater" & Diversity != "Shannon") %>%
  ggplot(aes(x = Time, y = mean.obs)) + 
  geom_quasirandom(data = alphadiversity_tab %>% filter(Diversity != "Shannon" & Subject != "Seawater"), 
                   aes(y = Obs, col = Phase),size = 3, dodge.width=0, alpha = 0.6, stroke = NA) + 
  geom_errorbar(aes(ymax = mean.obs+sd.obs, ymin = mean.obs-sd.obs, col = Phase), width = 0, linewidth = 0.7, alpha = 0.2)+
  geom_line(size = 1, alpha = 0.6, col = "black") +
  labs(x= "\nTime (days)", y = "\n ASV richness \n", subtitle = "Biofilm - AD 99% cluster" )+
  scale_x_continuous(breaks = c(3,7,10,15,23,29,44,57,71,85,99,113))+
  theme_bw(base_size = 10) +
  facet_wrap(.~Measure, scales = "free_y") +
  scale_color_manual(values = color_phase)+
  scale_fill_manual(values = color_phase)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(colour = "black", face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(colour = "black", face = "bold"),
        legend.position = "none",
        #legend.background = element_rect(linetype = 1, size = 0.5, colour = "lightgrey"),
        axis.text.x = element_text(angle = 45, vjust =1, hjust = 0.5))



#Investigate correlations

Observed_ADand16Sand18SandFeatures_filter_norm_tab %>% 
  ggplot(aes(y=Obs_16S, x=Obs_features, color = Phase)) +
  geom_point(alpha=0.6, size = 1) +
  #geom_smooth(method = lm, se=FALSE)+
  stat_cor(method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name="rho", size = 2)+
  scale_size(range = c(0.1, 10)) +
  labs(x = "\n\nmetabolic richness ", y = "\nAD richness", subtitle = "") +
  theme_bw(base_size = 8) +
  scale_fill_manual(values = c("#240785", "#f21395", "#e0b62b")) + 
  scale_color_manual(values = c("#240785", "#f21395", "#e0b62b")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        #legend.position = "none"
        )
#scale_y_continuous(limits = c(0,1000))

#Beta-diversity #####

#From phyloseq to vegan
Phylo_metabolome_filter_norm_dat <- as.data.frame((as(otu_table(Phylo_metabolome_filter_norm), "matrix"))) #filtered and normalized data
Phylo_metabolome_dat <- as.data.frame((as(otu_table(Phylo_metabolome), "matrix"))) # non filtered and normalized



#Add arbitrary value to matrix
Phylo_metabolome_filter_norm_dat <- Phylo_metabolome_filter_norm_dat + 1

#Log transformation
Phylo_metabolome_filter_norm_dat_log_transf <- log2(Phylo_metabolome_filter_norm_dat)

#Pareto Scaling 
Phylo_metabolome_filter_norm_dat_log_transf_scalled <- (pareto_scale(Phylo_metabolome_filter_norm_dat_log_transf))

#Distance matrix based on euclidean distances
#Phylo_metabolome_norm_dat_scalled_dist_euc = as.matrix((vegdist(Phylo_metabolome_norm_dat_scalled, "euclidean")))

#Metadata dataframe
Phylo_metabolome_meta_dat <- data.frame(sample_data(Phylo_metabolome))

## PCoA
#with filtered and normalized data
# Make distance matrix w. euclidean distance
dist.matrix = vegdist(t(as.matrix(Phylo_metabolome_filter_norm_dat_log_transf_scalled)), "euclidean")

#cmdscale - PCOA ordination with two dimensions and eigenvalues given, add = true to get positive eigenvalues
pcoa = cmdscale(dist.matrix, k = 2, eig = TRUE, add = TRUE) 
# Calculate variation explained for PCoA1 and 2
100*pcoa$eig / sum(pcoa$eig)
# PCoA1 = 24%
# PCoA2 = 16%

#Beta-dispr - Calculate multivariate dispersion from all 18 combinations
metabolome_dist_eu_bdisp_stages <- betadisper(dist.matrix, factor(paste(Phylo_metabolome_meta_dat$Phase,sep=" ")))

#Permutation test for F
set.seed(41)
permutest(metabolome_dist_eu_bdisp_stages, pairwise = TRUE, permutations = 999)

#PERMANOVA (Adonis)
adonis_metabolome <- adonis2(dist.matrix ~ Phase, Phylo_metabolome_meta_dat)

#Make ordination plot from betadisper centroids and points
centroids.metabolome = data.frame(metabolome_dist_eu_bdisp_stages$centroids[,1:2])
centroids.metabolome$Phase = rownames(centroids.metabolome)
points.metabolome = data.frame(metabolome_dist_eu_bdisp_stages$vectors[,1:2], Phase = Phylo_metabolome_meta_dat$Phase)

main_col <- c("#240785", "#f21395", "#e0b62b")

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
  labs(x = "PCoA1 [24.6%]", y = "\nPCoA2 [16.1%]", title = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = c("#240785", "#f21395", "#e0b62b")) + 
  scale_color_manual(values = c("#240785", "#f21395", "#e0b62b"))

PCoA.metabolome.plot
#ggsave(PCoA.metabolome.plot, file = 'PCoA.metabolome.plot.tiff', width = 5, height = 3)

# Envfit analysis (biplot)
set.seed(41)
envfit.annotations = envfit(pcoa, env = as.data.frame(t(otu_table(Phylo_metabolome_filter_norm))),
                            permutations = 1000, choices = c(1,2)) 

# Look at R values for all asvs, (Goodness of fit /squared correlation coeff)
plot(sort(envfit.annotations$vectors$r, decreasing = TRUE))

# Pick out only vectors with an r value above 0.5
sigIndx=which(envfit.annotations$vectors$r>0.2) 
              #& envfit.annotations$vectors$pvals<0.05)

# Find vector coordinates
arrows = envfit.annotations$vectors$arrows[sigIndx,]
scores.arrows = as.data.frame(scores(envfit.annotations, "vectors"))[sigIndx,] * 1.9
#en_coord_cont = as.data.frame(scores(ASVfit, "vectors")) * ordiArrowMul(ASVfit) # ordiArrowMul gives the factor for meaningfull scaling
colnames(scores.arrows) =c("PCoA1", "PCoA2")

# Add annotation  names
annotations_forfit = as.data.frame(tax_table(Phylo_metabolome_filter_norm))
annotations_forfit2 <- as.data.frame(annotations_forfit[row.names(annotations_forfit) %in% row.names(scores.arrows), ])
rownames_annotations = rownames(scores.arrows)
combined.envfit = cbind(scores.arrows, annotations_forfit2, rownames_annotations) %>%
  unite("Annotations", c("m.z", "NPC.superclass"), sep ="_")

#rownames_loadings <- with(tes, paste(family, ASV_genus, sep = "_"))
rownames_loadings <- with(combined.envfit, Annotations)
row.names(annotations_forfit2) <- make.unique((rownames_loadings))

# PCOA with loadings from envfit

ggplot(centroids.metabolome, aes(PCoA1, PCoA2, color = Phase)) +
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
  labs(x = "PCoA1 [24.6%]", y = "\nPCoA2 [16.1%]", title = "")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = c("#240785", "#f21395", "#e0b62b")) + 
  scale_color_manual(values = c("#240785", "#f21395", "#e0b62b")) + 
  geom_segment(aes(x = 0, y = 0, xend = PCoA1, yend = PCoA2), 
               data = scores.arrows, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = scores.arrows, aes(x = PCoA1, y = PCoA2), colour = "grey10", 
            fontface = "italic", size = 4, label = row.names(scores.arrows))



#### PCoA with Jaccard and eucledian##### 
par(mar = c(10, 4, 4, 2) + 0.1) # make more room on bottom margin
barplot(sort(taxa_sums(Phylo_metabolome_filter_norm), TRUE)[1:74]/nsamples(Phylo_metabolome_filter_norm), las=2, col = "#1B9E77")

GP.ord <- ordinate(Phylo_metabolome_filter, "PCoA", "jaccard")
p2 = plot_ordination(Phylo_metabolome_filter, GP.ord, type="samples", color="Phase", axes = c(1,2)) 
pcoa_jaccard_p =p2 + geom_point(size=7.5) + 
  theme_minimal() + scale_fill_manual(values = main_col) +
  scale_color_manual(values = main_col)

GP.ord <- ordinate(Phylo_metabolome_filter, "PCoA", "euclidean")
p2 = plot_ordination(Phylo_metabolome_filter, GP.ord, type="samples", color="Phase", axes = c(1,2)) 

pcoa_euclidean_p =p2 + geom_point(size=7.5) + 
  theme_minimal() + scale_fill_manual(values = main_col) +
  scale_color_manual(values = main_col)

pcoa_euclidean_p

#PCA #####
library(ggbiplot)
library(factoextra)
pca = prcomp(t(Phylo_metabolome_filter_norm_dat), scale = TRUE) 

fviz_eig(pca)

groups <- as.factor(Phylo_metabolome_meta_dat$Phase)
length(groups)
dim(pca)

fviz_pca_ind(pca,
             col.ind = groups, # color by groups
             palette = color_phase,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)


fviz_pca_biplot(pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)




### Procruster analysis #####

data(varespec)
# vare.dist16S <- vegdist((P16S_filter_norm_dat_clean), "bray")
# vare.distAD <- vegdist((AD_filter_norm_dat_clean), "bray")
# mds.null_AD <- monoMDS(vare.distAD, y = cmdscale(vare.distAD))
# mds.alt_AD <- monoMDS(vare.distAD)
# vare.proc_AD <- procrustes(mds.alt_AD, mds.null_AD, symmetric = TRUE)
# vare.proc_AD
# summary(vare.proc_AD)
# plot(mds.alt_AD)
# plot(vare.proc_AD, kind=2)
# residuals(vare.proc_AD)
# vare.pca <- prcomp(AD_filter_norm_dat_clean)
# scores(prcomp(AD_filter_norm_dat_clean))
PCoA.ord_metabolome <- metabolome_dist_eu_bdisp_stages$vectors
#load PCoA.ord_AD and PCoA.ord_16S from "CombinedAnalysis.R"
PCoA.ord_metabolome

PCoA.ord_16S
PCoA.ord_AD

#Clean up rownames so it match the same sample names for all ordinations
PCoA.ord_metabolome <- PCoA.ord_metabolome[rownames(PCoA.ord_metabolome) %in% rownames(PCoA.ord_AD),]
PCoA.ord_AD <- PCoA.ord_AD[rownames(PCoA.ord_AD) %in% rownames(PCoA.ord_metabolome),]
PCoA.ord_16S <- PCoA.ord_16S[rownames(PCoA.ord_16S) %in% rownames(PCoA.ord_metabolome),]
PCoA.ord_metabolome <- PCoA.ord_metabolome[rownames(PCoA.ord_metabolome) %in% rownames(PCoA.ord_16S),]

protest(X = PCoA.ord_metabolome, Y = PCoA.ord_16S, scores = "sites", permutations = 999) ###Significant but M2 is only -4.441e-16

protest(X = PCoA.ord_16S, Y = PCoA.ord_metabolome, scores = "sites", permutations = 999) ###Significant but M2 is only -4.441e-16
protest(X = PCoA.ord_AD, Y = PCoA.ord_metabolome, scores = "sites", permutations = 999) ###Significant but M2 is 0.661

