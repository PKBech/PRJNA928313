library(rRDPData)
library(debar)
library(tidyverse)
library(Rcpp)
library(plyr)
library("phyloseq")
library(ggplot2)
library(vegan)
library(ANCOMBC)
library(gridExtra)
library(microbiome)
library(stringr)
library(ggpubr)
library(cowplot)


#######
######ANCOM-BC differential analysis#######


###############
## Load data ##
###############
setwd("/Users/pernillekjersgaardbech/Documents/JH21_Biofilm/16S amplicon/")
ps_16S <- readRDS('AmpliconAnalysis/16S/PS-16S-decon.noEuk.noneg.rds')
ps_AD <- readRDS('AmpliconAnalysis/AD/ps_AD')
ps_18S <- readRDS('AmpliconAnalysis/18S/PS-18S-decon.noneg.rds')


refseq_ASV_look_up <- refseq(ps_18S)[names((refseq(ps_18S)))%in%c("ASV_1","ASV_2","ASV_134","ASV_394")]
writeXStringSet(refseq_ASV_look_up, "AmpliconAnalysis/18S/refseq_ASV_look_up.fa", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

###################
## Preprocessing ##
###################

## Remove Phaeobacter strain from AD ASVs and Seawater-T1-3 because it is an outlier ##
ps_AD <- subset_samples(ps_AD, Subject != "JH21" & sample_names(ps_AD) != "Seawater-T1-3")
# Remove outliers: "Succession-T1-6" "Succession-T1-8" "Succession-T2-9" "Succession-T2-7" from 16S phyloobject
# ps_16S_outliers <- c("Succession-T1-6","Succession-T1-8","Succession-T2-9","Succession-T2-7")
ps_16S <- subset_samples(ps_16S, sample_names(ps_16S) != "Succession-T1-6" & sample_names(ps_16S) != "Succession-T1-8" &
                           sample_names(ps_16S) != "Succession-T2-9" & sample_names(ps_16S) != "Succession-T2-7")

## Subset to only bio.element study ##
ps_AD_succession <- subset_samples(ps_AD, Subject == "Succession")
ps_16S_succession <- subset_samples(ps_16S, element.type == "bioelement")
ps_18S_succession <- subset_samples(ps_18S, element.type == "bioelement")

## Compare to groups (T5+6 and T7+8)
sample_data(ps_16S_succession)
ps_16S_succession_T5678 <- subset_samples(ps_16S_succession, timepoint %in% c("5","6","7","8") )
ps_18S_succession_T5678 <- subset_samples(ps_18S_succession, timepoint %in% c("5","6","7","8") )


#Add to metadata phyloseq objects
bryozoan_appearence <- as.data.frame(cbind(Sample = rownames(as.data.frame(sample_data(ps_18S_succession))), Bryozoans = c(rep("no", 6), rep("yes", 27), rep("no", 42), rep("yes", 27) ) ) )
rownames(bryozoan_appearence) <- bryozoan_appearence$Sample

#16S new sample_data for phyloseq object
ps_16S_succession_T5678_meta <- data.frame(sample_data(ps_16S_succession_T5678) )
str(ps_16S_succession_T5678_meta)
ps_16S_succession_T5678_meta <- ps_16S_succession_T5678_meta %>% mutate(Sample = rownames(ps_16S_succession_T5678_meta))

ps_16S_succession_T5678_meta_dat <- left_join(ps_16S_succession_T5678_meta,bryozoan_appearence, "Sample" )
rownames(ps_16S_succession_T5678_meta_dat) <- ps_16S_succession_T5678_meta_dat$Sample
str(ps_16S_succession_T5678_meta_dat)
str(ps_16S_succession_T5678_meta)
ps_16S_succession_T5678 <- phyloseq(otu_table(ps_16S_succession_T5678), tax_table(ps_16S_succession_T5678), sample_data(ps_16S_succession_T5678_meta_dat) )

#18S phylo new sample_data for phyloseq object
ps_18S_succession_T5678
ps_18S_succession_T5678_meta <- data.frame(sample_data(ps_18S_succession_T5678) )
str(ps_18S_succession_T5678_meta)
ps_18S_succession_T5678_meta <- ps_18S_succession_T5678_meta %>% mutate(Sample = rownames(ps_18S_succession_T5678_meta))

ps_18S_succession_T5678_meta_dat <- left_join(ps_18S_succession_T5678_meta,bryozoan_appearence, "Sample" )
rownames(ps_18S_succession_T5678_meta_dat) <- ps_18S_succession_T5678_meta_dat$Sample
str(ps_18S_succession_T5678_meta_dat)
str(ps_18S_succession_T5678_meta)
ps_18S_succession_T5678 <- phyloseq(otu_table(ps_18S_succession_T5678), tax_table(ps_18S_succession_T5678), sample_data(ps_18S_succession_T5678_meta_dat) )





###Run ANCOMBC differential analysis 16S no aggregation####
#Gloom/aggregate to class level

sample_data(ps_16S_succession_T5678)$Bryozoans <- factor(sample_data(ps_16S_succession_T5678)$Bryozoans, levels = c("no", "yes"))


#Run ANCOMBC differential analysis
out_ps_16S_succession_T5678 <- ancombc(phyloseq = ps_16S_succession_T5678, formula = "Bryozoans",
                                       p_adj_method = "holm", lib_cut = 1000,
                                       group = "Bryozoans", struc_zero = TRUE, neg_lb = FALSE,
                                       tol = 1e-5, max_iter = 100, conserve = TRUE,
                                       alpha = 0.01, global = TRUE)



out_ps_16S_succession_T5678$res$diff_abn
#Results 
res_out_ps_16S_succession_T5678 <- out_ps_16S_succession_T5678$res


#### Dataframe 16S between early and late stage for logFC
df_fig1_16S_succession_T5678 = data.frame(res_out_ps_16S_succession_T5678$beta * res_out_ps_16S_succession_T5678$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_16S_succession_T5678 = data.frame(res_out_ps_16S_succession_T5678$se * res_out_ps_16S_succession_T5678$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_16S_succession_T5678)[-1] = paste0(colnames(df_fig2_16S_succession_T5678)[-1], "SD")
colnames(df_fig1_16S_succession_T5678) <- sub("Bryozoans", "", colnames(df_fig1_16S_succession_T5678))
colnames(df_fig2_16S_succession_T5678) <- sub("Bryozoans", "", colnames(df_fig2_16S_succession_T5678))


df_fig_16S_succession_T5678 = df_fig1_16S_succession_T5678 %>% left_join(df_fig2_16S_succession_T5678, by = "taxon_id")  %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_16S_succession_T5678


df_fig_16S_succession_T5678 %>% filter(yes < 0) %>% arrange(desc(yes))
ps_16S_succession_T5678_tax_dat <- data.frame(tax_table(ps_16S_succession_T5678))
ps_16S_succession_T5678_tax_dat <- ps_16S_succession_T5678_tax_dat %>% mutate(taxon_id = rownames(ps_16S_succession_T5678_tax_dat) )

df_fig_16S_succession_T5678 <- left_join(df_fig_16S_succession_T5678, ps_16S_succession_T5678_tax_dat, "taxon_id")


df_fig_16S_succession_T5678 %>% filter(-3>yes) %>%
  ggplot(aes(x=yes, y=reorder(taxon_id, yes), label=taxon_id)) + geom_vline(xintercept = 0, colour="gray", size = 0.3) +
  geom_bar(stat = "identity", width = 0.7, aes(fill = (phylum), color = phylum),
           position = position_dodge(width = 0.7)) +
  labs(y = NULL, x="Log fold change",title = "ANCOM-BC diff. analysis\nbetween early and late biofilm microbiomes") +
  geom_errorbar(aes(xmin = yes - yesSD, xmax = yes + yesSD), width = 0.3,
                position = position_dodge(0.7), color = "black") +  
  theme_bw()

df_fig_16S_succession_T5678 %>% filter(2<yes) %>%
  ggplot(aes(x=yes, y=reorder(taxon_id, yes), label=taxon_id)) + geom_vline(xintercept = 0, colour="gray", size = 0.3) +
  geom_bar(stat = "identity", width = 0.7, aes(fill = (phylum), color = phylum),
           position = position_dodge(width = 0.7)) +
  labs(y = NULL, x="Log fold change",title = "ANCOM-BC diff. analysis\nbetween early and late biofilm microbiomes") +
  geom_errorbar(aes(xmin = yes - yesSD, xmax = yes + yesSD), width = 0.3,
                position = position_dodge(0.7), color = "black") +  
  theme_bw()



df_fig_16S_succession_T5678_no <- df_fig_16S_succession_T5678 %>% filter(-3>yes) %>%
  ggplot(aes(x=yes, y=reorder(taxon_id, yes), label=taxon_id)) + geom_vline(xintercept = 0, colour="gray", size = 0.3) +
  geom_bar(stat = "identity", width = 0.7, aes(fill = (phylum)),
           position = position_dodge(width = 0.7), color="black") +
  labs(y = NULL, x="Log fold change",title = "ANCOM-BC diff. analysis\nbetween early and late biofilm microbiomes") +
  geom_errorbar(aes(xmin = yes - yesSD, xmax = yes + yesSD), width = 0.3,
                position = position_dodge(0.7), color = "black") +  
  theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 


legend_df_fig_16S_succession_T5678_genus <- cowplot::get_legend(df_fig_16S_succession_T5678_genus_no + theme(legend.position="right") )


df_fig_16S_succession_T5678 %>% #filter(phylum %in% c("Acidobacteriota", "Myxococcota")) %>%
  #filter(-4>yes) %>%
  filter(taxon_id %in% c("ASV_129", "ASV_89", "ASV_146", "ASV_18") ) %>%
  ggplot(aes(x=yes, y=reorder(taxon_id, yes), label=taxon_id)) + geom_vline(xintercept = 0, colour="gray", size = 0.3) +
  geom_bar(stat = "identity", width = 0.7, aes(fill = (phylum)),
           position = position_dodge(width = 0.7), color = "black") +
  labs(y = NULL, x="Log fold change",title = "ANCOM-BC diff. analysis\non genera abundance in early and late biofilm microbiome") +
  geom_errorbar(aes(xmin = yes - yesSD, xmax = yes + yesSD), width = 0.3,
                position = position_dodge(0.7), color = "black") +
  theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank())


tiff("../ANCOM-BC_160522_16S_genus.tiff", units="in",  width=13.08661, height=7.93701, res=300)

grid.arrange( df_fig_16S_succession_T5678_genus_no, 
              legend_df_fig_16S_succession_T5678_genus,
              nrow=1, ncol = 2,
              layout_matrix = rbind(c(1,2)),
              widths = c(4, 0.8))

dev.off()


# Those that with positive log fold change is found more abundant together with bryozoans




###genus####
#Gloom/aggregate to genus level


ps_16S_succession_T5678_genus <- aggregate_taxa(ps_16S_succession_T5678, "genus")
sample_data(ps_16S_succession_T5678_genus)$Bryozoans <- factor(sample_data(ps_16S_succession_T5678_genus)$Bryozoans, levels = c("no", "yes"))


ps_18S_succession_T5678_genus <- aggregate_taxa(ps_18S_succession_T5678, "genus")
sample_data(ps_18S_succession_T5678_genus)$Bryozoans <- factor(sample_data(ps_18S_succession_T5678_genus)$Bryozoans, levels = c("no", "yes"))



#Run ANCOMBC differential analysis 16S genus #####
out_ps_16S_succession_T5678_genus <- ancombc(phyloseq = ps_16S_succession_T5678_genus, formula = "Bryozoans",
                                             p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                             group = "Bryozoans", struc_zero = TRUE, neg_lb = FALSE,
                                             tol = 1e-5, max_iter = 100, conserve = TRUE,
                                             alpha = 0.01, global = TRUE)



out_ps_16S_succession_T5678_genus$res$diff_abn
#Results 
res_out_ps_16S_succession_T5678_genus <- out_ps_16S_succession_T5678_genus$res

# #Add Relative abundances for color gradient in later plotting
# #Transform to TSS in percentage
# #Relative abundances 
# phylo_main_no_cont_NoOutliers_scalled <- transform_sample_counts(phylo_main_no_cont_NoOutliers, function(x) 100 * x/sum(x))
# phylo_main_no_cont_NoOutliers_scalled_SW <- subset_samples(phylo_main_no_cont_NoOutliers_scalled, Environment=="Planktonic suspension")
# phylo_main_no_cont_NoOutliers_scalled_SW_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_scalled_SW, Day=="4")
# phylo_main_no_cont_NoOutliers_scalled_SW_day4 <- subset_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4, Genus!="Phaeobacter_inhibens")
# 
# phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order <- aggregate_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4, "Order")
# # phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order_filter <- filter_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order, function (x) {sum(x > 0) > 0}, prune=TRUE)
# sample_data(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order)$Treatment, levels = c("WT", "Control", "dTDA"))
# 
# RelAbs_SW_Order_day4 = as.data.frame(otu_table(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order))
# RelAbs_SW_Order_day4$taxon_id <- rownames(RelAbs_SW_Order_day4)
# 
# 
# 
# colnames(RelAbs_SW_Order_day4)
# RelAbs_SW_Order_day4_long <- reshape(RelAbs_SW_Order_day4, direction = "long",
#                                      v.names = "Rel_Abs",
#                                      varying = 1:ncol(RelAbs_SW_Order_day4)-1,
#                                      times=c(rep('Control',9), rep('dTDA', 9), rep('WT',8)),
# )
# 
# #Mean Rel_Abs per treatment per tax
# RelAbs_SW_Order_day4_long_mean <- RelAbs_SW_Order_day4_long %>%
#   group_by(taxon_id, time) %>%
#   dplyr::summarise(Rel_Abs = mean(Rel_Abs)) 
# 


#dplyr::filter(RelAbs_SW_Order_day4_long_mean, grepl("Chloroplast",taxon_id))

#### Dataframe Day 4 for logFC
df_fig1_16S_succession_T5678_genus = data.frame(res_out_ps_16S_succession_T5678_genus$beta * res_out_ps_16S_succession_T5678_genus$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_16S_succession_T5678_genus = data.frame(res_out_ps_16S_succession_T5678_genus$se * res_out_ps_16S_succession_T5678_genus$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_16S_succession_T5678_genus)[-1] = paste0(colnames(df_fig2_16S_succession_T5678_genus)[-1], "SD")
colnames(df_fig1_16S_succession_T5678_genus) <- sub("Bryozoans", "", colnames(df_fig1_16S_succession_T5678_genus))
colnames(df_fig2_16S_succession_T5678_genus) <- sub("Bryozoans", "", colnames(df_fig2_16S_succession_T5678_genus))


df_fig_16S_succession_T5678_genus = df_fig1_16S_succession_T5678_genus %>% left_join(df_fig2_16S_succession_T5678_genus, by = "taxon_id")  %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_16S_succession_T5678_genus %>% filter(yes < 0) %>% arrange(desc(yes))
ps_16S_succession_T5678_genus_tax_dat <- data.frame(tax_table(ps_16S_succession_T5678_genus))
colnames(ps_16S_succession_T5678_genus_tax_dat)[7] <- "taxon_id"
df_fig_16S_succession_T5678_genus <- left_join(df_fig_16S_succession_T5678_genus, ps_16S_succession_T5678_genus_tax_dat, "taxon_id")

# Those that with positive log fold change is found more abundant together with bryozoans

df_fig_16S_succession_T5678_genus_clean <- rbind(df_fig_16S_succession_T5678_genus[df_fig_16S_succession_T5678_genus$yes>1,], df_fig_16S_succession_T5678_genus[df_fig_16S_succession_T5678_genus$yes<1,])

tiff("../ANCOM-BC_160522.tiff", units="in",  width=10.08661, height=10,.93701, res=300)

df_fig_16S_succession_T5678_genus %>% filter(phylum %in% c("Acidobacteriota", "Myxococcota")) %>% 
  ggplot(aes(x=yes, y=reorder(taxon_id, yes), label=taxon_id)) + geom_vline(xintercept = 0, colour="gray", size = 0.3) +
  geom_bar(stat = "identity", width = 0.7, aes(fill = (genus)),
           position = position_dodge(width = 0.7), color = "black") +
  labs(y = NULL, x="Log fold change",title = "ANCOM-BC diff. analysis\non genera abundance in early and late biofilm microbiome") +
  geom_errorbar(aes(xmin = yes - yesSD, xmax = yes + yesSD), width = 0.3,
                position = position_dodge(0.7), color = "black") +
  theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank())

dev.off()

# df_fig_16S_succession_T5678_genus_no <- df_fig_16S_succession_T5678_genus %>% filter(!(yes %in% (-1:1))) %>%
#   ggplot(aes(x=yes, y=reorder(taxon_id, yes), label=taxon_id)) + geom_vline(xintercept = 0, colour="gray", size = 0.3) +
#   geom_bar(stat = "identity", width = 0.7, aes(fill = (phylum)),
#            position = position_dodge(width = 0.7), color="black") + 
#   labs(y = NULL, x="Log fold change",title = "ANCOM-BC diff. analysis\nSignificant abundance genera in\nthe biofilm prior bryozoan bloom") +
#   geom_errorbar(aes(xmin = yes - yesSD, xmax = yes + yesSD), width = 0.3,
#                 position = position_dodge(0.7), color = "black") +  
#   theme_bw() +
#   theme(#axis.line = element_line(color='black'),
#     plot.background = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor = element_blank(), legend.position = "none") 
# 
# df_fig_16S_succession_T5678_genus_yes <- df_fig_16S_succession_T5678_genus %>% filter(0.5<yes) %>%
#   ggplot(aes(x=yes, y=reorder(taxon_id, yes), label=taxon_id)) + geom_vline(xintercept = 0, colour="gray", size = 0.3) +
#   geom_bar(stat = "identity", width = 0.7, aes(fill = (phylum)),
#            position = position_dodge(width = 0.7), color="black") + 
#   labs(y = NULL, x="Log fold change",title = "\nSignificant abundance genera present\ntogether with bryozoans") +
#   geom_errorbar(aes(xmin = yes - yesSD, xmax = yes + yesSD), width = 0.3,
#                 position = position_dodge(0.7), color = "black") +  
#   theme_bw() +
#   theme(#axis.line = element_line(color='black'),
#     plot.background = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor = element_blank(), legend.position = "none") 
# 
# 
# 
# legend_df_fig_16S_succession_T5678_genus_yes <- cowplot::get_legend(df_fig_16S_succession_T5678_genus_yes + theme(legend.position="right") )
# legend_df_fig_16S_succession_T5678_genus_no <- cowplot::get_legend(df_fig_16S_succession_T5678_genus_no + theme(legend.position="right") )
# 


tiff("../ANCOM-BC_160522.tiff", units="in",  width=10.08661, height=5.93701, res=300)

grid.arrange( df_fig_16S_succession_T5678_genus_no, legend_df_fig_16S_succession_T5678_genus_no,
              df_fig_16S_succession_T5678_genus_yes, legend_df_fig_16S_succession_T5678_genus_yes, 
              nrow=1, ncol = 4,
              layout_matrix = rbind(c(1,2,3,4)),
              widths = c(2.7,2, 2.7, 2))

dev.off()


###Run ANCOMBC differential analysis 16S Class####
#Gloom/aggregate to class level


ps_16S_succession_T5678 <- aggregate_taxa(ps_16S_succession_T5678, "class")
sample_data(ps_16S_succession_T5678)$Bryozoans <- factor(sample_data(ps_16S_succession_T5678)$Bryozoans, levels = c("no", "yes"))


############################################################
#Barplots of microbial composition
##########################################################


# agglomerate at Class level
ps_16S_succession_class <- aggregate_taxa(ps_16S_succession, "class")
# Transform to rel. abundance
ps_16S_succession_class_norm <- transform_sample_counts(ps_16S_succession_class, function(x) 100 * x/sum(x))
# Melt to long format
ps_16S_succession_class_norm_melt <- psmelt(ps_16S_succession_class_norm)
#Transform to percentage
ps_16S_succession_class_norm_melt <- aggregate(Abundance ~ OTU + Sample + timepoint, ps_16S_succession_class_norm_melt, sum)
ps_16S_succession_class_norm_melt_sum <- aggregate(Abundance ~ Sample + timepoint, ps_16S_succession_class_norm_melt, sum)
ps_16S_succession_class_norm_melt_sum <- left_join(ps_16S_succession_class_norm_melt, ps_16S_succession_class_norm_melt_sum, by=c("timepoint", "Sample"))
ps_16S_succession_class_norm_melt_sum_pct <- ps_16S_succession_class_norm_melt_sum %>% 
  mutate(Abundance_percentage = Abundance.x / Abundance.y *100)
ps_16S_succession_class_norm_melt_sum_pct$Abundance_percentage <- round(ps_16S_succession_class_norm_melt_sum_pct$Abundance_percentage, 2)
#Find the 11 most abundant classes for the whole dataset including Phaeobacter (Only 12 colors for plotting)
Top12_Class <- ps_16S_succession_class_norm_melt_sum_pct %>% group_by(OTU) %>%
  dplyr::summarise('Abundance_percentage' = mean(Abundance_percentage)) %>% arrange(desc(Abundance_percentage)) %>% head(11)

Top12_Class[[1]]

#Change Class names to "Others" if not among the 12 most abundant 
ps_16S_succession_class_norm_melt_sum_pct_1 <- mutate(ps_16S_succession_class_norm_melt_sum_pct, Class = ifelse(!OTU %in% Top12_Class[[1]], "Others", OTU))

#Clean up names
ps_16S_succession_class_norm_melt_sum_pct_1$Class <- ifelse(ps_16S_succession_class_norm_melt_sum_pct_1$Class=="Bacteria_Unclassified_Unclassified", "Unclassified",  ps_16S_succession_class_norm_melt_sum_pct_1$Class)
#Reorder and clean up table
ps_16S_succession_class_norm_melt_sum_pct_1$timepoint <- factor(ps_16S_succession_class_norm_melt_sum_pct_1$timepoint, levels=unique(ps_16S_succession_class_norm_melt_sum_pct_1$timepoint))
# phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Class <- factor(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Class , levels=c("Actinobacteria","Flavobacteriia","Planctomycetia","Sphingobacteriia","Alphaproteobacteria", "Betaproteobacteria", "dTDAproteobacteria","Epsilonproteobacteria","Gammaproteobacteria", "Others","Unclassified", "P. inhibens OTU"))
# phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Treatment <- factor(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Treatment , levels=c("Seawater", "Control", "WT", "dTDA"))
# phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1 <- arrange(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1, Day, Treatment)
# 
# #Arrange Sample order according to Days and Treatments
# Sample_ID_order <- distinct(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1, Sample, .keep_all = TRUE)
# Sample_ID_order <- unclass(arrange(Sample_ID_order, Day, Treatment)[5])$Sample
# 
# phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Sample <- factor(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Sample, levels=Sample_ID_order)
# 
# 


#Figure S2B -----
library(ggh4x)

# bar_plot_nolg <- ggplot(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1, aes(x=Sample, y = Abundance_percentage, fill = Class)) +
#   geom_bar(stat = "identity") + labs(y = "Relative abundance (%)", x="Day") +
#   scale_fill_brewer(palette = "Paired") + #facet_grid(cols=vars(Environment_1), scales = "free_x") +
#   facet_nested(. ~ Environment_1 + Day +  Treatment,  scales = "free_x") +
#   theme_bw() +
#   theme(#axis.line = element_line(color='black'),
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.title.x=element_blank(),
#     axis.text.x=element_blank(),
#     axis.ticks.x = element_blank(),
#     legend.position="none",
#     text = element_text(size=10))#,

# 
# #Order classes according to network figure legend
# phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Class <- factor(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1$Class, 
#                                                                         levels=c("Flavobacteriia","Gammaproteobacteria",
#                                                                                  "Alphaproteobacteria",
#                                                                                  "Actinobacteria", "Sphingobacteriia",
#                                                                                  "Planctomycetia",
#                                                                                  "Betaproteobacteria", "Epsilonproteobacteria", "Unclassified",
#                                                                                  "Others", "P. inhibens OTU"))

Paired_colors_expand_norm_shortcol <- c("#FFCE8F","#E6759E","#F89181","#343077","#f9c629","#f29d46","#ea3333","#f4722b","#fefd7d","#4e7677","#9C599E", "#f9c629","#2ace82")

bar_plot_nolg <- ggplot(phylo_main_no_cont_NoOutliers_Class_norm_melt_sum_pct_1, aes(x=Sample, y = Abundance_percentage, fill = Class)) +
  geom_bar(stat = "identity") + labs(y = "Relative abundance (%)", x="Day") +
  scale_fill_manual(
    values = Paired_colors_expand_norm_shortcol,
    aesthetics = c("fill") 
  ) +
  #scale_fill_(palette = "Paired") + #facet_grid(cols=vars(Environment_1), scales = "free_x") +
  facet_nested(. ~ Environment_1 + Day +  Treatment,  scales = "free_x") +
  theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank(),
    legend.position="none",
    text = element_text(size=10))#,




legend_bar_plot <- cowplot::get_legend(bar_plot_nolg + theme(legend.position="bottom") + guides(fill = guide_legend(title.position = "top", ncol = 5)))

#Figure S2 Combined ----
#Run ANCOMBC differential analysis
out_ps_16S_succession_T5678 <- ancombc(phyloseq = ps_16S_succession_T5678, formula = "Bryozoans",
                                                           p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                                           group = "Bryozoans", struc_zero = TRUE, neg_lb = FALSE,
                                                           tol = 1e-5, max_iter = 100, conserve = TRUE,
                                                           alpha = 0.01, global = TRUE)



out_ps_16S_succession_T5678$res$diff_abn
#Results 
res_out_ps_16S_succession_T5678 <- out_ps_16S_succession_T5678$res

# #Add Relative abundances for color gradient in later plotting
# #Transform to TSS in percentage
# #Relative abundances 
# phylo_main_no_cont_NoOutliers_scalled <- transform_sample_counts(phylo_main_no_cont_NoOutliers, function(x) 100 * x/sum(x))
# phylo_main_no_cont_NoOutliers_scalled_SW <- subset_samples(phylo_main_no_cont_NoOutliers_scalled, Environment=="Planktonic suspension")
# phylo_main_no_cont_NoOutliers_scalled_SW_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_scalled_SW, Day=="4")
# phylo_main_no_cont_NoOutliers_scalled_SW_day4 <- subset_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4, Genus!="Phaeobacter_inhibens")
# 
# phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order <- aggregate_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4, "Order")
# # phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order_filter <- filter_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order, function (x) {sum(x > 0) > 0}, prune=TRUE)
# sample_data(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order)$Treatment, levels = c("WT", "Control", "dTDA"))
# 
# RelAbs_SW_Order_day4 = as.data.frame(otu_table(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order))
# RelAbs_SW_Order_day4$taxon_id <- rownames(RelAbs_SW_Order_day4)
# 
# 
# 
# colnames(RelAbs_SW_Order_day4)
# RelAbs_SW_Order_day4_long <- reshape(RelAbs_SW_Order_day4, direction = "long",
#                                      v.names = "Rel_Abs",
#                                      varying = 1:ncol(RelAbs_SW_Order_day4)-1,
#                                      times=c(rep('Control',9), rep('dTDA', 9), rep('WT',8)),
# )
# 
# #Mean Rel_Abs per treatment per tax
# RelAbs_SW_Order_day4_long_mean <- RelAbs_SW_Order_day4_long %>%
#   group_by(taxon_id, time) %>%
#   dplyr::summarise(Rel_Abs = mean(Rel_Abs)) 
# 


#dplyr::filter(RelAbs_SW_Order_day4_long_mean, grepl("Chloroplast",taxon_id))

#### Dataframe Day 4 for logFC
df_fig1_16S_succession_T5678 = data.frame(res_out_ps_16S_succession_T5678$beta * res_out_ps_16S_succession_T5678$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_16S_succession_T5678 = data.frame(res_out_ps_16S_succession_T5678$se * res_out_ps_16S_succession_T5678$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_16S_succession_T5678)[-1] = paste0(colnames(df_fig2_16S_succession_T5678)[-1], "SD")
colnames(df_fig1_16S_succession_T5678) <- sub("Bryozoans", "", colnames(df_fig1_16S_succession_T5678))
colnames(df_fig2_16S_succession_T5678) <- sub("Bryozoans", "", colnames(df_fig2_16S_succession_T5678))


df_fig_16S_succession_T5678 = df_fig1_16S_succession_T5678 %>% left_join(df_fig2_16S_succession_T5678, by = "taxon_id")  %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_16S_succession_T5678

# Those that with positive log fold change is found more abundant together with bryozoans


###Run ANCOMBC differential analysis 16S Order####
#Gloom/aggregate to order level


ps_16S_succession_T5678_order <- aggregate_taxa(ps_16S_succession_T5678, "order")
sample_data(ps_16S_succession_T5678_order)$Bryozoans <- factor(sample_data(ps_16S_succession_T5678_order)$Bryozoans, levels = c("no", "yes"))


#Run ANCOMBC differential analysis
out_ps_16S_succession_T5678_order <- ancombc(phyloseq = ps_16S_succession_T5678_order, formula = "Bryozoans",
                                             p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                             group = "Bryozoans", struc_zero = TRUE, neg_lb = FALSE,
                                             tol = 1e-5, max_iter = 100, conserve = TRUE,
                                             alpha = 0.01, global = TRUE)



out_ps_16S_succession_T5678_order$res$diff_abn
#Results 
res_out_ps_16S_succession_T5678_order <- out_ps_16S_succession_T5678_order$res

# #Add Relative abundances for color gradient in later plotting
# #Transform to TSS in percentage
# #Relative abundances 
# phylo_main_no_cont_NoOutliers_scalled <- transform_sample_counts(phylo_main_no_cont_NoOutliers, function(x) 100 * x/sum(x))
# phylo_main_no_cont_NoOutliers_scalled_SW <- subset_samples(phylo_main_no_cont_NoOutliers_scalled, Environment=="Planktonic suspension")
# phylo_main_no_cont_NoOutliers_scalled_SW_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_scalled_SW, Day=="4")
# phylo_main_no_cont_NoOutliers_scalled_SW_day4 <- subset_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4, Genus!="Phaeobacter_inhibens")
# 
# phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order <- aggregate_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4, "Order")
# # phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order_filter <- filter_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order, function (x) {sum(x > 0) > 0}, prune=TRUE)
# sample_data(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order)$Treatment, levels = c("WT", "Control", "dTDA"))
# 
# RelAbs_SW_Order_day4 = as.data.frame(otu_table(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order))
# RelAbs_SW_Order_day4$taxon_id <- rownames(RelAbs_SW_Order_day4)
# 
# 
# 
# colnames(RelAbs_SW_Order_day4)
# RelAbs_SW_Order_day4_long <- reshape(RelAbs_SW_Order_day4, direction = "long",
#                                      v.names = "Rel_Abs",
#                                      varying = 1:ncol(RelAbs_SW_Order_day4)-1,
#                                      times=c(rep('Control',9), rep('dTDA', 9), rep('WT',8)),
# )
# 
# #Mean Rel_Abs per treatment per tax
# RelAbs_SW_Order_day4_long_mean <- RelAbs_SW_Order_day4_long %>%
#   group_by(taxon_id, time) %>%
#   dplyr::summarise(Rel_Abs = mean(Rel_Abs)) 
# 


#dplyr::filter(RelAbs_SW_Order_day4_long_mean, grepl("Chloroplast",taxon_id))

#### Dataframe Day 4 for logFC
df_fig1_16S_succession_T5678_order = data.frame(res_out_ps_16S_succession_T5678_order$beta * res_out_ps_16S_succession_T5678_order$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_16S_succession_T5678_order = data.frame(res_out_ps_16S_succession_T5678_order$se * res_out_ps_16S_succession_T5678_order$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_16S_succession_T5678_order)[-1] = paste0(colnames(df_fig2_16S_succession_T5678_order)[-1], "SD")
colnames(df_fig1_16S_succession_T5678_order) <- sub("Bryozoans", "", colnames(df_fig1_16S_succession_T5678_order))
colnames(df_fig2_16S_succession_T5678_order) <- sub("Bryozoans", "", colnames(df_fig2_16S_succession_T5678_order))


df_fig_16S_succession_T5678_order = df_fig1_16S_succession_T5678_order %>% left_join(df_fig2_16S_succession_T5678_order, by = "taxon_id")  %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 





###Run ANCOMBC differential analysis 16S phylum####
#Gloom/aggregate to phylum level


ps_16S_succession_T5678_phylum <- aggregate_taxa(ps_16S_succession_T5678, "phylum")
sample_data(ps_16S_succession_T5678_phylum)$Bryozoans <- factor(sample_data(ps_16S_succession_T5678_phylum)$Bryozoans, levels = c("no", "yes"))


#Run ANCOMBC differential analysis
out_ps_16S_succession_T5678_phylum <- ancombc(phyloseq = ps_16S_succession_T5678_phylum, formula = "Bryozoans",
                                             p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                             group = "Bryozoans", struc_zero = TRUE, neg_lb = FALSE,
                                             tol = 1e-5, max_iter = 100, conserve = TRUE,
                                             alpha = 0.01, global = TRUE)



out_ps_16S_succession_T5678_phylum$res$diff_abn
#Results 
res_out_ps_16S_succession_T5678_phylum <- out_ps_16S_succession_T5678_phylum$res

# #Add Relative abundances for color gradient in later plotting
# #Transform to TSS in percentage
# #Relative abundances 
# phylo_main_no_cont_NoOutliers_scalled <- transform_sample_counts(phylo_main_no_cont_NoOutliers, function(x) 100 * x/sum(x))
# phylo_main_no_cont_NoOutliers_scalled_SW <- subset_samples(phylo_main_no_cont_NoOutliers_scalled, Environment=="Planktonic suspension")
# phylo_main_no_cont_NoOutliers_scalled_SW_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_scalled_SW, Day=="4")
# phylo_main_no_cont_NoOutliers_scalled_SW_day4 <- subset_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4, Genus!="Phaeobacter_inhibens")
# 
# phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order <- aggregate_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4, "Order")
# # phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order_filter <- filter_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order, function (x) {sum(x > 0) > 0}, prune=TRUE)
# sample_data(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order)$Treatment, levels = c("WT", "Control", "dTDA"))
# 
# RelAbs_SW_Order_day4 = as.data.frame(otu_table(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order))
# RelAbs_SW_Order_day4$taxon_id <- rownames(RelAbs_SW_Order_day4)
# 
# 
# 
# colnames(RelAbs_SW_Order_day4)
# RelAbs_SW_Order_day4_long <- reshape(RelAbs_SW_Order_day4, direction = "long",
#                                      v.names = "Rel_Abs",
#                                      varying = 1:ncol(RelAbs_SW_Order_day4)-1,
#                                      times=c(rep('Control',9), rep('dTDA', 9), rep('WT',8)),
# )
# 
# #Mean Rel_Abs per treatment per tax
# RelAbs_SW_Order_day4_long_mean <- RelAbs_SW_Order_day4_long %>%
#   group_by(taxon_id, time) %>%
#   dplyr::summarise(Rel_Abs = mean(Rel_Abs)) 
# 


#dplyr::filter(RelAbs_SW_Order_day4_long_mean, grepl("Chloroplast",taxon_id))

#### Dataframe Day 4 for logFC
df_fig1_16S_succession_T5678_phylum = data.frame(res_out_ps_16S_succession_T5678_phylum$beta * res_out_ps_16S_succession_T5678_phylum$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_16S_succession_T5678_phylum = data.frame(res_out_ps_16S_succession_T5678_phylum$se * res_out_ps_16S_succession_T5678_phylum$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_16S_succession_T5678_phylum)[-1] = paste0(colnames(df_fig2_16S_succession_T5678_phylum)[-1], "SD")
colnames(df_fig1_16S_succession_T5678_phylum) <- sub("Bryozoans", "", colnames(df_fig1_16S_succession_T5678_phylum))
colnames(df_fig2_16S_succession_T5678_phylum) <- sub("Bryozoans", "", colnames(df_fig2_16S_succession_T5678_phylum))


df_fig_16S_succession_T5678_phylum = df_fig1_16S_succession_T5678_phylum %>% left_join(df_fig2_16S_succession_T5678_phylum, by = "taxon_id")  %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_16S_succession_T5678_phylum


# Those that with positive log fold change is found more abundant together with bryozoans


#Run ANCOMBC differential analysis 18S genus #####
out_ps_18S_succession_T5678_genus <- ancombc(phyloseq = ps_18S_succession_T5678_genus, formula = "Bryozoans",
                                             p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                             group = "Bryozoans", struc_zero = TRUE, neg_lb = FALSE,
                                             tol = 1e-5, max_iter = 100, conserve = TRUE,
                                             alpha = 0.01, global = TRUE)



out_ps_18S_succession_T5678_genus$res$diff_abn
#Results 
res_out_ps_18S_succession_T5678_genus <- out_ps_18S_succession_T5678_genus$res

# #Add Relative abundances for color gradient in later plotting
# #Transform to TSS in percentage
# #Relative abundances 
# phylo_main_no_cont_NoOutliers_scalled <- transform_sample_counts(phylo_main_no_cont_NoOutliers, function(x) 100 * x/sum(x))
# phylo_main_no_cont_NoOutliers_scalled_SW <- subset_samples(phylo_main_no_cont_NoOutliers_scalled, Environment=="Planktonic suspension")
# phylo_main_no_cont_NoOutliers_scalled_SW_day4 <- subset_samples(phylo_main_no_cont_NoOutliers_scalled_SW, Day=="4")
# phylo_main_no_cont_NoOutliers_scalled_SW_day4 <- subset_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4, Genus!="Phaeobacter_inhibens")
# 
# phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order <- aggregate_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4, "Order")
# # phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order_filter <- filter_taxa(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order, function (x) {sum(x > 0) > 0}, prune=TRUE)
# sample_data(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order)$Treatment <- factor(sample_data(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order)$Treatment, levels = c("WT", "Control", "dTDA"))
# 
# RelAbs_SW_Order_day4 = as.data.frame(otu_table(phylo_main_no_cont_NoOutliers_scalled_SW_day4_Order))
# RelAbs_SW_Order_day4$taxon_id <- rownames(RelAbs_SW_Order_day4)
# 
# 
# 
# colnames(RelAbs_SW_Order_day4)
# RelAbs_SW_Order_day4_long <- reshape(RelAbs_SW_Order_day4, direction = "long",
#                                      v.names = "Rel_Abs",
#                                      varying = 1:ncol(RelAbs_SW_Order_day4)-1,
#                                      times=c(rep('Control',9), rep('dTDA', 9), rep('WT',8)),
# )
# 
# #Mean Rel_Abs per treatment per tax
# RelAbs_SW_Order_day4_long_mean <- RelAbs_SW_Order_day4_long %>%
#   group_by(taxon_id, time) %>%
#   dplyr::summarise(Rel_Abs = mean(Rel_Abs)) 
# 


#dplyr::filter(RelAbs_SW_Order_day4_long_mean, grepl("Chloroplast",taxon_id))

#### Dataframe Day 4 for logFC
df_fig1_18S_succession_T5678_genus = data.frame(res_out_ps_18S_succession_T5678_genus$beta * res_out_ps_18S_succession_T5678_genus$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_18S_succession_T5678_genus = data.frame(res_out_ps_18S_succession_T5678_genus$se * res_out_ps_18S_succession_T5678_genus$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_18S_succession_T5678_genus)[-1] = paste0(colnames(df_fig2_18S_succession_T5678_genus)[-1], "SD")
colnames(df_fig1_18S_succession_T5678_genus) <- sub("Bryozoans", "", colnames(df_fig1_18S_succession_T5678_genus))
colnames(df_fig2_18S_succession_T5678_genus) <- sub("Bryozoans", "", colnames(df_fig2_18S_succession_T5678_genus))


df_fig_18S_succession_T5678_genus = df_fig1_18S_succession_T5678_genus %>% left_join(df_fig2_18S_succession_T5678_genus, by = "taxon_id")  %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_18S_succession_T5678_genus %>% filter(yes < 0) %>% arrange(desc(yes))
ps_18S_succession_T5678_genus_tax_dat <- data.frame(tax_table(ps_18S_succession_T5678_genus))
colnames(ps_18S_succession_T5678_genus_tax_dat)[7] <- "taxon_id"
df_fig_18S_succession_T5678_genus <- left_join(df_fig_18S_succession_T5678_genus, ps_18S_succession_T5678_genus_tax_dat, "taxon_id")

# Those that with positive log fold change is found more abundant together with bryozoans



df_fig_18S_succession_T5678_genus_no <- df_fig_18S_succession_T5678_genus %>% 
  ggplot(aes(x=yes, y=reorder(taxon_id, yes), label=taxon_id)) + geom_vline(xintercept = 0, colour="gray", size = 0.3) +
  geom_bar(stat = "identity", width = 0.7, aes(fill = (phylum)),
           position = position_dodge(width = 0.7), color="black") +
  labs(y = NULL, x="Log fold change",title = "ANCOM-BC diff. analysis\nbetween early and late biofilm microbiomes") +
  geom_errorbar(aes(xmin = yes - yesSD, xmax = yes + yesSD), width = 0.3,
                position = position_dodge(0.7), color = "black") +  
  theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 


legend_df_fig_18S_succession_T5678_genus <- cowplot::get_legend(df_fig_18S_succession_T5678_genus_no + theme(legend.position="right") )



tiff("../ANCOM-BC_160522_18S_genus.tiff", units="in",  width=13.08661, height=7.93701, res=300)

grid.arrange( df_fig_18S_succession_T5678_genus_no, 
              legend_df_fig_18S_succession_T5678_genus,
              nrow=1, ncol = 2,
              layout_matrix = rbind(c(1,2)),
              widths = c(4, 0.8))

dev.off()


###Run ANCOMBC differential analysis 18S no aggregation####
#Gloom/aggregate to class level

sample_data(ps_18S_succession_T5678)$Bryozoans <- factor(sample_data(ps_18S_succession_T5678)$Bryozoans, levels = c("no", "yes"))


#Run ANCOMBC differential analysis
out_ps_18S_succession_T5678 <- ancombc(phyloseq = ps_18S_succession_T5678, formula = "Bryozoans",
                                             p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                                             group = "Bryozoans", struc_zero = TRUE, neg_lb = FALSE,
                                             tol = 1e-5, max_iter = 100, conserve = TRUE,
                                             alpha = 0.01, global = TRUE)



out_ps_18S_succession_T5678$res$diff_abn
#Results 
res_out_ps_18S_succession_T5678 <- out_ps_18S_succession_T5678$res


#### Dataframe 18S between early and late stage for logFC
df_fig1_18S_succession_T5678 = data.frame(res_out_ps_18S_succession_T5678$beta * res_out_ps_18S_succession_T5678$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
df_fig2_18S_succession_T5678 = data.frame(res_out_ps_18S_succession_T5678$se * res_out_ps_18S_succession_T5678$diff_abn, check.names = FALSE) %>%
  rownames_to_column("taxon_id")
colnames(df_fig2_18S_succession_T5678)[-1] = paste0(colnames(df_fig2_18S_succession_T5678)[-1], "SD")
colnames(df_fig1_18S_succession_T5678) <- sub("Bryozoans", "", colnames(df_fig1_18S_succession_T5678))
colnames(df_fig2_18S_succession_T5678) <- sub("Bryozoans", "", colnames(df_fig2_18S_succession_T5678))


df_fig_18S_succession_T5678 = df_fig1_18S_succession_T5678 %>% left_join(df_fig2_18S_succession_T5678, by = "taxon_id")  %>%
  rowwise() %>%
  filter(any(across(is.numeric) != 0)) 

df_fig_18S_succession_T5678


df_fig_18S_succession_T5678 %>% filter(yes < 0) %>% arrange(desc(yes))
ps_18S_succession_T5678_tax_dat <- data.frame(tax_table(ps_18S_succession_T5678))
ps_18S_succession_T5678_tax_dat <- ps_18S_succession_T5678_tax_dat %>% mutate(taxon_id = rownames(ps_18S_succession_T5678_tax_dat) )

df_fig_18S_succession_T5678 <- left_join(df_fig_18S_succession_T5678, ps_18S_succession_T5678_tax_dat, "taxon_id")


df_fig_18S_succession_T5678 %>% filter(-3>yes) %>%
  ggplot(aes(x=yes, y=reorder(taxon_id, yes), label=taxon_id)) + geom_vline(xintercept = 0, colour="gray", size = 0.3) +
  geom_bar(stat = "identity", width = 0.7, aes(fill = (phylum), color = phylum),
           position = position_dodge(width = 0.7)) +
  labs(y = NULL, x="Log fold change",title = "ANCOM-BC diff. analysis\nbetween early and late biofilm microbiomes") +
  geom_errorbar(aes(xmin = yes - yesSD, xmax = yes + yesSD), width = 0.3,
                position = position_dodge(0.7), color = "black") +  
  theme_bw()

df_fig_18S_succession_T5678 %>% filter(2<yes) %>%
  ggplot(aes(x=yes, y=reorder(taxon_id, yes), label=taxon_id)) + geom_vline(xintercept = 0, colour="gray", size = 0.3) +
  geom_bar(stat = "identity", width = 0.7, aes(fill = (genus), color = phylum),
           position = position_dodge(width = 0.7)) +
  labs(y = NULL, x="Log fold change",title = "ANCOM-BC diff. analysis\nbetween early and late biofilm microbiomes") +
  geom_errorbar(aes(xmin = yes - yesSD, xmax = yes + yesSD), width = 0.3,
                position = position_dodge(0.7), color = "black") +  
  theme_bw()



df_fig_18S_succession_T5678_no <- df_fig_18S_succession_T5678 %>% filter(-3>yes) %>%
  ggplot(aes(x=yes, y=reorder(taxon_id, yes), label=taxon_id)) + geom_vline(xintercept = 0, colour="gray", size = 0.3) +
  geom_bar(stat = "identity", width = 0.7, aes(fill = (phylum)),
           position = position_dodge(width = 0.7), color="black") +
  labs(y = NULL, x="Log fold change",title = "ANCOM-BC diff. analysis\nbetween early and late biofilm microbiomes") +
  geom_errorbar(aes(xmin = yes - yesSD, xmax = yes + yesSD), width = 0.3,
                position = position_dodge(0.7), color = "black") +  
  theme_bw() +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 


legend_df_fig_18S_succession_T5678_genus <- cowplot::get_legend(df_fig_18S_succession_T5678_genus_no + theme(legend.position="right") )



tiff("../ANCOM-BC_160522_18S_genus.tiff", units="in",  width=13.08661, height=7.93701, res=300)

grid.arrange( df_fig_18S_succession_T5678_genus_no, 
              legend_df_fig_18S_succession_T5678_genus,
              nrow=1, ncol = 2,
              layout_matrix = rbind(c(1,2)),
              widths = c(4, 0.8))

dev.off()


# Those that with positive log fold change is found more abundant together with bryozoans
