rm(list=ls() )
setwd("/Users/pernillekjersgaardbech/Documents/JH21_Biofilm/AD_amplicons/")
library(rRDPData)
library(tidyverse)
library("devtools")
library("dplyr")
library("Biostrings")
library("phyloseq")
library(ggplot2)
library(vegan)
library(venn)
library(plotly)
library(RColorBrewer)
#install_github("Russel88/MicEco")
#library(MicEco)
library(ggbeeswarm)
library(gridExtra)

ps_16S <- readRDS('../16S amplicon/PS-16S-decon.noEuk.noneg.rds')
ps_AD <- readRDS('ps_AD')
### Evaluate rarefraction curves
# rcurves <- rcurve(ps_AD)
# 
# ggplot(data=rcurves, aes(x=Reads, y=Richness, group=Sample, col = ))+
#   geom_line()



#Remove Phaeobacter strain AD ASVs and Seawater-T1-3 because it is an outlier
ps_AD <- subset_samples(ps_AD, Subject != "JH21" & Sample != "Seawater-T1-3")
ps_AD_succession <- subset_samples(ps_AD, Subject == "Succession")
ps_16S_succession <- subset_samples(ps_16S, element.type == "bioelement")
sample_data(ps_16S)

#Filtering all taxa wiuth read sum < 100
ps_AD_filter <- filter_taxa(ps_AD, function (x) {sum(x > 100) > 0}, prune=TRUE)
ps_16S_filter <- filter_taxa(ps_16S, function (x) {sum(x > 100) > 0}, prune=TRUE)
ps_AD_succession_filter <- filter_taxa(ps_AD_succession, function (x) {sum(x > 100) > 0}, prune=TRUE)
ps_16S_succession_filter <- filter_taxa(ps_16S_succession, function (x) {sum(x > 100) > 0}, prune=TRUE)


#Scalling to total sum TSS
ps_AD_filter_norm = transform_sample_counts(ps_AD_filter, function(x) 100000 * x/sum(x+0.1))
ps_16S_filter_norm = transform_sample_counts(ps_16S_filter, function(x) 100000 * x/sum(x+0.1))
ps_AD_succession_filter_norm = transform_sample_counts(ps_AD_succession_filter, function(x) 100000 * x/sum(x+0.1))
ps_16S_succession_filter_norm = transform_sample_counts(ps_16S_succession_filter, function(x) 100000 * x/sum(x+0.1))


#Make phyloseq object to matrix - take out of Phyloseq and treat with normal vegan functions
AD_filter_norm_dat <- as.data.frame((as(otu_table(ps_AD_filter_norm), "matrix")))
AD_filter_norm_dat <- round(AD_filter_norm_dat)
str(AD_filter_norm_dat)


AD_succession_filter_norm_dat <- as.data.frame((as(otu_table(ps_AD_succession_filter_norm), "matrix")))
AD_succession_filter_norm_dat <- round(AD_succession_filter_norm_dat)
str(AD_filter_norm_dat)

P16S_filter_norm_dat <- as.data.frame((as(otu_table(ps_16S_filter_norm), "matrix")))
P16S_filter_norm_dat <- as.data.frame(round(t(P16S_filter_norm_dat)))
str(P16S_filter_norm_dat)

P16S_succession_filter_norm_dat <- as.data.frame((as(otu_table(ps_16S_succession_filter_norm), "matrix")))
P16S_succession_filter_norm_dat <- as.data.frame(round(t(P16S_succession_filter_norm_dat)))
str(P16S_succession_filter_norm_dat)

#remove empty rows
AD_filter_norm_dat_clean <- AD_filter_norm_dat[rowSums(AD_filter_norm_dat)!=0,]
AD_succession_filter_norm_dat_clean <- AD_succession_filter_norm_dat[rowSums(AD_succession_filter_norm_dat)!=0,]
P16S_filter_norm_dat_clean <- P16S_filter_norm_dat[rowSums(P16S_filter_norm_dat)!=0,]
P16S_succession_filter_norm_dat_clean <- P16S_succession_filter_norm_dat[rowSums(P16S_succession_filter_norm_dat)!=0,]


#Match metadata table with rownames of AD_filter_norm_dat_clean

AD_filter_norm_meta_dat <- data.frame(sample_data(ps_AD_filter_norm))
AD_filter_norm_meta_dat <- AD_filter_norm_meta_dat[rownames(AD_filter_norm_dat_clean), ]


AD_succession_filter_norm_meta_dat <- data.frame(sample_data(ps_AD_succession_filter_norm))
AD_succession_filter_norm_meta_dat <- AD_succession_filter_norm_meta_dat[rownames(AD_succession_filter_norm_dat_clean), ]


P16S_filter_norm_meta_dat <- data.frame(sample_data(ps_16S_filter_norm))
P16S_filter_norm_meta_dat <- P16S_filter_norm_meta_dat[rownames(P16S_filter_norm_dat_clean), ]


P16S_succession_filter_norm_meta_dat <- data.frame(sample_data(ps_16S_succession_filter_norm))
P16S_succession_filter_norm_meta_dat <- P16S_succession_filter_norm_meta_dat[rownames(P16S_succession_filter_norm_dat_clean), ]


# samples_to_move <-as.vector(which(is.na(rowSums(AD_filter_norm_dat))))

# AD_filter_norm_meta_dat <- AD_filter_norm_meta_dat[-samples_to_move,]
# 
# AD_filter_norm_dat <- round((AD_filter_norm_dat[-samples_to_move,]))

#AD_norm_dat_t <- t(AD_norm_dat)

# #Squar-root transform matrix
 sqrt_AD_filter_norm_dat_clean= sqrt(AD_filter_norm_dat_clean)
 sqrt_AD_succession_filter_norm_dat_clean= sqrt(AD_succession_filter_norm_dat_clean)
 sqrt_P16S_filter_norm_dat_clean= sqrt(P16S_filter_norm_dat_clean)
 sqrt_P16S_succession_filter_norm_dat_clean= sqrt(P16S_succession_filter_norm_dat_clean)
 
 
 # sqrt_noPhaeo_genus_norm_dat_t = sqrt(noPhaeo_genus_norm_dat_t)

#Bray-curtis distance matrix
AD_dist_bray = as.matrix((vegdist((sqrt_AD_filter_norm_dat_clean), "bray")))
str((AD_dist_bray))
AD_succession_dist_bray = as.matrix((vegdist((sqrt_AD_succession_filter_norm_dat_clean), "bray")))
str((AD_dist_bray))
P16S_dist_bray = as.matrix((vegdist((sqrt_P16S_filter_norm_dat_clean), "bray")))
str(P16S_dist_bray)
P16S_succession_dist_bray = as.matrix((vegdist((sqrt_P16S_succession_filter_norm_dat_clean), "bray")))
str(P16S_succession_dist_bray)

#betadisper #####
#bdisp_time <- betadisper(vegdist(sqrt_genus_norm_dat, "bray"), AD_filter_norm_meta_dat$Day)
# bdisp_treatment <- betadisper(vegdist(sqrt_noPhaeo_genus_norm_dat_t, "bray"), genus_meta_dat$Treatment)
# bdisp_environment <- betadisper(vegdist(sqrt_noPhaeo_genus_norm_dat_t, "bray"), genus_meta_dat$Environment)
# 
# AD_meta_dat_vec <- as.factor((paste( AD_filter_norm_meta_dat$Day, sep = "_")))
# bdisp_all <- betadisper(vegdist(AD_filter_norm_dat_clean, "bray"), AD_meta_dat_vec)


## Perform test
# anova(bdisp_time)
# anova(bdisp_treatment)
# anova(bdisp_environment)

# anova(bdisp_all)
# ## Permutation test for F
# permutest(bdisp_all, pairwise = TRUE, permutations = 999)
# # 
# # ## Tukey's Honest Significant Differences
# (bdisp.HSD <- TukeyHSD(bdisp_all))
# plot(bdisp.HSD)



#PERMANOVA (Adonis) ----
# dim(AD_dist_bray)
# dim(AD_dist_bray)
# adonis_AD <- adonis(AD_dist_bray ~ Subject, AD_filter_norm_meta_dat)
# adonis_AD_dat <- as.data.frame(adonis_AD$aov.tab)

# NDMS AD all subjects -----
NMDS_ord_bray = metaMDS(AD_dist_bray, k=3)
#NMDS_ord_bray_spec = metaMDS(sqrt_genus_norm_dat_t, distance = "bray", autotransform = FALSE, k=3)
str(NMDS_ord_bray)

stressplot(NMDS_ord_bray)

#build a data frame with NMDS coordinates and metadata
NMDS1 = NMDS_ord_bray$points[,1]
NMDS2 = NMDS_ord_bray$points[,2]
NMDS3 = NMDS_ord_bray$points[,3]


# str(NMDS_ord_bray)
# NMDS_ord_bray$species


NMDS = data.frame(NMDS1 = NMDS1, NMDS2 = NMDS2, NMDS2 = NMDS3, Subject = AD_filter_norm_meta_dat$Subject, Day = AD_filter_norm_meta_dat$Day)

#Set color theme
#fef0d9
#fdcc8a
#fc8d59
#d7301f
# cols_time <- c("0" = "#f7f7f7", "1" = "#fdcc8a", "4" = "#fc8d59", "10" = "#d7301f")
# NMDS$Day <- factor(NMDS$Day, levels = c("0", "1", "4","10"))

# cols_time <- c( "1" = "#fdcc8a", "4" = "#fc8d59", "10" = "#d7301f")
# NMDS$Day <- factor(NMDS$Day, levels = c("1", "4","10"))


p_subject_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Subject)) +
  #stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
   geom_point(shape = 21, alpha = 0.25, size = 2) 

p_time_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Day)) +
  #stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)
  # scale_colour_manual(
  #   values = cols_time,
  #   aesthetics = c("colour", "fill") 
  # )  + xlim(-0.58,0.58) + ylim(-0.4,0.4) + theme_bw() +
  # theme(#axis.line = element_line(color='black'),
  #   plot.background = element_blank(),
  #   panel.grid.major = element_blank(),
  #   panel.grid.minor = element_blank(), legend.position = "none") 

p_subject_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Subject)) +
  #stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2) 

p_time_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Day)) +
  #stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)

library(gridExtra)
grid.arrange(p_subject_NMDS12, p_time_NMDS12, p_subject_NMDS23, p_time_NMDS23, ncol=2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,4)),
             widths = c(2.7, 2.7), heights = c(2.7, 2.7))




# NDMS AD only succession -----
NMDS_ord_bray = metaMDS(AD_succession_dist_bray, k=3)
#NMDS_ord_bray_spec = metaMDS(sqrt_genus_norm_dat_t, distance = "bray", autotransform = FALSE, k=3)
str(NMDS_ord_bray)

stressplot(NMDS_ord_bray)

#build a data frame with NMDS coordinates and metadata
NMDS1 = NMDS_ord_bray$points[,1]
NMDS2 = NMDS_ord_bray$points[,2]
NMDS3 = NMDS_ord_bray$points[,3]


# str(NMDS_ord_bray)
# NMDS_ord_bray$species


NMDS = data.frame(NMDS1 = NMDS1, NMDS2 = NMDS2, NMDS2 = NMDS3, 
                  Chla_suc = AD_succession_filter_norm_meta_dat$Chla_suc, 
                  Chla_ass = AD_succession_filter_norm_meta_dat$Chla_ass, 
                  Day = AD_succession_filter_norm_meta_dat$Day)

#Set color theme
#fef0d9
#fdcc8a
#fc8d59
#d7301f
# cols_time <- c("0" = "#f7f7f7", "1" = "#fdcc8a", "4" = "#fc8d59", "10" = "#d7301f")
# NMDS$Day <- factor(NMDS$Day, levels = c("0", "1", "4","10"))

# cols_time <- c( "1" = "#fdcc8a", "4" = "#fc8d59", "10" = "#d7301f")
# NMDS$Day <- factor(NMDS$Day, levels = c("1", "4","10"))




p_time_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Day)) +
  #stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)

p_Chla_suc_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Chla_suc)) +
  #stat_ellipse() +
  scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)

p_Chla_ass_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Chla_ass)) +
  #stat_ellipse() +
  scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)
# scale_colour_manual(
#   values = cols_time,
#   aesthetics = c("colour", "fill") 
# )  + xlim(-0.58,0.58) + ylim(-0.4,0.4) + theme_bw() +
# theme(#axis.line = element_line(color='black'),
#   plot.background = element_blank(),
#   panel.grid.major = element_blank(),
#   panel.grid.minor = element_blank(), legend.position = "none") 


p_time_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Day)) +
  #stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)


p_Chla_suc_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Chla_suc)) +
  #stat_ellipse() +
  scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)

p_Chla_ass_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Chla_ass)) +
  #stat_ellipse() +
  scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)


grid.arrange(p_time_NMDS12, p_time_NMDS23, p_Chla_suc_NMDS12, p_Chla_suc_NMDS23, ncol=2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,4)),
             widths = c(2.7, 2.7), heights = c(2.7, 2.7))

grid.arrange(p_time_NMDS12, p_time_NMDS23, p_Chla_ass_NMDS12, p_Chla_ass_NMDS23, ncol=2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,4)),
             widths = c(2.7, 2.7), heights = c(2.7, 2.7))



# NDMS 16S  -----
NMDS_ord_bray = metaMDS(P16S_dist_bray, k=3)
#NMDS_ord_bray_spec = metaMDS(sqrt_genus_norm_dat_t, distance = "bray", autotransform = FALSE, k=3)
str(NMDS_ord_bray)

stressplot(NMDS_ord_bray)

#build a data frame with NMDS coordinates and metadata
NMDS1 = NMDS_ord_bray$points[,1]
NMDS2 = NMDS_ord_bray$points[,2]
NMDS3 = NMDS_ord_bray$points[,3]


# str(NMDS_ord_bray)
# NMDS_ord_bray$species


NMDS = data.frame(NMDS1 = NMDS1, NMDS2 = NMDS2, NMDS2 = NMDS3, 
                  Chla_suc = P16S_filter_norm_meta_dat$chla, 
                  Day = P16S_filter_norm_meta_dat$day,
                  Subject = P16S_filter_norm_meta_dat$element.type)

#Set color theme
#fef0d9
#fdcc8a
#fc8d59
#d7301f
# cols_time <- c("0" = "#f7f7f7", "1" = "#fdcc8a", "4" = "#fc8d59", "10" = "#d7301f")
# NMDS$Day <- factor(NMDS$Day, levels = c("0", "1", "4","10"))

# cols_time <- c( "1" = "#fdcc8a", "4" = "#fc8d59", "10" = "#d7301f")
# NMDS$Day <- factor(NMDS$Day, levels = c("1", "4","10"))




p_time_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Day)) +
  #stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)

p_Chla_suc_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Chla_suc)) +
  #stat_ellipse() +
  scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)

p_Chla_ass_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Subject)) +
  #stat_ellipse() +
  #scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)
# scale_colour_manual(
#   values = cols_time,
#   aesthetics = c("colour", "fill") 
# )  + xlim(-0.58,0.58) + ylim(-0.4,0.4) + theme_bw() +
# theme(#axis.line = element_line(color='black'),
#   plot.background = element_blank(),
#   panel.grid.major = element_blank(),
#   panel.grid.minor = element_blank(), legend.position = "none") 


p_time_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Day)) +
  #stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)


p_Chla_suc_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Chla_suc)) +
  #stat_ellipse() +
  scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)

p_Chla_ass_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Chla_ass)) +
  #stat_ellipse() +
  scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)


grid.arrange(p_time_NMDS12, p_time_NMDS23, p_Chla_suc_NMDS12, p_Chla_suc_NMDS23, ncol=2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,4)),
             widths = c(2.7, 2.7), heights = c(2.7, 2.7))








# NDMS 16S  succession-----
NMDS_ord_bray = metaMDS(P16S_succession_dist_bray, k=3)
#NMDS_ord_bray_spec = metaMDS(sqrt_genus_norm_dat_t, distance = "bray", autotransform = FALSE, k=3)
str(NMDS_ord_bray)

stressplot(NMDS_ord_bray)

#build a data frame with NMDS coordinates and metadata
NMDS1 = NMDS_ord_bray$points[,1]
NMDS2 = NMDS_ord_bray$points[,2]
NMDS3 = NMDS_ord_bray$points[,3]


# str(NMDS_ord_bray)
# NMDS_ord_bray$species


NMDS = data.frame(NMDS1 = NMDS1, NMDS2 = NMDS2, NMDS2 = NMDS3, 
                  Chla_suc = P16S_succession_filter_norm_meta_dat$chla, 
                  Day = P16S_succession_filter_norm_meta_dat$day,
                  Subject = P16S_succession_filter_norm_meta_dat$element.type)

#Set color theme
#fef0d9
#fdcc8a
#fc8d59
#d7301f
# cols_time <- c("0" = "#f7f7f7", "1" = "#fdcc8a", "4" = "#fc8d59", "10" = "#d7301f")
# NMDS$Day <- factor(NMDS$Day, levels = c("0", "1", "4","10"))

# cols_time <- c( "1" = "#fdcc8a", "4" = "#fc8d59", "10" = "#d7301f")
# NMDS$Day <- factor(NMDS$Day, levels = c("1", "4","10"))




p_time_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Day)) +
  #stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)

p_Chla_suc_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Chla_suc)) +
  #stat_ellipse() +
  scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)

p_Chla_ass_NMDS12 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, fill = Subject)) +
  #stat_ellipse() +
  #scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)
# scale_colour_manual(
#   values = cols_time,
#   aesthetics = c("colour", "fill") 
# )  + xlim(-0.58,0.58) + ylim(-0.4,0.4) + theme_bw() +
# theme(#axis.line = element_line(color='black'),
#   plot.background = element_blank(),
#   panel.grid.major = element_blank(),
#   panel.grid.minor = element_blank(), legend.position = "none") 


p_time_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Day)) +
  #stat_ellipse() +
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)


p_Chla_suc_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Chla_suc)) +
  #stat_ellipse() +
  scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)

p_Chla_ass_NMDS23 <- ggplot(NMDS, aes(x=NMDS3, y=NMDS2, fill = Chla_ass)) +
  #stat_ellipse() +
  scale_fill_gradient(low = "yellow", high = "#006400") + 
  theme_bw() +
  labs(title = "")  + 
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2)


grid.arrange(p_time_NMDS12, p_time_NMDS23, p_Chla_suc_NMDS12, p_Chla_suc_NMDS23, ncol=2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,4)),
             widths = c(2.7, 2.7), heights = c(2.7, 2.7))








#Richness ####
# sample_data(ps_AD_filter)
# plot_richness(ps_AD_filter, x="Day", measures=c("Shannon", "Observed"), color="Subject")


#Observed genera ####
Observed_16S_filter_norm <- data.frame(t(estimateR((round(t(otu_table(ps_16S_filter_norm)) )))))
Observed_16S_filter_norm_tab <- data.frame(cbind(Obs = Observed_16S_filter_norm$S.obs, Sample=rownames(Observed_16S_filter_norm) ))
colnames(Observed_16S_filter_norm_tab)[1] <- "Obs_16S"


Observed_AD_filter_norm <- data.frame(t(estimateR((round(otu_table(ps_AD_filter_norm)) ))))
Observed_AD_filter_norm_tab <- data.frame(cbind(Obs = Observed_AD_filter_norm$S.obs, Day=sample_data(ps_AD_filter_norm)$Day, Timepoint=sample_data(ps_AD_filter_norm)$Timepoint, Subject=sample_data(ps_AD_filter_norm)$element.type), Rep_tec=sample_data(ps_AD_filter_norm)$Rep, Sample=rownames(Observed_AD_filter_norm) )
colnames(Observed_AD_filter_norm_tab)[1] <- "Obs_AD"

#Merge 16S to AD tab to normalize the complexitiy or do correlations
Observed_ADand16S_filter_norm_tab <- left_join(Observed_AD_filter_norm_tab, Observed_16S_filter_norm_tab, "Sample")

Observed_ADand16S_filter_norm_tab$Obs_AD <- as.numeric(Observed_ADand16S_filter_norm_tab$Obs_AD)
Observed_ADand16S_filter_norm_tab$Obs_16S <- as.numeric(Observed_ADand16S_filter_norm_tab$Obs_16S)
Observed_ADand16S_filter_norm_tab <- left_join(Observed_ADand16S_filter_norm_tab, sample_data(ps_AD_filter_norm)[,c(1,6)], "Sample")



#First look for correlation pattern
Observed_ADand16S_filter_norm_tab %>%   ggplot(aes(y=Obs_AD, x=Obs_16S, color = Day)) +
  geom_point(alpha=0.5) +
  scale_size(range = c(0.1, 10))

#ITS NOT CORRELATED

#Plot AD complexity after normanilization with 16S complexity
Observed_ADand16S_filter_norm_tab <- Observed_ADand16S_filter_norm_tab %>% 
  mutate(
  AD_norm = Obs_AD /  Obs_16S
  ) 

Observed_ADand16S_filter_norm_tab$AD_norm[is.nan(Observed_ADand16S_filter_norm_tab$AD_norm)] <- 0
Observed_ADand16S_filter_norm_tab$AD_norm[is.infinite(Observed_ADand16S_filter_norm_tab$AD_norm)] <- 0
Observed_ADand16S_filter_norm_tab$Day <- as.numeric(Observed_ADand16S_filter_norm_tab$Day)
Observed_ADand16S_filter_norm_tab$Timepoint <- factor(Observed_ADand16S_filter_norm_tab$Timepoint, levels = c("0","T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12"))


str(Observed_ADand16S_filter_norm_tab)
#


#########Means and standard deviations#########
# Observed_AD_filter_norm_tab_stat <- Observed_AD_filter_norm_tab %>% group_by(Subject, Day) %>%
#   dplyr::summarise(n=n(),'Observed richness' = mean(Obs), sd = sd(Obs), max=max(Obs))


#Define colors
#colors <- c("Seawater at Day 0" ="grey", "Control" = "#80cdc1", "WT" = "#a6611a", "dTDA" = "#dfc27d")

## Beeswarm Plot -------

ggplot(Observed_ADand16S_filter_norm_tab, aes(x = Timepoint, y = AD_norm, col = Subject)) + 
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=2) + ylim(0,2)

p_Obs_16S <- ggplot(Observed_ADand16S_filter_norm_tab, aes(x = Timepoint, y = Obs_16S, col = Subject)) + 
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=2)

p_Obs_AD <- ggplot(Observed_ADand16S_filter_norm_tab, aes(x = Timepoint, y = Obs_AD, col = Subject)) + 
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=2)

grid.arrange(p_AD_norm, p_Obs_16S, p_Obs_AD, ncol=3)

  # stat_summary(aes(group = Treatment), fun.data=mean_sdl, fun.args = list(mult=1), 
  #              geom="errorbar", color="black", width=0.2, position=position_dodge(0.8)) +
  # stat_summary(aes(group = Treatment),fun=mean, geom="point", color="black", position=position_dodge(0.8)) +
  # facet_grid(cols = vars(Environment), scales="free_x") + 
  # xlab("Day") + ylab("Observed richness") + #guides(fill = FALSE) +
  # scale_colour_manual(
  #   values = colors,
  #   aesthetics = c("colour") 
  # ) + theme_bw() +
  # theme(#axis.line = element_line(color='black'),
  #   plot.background = element_blank(),
  #   panel.grid.major = element_blank(),
  #   panel.grid.minor = element_blank(), legend.position = "bottom",  text = element_text(size=10)) 
  # 
  # 


