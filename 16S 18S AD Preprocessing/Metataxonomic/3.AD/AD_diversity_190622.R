setwd("/Users/pernillekjersgaardbech/Documents/JH21_Biofilm/AD_amplicons/")


######################################
## PREPROCESSING OF PHYLOSEQOBJECTS ##
######################################


##############
## Packages ##
##############

library(tidyverse)
library("devtools")
library("dplyr")
library("Biostrings")
library("phyloseq")
library(ggplot2)
library(vegan)
library(venn)
library(plotly)
library(MicEco)
library(ggbeeswarm)
library(gridExtra)
library(ggpubr)

###############
## Load data ##
###############

ps_AD <- readRDS('../AD_amplicons/ps_AD_new_160622')


###################
## Preprocessing ##
###################

## Remove Phaeobacter strain from AD ASVs and Seawater-T1-3 because it is an outlier ##
ps_AD <- subset_samples(ps_AD, Subject != "JH21" & sample_names(ps_AD) != "Seawater-T1-3")
# Remove outliers: "Succession-T1-6" "Succession-T1-8" "Succession-T2-9" "Succession-T2-7" from 16S phyloobject
# ps_16S_outliers <- c("Succession-T1-6","Succession-T1-8","Succession-T2-9","Succession-T2-7")

## Subset to only bio.element study ##
ps_AD_succession <- subset_samples(ps_AD, Subject %in% c("Succession", "Bryozor"))
#try filter only ADs mapped to MAGs
#ADs_to_filter
#ps_AD_succession <- subset_taxa(ps_AD, AD_ASV %in% ADs_to_filter)
#colnames(tax_table(ps_AD_succession))

## #Remove ASV's below 0.01 %
ps_AD_pct = transform_sample_counts(ps_AD, function(x) 100 * x/sum(x))
ps_AD_succession_pct = transform_sample_counts(ps_AD_succession, function(x) 100 * x/sum(x))

ps_AD_filter <- filter_taxa(ps_AD_pct, function (x) {sum(x > 0.01) > 0}, prune=TRUE)
#Only bioelements PS
ps_AD_succession_filter <- filter_taxa(ps_AD_succession_pct, function (x)  {sum(x > 0.01) > 0}, prune=TRUE)
## Scalling to total sum TSS ##
#All sample types PS
ps_AD_filter_norm = transform_sample_counts(ps_AD_filter, function(x) 100000* x/sum(x+0.1))

#Only bioelements PS
ps_AD_succession_filter_norm = transform_sample_counts(ps_AD_succession_filter, function(x) 100000 * x/sum(x+0.1))




## Make phyloseq object to matrix ######
#AD - all sample types 
AD_filter_norm_dat <- as.data.frame((as(otu_table(ps_AD_filter_norm), "matrix")))
AD_filter_norm_dat <-  as.data.frame(round(t(AD_filter_norm_dat)))
AD_filter_norm_dat_clean <- AD_filter_norm_dat[rowSums(AD_filter_norm_dat)!=0,] # remove empty rows
str(AD_filter_norm_dat)
#AD - Only bioelements 
AD_succession_filter_norm_dat <- as.data.frame((as(otu_table(ps_AD_succession_filter_norm), "matrix")))
AD_succession_filter_norm_dat <-as.data.frame(round(t(AD_succession_filter_norm_dat)))
AD_succession_filter_norm_dat_clean <- AD_succession_filter_norm_dat[rowSums(AD_succession_filter_norm_dat)!=0,] #remove empty rows
str(AD_succession_filter_norm_dat_clean)


## Match metadata table with rownames of XX_filter_norm_dat_clean ##
#AD - all sample types PS
AD_filter_norm_meta_dat <- data.frame(sample_data(ps_AD_filter_norm))
AD_filter_norm_meta_dat <- AD_filter_norm_meta_dat[rownames(AD_filter_norm_dat_clean), ]
#AD - Only bioelements PS
AD_succession_filter_norm_meta_dat <- data.frame(sample_data(ps_AD_succession_filter_norm))
AD_succession_filter_norm_meta_dat <- AD_succession_filter_norm_meta_dat[rownames(AD_succession_filter_norm_dat_clean), ]


##############
## Richness ##
##############

## Estimate Richness ##
Observed_AD_filter_norm <- data.frame(t(estimateR((round(t(otu_table(ps_AD_filter_norm))) ))))
str(Observed_AD_filter_norm)

shannon_AD_filter_norm <- data.frame(diversity((round(t(otu_table(ps_AD_filter_norm)))), index = "shannon"))

str(shannon_AD_filter_norm)
Observed_AD_filter_norm_tab <- data.frame(cbind(Obs = Observed_AD_filter_norm$S.obs,shannon_AD_filter_norm$shannon,
                                                Day=sample_data(ps_AD_filter_norm)$Day, 
                                                Timepoint=sample_data(ps_AD_filter_norm)$Timepoint, 
                                                Subject=sample_data(ps_AD_filter_norm)$Subject), 
                                          Rep_tec=sample_data(ps_AD_filter_norm)$Rep, 
                                          Sample=rownames(Observed_AD_filter_norm) )
str(Observed_AD_filter_norm_tab)
colnames(Observed_AD_filter_norm_tab)[1] <- "Obs_AD"
colnames(Observed_AD_filter_norm_tab)[2] <- "Shannon_AD"
Observed_AD_filter_norm_tab$Obs_AD <- as.numeric(Observed_AD_filter_norm_tab$Obs_AD)
Observed_AD_filter_norm_tab$Shannon_AD <- as.numeric(Observed_AD_filter_norm_tab$Shannon_AD)

## Get Day to real days and not timepoints
Observed_AD_filter_norm_tab$Time = 0
Observed_AD_filter_norm_tab$Time[Observed_AD_filter_norm_tab$Day==0] <- 0
Observed_AD_filter_norm_tab$Time[Observed_AD_filter_norm_tab$Day==1] <- 3
Observed_AD_filter_norm_tab$Time[Observed_AD_filter_norm_tab$Day==2] <- 7
Observed_AD_filter_norm_tab$Time[Observed_AD_filter_norm_tab$Day==3] <- 10
Observed_AD_filter_norm_tab$Time[Observed_AD_filter_norm_tab$Day==4] <- 15
Observed_AD_filter_norm_tab$Time[Observed_AD_filter_norm_tab$Day==5] <- 23
Observed_AD_filter_norm_tab$Time[Observed_AD_filter_norm_tab$Day==6] <- 29
Observed_AD_filter_norm_tab$Time[Observed_AD_filter_norm_tab$Day==7] <- 44
Observed_AD_filter_norm_tab$Time[Observed_AD_filter_norm_tab$Day==8] <- 57
Observed_AD_filter_norm_tab$Time[Observed_AD_filter_norm_tab$Day==9] <- 71
Observed_AD_filter_norm_tab$Time[Observed_AD_filter_norm_tab$Day==10] <- 85
Observed_AD_filter_norm_tab$Time[Observed_AD_filter_norm_tab$Day==11] <- 99
Observed_AD_filter_norm_tab$Time[Observed_AD_filter_norm_tab$Day==12] <- 113

Observed_AD_filter_norm_tab$Time <- as.numeric(Observed_AD_filter_norm_tab$Time)

str(Observed_AD_filter_norm_tab)

Observed_AD_filter_norm_tab %>% filter(Subject != "Bryozor" & Subject != "JH21" & Shannon_AD > 0) %>% 
  ggplot(aes(x = Time, y = Shannon_AD, col = Subject, fill = Subject)) + 
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=1, alpha = 0.6) + 
  geom_smooth() + 
  labs(x= "\nTime (days)", y = "\nAD ", subtitle = "" )+
  theme_bw(base_size = 8)+
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none") + 
  scale_x_continuous(breaks = seq(0,120,by=10))
#scale_y_continuous(limits = c(0,2000))


Observed_AD_filter_norm_tab %>% filter(Subject != "Bryozor" & Subject != "JH21" & Obs_AD > 1) %>% 
  ggplot(aes(x = Time, y = Obs_AD, col = Subject, fill = Subject)) + 
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=1, alpha = 0.6)  +
  geom_smooth() + 
  labs(x= "\nTime (days)", y = "\nAD ", subtitle = "" )+
  theme_bw(base_size = 8)+
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none") + 
  scale_x_continuous(breaks = seq(0,120,by=10))
#scale_y_continuous(limits = c(0,2000))

###Beta-diversitet #############
# #Make phyloseq object to matrix - take out of Phyloseq and treat with normal vegan functions


#add real timepionts as numeric days to AD metadata frame

AD_filter_norm_meta_dat_joined <- left_join(AD_filter_norm_meta_dat, Observed_AD_filter_norm_tab, "Sample")
str(AD_filter_norm_meta_dat_joined)
AD_filter_succession_meta_dat_joined <- left_join(AD_succession_filter_norm_meta_dat, Observed_AD_filter_norm_tab, "Sample")
str(AD_filter_succession_meta_dat_joined)

AD_filter_norm_meta_dat <- AD_filter_norm_meta_dat_joined[,c(1:9,16)]
AD_filter_succession_meta_dat <- AD_filter_succession_meta_dat_joined[,c(1:9,16)]


## Add arbitary stages to all metadata
str(AD_filter_norm_meta_dat)
AD_filter_norm_meta_dat <- AD_filter_norm_meta_dat %>% mutate(stages = ifelse(Time>23, "late", "early"))

AD_filter_succession_meta_dat <- AD_filter_succession_meta_dat %>% mutate(stages = ifelse(Time>23, "late", "early"))
AD_filter_succession_meta_dat <- AD_filter_succession_meta_dat %>% mutate(stages = ifelse(Time==c(29, 44), "peak", stages))


#Beta-disper ----
## Calculate multivariate dispersions 
str(AD_filter_norm_dat_clean)
AD_dist_bray_bdisp_stages <- betadisper(vegdist(AD_filter_norm_dat_clean), "bray"), 
                                          factor(paste(AD_filter_norm_meta_dat$Subject.x,sep=" ")))

ADsuccession_dist_bray_bdisp_stages <- betadisper(vegdist(AD_filter_succession_norm_dat_clean, "bray"), 
                                        factor(paste(AD_filter_succession_meta_dat$stages,AD_filter_succession_meta_dat$Subject,sep=" ")))

plot(AD_dist_bray_bdisp_stages)
plot(ADsuccession_dist_bray_bdisp_stages)

adonis_all <- adonis(vegdist(AD_filter_norm_dat_clean, "bray") ~ Time * Subject.x * stages, AD_filter_norm_meta_dat)
adonis_all_dat <- as.data.frame(adonis_all$aov.tab)
adonis_all_dat

adonis_succession <- adonis(vegdist(AD_filter_succession_norm_dat_clean, "bray") ~ Time * stages, AD_filter_succession_meta_dat)
adonis_succession_dat <- as.data.frame(adonis_succession$aov.tab)
adonis_succession_dat


# NDMS AD for only bioelements-----
# 
# #Make phyloseq object to matrix - take out of Phyloseq and treat with normal vegan functions
# ps_AD_succession_filter_norm_dat <- as.data.frame((as(otu_table(ps_AD_succession_filter_norm), "matrix")))
# str(ps_AD_succession_filter_norm_dat)
# ps_AD_succession_filter_norm_dat <- round(ps_AD_succession_filter_norm_dat)
# #genus_norm_dat_t <- t(ps_AD_succession_filter_norm_dat)
# 
# 
# #Squar-root transform matrix
# sqrt_genus_norm_dat_t=sqrt(genus_norm_dat_t)
# 
# str(AD_succession_filter_norm_dat_clean)
# AD_succession_filter_norm_dat_clean_sqrt= sqrt(AD_succession_filter_norm_dat_clean)
str(AD_succession_filter_norm_dat_clean)
AD_filter_succession_norm_dist_bray = as.matrix((vegdist(AD_succession_filter_norm_dat_clean, "bray")))

NMDS_ord_bray_AD = metaMDS(AD_filter_succession_norm_dist_bray, k=2)
#NMDS_ord_bray_spec = metaMDS(sqrt_genus_norm_dat_t, distance = "bray", autotransform = FALSE, k=3)
str(NMDS_ord_bray_AD)

stressplot(NMDS_ord_bray_AD)

#build a data frame with NMDS coordinates and metadata
NMDS1_AD = NMDS_ord_bray_AD$points[,1]
NMDS2_AD = NMDS_ord_bray_AD$points[,2]
NMDS3_AD = NMDS_ord_bray_AD$points[,3]


# str(NMDS_ord_bray)
# NMDS_ord_bray$species
AD_filter_norm_meta_dat
P16S_filter_norm_meta_dat

NMDS_AD = data.frame(NMDS1 = NMDS1_AD, NMDS2 = NMDS2_AD, NMDS2 = NMDS3_AD, Time = AD_filter_norm_meta_dat_joined$Time, Chla = AD_filter_norm_meta_dat$Chla_suc, Stages = AD_filter_norm_meta_dat$Stages)


p_time_NMDS12_AD <- ggplot(NMDS_AD, aes(x=NMDS1, y=NMDS2, fill = Time)) +
  #stat_ellipse() +
  theme_bw(base_size = 8) +
  labs(subtitle = "AD") +
  geom_point(shape = 21,size = 2,colour = "black") +
  scale_fill_gradient(low = "#9111ab", high = "#e3b039") +
  geom_point(shape = 21, alpha = 0.25, size = 2) + theme_bw(base_size = 8) +
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 


p_bryo_NMDS12_AD <- ggplot(NMDS_AD, aes(x=NMDS1, y=NMDS2, fill = Stages)) +
  #stat_ellipse() +
  theme_bw(base_size = 8) +
  labs(title = "")  +
  geom_point(shape = 21,size = 2,colour = "black") +
  geom_point(shape = 21, alpha = 0.25, size = 2) + theme_bw(base_size = 8) + 
  scale_fill_manual(values = c("#4d44ab", "#e3b039")) + 
  theme(#axis.line = element_line(color='black'),
    plot.background = element_blank(),
    axis.title.y=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), legend.position = "none") 







