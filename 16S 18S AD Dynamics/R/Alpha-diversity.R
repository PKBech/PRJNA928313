
########################
## Richness dynamics ##
#######################

## Packages ##
library(tidyverse)
library("devtools")
library("Biostrings")
library("phyloseq")
library(vegan)
library(phyloseq)
library(venn)
library(plotly)
library(MicEco)
library(ggbeeswarm)
library(gridExtra)
library(ggpubr)
library(grid)
library(lme4)
library(emmeans)
library(biomeUtils) # For clustering
library(DECIPHER)
library(speedyseq)

#### Load data ####
#AD OBUs
load("16S 18S AD Preprocessing/Metataxonomic/3.AD/ps.AD.wTree_OBU99.RData") 
ps_AD_filtered_OBU99 

#16S
load("16S 18S AD Preprocessing/Metataxonomic/3.16S/ps_16S.reduced.wTree.20062022.RData") 
ps_16S_filtered <- ps.ASV.reduced.16S

#18S
load('16S 18S AD Preprocessing/Metataxonomic/3.18S/ps_18S.reduced.wTree.20062022.RData')
ps_18S_filtered <- ps_18S


#### Pre-processing, normalizing and cleaning-up tables ####

#Remove outliers
outliers <- c("Succession-T1-6","Succession-T1-8","Succession-T2-9","Succession-T2-7")
ps_AD_filtered_OBU99 <- subset_samples(ps_AD_filtered_OBU99, !sample_names(ps_AD_filtered_OBU99) %in% outliers)
ps_16S_filtered <- subset_samples(ps_16S_filtered, !sample_names(ps_16S_filtered) %in% outliers)
ps_18S_filtered <- subset_samples(ps_18S_filtered, !sample_names(ps_18S_filtered) %in% outliers)

#Scalling to total sum (100 000 reads)
#All sample types PS
ps_AD_filter_norm = transform_sample_counts(ps_AD_filtered_OBU99, function(x) 100000 * x/sum(x+0.1))
ps_16S_filter_norm = transform_sample_counts(ps_16S_filtered, function(x) 100000 * x/sum(x+0.1))
ps_18S_filter_norm = transform_sample_counts(ps_18S_filtered, function(x) 100000 * x/sum(x+0.1))

#Remove unwanted samples from AD and 16S
unwanted_AD <- c("Bryozor","JH21")
ps_AD_filter_norm <- subset_samples(ps_AD_filter_norm, !Subject %in% unwanted_AD)
unwanted_16S <- c("bryozoan","isolate")
ps_16S_filter_norm <- subset_samples(ps_16S_filter_norm, !element.type %in% unwanted_16S)


#### Richness/alpha-diversity ####

#16S#
# Observed richness
set.seed(41)
Observed_16S_filter_norm <- data.frame(t(estimateR((round(t(otu_table(ps_16S_filter_norm)))))))
#shannon diversity
shannon_16S_filter_norm <- diversity(round(t(otu_table(ps_16S_filter_norm))), index = "shannon")
#Bind together
Observed_16S_filter_norm_tab <- data.frame(cbind(Obs = Observed_16S_filter_norm$S.obs, shannon = shannon_16S_filter_norm, Sample=rownames(Observed_16S_filter_norm) ))
colnames(Observed_16S_filter_norm_tab)[1] <- "Obs_16S"
colnames(Observed_16S_filter_norm_tab)[2] <- "Shannon_16S"

#AD#
# Observed richness
set.seed(41)
Observed_AD_filter_norm <- data.frame(t(estimateR((round(t(otu_table(ps_AD_filter_norm))) ))))
# Shannon diversity
shannon_AD_filter_norm <- microbiome::diversity((round(t(otu_table(ps_AD_filter_norm)))), index = "shannon")
# Bind together
Observed_AD_filter_norm_tab <- data.frame(cbind(
  Obs = Observed_AD_filter_norm$S.obs, 
  shannon= shannon_AD_filter_norm, 
  
  Day=sample_data(ps_AD_filter_norm)$Day,
  Timepoint=sample_data(ps_AD_filter_norm)$Timepoint, 
  Subject=sample_data(ps_AD_filter_norm)$Subject), 
  Rep_tec=sample_data(ps_AD_filter_norm)$Rep, 
  Sample=rownames(Observed_AD_filter_norm) )

colnames(Observed_AD_filter_norm_tab)[1] <- "Obs_AD"
colnames(Observed_AD_filter_norm_tab)[2] <- "Shannon_AD"

#18S#
# Observed richness
set.seed(41)
Observed_18S_filter_norm <- data.frame(t(estimateR((round(t(otu_table(ps_18S_filter_norm)) )))))
# Shannon diversity
shannon_18S_filter_norm <- diversity(round(t(otu_table(ps_18S_filter_norm))), index = "shannon")
# Bind together
Observed_18S_filter_norm_tab <- data.frame(cbind(Obs = Observed_18S_filter_norm$S.obs, shannon = shannon_18S_filter_norm, Sample=rownames(Observed_18S_filter_norm) ))
colnames(Observed_18S_filter_norm_tab)[1] <- "Obs_18S"
colnames(Observed_18S_filter_norm_tab)[2] <- "Shannon_18S"

#Merge 16S, 18S to AD tabs - final diversity tabel#
#AD and 16S
Observed_ADand16S_filter_norm_tab <- left_join(Observed_AD_filter_norm_tab, Observed_16S_filter_norm_tab, "Sample")
#Polish columns
Observed_ADand16S_filter_norm_tab$Obs_AD <- as.numeric(Observed_ADand16S_filter_norm_tab$Obs_AD)
Observed_ADand16S_filter_norm_tab$Obs_16S <- as.numeric(Observed_ADand16S_filter_norm_tab$Obs_16S)
Observed_ADand16S_filter_norm_tab$Shannon_AD <- as.numeric(Observed_ADand16S_filter_norm_tab$Shannon_AD)
Observed_ADand16S_filter_norm_tab$Shannon_16S <- as.numeric(Observed_ADand16S_filter_norm_tab$Shannon_16S)
#Add 18S
Observed_ADand16Sand18S_filter_norm_tab <- left_join(Observed_ADand16S_filter_norm_tab, Observed_18S_filter_norm_tab, "Sample")
#Polish
Observed_ADand16Sand18S_filter_norm_tab$Obs_18S <- as.numeric(Observed_ADand16Sand18S_filter_norm_tab$Obs_18S)
Observed_ADand16Sand18S_filter_norm_tab$Shannon_18S <- as.numeric(Observed_ADand16Sand18S_filter_norm_tab$Shannon_18S)
Observed_ADand16Sand18S_filter_norm_tab$Day <- as.numeric(Observed_ADand16Sand18S_filter_norm_tab$Day)
# Change time to 'days'
Observed_ADand16Sand18S_filter_norm_tab$Time = 0
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==1] <- 3
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==2] <- 7
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==3] <- 10
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==4] <- 15
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==5] <- 23
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==6] <- 29
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==7] <- 44
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==8] <- 57
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==9] <- 71
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==10] <- 85
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==11] <- 99
Observed_ADand16Sand18S_filter_norm_tab$Time[Observed_ADand16Sand18S_filter_norm_tab$Day==12] <- 113

# Add factor: Community phase
Observed_ADand16Sand18S_filter_norm_tab$Phase = "Late"
Observed_ADand16Sand18S_filter_norm_tab$Phase[Observed_ADand16Sand18S_filter_norm_tab$Time %in% c(0,3,7,10,15,23)] <- "Early"
Observed_ADand16Sand18S_filter_norm_tab$Phase[Observed_ADand16Sand18S_filter_norm_tab$Time %in% c(29)] <- "Mid"
Observed_ADand16Sand18S_filter_norm_tab$Phase = factor(Observed_ADand16Sand18S_filter_norm_tab$Phase, levels= c("Early", "Mid", "Late"))

#transform to longformat
alphadiversity_tab = Observed_ADand16Sand18S_filter_norm_tab %>% 
  gather(key="Measure", value="Obs", "Shannon_AD", "Shannon_16S", "Shannon_18S", "Obs_18S", "Obs_16S", "Obs_AD") %>%
  select(-c(Day,Timepoint,Rep_tec,Sample)) %>%
  separate(col=Measure, into=c('Diversity', 'Measure'), sep='_')
  
#### Figures ####

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
Figure_1C = p_Obs.stats %>% filter(Subject != "Seawater" & Diversity != "Shannon") %>%
  ggplot(aes(x = Time, y = mean.obs)) + 
  geom_quasirandom(data = alphadiversity_tab %>% filter(Diversity != "Shannon" & Subject != "Seawater"), 
                   aes(y = Obs, col = Phase),size = 2, dodge.width=0, alpha = 0.6, stroke = NA) + 
  geom_errorbar(aes(ymax = mean.obs+sd.obs, ymin = mean.obs-sd.obs, col = Phase), width = 0, linewidth = 0.7, alpha = 0.2)+
  geom_line(size = 1, alpha = 0.6, col = "black") +
  labs(x= "\nTime (days)", y = "\n ASV/OBU richness \n")+
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

# ggsave(Figure_1C, file = "16S 18S AD Dynamics/Figures/Figure_1C.svg", height = 5.5, width =19, units = "cm")
# ggsave(Figure_1C, file = "16S 18S AD Dynamics/Figures/Figure_1C.png", height = 5.5, width =19, units = "cm")

## Tukey's test ####
unique(alphadiversity_tab$Diversity)

#Perfom tukeys test on observed richness for aov
aov_Obs_18S_succession <- alphadiversity_tab %>% filter(Diversity == "Obs" & 
                                                          Measure == "18S" &
                                                          Subject == "Succession" ) %>%
  aov(Obs ~ as.factor(Time), .)

TukeyHSD(aov_Obs_18S_succession)


# # Correlation pattern between amplicon richness ####

#AD vs 16S
corr_ADvs16S = Observed_ADand16Sand18S_filter_norm_tab %>% filter(!Subject %in% c("JH21", "Bryozor")) %>%
  ggplot(aes(y=Obs_AD, x=Obs_16S, color = Subject)) +
  geom_point(alpha=0.6, size = 1) +
  #geom_smooth(method = lm, se=FALSE)+
  stat_cor(method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name="rho", size = 2)+
  scale_size(range = c(0.1, 10)) +
  labs(x = "\n\n16S richness ", y = "\nAD richness ", subtitle = "") +
  theme_bw(base_size = 8) +
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        legend.position = "none")
#scale_y_continuous(limits = c(0,1000))

corr_ADvs16S

#AD vs 18S
corr_ADvs18S = Observed_ADand16Sand18S_filter_norm_tab %>% filter(!Subject %in% c("JH21", "Bryozor")) %>%
  ggplot(aes(y=Obs_AD, x=Obs_18S, color = Subject)) +
  geom_point(alpha=0.6, size = 1) +
  #geom_smooth(method = lm, se=FALSE)+
  stat_cor(method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name="rho", size = 2)+
  scale_size(range = c(0.1, 10)) +
  labs(x = "\n\n18S richness ", y = "\nAD richness ", subtitle ="")+
  theme_bw(base_size = 8) +
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none")
#scale_y_continuous(limits = c(0,1100)) +
#scale_x_continuous(limits = c(0,100))
corr_ADvs18S

#18S vs 16S
corr_16Svs18S <- Observed_ADand16Sand18S_filter_norm_tab %>% filter(!Subject %in% c("JH21", "Bryozor")) %>%
  ggplot(aes(y=Obs_18S, x=Obs_16S, color = Subject)) +
  geom_point(alpha=0.6, size = 1) +
  #geom_smooth(method = lm, se=FALSE)+
  stat_cor(method = "spearman", p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name="rho", size = 2)+
  scale_size(range = c(0.1, 10)) +
  labs(x = "\n\n16S richness ", y = "\n18S richness ", subtitle ="")+
  theme_bw(base_size = 8) +
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none")
#scale_y_continuous(limits = c(0,1200))

corr_16Svs18S

# Richness all environments ####

# 16S
p_Obs_16S.all <- Observed_ADand16Sand18S_filter_norm_tab %>%
  filter(Subject != "Bryozor" & Subject != "JH21" & Shannon_16S > 2) %>%
  ggplot(aes(x = Time, y = Obs_16S, col = Subject, fill = Subject)) +
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=1, alpha = 0.6) +
  geom_smooth() +
  labs(x= "\n\nTime (days)", y = "\n16S richness ", subtitle = "" )+
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
p_Obs_16S.all

p_Obs_AD.all <- Observed_ADand16Sand18S_filter_norm_tab %>%
  filter(Subject != "Bryozor" & Subject != "JH21" & Shannon_AD > 2) %>%
  ggplot(aes(x = Time, y = Obs_AD, col = Subject, fill = Subject)) +
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=1, alpha = 0.6)+
  geom_smooth() +
  labs(x= "\n\nTime (days)", y = "\nAD richness", subtitle = "")+
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme_bw(base_size = 8)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none")+
  scale_x_continuous(breaks = seq(0,120,by=10))
p_Obs_AD.all

p_Obs_18S.all = Observed_ADand16Sand18S_filter_norm_tab %>% filter(Subject != "Bryozor" & Subject != "JH21" & Shannon_18S > 2) %>%
  ggplot(aes(x = Time, y = Obs_18S, col = Subject, fill = Subject)) +
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=1, alpha = 0.6)+
  geom_smooth() +
  labs(x= "\n\nTime (days)", y = "\n18S richness", subtitle = "")+
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme_bw(base_size = 8) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none") +
  scale_x_continuous(breaks = seq(0,120,by=10))
#scale_y_continuous(limits = c(0,2000))
p_Obs_18S.all



# Shannon all environments ####
###############################
p_Shannon_16S = Observed_ADand16Sand18S_filter_norm_tab %>% filter(Subject != "Bryozor" & Subject != "JH21" & Shannon_16S > 0) %>%
  ggplot(aes(x = Time, y = Shannon_16S, col = Subject, fill = Subject)) +
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=1, alpha = 0.6) +
  geom_smooth() +
  labs(x= "\n\nTime (days)", y = " \n16S Shannon diversity\n", subtitle = "" )+
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
p_Shannon_16S

p_Shannon_AD = Observed_ADand16Sand18S_filter_norm_tab %>% filter(Subject != "Bryozor" & Subject != "JH21" & Shannon_AD > 0) %>%
  ggplot(aes(x = Time, y = Shannon_AD, col = Subject, fill = Subject)) +
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=1, alpha = 0.6)+
  geom_smooth() +
  labs(x= "\n\nTime (days)", y = " \nAD Shannon diversity\n", subtitle = "")+
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme_bw(base_size = 8)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none")+
  scale_x_continuous(breaks = seq(0,120,by=10))
#scale_y_continuous(limits = c(0,2000))
p_Shannon_AD

p_Shannon_18S = Observed_ADand16Sand18S_filter_norm_tab %>% filter(Subject != "Bryozor" & Subject != "JH21" & Shannon_18S > 0) %>%
  ggplot(aes(x = Time, y = Shannon_18S, col = Subject, fill = Subject)) +
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=1, alpha = 0.6)+
  geom_smooth() + ylim(0,6) +
  labs(x= "\n\nTime (days)", y = " \n18S Shannon diversity \n", subtitle = "")+
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme_bw(base_size = 8) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "none") +
  scale_x_continuous(breaks = seq(0,120,by=10))
#scale_y_continuous(limits = c(0,2000))
p_Shannon_18S

p_Shannon_18S_legend = Observed_ADand16Sand18S_filter_norm_tab %>% filter(Subject != "Bryozor" & Subject != "JH21" & Shannon_18S > 0) %>%
  ggplot(aes(x = Time, y = Shannon_18S, col = Subject, fill = Subject)) +
  geom_quasirandom(aes(col = Subject), dodge.width=.8, cex=1, alpha = 0.6)+
  geom_smooth() + ylim(0,6) +
  labs(x= "\n\nTime (days)", y = " \n18S Shannon diversity \n", subtitle = "")+
  scale_fill_manual(values = c("#f25a13", "#580a82"))+
  scale_color_manual(values = c("#f25a13", "#580a82"))+
  theme_bw(base_size = 8) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        legend.position = "bottom"
        ) +
  scale_x_continuous(breaks = seq(0,120,by=10))

legend.env = as_ggplot(ggpubr::get_legend(p_Shannon_18S_legend))

# Plot for supp figure
Figure_S2 = ggarrange(p_Shannon_16S,p_Shannon_18S, p_Shannon_AD,
                  p_Obs_16S.all, p_Obs_18S.all,p_Obs_AD.all,
                  corr_ADvs16S,corr_ADvs18S,corr_16Svs18S)


#ggsave(Figure_S2, file = "16S 18S AD Dynamics/Figures/Figure_S2.png", height = 16, width =18, units = "cm")

#ggsave(legend.env, file = "16S 18S AD Dynamics/Figures/legend_Figure_S2.png", height = 1, width =7, units = "cm")


#Linear model tested by emmeans ------
##Subset to succession and seawater
Observed_ADand16Sand18S_filter_norm_tab_filter <- Observed_ADand16Sand18S_filter_norm_tab %>% filter(Subject != "Bryozor" & Subject != "JH21" & Shannon_18S > 0 & Shannon_16S > 0 & Shannon_AD > 0)

str(Observed_ADand16Sand18S_filter_norm_tab_filter)
#Add replicate ID to dataframe
Observed_ADand16Sand18S_filter_norm_tab_filter$ID=gsub("(.*)(-T[0-9])","\\1",Observed_ADand16Sand18S_filter_norm_tab_filter$Sample)

#Make biorep as variable for lmer random effect
Observed_ADand16Sand18S_filter_norm_tab_filter = tidyr::separate(Observed_ADand16Sand18S_filter_norm_tab_filter,Sample,
                                                                 c("System", "Timepoint", "BioRep"), remove = FALSE)
Observed_ADand16Sand18S_filter_norm_tab_filter$Cage = paste(Observed_ADand16Sand18S_filter_norm_tab_filter$System,
                                                            Observed_ADand16Sand18S_filter_norm_tab_filter$BioRep, sep= "-")
Observed_ADand16Sand18S_filter_norm_tab_filter$Cage = as.factor(Observed_ADand16Sand18S_filter_norm_tab_filter$Cage)

#lmer
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_ObsAD.cage=lmer(Obs_AD~Time*Subject + (1|Cage), Observed_ADand16Sand18S_filter_norm_tab_filter)
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs16S.cage=lmer(Obs_16S~Time*Subject + (1|Cage), Observed_ADand16Sand18S_filter_norm_tab_filter)
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs18S.cage=lmer(Obs_18S~Time*Subject + (1|Cage), Observed_ADand16Sand18S_filter_norm_tab_filter)
car::Anova(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_ObsAD.cage)
car::Anova(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs18S.cage)
car::Anova(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs16S.cage)

#LMer for ADs 16S and 18S# Howver, this is only relevant for the succession not the seawater (no spatial effect from seawater samples)
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_ObsAD=lmer(Obs_AD~Time*Subject + (1|ID), Observed_ADand16Sand18S_filter_norm_tab_filter)
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs16S=lmer(Obs_16S~Time*Subject + (1|ID), Observed_ADand16Sand18S_filter_norm_tab_filter)
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs18S=lmer(Obs_18S~Time*Subject + (1|ID), Observed_ADand16Sand18S_filter_norm_tab_filter)

#LM only for observed and shannon
Observed_ADand16Sand18S_filter_norm_tab_filter_LM_ObsAD=lm(Obs_AD~Time*Subject , Observed_ADand16Sand18S_filter_norm_tab_filter)
Observed_ADand16Sand18S_filter_norm_tab_filter_LM_Obs16S=lm(Obs_16S~Time*Subject , Observed_ADand16Sand18S_filter_norm_tab_filter)
Observed_ADand16Sand18S_filter_norm_tab_filter_LM_Obs18S=lm(Obs_18S~Time*Subject , Observed_ADand16Sand18S_filter_norm_tab_filter)

# Observed_ADand16Sand18S_filter_norm_tab_filter_LM_ShaAD=lm(Shannon_AD~Time*Subject , Observed_ADand16Sand18S_filter_norm_tab_filter)
# Observed_ADand16Sand18S_filter_norm_tab_filter_LM_Sha16S=lm(Shannon_16S~Time*Subject , Observed_ADand16Sand18S_filter_norm_tab_filter)
# Observed_ADand16Sand18S_filter_norm_tab_filter_LM_Sha18S=lm(Shannon_18S~Time*Subject , Observed_ADand16Sand18S_filter_norm_tab_filter)
#

#significance for observed richness
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_ObsAD_EM=emmeans(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_ObsAD, ~Subject|Time)
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_ObsAD_EM_stat <- contrast(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_ObsAD_EM, interaction = "pairwise", adjust = "Bonferroni")
#Result is that over time AD OBS richness is significant greater in succession compared to seawater
# Time = 54.4:
#   Subject_pairwise      estimate   SE   df t.ratio p.value
# Seawater - Succession     -110 17.5 18.9  -6.295  <.0001
# 
# Degrees-of-freedom method: kenward-roger


Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs16S_EM=emmeans(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs16S, ~Subject|Time)
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs16S_EM_stat <- contrast(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs16S_EM, interaction = "pairwise", adjust = "Bonferroni")
#Result is that over time 16S OBS richness is significant greater in succession compared to seawater
# Time = 54.4:
#   Subject_pairwise      estimate  SE   df t.ratio p.value
# Seawater - Succession     -977 115 26.9  -8.486  <.0001
# 
# Degrees-of-freedom method: kenward-roger 


Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs18S_EM=emmeans(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs18S, ~Subject|Time)
Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs18S_EM_stat <- contrast(Observed_ADand16Sand18S_filter_norm_tab_filter_LMER_Obs18S_EM, interaction = "pairwise", adjust = "Bonferroni")
#Result is that over time 18S OBS richness is NOT significant diffrent between succession and seawater
# Time = 54.4:
#   Subject_pairwise      estimate   SE   df t.ratio p.value
# Seawater - Succession    -16.1 43.9 19.2  -0.367  0.7179
# 
# Degrees-of-freedom method: kenward-roger 



## Dunnett's test ----
#Make  Dunnest test with day 29 as "control"/ reference to all other points

#Make all combinations to individual groups
groups <- factor(paste(Observed_ADand16Sand18S_filter_norm_tab_filter$Subject,Observed_ADand16Sand18S_filter_norm_tab_filter$Time))



library(DescTools)
#ADs obs richness
DunnettTest_Day29_ObsAD <- DunnettTest(Observed_ADand16Sand18S_filter_norm_tab_filter$Obs_AD, groups, control = "Succession 29")
DunnettTest_Day29_ObsAD_dat <- as.data.frame(DunnettTest_Day29_ObsAD$`Succession 29`)
DunnettTest_Day29_ObsAD_dat$pval <- round(DunnettTest_Day29_ObsAD_dat$pval, digits = 3)
DunnettTest_Day29_ObsAD_dat$p.adj <- p.adjust(DunnettTest_Day29_ObsAD_dat$pval, method = "BH")
#Results ADs
# diff    lwr.ci     upr.ci  pval
# Seawater 10-Succession 29    -365.66667 -634.7135  -96.61981 0.002
# Seawater 113-Succession 29   -367.66667 -636.7135  -98.61981 0.001
# Seawater 15-Succession 29    -318.00000 -587.0469  -48.95314 0.010
# Seawater 23-Succession 29    -327.33333 -596.3802  -58.28648 0.007
# Seawater 29-Succession 29    -339.50000 -654.9854  -24.01460 0.026
# Seawater 44-Succession 29    -308.33333 -577.3802  -39.28648 0.013
# Seawater 57-Succession 29    -309.66667 -578.7135  -40.61981 0.013
# Seawater 71-Succession 29    -388.66667 -657.7135 -119.61981 0.001
# Seawater 85-Succession 29    -355.66667 -624.7135  -86.61981 0.002
# Seawater 99-Succession 29    -374.00000 -643.0469 -104.95314 0.001
# Succession 10-Succession 29  -224.66667 -414.9115  -34.42181 0.010
# Succession 113-Succession 29 -140.00000 -330.2449   50.24486 0.324
# Succession 15-Succession 29   -42.66667 -232.9115  147.57819 1.000
# Succession 23-Succession 29   -25.88889 -216.1337  164.35597 1.000
# Succession 44-Succession 29  -237.55556 -427.8004  -47.31070 0.005
# Succession 57-Succession 29  -250.85714 -454.2374  -47.47684 0.006
# Succession 71-Succession 29  -277.77778 -468.0226  -87.53292 0.000
# Succession 85-Succession 29  -318.55556 -508.8004 -128.31070 0.000
# Succession 99-Succession 29  -324.00000 -520.0999 -127.90009 0.000

#16S
DunnettTest_Day29_Obs16S <- DunnettTest(Observed_ADand16Sand18S_filter_norm_tab_filter$Obs_16S, groups, control = "Succession 29")
DunnettTest_Day29_Obs16S_dat <- as.data.frame(DunnettTest_Day29_Obs16S$`Succession 29`)
DunnettTest_Day29_Obs16S_dat$pval <- round(DunnettTest_Day29_Obs16S_dat$pval, digits = 3)

# diff      lwr.ci    upr.ci  pval
# Seawater 10-Succession 29    -1120.77778 -1728.76780 -512.7878 0.000
# Seawater 113-Succession 29   -1142.77778 -1750.76780 -534.7878 0.000
# Seawater 15-Succession 29    -1195.44444 -1803.43446 -587.4544 0.000
# Seawater 23-Succession 29    -1081.11111 -1689.10113 -473.1211 0.000
# Seawater 29-Succession 29    -1026.94444 -1739.87594 -314.0130 0.001
# Seawater 44-Succession 29    -1239.11111 -1847.10113 -631.1211 0.000
# Seawater 57-Succession 29     -872.44444 -1480.43446 -264.4544 0.001
# Seawater 71-Succession 29     -858.44444 -1466.43446 -250.4544 0.001
# Seawater 85-Succession 29    -1097.77778 -1705.76780 -489.7878 0.000
# Seawater 99-Succession 29    -1089.44444 -1697.43446 -481.4544 0.000
# Succession 10-Succession 29   -322.88889  -752.80276  107.0250 0.296
# Succession 113-Succession 29   450.11111    20.19724  880.0250 0.033
# Succession 15-Succession 29   -312.22222  -742.13609  117.6916 0.341
# Succession 23-Succession 29   -613.55556 -1043.46942 -183.6417 0.001
# Succession 44-Succession 29   -122.66667  -552.58053  307.2472 0.999
# Succession 57-Succession 29    232.12698  -227.47027  691.7242 0.832
# Succession 71-Succession 29    111.11111  -318.80276  541.0250 1.000
# Succession 85-Succession 29    -73.66667  -503.58053  356.2472 1.000
# Succession 99-Succession 29   -188.56944  -631.71451  254.5756 0.947

#And 18S
DunnettTest_Day23_Sha18S <- DunnettTest(Observed_ADand16Sand18S_filter_norm_tab_filter$Shannon_18S, groups, control = "Succession 23")
DunnettTest_Day23_Sha18S_dat <- as.data.frame(DunnettTest_Day23_Sha18S$`Succession 23`)
DunnettTest_Day23_Sha18S_dat$pval <- round(DunnettTest_Day23_Sha18S_dat$pval, digits = 3)
DunnettTest_Day23_Sha18S_dat

# diff    lwr.ci     upr.ci  pval
# Seawater 10-Succession 23    -278.111111 -617.6581   61.43587 0.194
# Seawater 113-Succession 23    -33.777778 -373.3248  305.76920 1.000
# Seawater 15-Succession 23    -302.777778 -642.3248   36.76920 0.117
# Seawater 23-Succession 23    -357.777778 -697.3248  -18.23080 0.032
# Seawater 29-Succession 23     -83.111111 -481.2652  315.04301 1.000
# Seawater 44-Succession 23    -268.444444 -607.9914   71.10253 0.233
# Seawater 57-Succession 23      37.888889 -301.6581  377.43587 1.000
# Seawater 71-Succession 23    -140.444444 -479.9914  199.10253 0.958
# Seawater 85-Succession 23       3.888889 -335.6581  343.43587 1.000
# Seawater 99-Succession 23    -117.444444 -456.9914  222.10253 0.992
# Succession 10-Succession 23   -61.333333 -301.4293  178.76264 1.000
# Succession 113-Succession 23   -8.888889 -248.9849  231.20708 1.000
# Succession 15-Succession 23    38.444444 -201.6515  278.54041 1.000
# Succession 29-Succession 23  -101.555556 -341.6515  138.54041 0.949
# Succession 44-Succession 23  -336.666667 -576.7626  -96.57070 0.001
# Succession 57-Succession 23  -123.682540 -380.3559  132.99085 0.873
# Succession 71-Succession 23  -317.444444 -557.5404  -77.34847 0.002
# Succession 85-Succession 23  -367.555556 -607.6515 -127.45959 0.000
# Succession 99-Succession 23  -145.861111 -393.3464  101.62415 0.647
