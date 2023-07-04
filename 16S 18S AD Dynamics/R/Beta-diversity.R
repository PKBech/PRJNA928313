
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

#Remove unwanted samples from AD and 16S
unwanted_AD <- c("Bryozor","JH21")
ps_AD_filtered_OBU99 <- subset_samples(ps_AD_filtered_OBU99, !Subject %in% unwanted_AD)
unwanted_16S <- c("bryozoan","isolate")
ps_16S_filtered <- subset_samples(ps_16S_filtered, !element.type %in% unwanted_16S)



#### beta-diversity Only for bioelements/succession ####
#Subset AD phyloseq to include only succession

ps_AD_succession_filter_norm <- subset_samples(ps_AD_filtered_OBU99, Subject == "Succession") %>%
  filter_taxa(., function (x) {sum(x > 0) > 0}, prune=TRUE) %>%
transform_sample_counts(., function(x) 100000 * x/sum(x+0.1))

ps_16S_succession_filter_norm <- subset_samples(ps_16S_filtered, element.type == "bioelement") %>%
  filter_taxa(., function (x) {sum(x > 0) > 0}, prune=TRUE) %>%
  transform_sample_counts(., function(x) 100000 * x/sum(x+0.1))

ps_18S_succession_filter_norm <- subset_samples(ps_18S_filtered, element.type == "bioelement") %>%
  filter_taxa(., function (x) {sum(x > 0) > 0}, prune=TRUE) %>%
  transform_sample_counts(., function(x) 100000 * x/sum(x+0.1))


#Make phyloseq object to matrix - take out of Phyloseq and treat with normal vegan functions
AD_filter_norm_dat <- as.data.frame(t(as(otu_table(ps_AD_succession_filter_norm), "matrix")))
AD_filter_norm_dat <- round(AD_filter_norm_dat)
str(AD_filter_norm_dat)

P16S_filter_norm_dat <- as.data.frame((as(otu_table(ps_16S_succession_filter_norm), "matrix")))
P16S_filter_norm_dat <- as.data.frame(round(t(P16S_filter_norm_dat)))
str(P16S_filter_norm_dat)

P18S_filter_norm_dat <- as.data.frame((as(otu_table(ps_18S_succession_filter_norm), "matrix")))
P18S_filter_norm_dat <- as.data.frame(round(t(P18S_filter_norm_dat)))
str(P16S_filter_norm_dat)


#remove empty rows
AD_filter_norm_dat_clean <- AD_filter_norm_dat[rowSums(AD_filter_norm_dat)!=0,]
P16S_filter_norm_dat_clean <- P16S_filter_norm_dat[rowSums(P16S_filter_norm_dat)!=0,]
P18S_filter_norm_dat_clean <- P18S_filter_norm_dat[rowSums(P18S_filter_norm_dat)!=0,]

#Match metadata table with rownames of X_filter_norm_dat_clean

AD_filter_norm_meta_dat <- data.frame(sample_data(ps_AD_succession_filter_norm))
#AD_filter_norm_meta_dat <- AD_succession_filter_norm_meta_dat
AD_filter_norm_meta_dat <- AD_filter_norm_meta_dat[rownames(AD_filter_norm_dat_clean), ]

#add real timepionts as numeric days to AD metadata frame
# Change time to 'days'
AD_filter_norm_meta_dat$Time = 0
AD_filter_norm_meta_dat$Time[AD_filter_norm_meta_dat$Day==1] <- 3
AD_filter_norm_meta_dat$Time[AD_filter_norm_meta_dat$Day==2] <- 7
AD_filter_norm_meta_dat$Time[AD_filter_norm_meta_dat$Day==3] <- 10
AD_filter_norm_meta_dat$Time[AD_filter_norm_meta_dat$Day==4] <- 15
AD_filter_norm_meta_dat$Time[AD_filter_norm_meta_dat$Day==5] <- 23
AD_filter_norm_meta_dat$Time[AD_filter_norm_meta_dat$Day==6] <- 29
AD_filter_norm_meta_dat$Time[AD_filter_norm_meta_dat$Day==7] <- 44
AD_filter_norm_meta_dat$Time[AD_filter_norm_meta_dat$Day==8] <- 57
AD_filter_norm_meta_dat$Time[AD_filter_norm_meta_dat$Day==9] <- 71
AD_filter_norm_meta_dat$Time[AD_filter_norm_meta_dat$Day==10] <- 85
AD_filter_norm_meta_dat$Time[AD_filter_norm_meta_dat$Day==11] <- 99
AD_filter_norm_meta_dat$Time[AD_filter_norm_meta_dat$Day==12] <- 113


## Add defined Phases to all metadata

phases <- as.data.frame(cbind(Sample = rownames(as.data.frame(sample_data(ps_16S_succession_filter_norm))),
                              phases = c(rep("Early", 4), rep("Late", 27),  rep("Early", 31), rep("Peak", 9),rep("Late", 27) ) ) )
phases
rownames(phases) <- phases$Sample
#Join to all dataframes
AD_filter_norm_meta_dat <- left_join(AD_filter_norm_meta_dat,phases, "Sample" )

P16S_filter_norm_meta_dat <- data.frame(sample_data(ps_16S_succession_filter_norm))
P16S_filter_norm_meta_dat <- P16S_filter_norm_meta_dat[rownames(P16S_filter_norm_dat_clean), ]
P16S_filter_norm_meta_dat <- P16S_filter_norm_meta_dat %>% mutate(Sample = rownames(P16S_filter_norm_meta_dat))
P16S_filter_norm_meta_dat <- left_join(P16S_filter_norm_meta_dat,phases, "Sample" )


P18S_filter_norm_meta_dat <- data.frame(sample_data(ps_18S_succession_filter_norm))
P18S_filter_norm_meta_dat <- P18S_filter_norm_meta_dat[rownames(P18S_filter_norm_dat_clean), ]
P18S_filter_norm_meta_dat <- P18S_filter_norm_meta_dat %>% mutate(Sample = rownames(P18S_filter_norm_meta_dat))
P18S_filter_norm_meta_dat <- left_join(P18S_filter_norm_meta_dat,phases, "Sample" )

#Sampling sizes:
AD_filter_norm_meta_dat %>% group_by(phases) %>% dplyr::summarise(n=n()) 
P16S_filter_norm_meta_dat %>% group_by(phases) %>% dplyr::summarise(n=n()) 
P18S_filter_norm_meta_dat %>% group_by(phases) %>% dplyr::summarise(n=n()) 


#Bray-curtis distance matrix
AD_dist_bray = as.matrix((vegdist((AD_filter_norm_dat_clean), "bray")))
str((AD_dist_bray))
P16S_dist_bray = as.matrix((vegdist((P16S_filter_norm_dat_clean), "bray")))
str(P16S_dist_bray)
P18S_dist_bray = as.matrix((vegdist((P18S_filter_norm_dat_clean), "bray")))
str(P18S_dist_bray)

#Beta-disper ----
## Calculate multivariate dispersions between phases
AD_dist_bray_bdisp_phases <- betadisper(vegdist((AD_filter_norm_dat_clean), "bray"),
                                        factor(paste(AD_filter_norm_meta_dat$phases,sep=" ")))
#PCoA
plot(AD_dist_bray_bdisp_phases)

#Test betadisper significans
## Permutation test for F
permutest(AD_dist_bray_bdisp_phases, pairwise = TRUE, permutations = 999)

## Tukey's Honest Significant Differences
(bdisp.HSD_AD <- TukeyHSD(AD_dist_bray_bdisp_phases))
plot(bdisp.HSD_AD)

## Calculate multivariate dispersions between phases

P16S_dist_bray_bdisp_phases <- betadisper(vegdist((P16S_filter_norm_dat_clean), "bray"),
                                          factor(paste(P16S_filter_norm_meta_dat$phases,sep=" "))) # use this file for PCoA plot (ggplot)


#Test betadisper significans
## Permutation test for F
permutest(P16S_dist_bray_bdisp_phases, pairwise = TRUE, permutations = 999)

## Tukey's Honest Significant Differences
(bdisp.HSD_16S <- TukeyHSD(P16S_dist_bray_bdisp_phases))
plot(bdisp.HSD_16S)




P18S_dist_bray_bdisp_phases <- betadisper(vegdist((P18S_filter_norm_dat_clean), "bray"),
                                          factor(paste(P18S_filter_norm_meta_dat$phases,sep=" ")))
plot(P18S_dist_bray_bdisp_phases)

#Test betadisper significans
## Permutation test for F
permutest(P18S_dist_bray_bdisp_phases, pairwise = TRUE, permutations = 999)

## Tukey's Honest Significant Differences
(bdisp.HSD_18S <- TukeyHSD(P18S_dist_bray_bdisp_phases))
plot(bdisp.HSD_18S)

#PERMANOVA (Adonis) ----
str(P16S_dist_bray)
str(P18S_dist_bray)
str(AD_dist_bray)



adonis_P16S_dist_bray <- adonis2(P16S_dist_bray ~ day * phases * biorep,
                                 P16S_filter_norm_meta_dat)
adonis_P16S_dist_bray

adonis_P18S_dist_bray <- adonis2(P18S_dist_bray ~ day + phases + biorep,
                                 P18S_filter_norm_meta_dat)
adonis_P18S_dist_bray

##Add temperature to AD metadata frame
colnames(AD_filter_norm_meta_dat)[3] <- "day"
AD_filter_norm_meta_dat <- left_join(AD_filter_norm_meta_dat, P18S_filter_norm_meta_dat, "Sample")

adonis_AD_dist_bray <-  adonis2(AD_dist_bray ~ day.x + phases.x * biorep,
                                AD_filter_norm_meta_dat)
adonis_AD_dist_bray


# write.table(adonis_Genus_dat, file = "adonis_Genus_dat.csv", sep = "\t", dec = ",")
# write.table(noPhao_adonis_Genus_dat, file = "noPhao_adonis_Genus_dat.csv", sep = "\t", dec = ",")
#



### Procruster analysis #####
data(varespec)

PCoA.ord_AD <- AD_dist_bray_bdisp_phases$vectors
PCoA.ord_16S <- P16S_dist_bray_bdisp_phases$vectors
PCoA.ord_18S <- P18S_dist_bray_bdisp_phases$vectors
#clean up rownames so it match the same sample names for all ordinations
PCoA.ord_16S <- PCoA.ord_16S[rownames(PCoA.ord_16S) %in% rownames(PCoA.ord_AD),]
PCoA.ord_AD <- PCoA.ord_AD[rownames(PCoA.ord_AD) %in% rownames(PCoA.ord_16S),]
PCoA.ord_18S <- PCoA.ord_18S[rownames(PCoA.ord_18S) %in% rownames(PCoA.ord_AD),]



cmdscale(AD_dist_bray)
str(PCoA.ord_AD)
cbind(PCoA.ord_16S, as.data.frame(AD_dist_bray_bdisp_phases$group[1:87]))

protest(X = PCoA.ord_16S, Y = PCoA.ord_AD, scores = "sites", permutations = 999)
#make data frame for procrustes residuals for 16S and AD procrustes
str(procrustes(X = PCoA.ord_16S, Y = PCoA.ord_AD, symmetric=TRUE), kind = 2)


residuals_dat_AD_16S <- as.data.frame(resid(procrustes(X = PCoA.ord_16S, Y = PCoA.ord_AD,symmetric=TRUE)))
residuals_dat_AD_16S <- cbind(residuals_dat_AD_16S, Sample = rownames(residuals_dat_AD_16S))
colnames(residuals_dat_AD_16S)[1] <- "Procrustes risidual"

residuals_dat_AD_16S <- left_join(residuals_dat_AD_16S, AD_filter_norm_meta_dat)
#Get the 25% (dashed), 50% (solid), and 75% (dashed) quantiles of the residuals.
quantiles_AD_16S <- as.data.frame(t(quantile(resid(procrustes(X = PCoA.ord_16S, Y = PCoA.ord_AD,symmetric=TRUE)))))
#Plot residuals
residuals_dat_AD_16S %>%
  ggplot(aes(x = as.factor(day.x), y = `Procrustes risidual`, col = phases.x, fill = phases.x)) +
  geom_quasirandom(aes(col = phases.x), dodge.width=.8, cex=1, alpha = 0.6) +
  #geom_smooth() +
  labs(x= "\nTime (days)", subtitle = "16S and AD ordinations association" )+
  theme_bw() +
  scale_fill_manual(values = c("#006400", "#f25a13", "#580a82")) +
  scale_color_manual(values = c("#006400", "#f25a13", "#580a82")) +
  geom_hline(yintercept=quantiles_AD_16S$`50%`) +
  geom_hline(yintercept=quantiles_AD_16S$`25%`, linetype="dashed") +
  geom_hline(yintercept=quantiles_AD_16S$`75%`, linetype="dashed") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        #axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        #legend.position = "none"
  )
protest(X = PCoA.ord_18S, Y = PCoA.ord_16S, scores = "sites", permutations = 999)
#make data frame for procrustes residuals for 16S and 18S procrustes
str(procrustes(X = PCoA.ord_18S, Y = PCoA.ord_16S,symmetric=TRUE), kind = 2)


residuals_dat_18S_16S <- as.data.frame(resid(procrustes(X = PCoA.ord_18S, Y = PCoA.ord_16S,symmetric=TRUE)))
residuals_dat_18S_16S <- cbind(residuals_dat_18S_16S, Sample = rownames(residuals_dat_18S_16S))
colnames(residuals_dat_18S_16S)[1] <- "Procrustes risidual"

residuals_dat_18S_16S <- left_join(residuals_dat_18S_16S, AD_filter_norm_meta_dat)
#Get the 25% (dashed), 50% (solid), and 75% (dashed) quantiles of the residuals.
quantiles_18S_16S <- as.data.frame(t(quantile(resid(procrustes(X = PCoA.ord_18S, Y = PCoA.ord_16S,symmetric=TRUE)))))
#Plot residuals
residuals_dat_18S_16S %>%
  ggplot(aes(x = as.factor(day.y), y = `Procrustes risidual`, col = phases.x, fill = phases.x)) +
  geom_quasirandom(aes(col = phases.x), dodge.width=.8, cex=1, alpha = 0.6) +
  #geom_smooth() +
  labs(x= "\nTime (days)", subtitle = "16S and 18S ordinations association" )+
  theme_bw() +
  scale_fill_manual(values = c("#006400", "#f25a13", "#580a82")) +
  scale_color_manual(values = c("#006400", "#f25a13", "#580a82")) +
  geom_hline(yintercept=quantiles_18S_16S$`50%`) +
  geom_hline(yintercept=quantiles_18S_16S$`25%`, linetype="dashed") +
  geom_hline(yintercept=quantiles_18S_16S$`75%`, linetype="dashed") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        #axis.title.x = element_text(margin = margin(t = -10)),
        axis.text = element_text(face = "bold"),
        #legend.position = "none"
  )


protest(X = PCoA.ord_18S, Y = PCoA.ord_AD, scores = "sites", permutations = 999)

# PCoA plot by Phases#####

# 16S (use data from betadisper)
centroids.16 = data.frame(P16S_dist_bray_bdisp_phases$centroids[,1:2])
centroids.16$Phase = as.factor(rownames(centroids.16))
points.16 = data.frame(P16S_dist_bray_bdisp_phases$vectors[,1:2], Phase = P16S_dist_bray_bdisp_phases$group)
centroids.16$Phase = factor(centroids.16$Phase, levels= c("Early", "Peak", "Late"))
points.16$Phase = factor(points.16$Phase, levels= c("Early", "Peak", "Late"))
#Centroids mean of first and second component
P16S_dist_bray_bdisp_phases$centroids[,1:2]


#Components
PCoA1 <- P16S_dist_bray_bdisp_phases$eig[1]/sum(P16S_dist_bray_bdisp_phases$eig)*100
PCoA1
PCoA2 <- P16S_dist_bray_bdisp_phases$eig[2]/sum(P16S_dist_bray_bdisp_phases$eig)*100
PCoA2
PCoA1+PCoA2



main_col <- c("#240785", "#e0b62b", "#f21395")

PCOA.16S = ggplot(centroids.16, aes(PCoA1, PCoA2, color = Phase)) +
  theme_bw(base_size = 12)+
  geom_point(data = points.16, aes(PCoA1, PCoA2, color = Phase), size = 1.5, alpha = 0.5) +
  geom_segment(aes(x = 0.3811575, y = -0.08426253, xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Phase == "Early",], size =1, alpha = 0.3, colour = "#240785") +
  geom_segment(aes(x =  0.1717421, y = 0.25524678, xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Phase == "Peak",], size =1, alpha = 0.3, colour = "#e0b62b") +
  geom_segment(aes(x =  -0.2435972, y = -0.00192419, xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Phase == "Late",], size =1, alpha = 0.3, colour = "#f21395" ) +
  geom_point(size = 4)+
  #facet_grid(Phase~.)+
  stat_ellipse(data = points.16, aes(group=Phase, color = Phase), size = 1, alpha = 0.3)+
  labs(x = "PCoA1 [36.1%]", y = "\nPCoA2 [12.0%]", title = "Prokaryotic")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = c("#240785", "#e0b62b", "#f21395")) +
  scale_color_manual(values = c("#240785", "#e0b62b", "#f21395"))
PCOA.16S

legend.pcoa = as_ggplot(ggpubr::get_legend(PCOA.16S))
#ggsave(legend.pcoa, file = "legend.pcoa.pdf", width = 1, height = 1)

# AD (use data from betadisper)
centroids.AD = data.frame(AD_dist_bray_bdisp_phases$centroids[,1:2])
points.AD = data.frame(AD_dist_bray_bdisp_phases$vectors[,1:2], Phase = AD_dist_bray_bdisp_phases$group)
#Components
PCoA1 <- AD_dist_bray_bdisp_phases$eig[1]/sum(AD_dist_bray_bdisp_phases$eig)*100
PCoA1
PCoA2 <- AD_dist_bray_bdisp_phases$eig[2]/sum(AD_dist_bray_bdisp_phases$eig)*100
PCoA2
PCoA1+PCoA2
#Centroids mean of first and second component
AD_dist_bray_bdisp_phases$centroids[,1:2]



PCOA.AD = ggplot(centroids.AD, aes(PCoA1, PCoA2, color = rownames(centroids.AD))) +
  theme_bw(base_size = 12)+
  geom_point(data = points.AD, aes(PCoA1, PCoA2, color = Phase), size = 1.5, alpha = 0.5) +
  geom_segment(aes(x = 0.34117984, y = -0.004116847, xend = PCoA1, yend = PCoA2),
               data = points.AD[points.AD$Phase == "Early",], size =1, alpha = 0.3, colour = "#240785") +
  geom_segment(aes(x =  -0.17626409, y = 0.002116555, xend = PCoA1, yend = PCoA2),
               data = points.AD[points.AD$Phase == "Late",], size =1, alpha = 0.3, colour = "#f21395" ) +
  geom_segment(aes(x = 0.02374626, y = 0.223536004, xend = PCoA1, yend = PCoA2),
               data = points.AD[points.AD$Phase == "Peak",], size =1, alpha = 0.3, colour = "#e0b62b") +
  geom_point(size = 4)+
  stat_ellipse(data = points.AD, aes(group=Phase, color = Phase), size = 1, alpha = 0.3)+
  labs(x = "PCoA1 [18.4%]", y = "\nPCoA2 [10.0%]", title = "AD")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = c("#240785", "#f21395", "#e0b62b")) +
  scale_color_manual(values = c("#240785", "#f21395", "#e0b62b"))

# 18S (use data from betadisper)
centroids.18S = data.frame(P18S_dist_bray_bdisp_phases$centroids[,1:2])
points.18S = data.frame(P18S_dist_bray_bdisp_phases$vectors[,1:2], Phase = P18S_dist_bray_bdisp_phases$group)
#Components
PCoA1 <- P18S_dist_bray_bdisp_phases$eig[1]/sum(P18S_dist_bray_bdisp_phases$eig)*100
PCoA1
PCoA2 <- P18S_dist_bray_bdisp_phases$eig[2]/sum(P18S_dist_bray_bdisp_phases$eig)*100
PCoA2
PCoA1+PCoA2
#Centroids mean of first and second component
P18S_dist_bray_bdisp_phases$centroids[,1:2]

PCOA.18S = ggplot(centroids.18S, aes(PCoA1, PCoA2, color = rownames(centroids.18S))) +
  theme_bw(base_size = 12)+
  geom_point(data = points.18S, aes(PCoA1, PCoA2, color = Phase), size = 1.5, alpha = 0.5) +
  geom_segment(aes(x = -0.3665717, y = -0.083343689, xend = PCoA1, yend = PCoA2),
               data = points.18S[points.18S$Phase == "Early",], size =1, alpha = 0.3, colour = "#240785") +
  geom_segment(aes(x =  0.2549065, y = -0.005344822, xend = PCoA1, yend = PCoA2),
               data = points.18S[points.18S$Phase == "Late",], size =1, alpha = 0.3, colour = "#f21395" ) +
  geom_segment(aes(x = -0.1308241, y = 0.319936332, xend = PCoA1, yend = PCoA2),
               data = points.18S[points.18S$Phase == "Peak",], size =1, alpha = 0.3, colour = "#e0b62b") +
  geom_point(size = 4)+
  stat_ellipse(data = points.18S, aes(group=Phase, color = Phase), size = 1, alpha = 0.3)+
  labs(x = "PCoA1 [27.2%]", y = "\nPCoA2 [12.9%]", title = "Eukaryotic")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = c("#240785", "#f21395", "#e0b62b")) +
  scale_color_manual(values = c("#240785", "#f21395", "#e0b62b"))

PCoA_Figure_1E = grid.arrange(PCOA.16S, PCOA.18S, PCOA.AD, nrow = 1)

ggsave(PCoA_Figure_1E, file = '16S 18S AD Dynamics/Figures/Figure_1E.png', width = 10.5, height = 3)


# PCoA plot by Time (categorical) #####
P16S_dist_bray_bdisp_time <- betadisper(vegdist((P16S_filter_norm_dat_clean), "bray"),
                                          factor(paste(P16S_filter_norm_meta_dat$day,sep=" "))) # use this file for PCoA plot (ggplot)
# 16S (use data from betadisper)
centroids.16 = data.frame(P16S_dist_bray_bdisp_time$centroids[,1:2])
centroids.16$Time = as.factor(rownames(centroids.16))
points.16 = data.frame(P16S_dist_bray_bdisp_time$vectors[,1:2], Time = P16S_dist_bray_bdisp_time$group)
unique(centroids.16$Time)
centroids.16$Time = factor(centroids.16$Time, levels= c("3","7", "10", "15","23", 
                                                        "29","44", "57", "71",
                                                        "85", "99", "113"))

points.16$Time = factor(points.16$Time, levels= c("3", "7","10", "15","23", 
                                                    "29","44", "57", "71",
                                                    "85", "99", "113"))
#Centroids mean of first and second component
P16S_dist_bray_bdisp_time$centroids[,1:2]

#Components
PCoA1 <- P16S_dist_bray_bdisp_time$eig[1]/sum(P16S_dist_bray_bdisp_time$eig)*100
PCoA1
PCoA2 <- P16S_dist_bray_bdisp_time$eig[2]/sum(P16S_dist_bray_bdisp_time$eig)*100
PCoA2
PCoA1+PCoA2


# time_colors = c("#f8df25", "#f8df25", "#fdc627","#f99a3e", "#e97257", "#de6164",
#                          "#c43e7f", "#b42e8d", "#a21d9a", 
#                          "#6100a7","#4903a0",
#                          "#2f0596","#2f0596", "#0d0887")
#                          
time_colors = c("#f8df25", "#f8df25", "#fdc627","#f99a3e", "#e97257", "#de6164",
                         "#b42e8d", "#a21d9a", 
                         "#6100a7",
                         "#2f0596","#2f0596", "#0d0887")
                         

PCOA.16S <- ggplot(centroids.16, aes(PCoA1, PCoA2, color = Time)) +
  theme_bw(base_size = 12)+
  geom_point(data = points.16, aes(PCoA1, PCoA2, color = Time), size = 1.5, alpha = 0.5) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[6,1], y = P16S_dist_bray_bdisp_time$centroids[6,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "3",], size =1, alpha = 0.3, colour = time_colors[1]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[9,1], y = P16S_dist_bray_bdisp_time$centroids[9,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "7",], size =1, alpha = 0.3, colour = time_colors[2]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[1,1], y = P16S_dist_bray_bdisp_time$centroids[1,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "10",], size =1, alpha = 0.3, colour = time_colors[3]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[3,1], y = P16S_dist_bray_bdisp_time$centroids[3,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "15",], size =1, alpha = 0.3, colour = time_colors[4]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[4,1], y = P16S_dist_bray_bdisp_time$centroids[4,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "23",], size =1, alpha = 0.3, colour = time_colors[5]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[5,1], y = P16S_dist_bray_bdisp_time$centroids[5,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "29",], size =1, alpha = 0.3, colour = time_colors[6]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[7,1], y = P16S_dist_bray_bdisp_time$centroids[7,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "44",], size =1, alpha = 0.3, colour = time_colors[7]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[8,1], y = P16S_dist_bray_bdisp_time$centroids[8,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "57",], size =1, alpha = 0.3, colour = time_colors[8]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[10,1], y = P16S_dist_bray_bdisp_time$centroids[10,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "71",], size =1, alpha = 0.3, colour = time_colors[9]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[11,1], y = P16S_dist_bray_bdisp_time$centroids[11,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "85",], size =1, alpha = 0.3, colour = time_colors[10]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[12,1], y = P16S_dist_bray_bdisp_time$centroids[12,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "99",], size =1, alpha = 0.3, colour = time_colors[11]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[2,1], y = P16S_dist_bray_bdisp_time$centroids[2,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "113",], size =1, alpha = 0.3, colour = time_colors[12]) +
  geom_point(size = 4) +
  #facet_grid(Time~.)+
  stat_ellipse(data = points.16, aes(group=Time, color = Time), size = 1, alpha = 0.3)+
  labs(x = "PCoA1 [18.4%]", y = "\nPCoA2 [11.9%]", title = "Prokaryotic")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = time_colors) +
  scale_color_manual(values = time_colors)
PCOA.16S


PCOA.16S_legend <- ggplot(centroids.16, aes(PCoA1, PCoA2, color = Time)) +
  theme_bw(base_size = 12)+
  geom_point(data = points.16, aes(PCoA1, PCoA2, color = Time), size = 1.5, alpha = 0.5) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[6,1], y = P16S_dist_bray_bdisp_time$centroids[6,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "3",], size =1, alpha = 0.3, colour = time_colors[1]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[9,1], y = P16S_dist_bray_bdisp_time$centroids[9,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "7",], size =1, alpha = 0.3, colour = time_colors[2]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[1,1], y = P16S_dist_bray_bdisp_time$centroids[1,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "10",], size =1, alpha = 0.3, colour = time_colors[3]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[3,1], y = P16S_dist_bray_bdisp_time$centroids[3,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "15",], size =1, alpha = 0.3, colour = time_colors[4]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[4,1], y = P16S_dist_bray_bdisp_time$centroids[4,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "23",], size =1, alpha = 0.3, colour = time_colors[5]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[5,1], y = P16S_dist_bray_bdisp_time$centroids[5,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "29",], size =1, alpha = 0.3, colour = time_colors[6]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[7,1], y = P16S_dist_bray_bdisp_time$centroids[7,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "44",], size =1, alpha = 0.3, colour = time_colors[7]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[8,1], y = P16S_dist_bray_bdisp_time$centroids[8,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "57",], size =1, alpha = 0.3, colour = time_colors[8]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[10,1], y = P16S_dist_bray_bdisp_time$centroids[10,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "71",], size =1, alpha = 0.3, colour = time_colors[9]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[11,1], y = P16S_dist_bray_bdisp_time$centroids[11,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "85",], size =1, alpha = 0.3, colour = time_colors[10]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[12,1], y = P16S_dist_bray_bdisp_time$centroids[12,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "99",], size =1, alpha = 0.3, colour = time_colors[11]) +
  geom_segment(aes(x = P16S_dist_bray_bdisp_time$centroids[2,1], y = P16S_dist_bray_bdisp_time$centroids[2,2], xend = PCoA1, yend = PCoA2),
               data = points.16[points.16$Time == "113",], size =1, alpha = 0.3, colour = time_colors[12]) +
  geom_point(size = 4) +
  #facet_grid(Time~.)+
  stat_ellipse(data = points.16, aes(group=Time, color = Time), size = 1, alpha = 0.3)+
  labs(x = "PCoA1 [18.4%]", y = "\nPCoA2 [11.9%]", title = "Prokaryotic",color='Day' ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        #legend.position = "none", 
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = time_colors) +
  scale_color_manual(values = time_colors)

#ggsave(legend.pcoa, file = "legend.pcoa.pdf", width = 1, height = 1)

# AD (use data from betadisper)
AD_dist_bray_bdisp_time <- betadisper(vegdist((AD_filter_norm_dat_clean), "bray"),
                                      factor(paste(AD_filter_norm_meta_dat$Time,sep=" ")))

# 16S (use data from betadisper)
centroids.AD = data.frame(AD_dist_bray_bdisp_time$centroids[,1:2])
centroids.AD$Time = as.factor(rownames(centroids.AD))
points.AD = data.frame(AD_dist_bray_bdisp_time$vectors[,1:2], Time = AD_dist_bray_bdisp_time$group)
unique(centroids.AD$Time)
centroids.AD$Time = factor(centroids.AD$Time, levels= c( "10", "15","23", 
                                                        "29","44", "57", "71",
                                                        "85", "99", "113"))

points.AD$Time = factor(points.AD$Time, levels= c("10", "15","23", 
                                                  "29","44", "57", "71",
                                                  "85", "99", "113"))
#Centroids mean of first and second component
AD_dist_bray_bdisp_time$centroids[,1:2]

#Components
PCoA1 <- AD_dist_bray_bdisp_time$eig[1]/sum(AD_dist_bray_bdisp_time$eig)*100
PCoA1
PCoA2 <- AD_dist_bray_bdisp_time$eig[2]/sum(AD_dist_bray_bdisp_time$eig)*100
PCoA2
PCoA1+PCoA2

time_colors_10 = c("#fdc627","#f99a3e", "#e97257", "#de6164",
                         "#b42e8d", "#a21d9a", 
                         "#6100a7",
                         "#2f0596","#2f0596", "#0d0887")
                         

PCOA.AD <- ggplot(centroids.AD, aes(PCoA1, PCoA2, color = Time)) +
  theme_bw(base_size = 12)+
  geom_point(data = points.AD, aes(PCoA1, PCoA2, color = Time), size = 1.5, alpha = 0.5) +
  # geom_segment(aes(x = AD_dist_bray_bdisp_time$centroids[6,1], y = AD_dist_bray_bdisp_time$centroids[6,2], xend = PCoA1, yend = PCoA2),
  #              data = points.AD[points.AD$Time == "3",], size =1, alpha = 0.3, colour = time_colors_10[1]) +
  # geom_segment(aes(x = AD_dist_bray_bdisp_time$centroids[9,1], y = AD_dist_bray_bdisp_time$centroids[9,2], xend = PCoA1, yend = PCoA2),
  #              data = points.AD[points.AD$Time == "7",], size =1, alpha = 0.3, colour = time_colors_10[2]) +
  geom_segment(aes(x = AD_dist_bray_bdisp_time$centroids[1,1], y = AD_dist_bray_bdisp_time$centroids[1,2], xend = PCoA1, yend = PCoA2),
               data = points.AD[points.AD$Time == "10",], size =1, alpha = 0.3, colour = time_colors_10[1]) +
  geom_segment(aes(x = AD_dist_bray_bdisp_time$centroids[3,1], y = AD_dist_bray_bdisp_time$centroids[3,2], xend = PCoA1, yend = PCoA2),
               data = points.AD[points.AD$Time == "15",], size =1, alpha = 0.3, colour = time_colors_10[2]) +
  geom_segment(aes(x = AD_dist_bray_bdisp_time$centroids[4,1], y = AD_dist_bray_bdisp_time$centroids[4,2], xend = PCoA1, yend = PCoA2),
               data = points.AD[points.AD$Time == "23",], size =1, alpha = 0.3, colour = time_colors_10[3]) +
  geom_segment(aes(x = AD_dist_bray_bdisp_time$centroids[5,1], y = AD_dist_bray_bdisp_time$centroids[5,2], xend = PCoA1, yend = PCoA2),
               data = points.AD[points.AD$Time == "29",], size =1, alpha = 0.3, colour = time_colors_10[4]) +
  geom_segment(aes(x = AD_dist_bray_bdisp_time$centroids[6,1], y = AD_dist_bray_bdisp_time$centroids[6,2], xend = PCoA1, yend = PCoA2),
               data = points.AD[points.AD$Time == "44",], size =1, alpha = 0.3, colour = time_colors_10[5]) +
  geom_segment(aes(x = AD_dist_bray_bdisp_time$centroids[7,1], y = AD_dist_bray_bdisp_time$centroids[7,2], xend = PCoA1, yend = PCoA2),
               data = points.AD[points.AD$Time == "57",], size =1, alpha = 0.3, colour = time_colors_10[6]) +
  geom_segment(aes(x = AD_dist_bray_bdisp_time$centroids[8,1], y = AD_dist_bray_bdisp_time$centroids[8,2], xend = PCoA1, yend = PCoA2),
               data = points.AD[points.AD$Time == "71",], size =1, alpha = 0.3, colour = time_colors_10[7]) +
  geom_segment(aes(x = AD_dist_bray_bdisp_time$centroids[9,1], y = AD_dist_bray_bdisp_time$centroids[9,2], xend = PCoA1, yend = PCoA2),
               data = points.AD[points.AD$Time == "85",], size =1, alpha = 0.3, colour = time_colors_10[8]) +
  geom_segment(aes(x = AD_dist_bray_bdisp_time$centroids[10,1], y = AD_dist_bray_bdisp_time$centroids[10,2], xend = PCoA1, yend = PCoA2),
               data = points.AD[points.AD$Time == "99",], size =1, alpha = 0.3, colour = time_colors_10[9]) +
  geom_segment(aes(x = AD_dist_bray_bdisp_time$centroids[2,1], y = AD_dist_bray_bdisp_time$centroids[2,2], xend = PCoA1, yend = PCoA2),
               data = points.AD[points.AD$Time == "113",], size =1, alpha = 0.3, colour = time_colors_10[10]) +
  geom_point(size = 4) +
  #facet_grid(Time~.)+
  stat_ellipse(data = points.AD, aes(group=Time, color = Time), size = 1, alpha = 0.3)+
  labs(x = "PCoA1 [18.4%]", y = "\nPCoA2 [9.9%]", title = "AD")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = time_colors_10) +
  scale_color_manual(values = time_colors_10)
PCOA.AD

# 18S (use data from betadisper)
P18S_dist_bray_bdisp_time <- betadisper(vegdist((P18S_filter_norm_dat_clean), "bray"),
                                        factor(paste(P18S_filter_norm_meta_dat$day,sep=" "))) # use this file for PCoA plot (ggplot)
# 16S (use data from betadisper)
centroids.18 = data.frame(P18S_dist_bray_bdisp_time$centroids[,1:2])
centroids.18$Time = as.factor(rownames(centroids.18))
points.18 = data.frame(P18S_dist_bray_bdisp_time$vectors[,1:2], Time = P18S_dist_bray_bdisp_time$group)
unique(centroids.18$Time)
centroids.18$Time = factor(centroids.18$Time, levels= c("3","7", "10", "15","23", 
                                                        "29","44", "57", "71",
                                                        "85", "99", "113"))

points.18$Time = factor(points.18$Time, levels= c("3", "7","10", "15","23", 
                                                  "29","44", "57", "71",
                                                  "85", "99", "113"))
#Centroids mean of first and second component
P18S_dist_bray_bdisp_time$centroids[,1:2]

#Components
PCoA1 <- P18S_dist_bray_bdisp_time$eig[1]/sum(P18S_dist_bray_bdisp_time$eig)*100
PCoA1
PCoA2 <- P18S_dist_bray_bdisp_time$eig[2]/sum(P18S_dist_bray_bdisp_time$eig)*100
PCoA2
PCoA1+PCoA2

PCOA.18S <- ggplot(centroids.18, aes(PCoA1, PCoA2, color = Time)) +
  theme_bw(base_size = 12)+
  geom_point(data = points.18, aes(PCoA1, PCoA2, color = Time), size = 1.5, alpha = 0.5) +
  geom_segment(aes(x = P18S_dist_bray_bdisp_time$centroids[6,1], y = P18S_dist_bray_bdisp_time$centroids[6,2], xend = PCoA1, yend = PCoA2),
               data = points.18[points.18$Time == "3",], size =1, alpha = 0.3, colour = time_colors[1]) +
  geom_segment(aes(x = P18S_dist_bray_bdisp_time$centroids[9,1], y = P18S_dist_bray_bdisp_time$centroids[9,2], xend = PCoA1, yend = PCoA2),
               data = points.18[points.18$Time == "7",], size =1, alpha = 0.3, colour = time_colors[2]) +
  geom_segment(aes(x = P18S_dist_bray_bdisp_time$centroids[1,1], y = P18S_dist_bray_bdisp_time$centroids[1,2], xend = PCoA1, yend = PCoA2),
               data = points.18[points.18$Time == "10",], size =1, alpha = 0.3, colour = time_colors[3]) +
  geom_segment(aes(x = P18S_dist_bray_bdisp_time$centroids[3,1], y = P18S_dist_bray_bdisp_time$centroids[3,2], xend = PCoA1, yend = PCoA2),
               data = points.18[points.18$Time == "15",], size =1, alpha = 0.3, colour = time_colors[4]) +
  geom_segment(aes(x = P18S_dist_bray_bdisp_time$centroids[4,1], y = P18S_dist_bray_bdisp_time$centroids[4,2], xend = PCoA1, yend = PCoA2),
               data = points.18[points.18$Time == "23",], size =1, alpha = 0.3, colour = time_colors[5]) +
  geom_segment(aes(x = P18S_dist_bray_bdisp_time$centroids[5,1], y = P18S_dist_bray_bdisp_time$centroids[5,2], xend = PCoA1, yend = PCoA2),
               data = points.18[points.18$Time == "29",], size =1, alpha = 0.3, colour = time_colors[6]) +
  geom_segment(aes(x = P18S_dist_bray_bdisp_time$centroids[7,1], y = P18S_dist_bray_bdisp_time$centroids[7,2], xend = PCoA1, yend = PCoA2),
               data = points.18[points.18$Time == "44",], size =1, alpha = 0.3, colour = time_colors[7]) +
  geom_segment(aes(x = P18S_dist_bray_bdisp_time$centroids[8,1], y = P18S_dist_bray_bdisp_time$centroids[8,2], xend = PCoA1, yend = PCoA2),
               data = points.18[points.18$Time == "57",], size =1, alpha = 0.3, colour = time_colors[8]) +
  geom_segment(aes(x = P18S_dist_bray_bdisp_time$centroids[10,1], y = P18S_dist_bray_bdisp_time$centroids[10,2], xend = PCoA1, yend = PCoA2),
               data = points.18[points.18$Time == "71",], size =1, alpha = 0.3, colour = time_colors[9]) +
  geom_segment(aes(x = P18S_dist_bray_bdisp_time$centroids[11,1], y = P18S_dist_bray_bdisp_time$centroids[11,2], xend = PCoA1, yend = PCoA2),
               data = points.18[points.18$Time == "85",], size =1, alpha = 0.3, colour = time_colors[10]) +
  geom_segment(aes(x = P18S_dist_bray_bdisp_time$centroids[12,1], y = P18S_dist_bray_bdisp_time$centroids[12,2], xend = PCoA1, yend = PCoA2),
               data = points.18[points.18$Time == "99",], size =1, alpha = 0.3, colour = time_colors[11]) +
  geom_segment(aes(x = P18S_dist_bray_bdisp_time$centroids[2,1], y = P18S_dist_bray_bdisp_time$centroids[2,2], xend = PCoA1, yend = PCoA2),
               data = points.18[points.18$Time == "113",], size =1, alpha = 0.3, colour = time_colors[12]) +
  geom_point(size = 4) +
  #facet_grid(Time~.)+
  stat_ellipse(data = points.18, aes(group=Time, color = Time), size = 1, alpha = 0.3)+
  labs(x = "PCoA1 [27.2%]", y = "\nPCoA2 [12.9%]", title = "Eukaryotes")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "none",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = time_colors) +
  scale_color_manual(values = time_colors)
PCOA.18S



PCoA_Figure_S3 = grid.arrange(PCOA.16S, PCOA.18S, PCOA.AD, nrow = 1)

ggsave(PCoA_Figure_S3, file = '16S 18S AD Dynamics/Figures/Figure_S3.png', width = 9, height = 4)


legend.pcoa = as_ggplot(ggpubr::get_legend(PCOA.16S_legend)
ggsave(legend.pcoa, file = '16S 18S AD Dynamics/Figures/Figure_S3_legend.svg', width = 3, height = 5)










