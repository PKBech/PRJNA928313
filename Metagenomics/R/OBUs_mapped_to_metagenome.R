#### OBUs mapped to MAGs results from nBLAST analysis ####

library("devtools")
library("dplyr")
library(ggplot2)


###############
## Load data ##
###############

ps_AD_filtered_OBU99 <- readRDS(file = "Metagenomics/Data/ps_AD_filtered_OBU99.rds")


## Subset to only bio.element study ##
ps_AD_filtered_OBU99_succession <- subset_samples(ps_AD_filtered_OBU99, Subject == "Succession")

##filter low abundant OBUs out
#ps_AD_filtered_OBU99_succession <- filter_taxa(ps_AD_filtered_OBU99_succession, function (x) {sum(x > 1000) > 3}, prune=TRUE)


#Extract fasta file with OBUs for succession and do nblast
# ps_AD_filtered_OBU99_succession %>%
#   refseq() %>%
#   Biostrings::writeXStringSet("Metagenomics/Data/ps_AD_filtered_OBU99_succession.fa", append=FALSE,
#                               compress=FALSE, compression_level=NA, format="fasta")
# 
#make nblast against co-assembly and MAGs in shell --> -->


ps_AD_filtered_OBU99_succession_dat <- as.data.frame((as(otu_table(ps_AD_filtered_OBU99_succession), "matrix")))
ps_AD_filtered_OBU99_succession_dat <-as.data.frame(round(t(ps_AD_filtered_OBU99_succession_dat)))
AD_succession_filter_norm_dat_clean <- ps_AD_filtered_OBU99_succession_dat[,colSums(ps_AD_filtered_OBU99_succession_dat)!=0] 
str(AD_succession_filter_norm_dat_clean)

#Get relative abundances per OBU
ps_AD_succession_filter_norm = transform_sample_counts(ps_AD_filtered_OBU99_succession, function(x) 100000 * x/sum(x+0.1))

# Melt to long format
ps_AD_succession_filter_norm_melt <- psmelt(ps_AD_succession_filter_norm)
#Transform to percentage
ps_AD_succession_filter_norm_melt <- aggregate(Abundance ~ OTU + Timepoint, ps_AD_succession_filter_norm_melt, sum)
ps_AD_succession_filter_norm_melt_pct <- ps_AD_succession_filter_norm_melt %>% 
  mutate(Abundance_percentage = Abundance/100000 *100)

# mapped_16S_OBUs <- read.csv("16S_mapped.out.csv", sep = "\t", header = FALSE)
# str(mapped_16S_OBUs)
# colnames(mapped_16S_OBUs) <- c("qseqid", "sseqid", "pident", "length", "mismatch",
#                                "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
# 
# sseqid <- mapped_16S_OBUs$sseqid
# strain <- substring(sseqid,16)
# 
# mapped_16S_OBUs$strain <- strain

mapped_AD_OBUs <- read.csv("Metagenomics/Data/AD_mapped_MAGs.csv", sep = "\t", header = FALSE)
mapped_AD_OBUs_co_assembly <- read.csv("Metagenomics/Data/AD_mapped_final.contigs.csv", sep = "\t", header = FALSE)
#mapped_AD_OBUs_JH21_MAG_00132 <- read.csv("AD_mapped_JH21_MAG_00132_070722.csv", sep = "\t", header = FALSE)


str(mapped_AD_OBUs)
colnames(mapped_AD_OBUs) <- c("qseqid", "sseqid", "pident", "length", "mismatch",
                               "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

colnames(mapped_AD_OBUs_co_assembly) <- c("qseqid", "sseqid", "pident", "length", "mismatch",
                              "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")


# colnames(mapped_AD_OBUs_JH21_MAG_00132) <- c("qseqid", "sseqid", "pident", "length", "mismatch",
#                                           "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
sseqid <- mapped_AD_OBUs$sseqid
strain <- substring(sseqid,16)

mapped_AD_OBUs$strain <- strain

### remove duplicates from each dataframe based on highest procentage and with mapping quary length > 200

#mapped_AD_OBUs_clean <- mapped_AD_OBUs %>% filter(length>200, pident>97) %>% group_by(pident,length, sseqid) %>% summarise(qseqid = max(qseqid))
mapped_AD_OBUs_clean <- mapped_AD_OBUs %>%filter(length>200, pident>95) %>% group_by(pident,length, sseqid) %>% dplyr::summarise(qseqid = max(qseqid))
mapped_AD_OBUs_clean <- mapped_AD_OBUs_clean[!duplicated(mapped_AD_OBUs_clean$qseqid),]

mapped_AD_OBUs_co_assembly_clean <- mapped_AD_OBUs_co_assembly %>%filter(length>200, pident>95) %>% group_by(pident,length, sseqid) %>% dplyr::summarise(qseqid = max(qseqid))
mapped_AD_OBUs_co_assembly_clean <- mapped_AD_OBUs_co_assembly_clean[!duplicated(mapped_AD_OBUs_co_assembly_clean$qseqid),]

# mapped_AD_OBUs_JH21_MAG_00132_clean <- mapped_AD_OBUs_JH21_MAG_00132 %>%filter(length>200, pident>95) %>% group_by(pident,length, sseqid) %>% dplyr::summarise(qseqid = max(qseqid))
# mapped_AD_OBUs_JH21_MAG_00132_clean <- mapped_AD_OBUs_JH21_MAG_00132_clean[!duplicated(mapped_AD_OBUs_JH21_MAG_00132_clean$qseqid),]


ASV_AD_names_total <- colnames(AD_succession_filter_norm_dat_clean)
length(ASV_AD_names_total)

length(mapped_AD_OBUs_clean$qseqid)

length(mapped_AD_OBUs_clean$qseqid)/
length(ASV_AD_names_total)*100

length(mapped_AD_OBUs_co_assembly_clean$qseqid)/
  length(ASV_AD_names_total)*100

length(ASV_AD_names_total) - length(mapped_AD_OBUs_co_assembly_clean$qseqid)

100 - 67.20602


mapped_AD_OBUs_co_assembly_clean_dat <- data.frame(ASV = mapped_AD_OBUs_co_assembly_clean$qseqid, mapping = rep("OBUs mapped to Co-assembly"))
mapped_AD_OBUs_clean_dat <- data.frame(ASV = mapped_AD_OBUs_clean$qseqid, mapping = rep("OBUs mapped to MAGs"))

unmapped_AD_OBUs <- ASV_AD_names_total[!ASV_AD_names_total %in% mapped_AD_OBUs_co_assembly_clean$qseqid]
length(unmapped_AD_OBUs)


unmapped_AD_OBUs_clean_dat <- data.frame(ASV = unmapped_AD_OBUs, mapping = rep("Unmapped OBUs"))

length(ASV_AD_names_total) - length(mapped_AD_OBUs_co_assembly_clean$qseqid) 
length(unmapped_AD_OBUs_clean_dat$ASV)
length(mapped_AD_OBUs_clean_dat$ASV)
length(mapped_AD_OBUs_co_assembly_clean_dat$ASV)

##Merge all 
mapping_rate_AD_OBUs_succession <-rbind(mapped_AD_OBUs_co_assembly_clean_dat, mapped_AD_OBUs_clean_dat, unmapped_AD_OBUs_clean_dat)


colnames(ps_AD_succession_filter_norm_melt_pct)[1] <- "ASV"
mapping_rate_AD_OBUs_succession_abundance <- left_join(ps_AD_succession_filter_norm_melt_pct, mapping_rate_AD_OBUs_succession)
#na.omit(mapping_rate_AD_OBUs_succession_abundance) %>% group_by(mapping) %>% dplyr::count(ASV) %>% dplyr::count(mapping)

mapping_rate_AD_OBUs_succession_abundance$Abundance_percentage_log <- log2(mapping_rate_AD_OBUs_succession_abundance$Abundance_percentage)

mapping_rate_AD_OBUs_succession_abundance$Abundance_percentage_log <- ifelse(abs(mapping_rate_AD_OBUs_succession_abundance$Abundance_percentage_log)==Inf, NA, mapping_rate_AD_OBUs_succession_abundance$Abundance_percentage_log)


mapping_rate_AD_OBUs_succession_abundance$Time = 0
mapping_rate_AD_OBUs_succession_abundance$Time[mapping_rate_AD_OBUs_succession_abundance$Timepoint=="T1"] <- 0
mapping_rate_AD_OBUs_succession_abundance$Time[mapping_rate_AD_OBUs_succession_abundance$Timepoint=="T2"] <- 7
mapping_rate_AD_OBUs_succession_abundance$Time[mapping_rate_AD_OBUs_succession_abundance$Timepoint=="T3"] <- 10
mapping_rate_AD_OBUs_succession_abundance$Time[mapping_rate_AD_OBUs_succession_abundance$Timepoint=="T4"] <- 15
mapping_rate_AD_OBUs_succession_abundance$Time[mapping_rate_AD_OBUs_succession_abundance$Timepoint=="T5"] <- 23
mapping_rate_AD_OBUs_succession_abundance$Time[mapping_rate_AD_OBUs_succession_abundance$Timepoint=="T6"] <- 29
mapping_rate_AD_OBUs_succession_abundance$Time[mapping_rate_AD_OBUs_succession_abundance$Timepoint=="T7"] <- 44
mapping_rate_AD_OBUs_succession_abundance$Time[mapping_rate_AD_OBUs_succession_abundance$Timepoint=="T8"] <- 58
mapping_rate_AD_OBUs_succession_abundance$Time[mapping_rate_AD_OBUs_succession_abundance$Timepoint=="T9"] <- 71
mapping_rate_AD_OBUs_succession_abundance$Time[mapping_rate_AD_OBUs_succession_abundance$Timepoint=="T10"] <- 85
mapping_rate_AD_OBUs_succession_abundance$Time[mapping_rate_AD_OBUs_succession_abundance$Timepoint=="T11"] <- 99
mapping_rate_AD_OBUs_succession_abundance$Time[mapping_rate_AD_OBUs_succession_abundance$Timepoint=="T12"] <- 113

mapping_rate_AD_OBUs_succession_abundance$Time <- factor(mapping_rate_AD_OBUs_succession_abundance$Time)
summary(unique(mapping_rate_AD_OBUs_succession_abundance$ASV))
mapping_rate_AD_OBUs_succession_abundance$mapping <- factor(mapping_rate_AD_OBUs_succession_abundance$mapping, levels = c("OBUs mapped to MAGs","OBUs mapped to Co-assembly", "Unmapped OBUs"))


ASV_AD_names_total <- colnames(AD_succession_filter_norm_dat_clean)
length(ASV_AD_names_total)

length(mapped_AD_OBUs_clean$qseqid)

length(mapped_AD_OBUs_clean$qseqid)/
  length(ASV_AD_names_total)*100

length(mapped_AD_OBUs_co_assembly_clean$qseqid)/
  length(ASV_AD_names_total)*100

length(ASV_AD_names_total) - length(mapped_AD_OBUs_co_assembly_clean$qseqid)


Mapped_OBUs <- mapping_rate_AD_OBUs_succession_abundance %>% filter(Time != 0 & Time != 7,Abundance_percentage > 0 ) %>% 
  ggplot(aes(x = mapping, y = Abundance_percentage_log)) + 
  geom_violin(col="black") + #ylim(0.01,1) +
  geom_violin() + #ylim(0.01,1) +
  theme_bw() +
  #scale_fill_manual(values = main_col) + 
  ylab("Log(Relative Abundance)") + xlab("Days") +
  #scale_color_discrete(guide = FALSE) +  # turn legend on or off
  #geom_line() +  #facet_nested(. ~ Environment_1 + Day +  Treatment,  scales = "free_x") +
  theme_bw(base_size = 12) + #geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme(#axis.line = element_line(color='black'),
    #axis.text.x = element_text(angle = 45),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.title.x=element_blank(),
    #axis.ticks.x=element_blank(),
    #axis.text.x = element_blank(),
  ) + 
  #guides(fill=guide_legend(title="")) +
 #geom_text(data=pct, aes(x = mapping, label=pct)
  annotate("text", label = "166/934\n~ 17.8%", x =1, size = 3, y = -0.4,  colour = "black") + 
  annotate("text", label = "304/934\n~ 32.5%", x =2, size = 3, y = -0.4,  colour = "black") +
  annotate("text", label = "603/934\n~ 65.5%", x =3, size = 3, y = -0.4,  colour = "black") 
  

ggsave(Mapped_OBUs, file = "Metagenomics/Figures/Mapped_OBUs.svg", width=4, height=4)


