setwd("/Users/pernillekjersgaardbech/Documents/JH21_Biofilm/AD_amplicons/")
library("devtools")
library("SpiecEasi")
library("phyloseq")
library("NetCoMi")
library("metagMisc") #Split phyloseq object into different groups
library(remotes)
library(microbiome)
library("dplyr")
library(phyloseq)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(magick)
library(cowplot)

ps_18S <- readRDS("18S/phys.18s.rds")
sample_data(ps_18S_succession)
#Add ASV name as column to tax table
ps_18S_tax <- data.frame(tax_table(ps_18S))
ps_18S_tax <- ps_18S_tax %>% mutate(ASV = rownames(ps_18S_tax))
ps_18S <- phyloseq(otu_table(ps_18S), sample_data(ps_18S), tax_table(as.matrix(ps_18S_tax)) )
tax_table(ps_18S)

#Select only succession 
ps_18S_succession <- subset_samples(ps_18S, element.type == "bioelement" & sample.type!="negative")

#Remove outliers
ps_18S_succession <- subset_samples(ps_18S_succession, sample_names(ps_18S_succession) != "Succession-T1-6" & sample_names(ps_18S_succession) != "Succession-T1-8" &
                                      sample_names(ps_18S_succession) != "Succession-T2-9" & sample_names(ps_18S_succession) != "Succession-T2-7")


#Subset to each timepoint and make network for each time point
#ps_18S_succession_T2 <- subset_samples(ps_18S_succession, timepoint == "2") #only 4 sample - not enough for correlations
ps_18S_succession_T3 <- subset_samples(ps_18S_succession, timepoint == "3") # n=9
ps_18S_succession_T4 <- subset_samples(ps_18S_succession, timepoint == "4") # n=9
ps_18S_succession_T5 <- subset_samples(ps_18S_succession, timepoint == "5") # n=9
ps_18S_succession_T6 <- subset_samples(ps_18S_succession, timepoint == "6") # n=9
ps_18S_succession_T7 <- subset_samples(ps_18S_succession, timepoint == "7") # n=9
ps_18S_succession_T8 <- subset_samples(ps_18S_succession, timepoint == "8") # n=9
ps_18S_succession_T9 <- subset_samples(ps_18S_succession, timepoint == "9") # n=9
ps_18S_succession_T10 <- subset_samples(ps_18S_succession, timepoint == "10") # n=9
ps_18S_succession_T11 <- subset_samples(ps_18S_succession, timepoint == "11") # n=9
ps_18S_succession_T12 <- subset_samples(ps_18S_succession, timepoint == "12") # n=9

####_T3######
net_single_T3 <- netConstruct(data = ps_18S_succession_T3, #Samples from the WT treatment at day 1 for the Biofilm
                           measure = "sparcc",  #Correlation method
                           filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 100),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           # filtTax = "highestVar",
                           # filtTaxPar = list(highestVar = 100),
                           # filtSamp = "totalReads",
                           # filtSampPar = list(totalReads = 1000),
                           # filtTax = c("totalReads", "numbSamp"), #Filter on minimum total reads per genus and number of samples the genera should be present
                           # filtTaxPar = list(totalReads = 1000, numbSamp = 20),#Filter on minimum total reads per genus = 100 and number of samples = 5 the genera should be present
                           zeroMethod = "none", normMethod = "none", #No normalization, Zeros not important
                           sparsMethod = "threshold", thresh = 0.7, #Filter all edges with strength below 0.6
                           dissFunc = "signed")

props_single_T3 <- netAnalyze(net_single_T3, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = c("eigenvector", "degree", "closeness"),
                           weightDeg = FALSE, normDeg = FALSE)

#?summary.microNetProps
summary(props_single_T3, numbNodes = 2L)

# p <- plot(props_single, 
#           nodeColor = "cluster", 
#           nodeSize = "eigenvector",
#           title1 = "Network on genus level with SparCC associations", 
#           showTitle = TRUE,
#           cexTitle = 2.3)

plot_props_single_T3 <- plot(props_single_T3,
                                  # nodeFilter = "highestDegree",
                                  # nodeFilterPar = 60,
                                  nodeColor = "cluster",
                                  nodeSize = "eigenvector",
                                  edgeWidth = 0.3,
                                  # edgeFilter = "highestWeight",
                                  # edgeFilterPar = 300,
                                  borderCol = "gray40",
                                  title1 = "Network on ASV level with SparCC measures",
                                  showTitle = TRUE,
                                  cexTitle = 1,
                                  cexLabels = 0.8,
                                  #labelLength = 15,
                                  #labelPattern = c(120,"'",4),
                                  labelScale = FALSE, repulsion = 0.9,
                                  shortenLabels = "none",
                                  #charToRm = "Bacteria_",
                                  nodeSizeSpread = 4, nodeTransp = 65, hubTransp = 50,
                                  hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
                                  edgeTranspHigh = 40,
                                  mar = c(2,2,3,1), set.seed(1000))


labels1 <- plot_props_single_T3$labels$labels1
#labels2 <- plot_ps_AD_16S_succession_2_5$labels$labels2

nodeNames=list(labels1=labels1)

names(nodeNames[[1]])=labels1


# nodeNames[[1]]=ifelse(nodeNames[[1]]!="Phaeobacter_inhibens", "",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]!="Phaeobacter_inhibens", "",nodeNames[[2]])
# 
# nodeNames[[1]]=ifelse(nodeNames[[1]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[2]])


nodecol1 <- plot_props_single_T3$labels$labels1
#nodecol2 <- plot_WTdTDA_SparCC_Biofilm_day4_100$labels$labels2
nodecol1_dat <- as.data.frame(nodecol1)
colnames(nodecol1_dat) <- "ASV"

dim(nodecol1_dat)
length(nodecol1)


tax_18S <- as.data.frame(tax_table(ps_18S_succession_T3))

str(nodecol1_dat)
nodecol1_dat_all <- left_join(nodecol1_dat, tax_18S, by ='ASV')
dim(nodecol1_dat_all)

nodecol_order <- nodecol1_dat_all$class

names(nodecol_order) <- nodecol1

#nodecol_order=ifelse(is.na(nodecol_order), "Unclassified", nodecol_order)
unique(nodecol_order)

Paired_colors <- brewer.pal(n = 12, name = "Paired")
# # OrRd_colors <- brewer.pal(n = 9, name = "OrRd")
# # PuRd_colors <- brewer.pal(n = 6, name = "PuRd")[3:6]
# # Blues_colors <- brewer.pal(n = 6, name = "Blues")[c(1,3)]
# # Spectral_colors <- brewer.pal(n = 11, name = "Spectral")[1:5]
# # PRGn_colors <- brewer.pal(10, name = "PRGn")[1:3]
# # yellow_colors <- c("#FFD300", "#FFDDAF", "#EFFD5F", "#CEB180", "#FFFDD0")
# # # Sea_Sunset_colors <- c("#343077","#9C599E", "#E6759E", "#F89181", "#FFCE8F")
# # # Tropical_colors <- c("#a7e351","#2ace82", "#32e7c8", "#f9c629", "#f29d46") 
# # Sea_Sunset_colors <- rev(c("#343077","#F89181","#9C599E", "#E6759E", "#FFCE8F"))
# # Tropical_colors <- c("#f9c629", "#f29d46", "#a7e351","#2ace82", "#32e7c8") 
# # 
# Peep_beach_colors <- c("#f1bb4e", "#fefd7d", "#9d976f", "#4e7677") 
# Mountain_colors <- c("#ea3333","#f4722b", "#b20f59", "#dd1885", "#9f0772") 
# # 
# PurpleAwesome <- c("#e0b62b","#f25a13", "#2b55e0", "#240785", "#9111ab", "#e33963", "#e35b39", "#e3b039", 
#                    "#580a82", "#060896", "#0d4582", "#c439e3", "#4d44ab", "#580a82") 
# 
# show_col(PurpleAwesome)
# # # 
# # # Paired_colors_expand_other <- sample(c(yellow_colors, PuRd_colors, Blues_colors, Spectral_colors, PRGn_colors))
# # # Paired_colors_expand <- sample(c(Paired_colors, brewer.pal(n = 6, name = "PuRd")))
# # 
# Paired_colors_expand_norm <- (c(Paired_colors,Peep_beach_colors, Mountain_colors))
# # 
# #show_col(Paired_colors)
# # #show_col(brewer.pal(n = 12, name = "Paired"))
# # 
# 
nodecol_order <- factor(nodecol_order, levels=unique(nodecol_order))
# nodecol_order <- factor(nodecol_order, levels=c("Flavobacteriia","Gammaproteobacteria",
#                                                 "Cytophagia","Alphaproteobacteria",
#                                                 "Actinobacteria", "Sphingobacteriia",
#                                                 "Caldilineae","Deltaproteobacteria",
#                                                 "Opitutae", "Planctomycetia",
#                                                 "Deinococci",
#                                                 "Holophagae", "Chlamydiia", "Phycisphaerae",
#                                                 "Betaproteobacteria", "Unclassified", "Epsilonproteobacteria",
#                                                 "Verrucomicrobiae"))
# #Paired_colors_expand_other <- sample(c(Tropical_colors,Peep_beach_colors,Sea_Sunset_colors,Mountain_colors))
# 
# Paired_colors_expand_norm <- (c(Tropical_colors,Sea_Sunset_colors,Mountain_colors, Peep_beach_colors))
# 
# 
# show_col(Paired_colors_expand_norm)
# 

tiff("../Network_props_single_T3.tiff", units="in", width=5, height=5, res=500)

plot(props_single_T3, 
     sameLayout = TRUE, 
     labels = nodeNames,
     nodeColor = "cluster", featVecCol = nodecol_order, colorVec = Paired_colors,
     nodeSize = "eigenvector", 
     edgeWidth = 0.3,
     # edgeFilter = "highestWeight",
     # edgeFilterPar = 300,
     borderCol = "gray40",
     title1 = "Network on ASV level with SparCC measures; thresh = 0.7; HighFreq = 100",
     showTitle = TRUE,
     cexTitle = 1,
     cexLabels = 0.8,
     #labelLength = 15,
     #labelPattern = c(120,"'",4),
     labelScale = FALSE, repulsion = 0.9,
     shortenLabels = "none",
     #charToRm = "Bacteria_",
     nodeSizeSpread = 4, nodeTransp = 20, hubTransp = 20,
     hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
     edgeTranspHigh = 40, rmSingles = TRUE,
     mar = c(2,2,3,1), set.seed(1000))

dev.off()

####_T4######
net_single_T4 <- netConstruct(data = ps_18S_succession_T4, #Samples from the WT treatment at day 1 for the Biofilm
                              measure = "sparcc",  #Correlation method
                              filtTax = "highestFreq",
                              filtTaxPar = list(highestFreq = 100),
                              filtSamp = "totalReads",
                              filtSampPar = list(totalReads = 1000),
                              # filtTax = "highestVar",
                              # filtTaxPar = list(highestVar = 100),
                              # filtSamp = "totalReads",
                              # filtSampPar = list(totalReads = 1000),
                              # filtTax = c("totalReads", "numbSamp"), #Filter on minimum total reads per genus and number of samples the genera should be present
                              # filtTaxPar = list(totalReads = 1000, numbSamp = 20),#Filter on minimum total reads per genus = 100 and number of samples = 5 the genera should be present
                              zeroMethod = "none", normMethod = "none", #No normalization, Zeros not important
                              sparsMethod = "threshold", thresh = 0.7, #Filter all edges with strength below 0.6
                              dissFunc = "signed")

props_single_T4 <- netAnalyze(net_single_T4, 
                              centrLCC = TRUE,
                              clustMethod = "cluster_fast_greedy",
                              hubPar = c("eigenvector", "degree", "closeness"),
                              weightDeg = FALSE, normDeg = FALSE)

#?summary.microNetProps
summary(props_single_T4, numbNodes = 2L)

# p <- plot(props_single, 
#           nodeColor = "cluster", 
#           nodeSize = "eigenvector",
#           title1 = "Network on genus level with SparCC associations", 
#           showTitle = TRUE,
#           cexTitle = 2.3)

plot_props_single_T4 <- plot(props_single_T4,
                             # nodeFilter = "highestDegree",
                             # nodeFilterPar = 60,
                             nodeColor = "cluster",
                             nodeSize = "eigenvector",
                             edgeWidth = 0.3,
                             # edgeFilter = "highestWeight",
                             # edgeFilterPar = 300,
                             borderCol = "gray40",
                             title1 = "Network on ASV level with SparCC measures",
                             showTitle = TRUE,
                             cexTitle = 1,
                             cexLabels = 0.8,
                             #labelLength = 15,
                             #labelPattern = c(120,"'",4),
                             labelScale = FALSE, repulsion = 0.9,
                             shortenLabels = "none",
                             #charToRm = "Bacteria_",
                             nodeSizeSpread = 4, nodeTransp = 65, hubTransp = 50,
                             hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
                             edgeTranspHigh = 40,
                             mar = c(2,2,3,1), set.seed(1000))


labels1 <- plot_props_single_T4$labels$labels1
#labels2 <- plot_ps_AD_16S_succession_2_5$labels$labels2

nodeNames=list(labels1=labels1)

names(nodeNames[[1]])=labels1


# nodeNames[[1]]=ifelse(nodeNames[[1]]!="Phaeobacter_inhibens", "",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]!="Phaeobacter_inhibens", "",nodeNames[[2]])
# 
# nodeNames[[1]]=ifelse(nodeNames[[1]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[2]])


nodecol1 <- plot_props_single_T4$labels$labels1
#nodecol2 <- plot_WTdTDA_SparCC_Biofilm_day4_100$labels$labels2
nodecol1_dat <- as.data.frame(nodecol1)
colnames(nodecol1_dat) <- "ASV"

dim(nodecol1_dat)
length(nodecol1)


tax_18S <- as.data.frame(tax_table(ps_18S_succession_T4))

str(nodecol1_dat)
nodecol1_dat_all <- left_join(nodecol1_dat, tax_18S, by ='ASV')
dim(nodecol1_dat_all)

nodecol_order <- nodecol1_dat_all$class

names(nodecol_order) <- nodecol1

#nodecol_order=ifelse(is.na(nodecol_order), "Unclassified", nodecol_order)
unique(nodecol_order)

Paired_colors <- brewer.pal(n = 12, name = "Paired")
# # OrRd_colors <- brewer.pal(n = 9, name = "OrRd")
# # PuRd_colors <- brewer.pal(n = 6, name = "PuRd")[3:6]
# # Blues_colors <- brewer.pal(n = 6, name = "Blues")[c(1,3)]
# # Spectral_colors <- brewer.pal(n = 11, name = "Spectral")[1:5]
# # PRGn_colors <- brewer.pal(10, name = "PRGn")[1:3]
# # yellow_colors <- c("#FFD300", "#FFDDAF", "#EFFD5F", "#CEB180", "#FFFDD0")
# # # Sea_Sunset_colors <- c("#343077","#9C599E", "#E6759E", "#F89181", "#FFCE8F")
# # # Tropical_colors <- c("#a7e351","#2ace82", "#32e7c8", "#f9c629", "#f29d46") 
# # Sea_Sunset_colors <- rev(c("#343077","#F89181","#9C599E", "#E6759E", "#FFCE8F"))
# # Tropical_colors <- c("#f9c629", "#f29d46", "#a7e351","#2ace82", "#32e7c8") 
# # 
# Peep_beach_colors <- c("#f1bb4e", "#fefd7d", "#9d976f", "#4e7677") 
# Mountain_colors <- c("#ea3333","#f4722b", "#b20f59", "#dd1885", "#9f0772") 
# # 
# PurpleAwesome <- c("#e0b62b","#f25a13", "#2b55e0", "#240785", "#9111ab", "#e33963", "#e35b39", "#e3b039", 
#                    "#580a82", "#060896", "#0d4582", "#c439e3", "#4d44ab", "#580a82") 
# 
# show_col(PurpleAwesome)
# # # 
# # # Paired_colors_expand_other <- sample(c(yellow_colors, PuRd_colors, Blues_colors, Spectral_colors, PRGn_colors))
# # # Paired_colors_expand <- sample(c(Paired_colors, brewer.pal(n = 6, name = "PuRd")))
# # 
# Paired_colors_expand_norm <- (c(Paired_colors,Peep_beach_colors, Mountain_colors))
# # 
# #show_col(Paired_colors)
# # #show_col(brewer.pal(n = 12, name = "Paired"))
# # 
# 
nodecol_order <- factor(nodecol_order, levels=unique(nodecol_order))
# nodecol_order <- factor(nodecol_order, levels=c("Flavobacteriia","Gammaproteobacteria",
#                                                 "Cytophagia","Alphaproteobacteria",
#                                                 "Actinobacteria", "Sphingobacteriia",
#                                                 "Caldilineae","Deltaproteobacteria",
#                                                 "Opitutae", "Planctomycetia",
#                                                 "Deinococci",
#                                                 "Holophagae", "Chlamydiia", "Phycisphaerae",
#                                                 "Betaproteobacteria", "Unclassified", "Epsilonproteobacteria",
#                                                 "Verrucomicrobiae"))
# #Paired_colors_expand_other <- sample(c(Tropical_colors,Peep_beach_colors,Sea_Sunset_colors,Mountain_colors))
# 
# Paired_colors_expand_norm <- (c(Tropical_colors,Sea_Sunset_colors,Mountain_colors, Peep_beach_colors))
# 
# 
# show_col(Paired_colors_expand_norm)
# 

tiff("../Network_props_single_T4.tiff", units="in", width=5, height=5, res=500)

plot(props_single_T4, 
     sameLayout = TRUE, 
     labels = nodeNames,
     nodeColor = "cluster", featVecCol = nodecol_order, colorVec = Paired_colors,
     nodeSize = "eigenvector", 
     edgeWidth = 0.3,
     # edgeFilter = "highestWeight",
     # edgeFilterPar = 300,
     borderCol = "gray40",
     title1 = "Network on ASV level with SparCC measures; thresh = 0.7; HighFreq = 100",
     showTitle = TRUE,
     cexTitle = 1,
     cexLabels = 0.8,
     #labelLength = 15,
     #labelPattern = c(120,"'",4),
     labelScale = FALSE, repulsion = 0.9,
     shortenLabels = "none",
     #charToRm = "Bacteria_",
     nodeSizeSpread = 4, nodeTransp = 20, hubTransp = 20,
     hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
     edgeTranspHigh = 40, rmSingles = TRUE,
     mar = c(2,2,3,1), set.seed(1000))

dev.off()









####_T6######
net_single_T6 <- netConstruct(data = ps_18S_succession_T6, #Samples from the WT treatment at day 1 for the Biofilm
                              measure = "sparcc",  #Correlation method
                              filtTax = "highestFreq",
                              filtTaxPar = list(highestFreq = 100),
                              filtSamp = "totalReads",
                              filtSampPar = list(totalReads = 1000),
                              # filtTax = "highestVar",
                              # filtTaxPar = list(highestVar = 100),
                              # filtSamp = "totalReads",
                              # filtSampPar = list(totalReads = 1000),
                              # filtTax = c("totalReads", "numbSamp"), #Filter on minimum total reads per genus and number of samples the genera should be present
                              # filtTaxPar = list(totalReads = 1000, numbSamp = 20),#Filter on minimum total reads per genus = 100 and number of samples = 5 the genera should be present
                              zeroMethod = "none", normMethod = "none", #No normalization, Zeros not important
                              sparsMethod = "threshold", thresh = 0.7, #Filter all edges with strength below 0.6
                              dissFunc = "signed")

props_single_T6 <- netAnalyze(net_single_T6, 
                              centrLCC = TRUE,
                              clustMethod = "cluster_fast_greedy",
                              hubPar = c("eigenvector", "degree", "closeness"),
                              weightDeg = FALSE, normDeg = FALSE)

#?summary.microNetProps
summary(props_single_T6, numbNodes = 2L)

# p <- plot(props_single, 
#           nodeColor = "cluster", 
#           nodeSize = "eigenvector",
#           title1 = "Network on genus level with SparCC associations", 
#           showTitle = TRUE,
#           cexTitle = 2.3)

plot_props_single_T6 <- plot(props_single_T6,
                             # nodeFilter = "highestDegree",
                             # nodeFilterPar = 60,
                             nodeColor = "cluster",
                             nodeSize = "eigenvector",
                             edgeWidth = 0.3,
                             # edgeFilter = "highestWeight",
                             # edgeFilterPar = 300,
                             borderCol = "gray40",
                             title1 = "Network on ASV level with SparCC measures",
                             showTitle = TRUE,
                             cexTitle = 1,
                             cexLabels = 0.8,
                             #labelLength = 15,
                             #labelPattern = c(120,"'",4),
                             labelScale = FALSE, repulsion = 0.9,
                             shortenLabels = "none",
                             #charToRm = "Bacteria_",
                             nodeSizeSpread = 4, nodeTransp = 65, hubTransp = 50,
                             hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
                             edgeTranspHigh = 40,
                             mar = c(2,2,3,1), set.seed(1000))


labels1 <- plot_props_single_T6$labels$labels1
#labels2 <- plot_ps_AD_16S_succession_2_5$labels$labels2

nodeNames=list(labels1=labels1)

names(nodeNames[[1]])=labels1


# nodeNames[[1]]=ifelse(nodeNames[[1]]!="Phaeobacter_inhibens", "",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]!="Phaeobacter_inhibens", "",nodeNames[[2]])
# 
# nodeNames[[1]]=ifelse(nodeNames[[1]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[2]])


nodecol1 <- plot_props_single_T6$labels$labels1
#nodecol2 <- plot_WTdTDA_SparCC_Biofilm_day4_100$labels$labels2
nodecol1_dat <- as.data.frame(nodecol1)
colnames(nodecol1_dat) <- "ASV"

dim(nodecol1_dat)
length(nodecol1)


tax_18S <- as.data.frame(tax_table(ps_18S_succession_T6))

str(nodecol1_dat)
nodecol1_dat_all <- left_join(nodecol1_dat, tax_18S, by ='ASV')
dim(nodecol1_dat_all)

nodecol_order <- nodecol1_dat_all$class

names(nodecol_order) <- nodecol1

#nodecol_order=ifelse(is.na(nodecol_order), "Unclassified", nodecol_order)
unique(nodecol_order)

Paired_colors <- brewer.pal(n = 12, name = "Paired")
# # OrRd_colors <- brewer.pal(n = 9, name = "OrRd")
# # PuRd_colors <- brewer.pal(n = 6, name = "PuRd")[3:6]
# # Blues_colors <- brewer.pal(n = 6, name = "Blues")[c(1,3)]
# # Spectral_colors <- brewer.pal(n = 11, name = "Spectral")[1:5]
# # PRGn_colors <- brewer.pal(10, name = "PRGn")[1:3]
# # yellow_colors <- c("#FFD300", "#FFDDAF", "#EFFD5F", "#CEB180", "#FFFDD0")
# # # Sea_Sunset_colors <- c("#343077","#9C599E", "#E6759E", "#F89181", "#FFCE8F")
# # # Tropical_colors <- c("#a7e351","#2ace82", "#32e7c8", "#f9c629", "#f29d46") 
# # Sea_Sunset_colors <- rev(c("#343077","#F89181","#9C599E", "#E6759E", "#FFCE8F"))
# # Tropical_colors <- c("#f9c629", "#f29d46", "#a7e351","#2ace82", "#32e7c8") 
# # 
# Peep_beach_colors <- c("#f1bb4e", "#fefd7d", "#9d976f", "#4e7677") 
# Mountain_colors <- c("#ea3333","#f4722b", "#b20f59", "#dd1885", "#9f0772") 
# # 
# PurpleAwesome <- c("#e0b62b","#f25a13", "#2b55e0", "#240785", "#9111ab", "#e33963", "#e35b39", "#e3b039", 
#                    "#580a82", "#060896", "#0d4582", "#c439e3", "#4d44ab", "#580a82") 
# 
# show_col(PurpleAwesome)
# # # 
# # # Paired_colors_expand_other <- sample(c(yellow_colors, PuRd_colors, Blues_colors, Spectral_colors, PRGn_colors))
# # # Paired_colors_expand <- sample(c(Paired_colors, brewer.pal(n = 6, name = "PuRd")))
# # 
# Paired_colors_expand_norm <- (c(Paired_colors,Peep_beach_colors, Mountain_colors))
# # 
# #show_col(Paired_colors)
# # #show_col(brewer.pal(n = 12, name = "Paired"))
# # 
# 
nodecol_order <- factor(nodecol_order, levels=unique(nodecol_order))
# nodecol_order <- factor(nodecol_order, levels=c("Flavobacteriia","Gammaproteobacteria",
#                                                 "Cytophagia","Alphaproteobacteria",
#                                                 "Actinobacteria", "Sphingobacteriia",
#                                                 "Caldilineae","Deltaproteobacteria",
#                                                 "Opitutae", "Planctomycetia",
#                                                 "Deinococci",
#                                                 "Holophagae", "Chlamydiia", "Phycisphaerae",
#                                                 "Betaproteobacteria", "Unclassified", "Epsilonproteobacteria",
#                                                 "Verrucomicrobiae"))
# #Paired_colors_expand_other <- sample(c(Tropical_colors,Peep_beach_colors,Sea_Sunset_colors,Mountain_colors))
# 
# Paired_colors_expand_norm <- (c(Tropical_colors,Sea_Sunset_colors,Mountain_colors, Peep_beach_colors))
# 
# 
# show_col(Paired_colors_expand_norm)
# 

tiff("../Network_props_single_T6.tiff", units="in", width=5, height=5, res=500)

plot(props_single_T6, 
     sameLayout = TRUE, 
     labels = nodeNames,
     nodeColor = "cluster", featVecCol = nodecol_order, colorVec = Paired_colors,
     nodeSize = "eigenvector", 
     edgeWidth = 0.3,
     # edgeFilter = "highestWeight",
     # edgeFilterPar = 300,
     borderCol = "gray40",
     title1 = "Network on ASV level with SparCC measures; thresh = 0.7; HighFreq = 100",
     showTitle = TRUE,
     cexTitle = 1,
     cexLabels = 0.8,
     #labelLength = 15,
     #labelPattern = c(120,"'",4),
     labelScale = FALSE, repulsion = 0.9,
     shortenLabels = "none",
     #charToRm = "Bacteria_",
     nodeSizeSpread = 4, nodeTransp = 20, hubTransp = 20,
     hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
     edgeTranspHigh = 40, rmSingles = TRUE,
     mar = c(2,2,3,1), set.seed(1000))


dev.off()















####_T7######
net_single_T7 <- netConstruct(data = ps_18S_succession_T7, #Samples from the WT treatment at day 1 for the Biofilm
                              measure = "sparcc",  #Correlation method
                              filtTax = "highestFreq",
                              filtTaxPar = list(highestFreq = 100),
                              filtSamp = "totalReads",
                              filtSampPar = list(totalReads = 1000),
                              # filtTax = "highestVar",
                              # filtTaxPar = list(highestVar = 100),
                              # filtSamp = "totalReads",
                              # filtSampPar = list(totalReads = 1000),
                              # filtTax = c("totalReads", "numbSamp"), #Filter on minimum total reads per genus and number of samples the genera should be present
                              # filtTaxPar = list(totalReads = 1000, numbSamp = 20),#Filter on minimum total reads per genus = 100 and number of samples = 5 the genera should be present
                              zeroMethod = "none", normMethod = "none", #No normalization, Zeros not important
                              sparsMethod = "threshold", thresh = 0.7, #Filter all edges with strength below 0.6
                              dissFunc = "signed")

props_single_T7 <- netAnalyze(net_single_T7, 
                              centrLCC = TRUE,
                              clustMethod = "cluster_fast_greedy",
                              hubPar = c("eigenvector", "degree", "closeness"),
                              weightDeg = FALSE, normDeg = FALSE)

#?summary.microNetProps
summary(props_single_T7, numbNodes = 2L)

# p <- plot(props_single, 
#           nodeColor = "cluster", 
#           nodeSize = "eigenvector",
#           title1 = "Network on genus level with SparCC associations", 
#           showTitle = TRUE,
#           cexTitle = 2.3)

plot_props_single_T7 <- plot(props_single_T7,
                             # nodeFilter = "highestDegree",
                             # nodeFilterPar = 60,
                             nodeColor = "cluster",
                             nodeSize = "eigenvector",
                             edgeWidth = 0.3,
                             # edgeFilter = "highestWeight",
                             # edgeFilterPar = 300,
                             borderCol = "gray40",
                             title1 = "Network on ASV level with SparCC measures",
                             showTitle = TRUE,
                             cexTitle = 1,
                             cexLabels = 0.8,
                             #labelLength = 15,
                             #labelPattern = c(120,"'",4),
                             labelScale = FALSE, repulsion = 0.9,
                             shortenLabels = "none",
                             #charToRm = "Bacteria_",
                             nodeSizeSpread = 4, nodeTransp = 65, hubTransp = 50,
                             hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
                             edgeTranspHigh = 40,
                             mar = c(2,2,3,1), set.seed(1000))


labels1 <- plot_props_single_T7$labels$labels1
#labels2 <- plot_ps_AD_16S_succession_2_5$labels$labels2

nodeNames=list(labels1=labels1)

names(nodeNames[[1]])=labels1


# nodeNames[[1]]=ifelse(nodeNames[[1]]!="Phaeobacter_inhibens", "",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]!="Phaeobacter_inhibens", "",nodeNames[[2]])
# 
# nodeNames[[1]]=ifelse(nodeNames[[1]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[2]])


nodecol1 <- plot_props_single_T7$labels$labels1
#nodecol2 <- plot_WTdTDA_SparCC_Biofilm_day4_100$labels$labels2
nodecol1_dat <- as.data.frame(nodecol1)
colnames(nodecol1_dat) <- "ASV"

dim(nodecol1_dat)
length(nodecol1)


tax_18S <- as.data.frame(tax_table(ps_18S_succession_T7))

str(nodecol1_dat)
nodecol1_dat_all <- left_join(nodecol1_dat, tax_18S, by ='ASV')
dim(nodecol1_dat_all)

nodecol_order <- nodecol1_dat_all$class

names(nodecol_order) <- nodecol1

#nodecol_order=ifelse(is.na(nodecol_order), "Unclassified", nodecol_order)
unique(nodecol_order)

Paired_colors <- brewer.pal(n = 12, name = "Paired")
# # OrRd_colors <- brewer.pal(n = 9, name = "OrRd")
# # PuRd_colors <- brewer.pal(n = 6, name = "PuRd")[3:6]
# # Blues_colors <- brewer.pal(n = 6, name = "Blues")[c(1,3)]
# # Spectral_colors <- brewer.pal(n = 11, name = "Spectral")[1:5]
# # PRGn_colors <- brewer.pal(10, name = "PRGn")[1:3]
# # yellow_colors <- c("#FFD300", "#FFDDAF", "#EFFD5F", "#CEB180", "#FFFDD0")
# # # Sea_Sunset_colors <- c("#343077","#9C599E", "#E6759E", "#F89181", "#FFCE8F")
# # # Tropical_colors <- c("#a7e351","#2ace82", "#32e7c8", "#f9c629", "#f29d46") 
# # Sea_Sunset_colors <- rev(c("#343077","#F89181","#9C599E", "#E6759E", "#FFCE8F"))
# # Tropical_colors <- c("#f9c629", "#f29d46", "#a7e351","#2ace82", "#32e7c8") 
# # 
# Peep_beach_colors <- c("#f1bb4e", "#fefd7d", "#9d976f", "#4e7677") 
# Mountain_colors <- c("#ea3333","#f4722b", "#b20f59", "#dd1885", "#9f0772") 
# # 
# PurpleAwesome <- c("#e0b62b","#f25a13", "#2b55e0", "#240785", "#9111ab", "#e33963", "#e35b39", "#e3b039", 
#                    "#580a82", "#060896", "#0d4582", "#c439e3", "#4d44ab", "#580a82") 
# 
# show_col(PurpleAwesome)
# # # 
# # # Paired_colors_expand_other <- sample(c(yellow_colors, PuRd_colors, Blues_colors, Spectral_colors, PRGn_colors))
# # # Paired_colors_expand <- sample(c(Paired_colors, brewer.pal(n = 6, name = "PuRd")))
# # 
# Paired_colors_expand_norm <- (c(Paired_colors,Peep_beach_colors, Mountain_colors))
# # 
# #show_col(Paired_colors)
# # #show_col(brewer.pal(n = 12, name = "Paired"))
# # 
# 
nodecol_order <- factor(nodecol_order, levels=unique(nodecol_order))
# nodecol_order <- factor(nodecol_order, levels=c("Flavobacteriia","Gammaproteobacteria",
#                                                 "Cytophagia","Alphaproteobacteria",
#                                                 "Actinobacteria", "Sphingobacteriia",
#                                                 "Caldilineae","Deltaproteobacteria",
#                                                 "Opitutae", "Planctomycetia",
#                                                 "Deinococci",
#                                                 "Holophagae", "Chlamydiia", "Phycisphaerae",
#                                                 "Betaproteobacteria", "Unclassified", "Epsilonproteobacteria",
#                                                 "Verrucomicrobiae"))
# #Paired_colors_expand_other <- sample(c(Tropical_colors,Peep_beach_colors,Sea_Sunset_colors,Mountain_colors))
# 
# Paired_colors_expand_norm <- (c(Tropical_colors,Sea_Sunset_colors,Mountain_colors, Peep_beach_colors))
# 
# 
# show_col(Paired_colors_expand_norm)
# 

tiff("../Network_props_single_T7.tiff", units="in", width=5, height=5, res=500)

plot(props_single_T7, 
     sameLayout = TRUE, 
     labels = nodeNames,
     nodeColor = "cluster", featVecCol = nodecol_order, colorVec = Paired_colors,
     nodeSize = "eigenvector", 
     edgeWidth = 0.3,
     # edgeFilter = "highestWeight",
     # edgeFilterPar = 300,
     borderCol = "gray40",
     title1 = "Network on ASV level with SparCC measures; thresh = 0.7; HighFreq = 100",
     showTitle = TRUE,
     cexTitle = 1,
     cexLabels = 0.8,
     #labelLength = 15,
     #labelPattern = c(120,"'",4),
     labelScale = FALSE, repulsion = 0.9,
     shortenLabels = "none",
     #charToRm = "Bacteria_",
     nodeSizeSpread = 4, nodeTransp = 20, hubTransp = 20,
     hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
     edgeTranspHigh = 40, rmSingles = TRUE,
     mar = c(2,2,3,1), set.seed(1000))


dev.off()













####_T8######
net_single_T8 <- netConstruct(data = ps_18S_succession_T8, #Samples from the WT treatment at day 1 for the Biofilm
                              measure = "sparcc",  #Correlation method
                              filtTax = "highestFreq",
                              filtTaxPar = list(highestFreq = 100),
                              filtSamp = "totalReads",
                              filtSampPar = list(totalReads = 1000),
                              # filtTax = "highestVar",
                              # filtTaxPar = list(highestVar = 100),
                              # filtSamp = "totalReads",
                              # filtSampPar = list(totalReads = 1000),
                              # filtTax = c("totalReads", "numbSamp"), #Filter on minimum total reads per genus and number of samples the genera should be present
                              # filtTaxPar = list(totalReads = 1000, numbSamp = 20),#Filter on minimum total reads per genus = 100 and number of samples = 5 the genera should be present
                              zeroMethod = "none", normMethod = "none", #No normalization, Zeros not important
                              sparsMethod = "threshold", thresh = 0.7, #Filter all edges with strength below 0.6
                              dissFunc = "signed")

props_single_T8 <- netAnalyze(net_single_T8, 
                              centrLCC = TRUE,
                              clustMethod = "cluster_fast_greedy",
                              hubPar = c("eigenvector", "degree", "closeness"),
                              weightDeg = FALSE, normDeg = FALSE)

#?summary.microNetProps
summary(props_single_T8, numbNodes = 2L)

# p <- plot(props_single, 
#           nodeColor = "cluster", 
#           nodeSize = "eigenvector",
#           title1 = "Network on genus level with SparCC associations", 
#           showTitle = TRUE,
#           cexTitle = 2.3)

plot_props_single_T8 <- plot(props_single_T8,
                             # nodeFilter = "highestDegree",
                             # nodeFilterPar = 60,
                             nodeColor = "cluster",
                             nodeSize = "eigenvector",
                             edgeWidth = 0.3,
                             # edgeFilter = "highestWeight",
                             # edgeFilterPar = 300,
                             borderCol = "gray40",
                             title1 = "Network on ASV level with SparCC measures",
                             showTitle = TRUE,
                             cexTitle = 1,
                             cexLabels = 0.8,
                             #labelLength = 15,
                             #labelPattern = c(120,"'",4),
                             labelScale = FALSE, repulsion = 0.9,
                             shortenLabels = "none",
                             #charToRm = "Bacteria_",
                             nodeSizeSpread = 4, nodeTransp = 65, hubTransp = 50,
                             hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
                             edgeTranspHigh = 40,
                             mar = c(2,2,3,1), set.seed(1000))


labels1 <- plot_props_single_T8$labels$labels1
#labels2 <- plot_ps_AD_16S_succession_2_5$labels$labels2

nodeNames=list(labels1=labels1)

names(nodeNames[[1]])=labels1


# nodeNames[[1]]=ifelse(nodeNames[[1]]!="Phaeobacter_inhibens", "",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]!="Phaeobacter_inhibens", "",nodeNames[[2]])
# 
# nodeNames[[1]]=ifelse(nodeNames[[1]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[2]])


nodecol1 <- plot_props_single_T8$labels$labels1
#nodecol2 <- plot_WTdTDA_SparCC_Biofilm_day4_100$labels$labels2
nodecol1_dat <- as.data.frame(nodecol1)
colnames(nodecol1_dat) <- "ASV"

dim(nodecol1_dat)
length(nodecol1)


tax_18S <- as.data.frame(tax_table(ps_18S_succession_T8))

str(nodecol1_dat)
nodecol1_dat_all <- left_join(nodecol1_dat, tax_18S, by ='ASV')
dim(nodecol1_dat_all)

nodecol_order <- nodecol1_dat_all$class

names(nodecol_order) <- nodecol1

#nodecol_order=ifelse(is.na(nodecol_order), "Unclassified", nodecol_order)
unique(nodecol_order)

Paired_colors <- brewer.pal(n = 12, name = "Paired")
# # OrRd_colors <- brewer.pal(n = 9, name = "OrRd")
# # PuRd_colors <- brewer.pal(n = 6, name = "PuRd")[3:6]
# # Blues_colors <- brewer.pal(n = 6, name = "Blues")[c(1,3)]
# # Spectral_colors <- brewer.pal(n = 11, name = "Spectral")[1:5]
# # PRGn_colors <- brewer.pal(10, name = "PRGn")[1:3]
# # yellow_colors <- c("#FFD300", "#FFDDAF", "#EFFD5F", "#CEB180", "#FFFDD0")
# # # Sea_Sunset_colors <- c("#343077","#9C599E", "#E6759E", "#F89181", "#FFCE8F")
# # # Tropical_colors <- c("#a7e351","#2ace82", "#32e7c8", "#f9c629", "#f29d46") 
# # Sea_Sunset_colors <- rev(c("#343077","#F89181","#9C599E", "#E6759E", "#FFCE8F"))
# # Tropical_colors <- c("#f9c629", "#f29d46", "#a7e351","#2ace82", "#32e7c8") 
# # 
# Peep_beach_colors <- c("#f1bb4e", "#fefd7d", "#9d976f", "#4e7677") 
# Mountain_colors <- c("#ea3333","#f4722b", "#b20f59", "#dd1885", "#9f0772") 
# # 
# PurpleAwesome <- c("#e0b62b","#f25a13", "#2b55e0", "#240785", "#9111ab", "#e33963", "#e35b39", "#e3b039", 
#                    "#580a82", "#060896", "#0d4582", "#c439e3", "#4d44ab", "#580a82") 
# 
# show_col(PurpleAwesome)
# # # 
# # # Paired_colors_expand_other <- sample(c(yellow_colors, PuRd_colors, Blues_colors, Spectral_colors, PRGn_colors))
# # # Paired_colors_expand <- sample(c(Paired_colors, brewer.pal(n = 6, name = "PuRd")))
# # 
# Paired_colors_expand_norm <- (c(Paired_colors,Peep_beach_colors, Mountain_colors))
# # 
# #show_col(Paired_colors)
# # #show_col(brewer.pal(n = 12, name = "Paired"))
# # 
# 
nodecol_order <- factor(nodecol_order, levels=unique(nodecol_order))
# nodecol_order <- factor(nodecol_order, levels=c("Flavobacteriia","Gammaproteobacteria",
#                                                 "Cytophagia","Alphaproteobacteria",
#                                                 "Actinobacteria", "Sphingobacteriia",
#                                                 "Caldilineae","Deltaproteobacteria",
#                                                 "Opitutae", "Planctomycetia",
#                                                 "Deinococci",
#                                                 "Holophagae", "Chlamydiia", "Phycisphaerae",
#                                                 "Betaproteobacteria", "Unclassified", "Epsilonproteobacteria",
#                                                 "Verrucomicrobiae"))
# #Paired_colors_expand_other <- sample(c(Tropical_colors,Peep_beach_colors,Sea_Sunset_colors,Mountain_colors))
# 
# Paired_colors_expand_norm <- (c(Tropical_colors,Sea_Sunset_colors,Mountain_colors, Peep_beach_colors))
# 
# 
# show_col(Paired_colors_expand_norm)
# 

tiff("../Network_props_single_T8.tiff", units="in", width=5, height=5, res=500)

plot(props_single_T8, 
     sameLayout = TRUE, 
     labels = nodeNames,
     nodeColor = "cluster", featVecCol = nodecol_order, colorVec = Paired_colors,
     nodeSize = "eigenvector", 
     edgeWidth = 0.3,
     # edgeFilter = "highestWeight",
     # edgeFilterPar = 300,
     borderCol = "gray40",
     title1 = "Network on ASV level with SparCC measures; thresh = 0.7; HighFreq = 100",
     showTitle = TRUE,
     cexTitle = 1,
     cexLabels = 0.8,
     #labelLength = 15,
     #labelPattern = c(120,"'",4),
     labelScale = FALSE, repulsion = 0.9,
     shortenLabels = "none",
     #charToRm = "Bacteria_",
     nodeSizeSpread = 4, nodeTransp = 20, hubTransp = 20,
     hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
     edgeTranspHigh = 40, rmSingles = TRUE,
     mar = c(2,2,3,1), set.seed(1000))


dev.off()













####_T9######
net_single_T9 <- netConstruct(data = ps_18S_succession_T9, #Samples from the WT treatment at day 1 for the Biofilm
                              measure = "sparcc",  #Correlation method
                              filtTax = "highestFreq",
                              filtTaxPar = list(highestFreq = 100),
                              filtSamp = "totalReads",
                              filtSampPar = list(totalReads = 1000),
                              # filtTax = "highestVar",
                              # filtTaxPar = list(highestVar = 100),
                              # filtSamp = "totalReads",
                              # filtSampPar = list(totalReads = 1000),
                              # filtTax = c("totalReads", "numbSamp"), #Filter on minimum total reads per genus and number of samples the genera should be present
                              # filtTaxPar = list(totalReads = 1000, numbSamp = 20),#Filter on minimum total reads per genus = 100 and number of samples = 5 the genera should be present
                              zeroMethod = "none", normMethod = "none", #No normalization, Zeros not important
                              sparsMethod = "threshold", thresh = 0.7, #Filter all edges with strength below 0.6
                              dissFunc = "signed")

props_single_T9 <- netAnalyze(net_single_T9, 
                              centrLCC = TRUE,
                              clustMethod = "cluster_fast_greedy",
                              hubPar = c("eigenvector", "degree", "closeness"),
                              weightDeg = FALSE, normDeg = FALSE)

#?summary.microNetProps
summary(props_single_T9, numbNodes = 2L)

# p <- plot(props_single, 
#           nodeColor = "cluster", 
#           nodeSize = "eigenvector",
#           title1 = "Network on genus level with SparCC associations", 
#           showTitle = TRUE,
#           cexTitle = 2.3)

plot_props_single_T9 <- plot(props_single_T9,
                             # nodeFilter = "highestDegree",
                             # nodeFilterPar = 60,
                             nodeColor = "cluster",
                             nodeSize = "eigenvector",
                             edgeWidth = 0.3,
                             # edgeFilter = "highestWeight",
                             # edgeFilterPar = 300,
                             borderCol = "gray40",
                             title1 = "Network on ASV level with SparCC measures",
                             showTitle = TRUE,
                             cexTitle = 1,
                             cexLabels = 0.8,
                             #labelLength = 15,
                             #labelPattern = c(120,"'",4),
                             labelScale = FALSE, repulsion = 0.9,
                             shortenLabels = "none",
                             #charToRm = "Bacteria_",
                             nodeSizeSpread = 4, nodeTransp = 65, hubTransp = 50,
                             hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
                             edgeTranspHigh = 40,
                             mar = c(2,2,3,1), set.seed(1000))


labels1 <- plot_props_single_T9$labels$labels1
#labels2 <- plot_ps_AD_16S_succession_2_5$labels$labels2

nodeNames=list(labels1=labels1)

names(nodeNames[[1]])=labels1


# nodeNames[[1]]=ifelse(nodeNames[[1]]!="Phaeobacter_inhibens", "",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]!="Phaeobacter_inhibens", "",nodeNames[[2]])
# 
# nodeNames[[1]]=ifelse(nodeNames[[1]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[2]])


nodecol1 <- plot_props_single_T9$labels$labels1
#nodecol2 <- plot_WTdTDA_SparCC_Biofilm_day4_100$labels$labels2
nodecol1_dat <- as.data.frame(nodecol1)
colnames(nodecol1_dat) <- "ASV"

dim(nodecol1_dat)
length(nodecol1)


tax_18S <- as.data.frame(tax_table(ps_18S_succession_T9))

str(nodecol1_dat)
nodecol1_dat_all <- left_join(nodecol1_dat, tax_18S, by ='ASV')
dim(nodecol1_dat_all)

nodecol_order <- nodecol1_dat_all$class

names(nodecol_order) <- nodecol1

#nodecol_order=ifelse(is.na(nodecol_order), "Unclassified", nodecol_order)
unique(nodecol_order)

Paired_colors <- brewer.pal(n = 12, name = "Paired")
# # OrRd_colors <- brewer.pal(n = 9, name = "OrRd")
# # PuRd_colors <- brewer.pal(n = 6, name = "PuRd")[3:6]
# # Blues_colors <- brewer.pal(n = 6, name = "Blues")[c(1,3)]
# # Spectral_colors <- brewer.pal(n = 11, name = "Spectral")[1:5]
# # PRGn_colors <- brewer.pal(10, name = "PRGn")[1:3]
# # yellow_colors <- c("#FFD300", "#FFDDAF", "#EFFD5F", "#CEB180", "#FFFDD0")
# # # Sea_Sunset_colors <- c("#343077","#9C599E", "#E6759E", "#F89181", "#FFCE8F")
# # # Tropical_colors <- c("#a7e351","#2ace82", "#32e7c8", "#f9c629", "#f29d46") 
# # Sea_Sunset_colors <- rev(c("#343077","#F89181","#9C599E", "#E6759E", "#FFCE8F"))
# # Tropical_colors <- c("#f9c629", "#f29d46", "#a7e351","#2ace82", "#32e7c8") 
# # 
# Peep_beach_colors <- c("#f1bb4e", "#fefd7d", "#9d976f", "#4e7677") 
# Mountain_colors <- c("#ea3333","#f4722b", "#b20f59", "#dd1885", "#9f0772") 
# # 
# PurpleAwesome <- c("#e0b62b","#f25a13", "#2b55e0", "#240785", "#9111ab", "#e33963", "#e35b39", "#e3b039", 
#                    "#580a82", "#060896", "#0d4582", "#c439e3", "#4d44ab", "#580a82") 
# 
# show_col(PurpleAwesome)
# # # 
# # # Paired_colors_expand_other <- sample(c(yellow_colors, PuRd_colors, Blues_colors, Spectral_colors, PRGn_colors))
# # # Paired_colors_expand <- sample(c(Paired_colors, brewer.pal(n = 6, name = "PuRd")))
# # 
# Paired_colors_expand_norm <- (c(Paired_colors,Peep_beach_colors, Mountain_colors))
# # 
# #show_col(Paired_colors)
# # #show_col(brewer.pal(n = 12, name = "Paired"))
# # 
# 
nodecol_order <- factor(nodecol_order, levels=unique(nodecol_order))
# nodecol_order <- factor(nodecol_order, levels=c("Flavobacteriia","Gammaproteobacteria",
#                                                 "Cytophagia","Alphaproteobacteria",
#                                                 "Actinobacteria", "Sphingobacteriia",
#                                                 "Caldilineae","Deltaproteobacteria",
#                                                 "Opitutae", "Planctomycetia",
#                                                 "Deinococci",
#                                                 "Holophagae", "Chlamydiia", "Phycisphaerae",
#                                                 "Betaproteobacteria", "Unclassified", "Epsilonproteobacteria",
#                                                 "Verrucomicrobiae"))
# #Paired_colors_expand_other <- sample(c(Tropical_colors,Peep_beach_colors,Sea_Sunset_colors,Mountain_colors))
# 
# Paired_colors_expand_norm <- (c(Tropical_colors,Sea_Sunset_colors,Mountain_colors, Peep_beach_colors))
# 
# 
# show_col(Paired_colors_expand_norm)
# 

tiff("../Network_props_single_T9.tiff", units="in", width=5, height=5, res=500)

plot(props_single_T9, 
     sameLayout = TRUE, 
     labels = nodeNames,
     nodeColor = "cluster", featVecCol = nodecol_order, colorVec = Paired_colors,
     nodeSize = "eigenvector", 
     edgeWidth = 0.3,
     # edgeFilter = "highestWeight",
     # edgeFilterPar = 300,
     borderCol = "gray40",
     title1 = "Network on ASV level with SparCC measures; thresh = 0.7; HighFreq = 100",
     showTitle = TRUE,
     cexTitle = 1,
     cexLabels = 0.8,
     #labelLength = 15,
     #labelPattern = c(120,"'",4),
     labelScale = FALSE, repulsion = 0.9,
     shortenLabels = "none",
     #charToRm = "Bacteria_",
     nodeSizeSpread = 4, nodeTransp = 20, hubTransp = 20,
     hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
     edgeTranspHigh = 40, rmSingles = TRUE,
     mar = c(2,2,3,1), set.seed(1000))


dev.off()













####_T10######
net_single_T10 <- netConstruct(data = ps_18S_succession_T10, #Samples from the WT treatment at day 1 for the Biofilm
                              measure = "sparcc",  #Correlation method
                              filtTax = "highestFreq",
                              filtTaxPar = list(highestFreq = 100),
                              filtSamp = "totalReads",
                              filtSampPar = list(totalReads = 1000),
                              # filtTax = "highestVar",
                              # filtTaxPar = list(highestVar = 100),
                              # filtSamp = "totalReads",
                              # filtSampPar = list(totalReads = 1000),
                              # filtTax = c("totalReads", "numbSamp"), #Filter on minimum total reads per genus and number of samples the genera should be present
                              # filtTaxPar = list(totalReads = 1000, numbSamp = 20),#Filter on minimum total reads per genus = 100 and number of samples = 5 the genera should be present
                              zeroMethod = "none", normMethod = "none", #No normalization, Zeros not important
                              sparsMethod = "threshold", thresh = 0.7, #Filter all edges with strength below 0.6
                              dissFunc = "signed")

props_single_T10 <- netAnalyze(net_single_T10, 
                              centrLCC = TRUE,
                              clustMethod = "cluster_fast_greedy",
                              hubPar = c("eigenvector", "degree", "closeness"),
                              weightDeg = FALSE, normDeg = FALSE)

#?summary.microNetProps
summary(props_single_T10, numbNodes = 2L)

# p <- plot(props_single, 
#           nodeColor = "cluster", 
#           nodeSize = "eigenvector",
#           title1 = "Network on genus level with SparCC associations", 
#           showTitle = TRUE,
#           cexTitle = 2.3)

plot_props_single_T10 <- plot(props_single_T10,
                             # nodeFilter = "highestDegree",
                             # nodeFilterPar = 60,
                             nodeColor = "cluster",
                             nodeSize = "eigenvector",
                             edgeWidth = 0.3,
                             # edgeFilter = "highestWeight",
                             # edgeFilterPar = 300,
                             borderCol = "gray40",
                             title1 = "Network on ASV level with SparCC measures",
                             showTitle = TRUE,
                             cexTitle = 1,
                             cexLabels = 0.8,
                             #labelLength = 15,
                             #labelPattern = c(120,"'",4),
                             labelScale = FALSE, repulsion = 0.9,
                             shortenLabels = "none",
                             #charToRm = "Bacteria_",
                             nodeSizeSpread = 4, nodeTransp = 65, hubTransp = 50,
                             hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
                             edgeTranspHigh = 40,
                             mar = c(2,2,3,1), set.seed(1000))


labels1 <- plot_props_single_T10$labels$labels1
#labels2 <- plot_ps_AD_16S_succession_2_5$labels$labels2

nodeNames=list(labels1=labels1)

names(nodeNames[[1]])=labels1


# nodeNames[[1]]=ifelse(nodeNames[[1]]!="Phaeobacter_inhibens", "",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]!="Phaeobacter_inhibens", "",nodeNames[[2]])
# 
# nodeNames[[1]]=ifelse(nodeNames[[1]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[2]])


nodecol1 <- plot_props_single_T10$labels$labels1
#nodecol2 <- plot_WTdTDA_SparCC_Biofilm_day4_100$labels$labels2
nodecol1_dat <- as.data.frame(nodecol1)
colnames(nodecol1_dat) <- "ASV"

dim(nodecol1_dat)
length(nodecol1)


tax_18S <- as.data.frame(tax_table(ps_18S_succession_T10))

str(nodecol1_dat)
nodecol1_dat_all <- left_join(nodecol1_dat, tax_18S, by ='ASV')
dim(nodecol1_dat_all)

nodecol_order <- nodecol1_dat_all$class

names(nodecol_order) <- nodecol1

#nodecol_order=ifelse(is.na(nodecol_order), "Unclassified", nodecol_order)
unique(nodecol_order)

Paired_colors <- brewer.pal(n = 12, name = "Paired")
# # OrRd_colors <- brewer.pal(n = 9, name = "OrRd")
# # PuRd_colors <- brewer.pal(n = 6, name = "PuRd")[3:6]
# # Blues_colors <- brewer.pal(n = 6, name = "Blues")[c(1,3)]
# # Spectral_colors <- brewer.pal(n = 11, name = "Spectral")[1:5]
# # PRGn_colors <- brewer.pal(10, name = "PRGn")[1:3]
# # yellow_colors <- c("#FFD300", "#FFDDAF", "#EFFD5F", "#CEB180", "#FFFDD0")
# # # Sea_Sunset_colors <- c("#343077","#9C599E", "#E6759E", "#F89181", "#FFCE8F")
# # # Tropical_colors <- c("#a7e351","#2ace82", "#32e7c8", "#f9c629", "#f29d46") 
# # Sea_Sunset_colors <- rev(c("#343077","#F89181","#9C599E", "#E6759E", "#FFCE8F"))
# # Tropical_colors <- c("#f9c629", "#f29d46", "#a7e351","#2ace82", "#32e7c8") 
# # 
# Peep_beach_colors <- c("#f1bb4e", "#fefd7d", "#9d976f", "#4e7677") 
# Mountain_colors <- c("#ea3333","#f4722b", "#b20f59", "#dd1885", "#9f0772") 
# # 
# PurpleAwesome <- c("#e0b62b","#f25a13", "#2b55e0", "#240785", "#9111ab", "#e33963", "#e35b39", "#e3b039", 
#                    "#580a82", "#060896", "#0d4582", "#c439e3", "#4d44ab", "#580a82") 
# 
# show_col(PurpleAwesome)
# # # 
# # # Paired_colors_expand_other <- sample(c(yellow_colors, PuRd_colors, Blues_colors, Spectral_colors, PRGn_colors))
# # # Paired_colors_expand <- sample(c(Paired_colors, brewer.pal(n = 6, name = "PuRd")))
# # 
# Paired_colors_expand_norm <- (c(Paired_colors,Peep_beach_colors, Mountain_colors))
# # 
# #show_col(Paired_colors)
# # #show_col(brewer.pal(n = 12, name = "Paired"))
# # 
# 
nodecol_order <- factor(nodecol_order, levels=unique(nodecol_order))
# nodecol_order <- factor(nodecol_order, levels=c("Flavobacteriia","Gammaproteobacteria",
#                                                 "Cytophagia","Alphaproteobacteria",
#                                                 "Actinobacteria", "Sphingobacteriia",
#                                                 "Caldilineae","Deltaproteobacteria",
#                                                 "Opitutae", "Planctomycetia",
#                                                 "Deinococci",
#                                                 "Holophagae", "Chlamydiia", "Phycisphaerae",
#                                                 "Betaproteobacteria", "Unclassified", "Epsilonproteobacteria",
#                                                 "Verrucomicrobiae"))
# #Paired_colors_expand_other <- sample(c(Tropical_colors,Peep_beach_colors,Sea_Sunset_colors,Mountain_colors))
# 
# Paired_colors_expand_norm <- (c(Tropical_colors,Sea_Sunset_colors,Mountain_colors, Peep_beach_colors))
# 
# 
# show_col(Paired_colors_expand_norm)
# 

tiff("../Network_props_single_T10.tiff", units="in", width=5, height=5, res=500)

plot(props_single_T10, 
     sameLayout = TRUE, 
     labels = nodeNames,
     nodeColor = "cluster", featVecCol = nodecol_order, colorVec = Paired_colors,
     nodeSize = "eigenvector", 
     edgeWidth = 0.3,
     # edgeFilter = "highestWeight",
     # edgeFilterPar = 300,
     borderCol = "gray40",
     title1 = "Network on ASV level with SparCC measures; thresh = 0.7; HighFreq = 100",
     showTitle = TRUE,
     cexTitle = 1,
     cexLabels = 0.8,
     #labelLength = 15,
     #labelPattern = c(120,"'",4),
     labelScale = FALSE, repulsion = 0.9,
     shortenLabels = "none",
     #charToRm = "Bacteria_",
     nodeSizeSpread = 4, nodeTransp = 20, hubTransp = 20,
     hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
     edgeTranspHigh = 40, rmSingles = TRUE,
     mar = c(2,2,3,1), set.seed(1000))


dev.off()













####_T11######
net_single_T11 <- netConstruct(data = ps_18S_succession_T11, #Samples from the WT treatment at day 1 for the Biofilm
                              measure = "sparcc",  #Correlation method
                              filtTax = "highestFreq",
                              filtTaxPar = list(highestFreq = 100),
                              filtSamp = "totalReads",
                              filtSampPar = list(totalReads = 1000),
                              # filtTax = "highestVar",
                              # filtTaxPar = list(highestVar = 100),
                              # filtSamp = "totalReads",
                              # filtSampPar = list(totalReads = 1000),
                              # filtTax = c("totalReads", "numbSamp"), #Filter on minimum total reads per genus and number of samples the genera should be present
                              # filtTaxPar = list(totalReads = 1000, numbSamp = 20),#Filter on minimum total reads per genus = 100 and number of samples = 5 the genera should be present
                              zeroMethod = "none", normMethod = "none", #No normalization, Zeros not important
                              sparsMethod = "threshold", thresh = 0.7, #Filter all edges with strength below 0.6
                              dissFunc = "signed")

props_single_T11 <- netAnalyze(net_single_T11, 
                              centrLCC = TRUE,
                              clustMethod = "cluster_fast_greedy",
                              hubPar = c("eigenvector", "degree", "closeness"),
                              weightDeg = FALSE, normDeg = FALSE)

#?summary.microNetProps
summary(props_single_T11, numbNodes = 2L)

# p <- plot(props_single, 
#           nodeColor = "cluster", 
#           nodeSize = "eigenvector",
#           title1 = "Network on genus level with SparCC associations", 
#           showTitle = TRUE,
#           cexTitle = 2.3)

plot_props_single_T11 <- plot(props_single_T11,
                             # nodeFilter = "highestDegree",
                             # nodeFilterPar = 60,
                             nodeColor = "cluster",
                             nodeSize = "eigenvector",
                             edgeWidth = 0.3,
                             # edgeFilter = "highestWeight",
                             # edgeFilterPar = 300,
                             borderCol = "gray40",
                             title1 = "Network on ASV level with SparCC measures",
                             showTitle = TRUE,
                             cexTitle = 1,
                             cexLabels = 0.8,
                             #labelLength = 15,
                             #labelPattern = c(120,"'",4),
                             labelScale = FALSE, repulsion = 0.9,
                             shortenLabels = "none",
                             #charToRm = "Bacteria_",
                             nodeSizeSpread = 4, nodeTransp = 65, hubTransp = 50,
                             hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
                             edgeTranspHigh = 40,
                             mar = c(2,2,3,1), set.seed(1000))


labels1 <- plot_props_single_T11$labels$labels1
#labels2 <- plot_ps_AD_16S_succession_2_5$labels$labels2

nodeNames=list(labels1=labels1)

names(nodeNames[[1]])=labels1


# nodeNames[[1]]=ifelse(nodeNames[[1]]!="Phaeobacter_inhibens", "",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]!="Phaeobacter_inhibens", "",nodeNames[[2]])
# 
# nodeNames[[1]]=ifelse(nodeNames[[1]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[2]])


nodecol1 <- plot_props_single_T11$labels$labels1
#nodecol2 <- plot_WTdTDA_SparCC_Biofilm_day4_100$labels$labels2
nodecol1_dat <- as.data.frame(nodecol1)
colnames(nodecol1_dat) <- "ASV"

dim(nodecol1_dat)
length(nodecol1)


tax_18S <- as.data.frame(tax_table(ps_18S_succession_T11))

str(nodecol1_dat)
nodecol1_dat_all <- left_join(nodecol1_dat, tax_18S, by ='ASV')
dim(nodecol1_dat_all)

nodecol_order <- nodecol1_dat_all$class

names(nodecol_order) <- nodecol1

#nodecol_order=ifelse(is.na(nodecol_order), "Unclassified", nodecol_order)
unique(nodecol_order)

Paired_colors <- brewer.pal(n = 12, name = "Paired")
# # OrRd_colors <- brewer.pal(n = 9, name = "OrRd")
# # PuRd_colors <- brewer.pal(n = 6, name = "PuRd")[3:6]
# # Blues_colors <- brewer.pal(n = 6, name = "Blues")[c(1,3)]
# # Spectral_colors <- brewer.pal(n = 11, name = "Spectral")[1:5]
# # PRGn_colors <- brewer.pal(10, name = "PRGn")[1:3]
# # yellow_colors <- c("#FFD300", "#FFDDAF", "#EFFD5F", "#CEB180", "#FFFDD0")
# # # Sea_Sunset_colors <- c("#343077","#9C599E", "#E6759E", "#F89181", "#FFCE8F")
# # # Tropical_colors <- c("#a7e351","#2ace82", "#32e7c8", "#f9c629", "#f29d46") 
# # Sea_Sunset_colors <- rev(c("#343077","#F89181","#9C599E", "#E6759E", "#FFCE8F"))
# # Tropical_colors <- c("#f9c629", "#f29d46", "#a7e351","#2ace82", "#32e7c8") 
# # 
# Peep_beach_colors <- c("#f1bb4e", "#fefd7d", "#9d976f", "#4e7677") 
# Mountain_colors <- c("#ea3333","#f4722b", "#b20f59", "#dd1885", "#9f0772") 
# # 
# PurpleAwesome <- c("#e0b62b","#f25a13", "#2b55e0", "#240785", "#9111ab", "#e33963", "#e35b39", "#e3b039", 
#                    "#580a82", "#060896", "#0d4582", "#c439e3", "#4d44ab", "#580a82") 
# 
# show_col(PurpleAwesome)
# # # 
# # # Paired_colors_expand_other <- sample(c(yellow_colors, PuRd_colors, Blues_colors, Spectral_colors, PRGn_colors))
# # # Paired_colors_expand <- sample(c(Paired_colors, brewer.pal(n = 6, name = "PuRd")))
# # 
# Paired_colors_expand_norm <- (c(Paired_colors,Peep_beach_colors, Mountain_colors))
# # 
# #show_col(Paired_colors)
# # #show_col(brewer.pal(n = 12, name = "Paired"))
# # 
# 
nodecol_order <- factor(nodecol_order, levels=unique(nodecol_order))
# nodecol_order <- factor(nodecol_order, levels=c("Flavobacteriia","Gammaproteobacteria",
#                                                 "Cytophagia","Alphaproteobacteria",
#                                                 "Actinobacteria", "Sphingobacteriia",
#                                                 "Caldilineae","Deltaproteobacteria",
#                                                 "Opitutae", "Planctomycetia",
#                                                 "Deinococci",
#                                                 "Holophagae", "Chlamydiia", "Phycisphaerae",
#                                                 "Betaproteobacteria", "Unclassified", "Epsilonproteobacteria",
#                                                 "Verrucomicrobiae"))
# #Paired_colors_expand_other <- sample(c(Tropical_colors,Peep_beach_colors,Sea_Sunset_colors,Mountain_colors))
# 
# Paired_colors_expand_norm <- (c(Tropical_colors,Sea_Sunset_colors,Mountain_colors, Peep_beach_colors))
# 
# 
# show_col(Paired_colors_expand_norm)
# 

tiff("../Network_props_single_T11.tiff", units="in", width=5, height=5, res=500)

plot(props_single_T11, 
     sameLayout = TRUE, 
     labels = nodeNames,
     nodeColor = "cluster", featVecCol = nodecol_order, colorVec = Paired_colors,
     nodeSize = "eigenvector", 
     edgeWidth = 0.3,
     # edgeFilter = "highestWeight",
     # edgeFilterPar = 300,
     borderCol = "gray40",
     title1 = "Network on ASV level with SparCC measures; thresh = 0.7; HighFreq = 100",
     showTitle = TRUE,
     cexTitle = 1,
     cexLabels = 0.8,
     #labelLength = 15,
     #labelPattern = c(120,"'",4),
     labelScale = FALSE, repulsion = 0.9,
     shortenLabels = "none",
     #charToRm = "Bacteria_",
     nodeSizeSpread = 4, nodeTransp = 20, hubTransp = 20,
     hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
     edgeTranspHigh = 40, rmSingles = TRUE,
     mar = c(2,2,3,1), set.seed(1000))


dev.off()













####_T12######
net_single_T12 <- netConstruct(data = ps_18S_succession_T12, #Samples from the WT treatment at day 1 for the Biofilm
                              measure = "sparcc",  #Correlation method
                              filtTax = "highestFreq",
                              filtTaxPar = list(highestFreq = 100),
                              filtSamp = "totalReads",
                              filtSampPar = list(totalReads = 1000),
                              # filtTax = "highestVar",
                              # filtTaxPar = list(highestVar = 100),
                              # filtSamp = "totalReads",
                              # filtSampPar = list(totalReads = 1000),
                              # filtTax = c("totalReads", "numbSamp"), #Filter on minimum total reads per genus and number of samples the genera should be present
                              # filtTaxPar = list(totalReads = 1000, numbSamp = 20),#Filter on minimum total reads per genus = 100 and number of samples = 5 the genera should be present
                              zeroMethod = "none", normMethod = "none", #No normalization, Zeros not important
                              sparsMethod = "threshold", thresh = 0.7, #Filter all edges with strength below 0.6
                              dissFunc = "signed")

props_single_T12 <- netAnalyze(net_single_T12, 
                              centrLCC = TRUE,
                              clustMethod = "cluster_fast_greedy",
                              hubPar = c("eigenvector", "degree", "closeness"),
                              weightDeg = FALSE, normDeg = FALSE)

#?summary.microNetProps
summary(props_single_T12, numbNodes = 2L)

# p <- plot(props_single, 
#           nodeColor = "cluster", 
#           nodeSize = "eigenvector",
#           title1 = "Network on genus level with SparCC associations", 
#           showTitle = TRUE,
#           cexTitle = 2.3)

plot_props_single_T12 <- plot(props_single_T12,
                             # nodeFilter = "highestDegree",
                             # nodeFilterPar = 60,
                             nodeColor = "cluster",
                             nodeSize = "eigenvector",
                             edgeWidth = 0.3,
                             # edgeFilter = "highestWeight",
                             # edgeFilterPar = 300,
                             borderCol = "gray40",
                             title1 = "Network on ASV level with SparCC measures",
                             showTitle = TRUE,
                             cexTitle = 1,
                             cexLabels = 0.8,
                             #labelLength = 15,
                             #labelPattern = c(120,"'",4),
                             labelScale = FALSE, repulsion = 0.9,
                             shortenLabels = "none",
                             #charToRm = "Bacteria_",
                             nodeSizeSpread = 4, nodeTransp = 65, hubTransp = 50,
                             hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
                             edgeTranspHigh = 40,
                             mar = c(2,2,3,1), set.seed(1000))


labels1 <- plot_props_single_T12$labels$labels1
#labels2 <- plot_ps_AD_16S_succession_2_5$labels$labels2

nodeNames=list(labels1=labels1)

names(nodeNames[[1]])=labels1


# nodeNames[[1]]=ifelse(nodeNames[[1]]!="Phaeobacter_inhibens", "",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]!="Phaeobacter_inhibens", "",nodeNames[[2]])
# 
# nodeNames[[1]]=ifelse(nodeNames[[1]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[1]])
# nodeNames[[2]]=ifelse(nodeNames[[2]]=="Phaeobacter_inhibens", "P. inhibens OTU",nodeNames[[2]])


nodecol1 <- plot_props_single_T12$labels$labels1
#nodecol2 <- plot_WTdTDA_SparCC_Biofilm_day4_100$labels$labels2
nodecol1_dat <- as.data.frame(nodecol1)
colnames(nodecol1_dat) <- "ASV"

dim(nodecol1_dat)
length(nodecol1)


tax_18S <- as.data.frame(tax_table(ps_18S_succession_T12))

str(nodecol1_dat)
nodecol1_dat_all <- left_join(nodecol1_dat, tax_18S, by ='ASV')
dim(nodecol1_dat_all)

nodecol_order <- nodecol1_dat_all$class

names(nodecol_order) <- nodecol1

#nodecol_order=ifelse(is.na(nodecol_order), "Unclassified", nodecol_order)
unique(nodecol_order)

Paired_colors <- brewer.pal(n = 12, name = "Paired")
# # OrRd_colors <- brewer.pal(n = 9, name = "OrRd")
# # PuRd_colors <- brewer.pal(n = 6, name = "PuRd")[3:6]
# # Blues_colors <- brewer.pal(n = 6, name = "Blues")[c(1,3)]
# # Spectral_colors <- brewer.pal(n = 11, name = "Spectral")[1:5]
# # PRGn_colors <- brewer.pal(10, name = "PRGn")[1:3]
# # yellow_colors <- c("#FFD300", "#FFDDAF", "#EFFD5F", "#CEB180", "#FFFDD0")
# # # Sea_Sunset_colors <- c("#343077","#9C599E", "#E6759E", "#F89181", "#FFCE8F")
# # # Tropical_colors <- c("#a7e351","#2ace82", "#32e7c8", "#f9c629", "#f29d46") 
# # Sea_Sunset_colors <- rev(c("#343077","#F89181","#9C599E", "#E6759E", "#FFCE8F"))
# # Tropical_colors <- c("#f9c629", "#f29d46", "#a7e351","#2ace82", "#32e7c8") 
# # 
# Peep_beach_colors <- c("#f1bb4e", "#fefd7d", "#9d976f", "#4e7677") 
# Mountain_colors <- c("#ea3333","#f4722b", "#b20f59", "#dd1885", "#9f0772") 
# # 
# PurpleAwesome <- c("#e0b62b","#f25a13", "#2b55e0", "#240785", "#9111ab", "#e33963", "#e35b39", "#e3b039", 
#                    "#580a82", "#060896", "#0d4582", "#c439e3", "#4d44ab", "#580a82") 
# 
# show_col(PurpleAwesome)
# # # 
# # # Paired_colors_expand_other <- sample(c(yellow_colors, PuRd_colors, Blues_colors, Spectral_colors, PRGn_colors))
# # # Paired_colors_expand <- sample(c(Paired_colors, brewer.pal(n = 6, name = "PuRd")))
# # 
# Paired_colors_expand_norm <- (c(Paired_colors,Peep_beach_colors, Mountain_colors))
# # 
# #show_col(Paired_colors)
# # #show_col(brewer.pal(n = 12, name = "Paired"))
# # 
# 
nodecol_order <- factor(nodecol_order, levels=unique(nodecol_order))
# nodecol_order <- factor(nodecol_order, levels=c("Flavobacteriia","Gammaproteobacteria",
#                                                 "Cytophagia","Alphaproteobacteria",
#                                                 "Actinobacteria", "Sphingobacteriia",
#                                                 "Caldilineae","Deltaproteobacteria",
#                                                 "Opitutae", "Planctomycetia",
#                                                 "Deinococci",
#                                                 "Holophagae", "Chlamydiia", "Phycisphaerae",
#                                                 "Betaproteobacteria", "Unclassified", "Epsilonproteobacteria",
#                                                 "Verrucomicrobiae"))
# #Paired_colors_expand_other <- sample(c(Tropical_colors,Peep_beach_colors,Sea_Sunset_colors,Mountain_colors))
# 
# Paired_colors_expand_norm <- (c(Tropical_colors,Sea_Sunset_colors,Mountain_colors, Peep_beach_colors))
# 
# 
# show_col(Paired_colors_expand_norm)
# 

tiff("../Network_props_single_T12.tiff", units="in", width=5, height=5, res=500)

plot(props_single_T12, 
     sameLayout = TRUE, 
     labels = nodeNames,
     nodeColor = "cluster", featVecCol = nodecol_order, colorVec = Paired_colors,
     nodeSize = "eigenvector", 
     edgeWidth = 0.3,
     # edgeFilter = "highestWeight",
     # edgeFilterPar = 300,
     borderCol = "gray40",
     title1 = "Network on ASV level with SparCC measures; thresh = 0.7; HighFreq = 100",
     showTitle = TRUE,
     cexTitle = 1,
     cexLabels = 0.8,
     #labelLength = 15,
     #labelPattern = c(120,"'",4),
     labelScale = FALSE, repulsion = 0.9,
     shortenLabels = "none",
     #charToRm = "Bacteria_",
     nodeSizeSpread = 4, nodeTransp = 20, hubTransp = 20,
     hubBorderWidth = 2, cexNodes = 1, edgeTranspLow = 0,
     edgeTranspHigh = 40, rmSingles = TRUE,
     mar = c(2,2,3,1), set.seed(1000))


dev.off()












