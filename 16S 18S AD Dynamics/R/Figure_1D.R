########################################

### Top three most abundant 18S ASV ###

#######################################

library("ggplot2")
library("phyloseq")
library("dplyr")
library("tidyr")
#library("wesanderson")

#################
### Load data ###
#################

load('16S 18S AD Preprocessing/Metataxonomic/3.18S/ps_18S.reduced.wTree.20062022.RData')
ps_18S_filtered <- ps_18S


#########################
### Filter out ASV1-3 ###
#########################

#Get ASV as column
ps_18S_filtered.tax <- data.frame(tax_table(ps_18S_filtered))
ps_18S_filtered.tax<-ps_18S_filtered.tax %>% mutate(ASV = rownames(ps_18S_filtered.tax))
tax_table(ps_18S_filtered) <- as.matrix(ps_18S_filtered.tax)

#Make df with top abundant ASV
df.18S = ps_18S_filtered %>%
  #tax_glom("family") %>%
  subset_samples(element.type == "bioelement") %>%
  transform_sample_counts(function(x) 100*{x/sum(x)} ) %>%
  subset_taxa(ASV == "18S_ASV_1" | ASV == "18S_ASV_2" | ASV == "18S_ASV_3") %>%
  psmelt()                                          

#Change tax names
df.18S$ASV[df.18S$ASV == "18S_ASV_1"] <- "ASV1 (Copepoda)"
df.18S$ASV[df.18S$ASV == "18S_ASV_2"] <- "ASV2 (Cheilostomatida)"
df.18S$ASV[df.18S$ASV == "18S_ASV_3"] <- "ASV3 (Copepoda)"

#Colors (ofc...)
df.18S$ASV = factor(df.18S$ASV, levels= c("ASV2 (Cheilostomatida)","ASV1 (Copepoda)", "ASV3 (Copepoda)"))
Purplerain <- c("#e0b62b","#f25a13", "#2b55e0", "#240785", "#9111ab", "#e33963", "#e35b39", "#e3b039",
                "#580a82", "#060896", "#0d4582", "#c439e3", "#4d44ab")
scales::show_col(Purplerain)


#Plot
Top3.18S = ggplot(df.18S, aes(x = as.factor(biorep), y = Abundance, fill = ASV)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance (%) \n") +
  xlab(" ") +
  theme_bw(base_size = 12)+
  scale_fill_manual(values = c("#e0b62b","#c439e3","#580a82"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"), 
        axis.title.x = element_text(margin = margin(t = -10)),
        #axis.text.x = element_blank(),
        legend.position = "right",
        strip.background = element_rect(fill = "white"),
        strip.placement = "inside",
        strip.text.x = element_text(face = "bold", size = 12)) +
  facet_wrap(.~day, strip.position="bottom", ncol =  12)

ggsave(Top3.18S, file = "16S 18S AD Dynamics/Figures/Figure_1D.png", width = 13, height = 3)

ggsave(Top3.18S, file = "16S 18S AD Dynamics/Figures/Figure_1D.svg", width = 13, height = 3)

#Extract fastas 

sequences <-  ps_18S_filtered %>%
  #tax_glom("family") %>%
  subset_samples(element.type == "bioelement") %>%
  transform_sample_counts(function(x) 100*{x/sum(x)} ) %>%
  subset_taxa(ASV == "18S_ASV_1" | ASV == "18S_ASV_2" | ASV == "18S_ASV_3") %>%
  refseq(.)

library(ShortRead)
library(Biostrings)
writeFasta(sequences, "16S 18S AD Dynamics/R/18S_top3_ASVs.fa", mode = "a")


###  Diatom composition ######
ps_18S_filtered %>%
  #tax_glom("family") %>%
  subset_samples(element.type == "bioelement") %>%
  transform_sample_counts(function(x) 100*{x/sum(x)} ) %>%
  tax_table(.) %>% as.data.frame(.)  %>% subset(., apply(., 1, function(x) any(grepl("SAR", x)))) %>% 
select(order, class, family) %>% unique(.)

#Cercozoa size 2-8um
#Alveolates (dinoflagellates, ciliates and apicomplexans) and Rhizarians are the most common microbial eukaryotes in temperate Appalachian karst caves


sample_data(ps_18S_filtered)
ps_18S_filtered %>%
  #tax_glom("family") %>%
  subset_samples(element.type == "bioelement") %>%
  #transform_sample_counts(function(x) 100*{x/sum(x)} ) %>% 
  subset_taxa(., class == "Alveolata") %>% 
  #plot_richness(., x="day",measures=c("Observed"))
  plot_bar(., fill="order", x="day", facet_grid=~order) 
  
  ps_18S_filtered %>%
    #tax_glom("family") %>%
    subset_samples(element.type == "bioelement") %>%
    #transform_sample_counts(function(x) 100*{x/sum(x)} ) %>% 
    subset_taxa(., family == "Bacillariophyceae") %>% 
    plot_richness(., x="day",measures=c("Observed"))

  
  ps_18S_filtered %>%
    #tax_glom("family") %>%
    subset_samples(element.type == "bioelement") %>%
    #transform_sample_counts(function(x) 100*{x/sum(x)} ) %>% 
    subset_taxa(., order == "Ciliophora") %>% 
    plot_richness(., x="day",measures=c("Observed"))


  library(gmp)
  library(tidyverse)
  library(ANCOMBC)
  ps_18S_filtered_bioelements <- ps_18S_filtered %>%
    #tax_glom("family") %>%
    subset_samples(element.type == "bioelement") %>%
    filter_taxa(., function (x) {sum(x > 10) > 3}, prune=TRUE)
  
  sample_data(ps_18S_filtered_bioelements) <- data.frame(sample_data(ps_18S_filtered_bioelements)) %>% 
  mutate(., Phase_major = ifelse(day <= 29, "Early", "Late"))
  
  #remove two first timepoints
  
  ps_18S_filtered_bioelements <- ps_18S_filtered_bioelements %>% subset_samples(!day < 10)
  
  sample_data(ps_18S_filtered_bioelements) <- data.frame(sample_data(ps_18S_filtered_bioelements)) %>% mutate(sample_name = rownames(.))
  ps_18S_filtered_ANCOM <- ancombc(phyloseq = ps_18S_filtered_bioelements, formula = "Phase_major",
                                           p_adj_method = "holm",  lib_cut = 1000, tax_level = "family",
                                           group = "Phase_major", struc_zero = TRUE, neg_lb = FALSE,
                                           tol = 1e-5, max_iter = 100, conserve = TRUE,
                                           alpha = 0.001, global = TRUE)
  
  
  ps_18S_filtered_ANCOM_lfc <- ps_18S_filtered_ANCOM$res$lfc 
  ps_18S_filtered_ANCOM_lfc$sig <- ps_18S_filtered_ANCOM$res$diff_abn$Phase_majorLate
  
  ps_18S_filtered_ANCOM_lfc_sig <- ps_18S_filtered_ANCOM_lfc %>% filter(sig == TRUE)
  ps_18S_filtered_ANCOM_lfc_sig$taxon <- gsub("family:","",ps_18S_filtered_ANCOM_lfc_sig$taxon)
  
  #Significant 18S families in Early vs. Late
  as.data.frame(tax_table(ps_18S_filtered_bioelements))
  subset_taxa(ps_18S_filtered_bioelements, family %in% ps_18S_filtered_ANCOM_lfc_sig$taxon) %>%
    #subset_samples(., Phase == "Early") %>%
    #filter_taxa(., function (x) {sum(x > 1000) > 10}, prune=TRUE) %>%
    subset_taxa(., phylum != "SAR") %>%
    transform_sample_counts(., function(x) x/sum(x)) %>%  
    #na.omit(.) %>%  
    plot_bar(., "sample_name", fill="genus") + 
    facet_grid(~day, scales="free") + 
    theme(
      #legend.position = "none", 
      axis.text.x=element_blank())
  
  #Significant 18S families in Early vs. Late
  subset_taxa(ps_18S_filtered_bioelements, family %in% ps_18S_filtered_ANCOM_lfc_sig$taxon) %>%
    #subset_samples(., Phase == "Early") %>%
    #filter_taxa(., function (x) {sum(x > 0) > 10}, prune=TRUE) %>%
    subset_taxa(., family == "Ochromonadales") %>%
    transform_sample_counts(., function(x) x/sum(x)) %>%  
    plot_bar(., "sample_name", fill="genus") + 
    facet_grid(~day, scales="free") + 
    theme(
      #legend.position = "none", 
      axis.text.x=element_blank())
  
