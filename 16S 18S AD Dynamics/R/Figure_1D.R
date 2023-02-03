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

ggsave(Top3.18S, file = "16S 18S AD Dynamics/Figures/Figure_1D.png", width = 13, height = 4)

