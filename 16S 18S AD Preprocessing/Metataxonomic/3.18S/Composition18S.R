########################################

### Top three most abundant 18S ASV ###

#######################################

library("ggplot2")
library("phyloseq")
library("dplyr")
library("tidyr")
library("wesanderson")

#################
### Load data ###
#################

load('AmpliconAnalysis/18S/ps_18S.reduced.wTree.20062022.RData') #nasuh
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

ggsave(Top3.18S, file = "Figures/top3_18S_barplot.tiff", width = 13, height = 4)







# Phaeobacter only
phaeobacter = phys.decon.noEuk.noneg %>%
  subset_samples(element.type == "bioelement") %>%
  subset_samples(day != "3" & day != "7", prune = TRUE)%>%
  filter_taxa(function (x) {sum(x > 100) > 0}, prune=TRUE) %>%
  transform_sample_counts(function(x) 100*{x/sum(x)} ) %>%
  subset_taxa(genus == "Phaeobacter")%>%
  psmelt()

phaeobacter$biorep <- as.factor(phaeobacter$biorep)
library(viridis)
ggplot(phaeobacter, aes(x = day, y = Abundance, color = OTU)) +
  geom_point(alpha = 0.7, size = 3) +
  ylab("Relative abundance (%)\n") +
  xlab("\n Time (days)")+
  geom_smooth(aes(color = OTU, fill = OTU), method = "lm")+
  #scale_color_manual(values = )+
  theme_bw()+
  theme(legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = c(0.1,0.9))+
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"))+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(0,120,by=10))+
  scale_color_manual(values = c("#E69F00", "#56B4E9"))+
  scale_fill_manual(values = c("#E69F00", "#56B4E9"))








sum.phaeobacter = phaeobacter %>%
  group_by(day) %>%
  summarise(mean.RA = mean(Abundance),
            sd.RA = sd(Abundance))

# Plot with no division between lineages
ggplot(df.genus, aes(x = biorep, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity") +
  ylab("Relative abundance (%) \n") +
  xlab("\n biorep")+
  #scale_fill_manual(values = col_vector)+
  facet_grid(.~day) +
  ggtitle("phylum")

# Plot with division into lineages
comp.lineage.RA = ggplot(df.genus, aes(x = day, y = Abundance, fill = genus)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("Relative abundance (%)\n") +
  xlab("\n Time (day)")+
  scale_fill_manual(values = col_vector)+
  facet_grid(bio.rep~exp.condition) +
  guides(fill=guide_legend(title="Genus", ncol = 4, byrow = TRUE))+
  #scale_x_continuous(limits=c(0,70), breaks=seq(0,70,by =7))+
  theme(strip.background = element_rect(color="black", fill="#FFFFFF"),
        strip.text = element_text(face = "bold"),
        strip.text.y = element_text(angle=0)) +
  theme(legend.text = element_text(face = "bold.italic"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")+
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"))+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
comp.lineage.RA
