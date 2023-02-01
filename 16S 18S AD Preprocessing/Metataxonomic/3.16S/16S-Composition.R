#################

### Find Myxo ###

################

library("ggplot2")
library("phyloseq")
library("dplyr")
library("tidyr")
library(MicEco)

#################
### Load data ###
#################

load("AmpliconAnalysis/16S/ps_16S.reduced.wTree.20062022.RData") # nasuh
ps_16S_filtered <- ps.ASV.reduced.16S

#Get ASV as column
ps_16S_filtered.tax <- data.frame(tax_table(ps_16S_filtered))
ps_16S_filtered.tax<-ps_16S_filtered.tax %>% mutate(ASV = rownames(ps_16S_filtered.tax))
tax_table(ps_16S_filtered) <- as.matrix(ps_16S_filtered.tax)

#Make df with top abundant ASV
df.16S = ps_16S_filtered %>%
  subset_samples(element.type == "bioelement") %>%
  transform_sample_counts(function(x) 100*{x/sum(x)} ) %>%
  #subset_taxa(phylum == "Myxococcota") %>%
  subset_taxa(genus == "Pseudoalteromonas") %>%
  #tax_glom("genus") %>%
  psmelt()     

#Colors
Purplerain <- c("#9111ab","#f25a13","#240785", "#e0b62b", "#580a82", "#f21395","#2b55e0",  "#e33963", "#f21395", "#7a202c", "#d1bb8a", "#e06f55", "#b78bd6", "#5674a8")
scales::show_col(Purplerain)

#Plot
barplot.myxo = ggplot(df.16S, aes(x = as.factor(biorep), y = Abundance, fill = family)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  ylab("Relative abundance (%) \n") +
  xlab(" ") +
  theme_bw(base_size = 12)+
  scale_fill_manual(values = Purplerain)+
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

ggsave(barplot.myxo, file = "barplot.myxo.svg", width = 18, height = 12, units = "cm")


#
ggplot(df.16S, aes(x = as.factor(biorep), y = Abundance, fill = species)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  ylab("Relative abundance (%) \n") +
  xlab("") +
  ggtitle("Pseudoalteromonas species")+
  theme_bw(base_size = 12)+
  scale_fill_manual(values = Purplerain)+
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

# Pseudomonas

#Make df with top abundant ASV
df.16S.p = ps_16S_filtered %>%
  subset_samples(element.type == "bioelement") %>%
  transform_sample_counts(function(x) 100*{x/sum(x)} ) %>%
  #subset_taxa(phylum == "Myxococcota") %>%
  subset_taxa(genus == "Pseudomonas") %>%
  #tax_glom("genus") %>%
  psmelt()     

#Colors
Purplerain <- c("#9111ab","#f25a13","#240785", "#e0b62b", "#580a82", "#f21395","#2b55e0",  "#e33963", "#f21395", "#7a202c", "#d1bb8a", "#e06f55", "#b78bd6", "#5674a8")
scales::show_col(Purplerain)


#
ggplot(df.16S.p, aes(x = as.factor(biorep), y = Abundance, fill = species)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  ylab("Relative abundance (%) \n") +
  xlab("") +
  ggtitle("Pseudomonas species")+
  theme_bw(base_size = 12)+
  scale_fill_manual(values = Purplerain)+
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

# ANCOM

library(rRDPData)
library(debar)
library(Rcpp)
library(vegan)
library(ANCOMBC)
library(gridExtra)
library(microbiome)
library(stringr)
library(ggpubr)
library(cowplot)

#Only bioelements
ps_16S_succession <- subset_samples(ps_16S_filtered, element.type == "bioelement")
#Subset to peak and late
ps_16S_succession_T5678 <- subset_samples(ps_16S_succession, timepoint %in% c("5","6","7","8") )

#Modify metadata to include phases
ps_16S_succession_T5678_meta <- data.frame(sample_data(ps_16S_succession_T5678))
#Add levels
ps_16S_succession_T5678_meta$Phase[ps_16S_succession_T5678_meta$timepoint < 7] <- "Peak"
ps_16S_succession_T5678_meta$Phase[ps_16S_succession_T5678_meta$timepoint > 6] <- "Late"
#Back in with metadata
ps_16S_succession_T5678 = phyloseq(otu_table(ps_16S_succession_T5678), 
                                   tax_table(ps_16S_succession_T5678), 
                                   refseq(ps_16S_succession_T5678), 
                                   phy_tree(ps_16S_succession_T5678),
                                   sample_data(ps_16S_succession_T5678_meta))

#Make factor
sample_data(ps_16S_succession_T5678)$Phase <- factor(sample_data(ps_16S_succession_T5678)$Phase, levels = c("Peak", "Late"))

##Make df with top abundant ASV
ps_16S_succession_T5678.gen = aggregate_taxa(ps_16S_succession_T5678, "genus")


#Run ANCOMBC differential analysis
out_ps_16S_succession_T5678 <- ancombc(phyloseq = ps_16S_succession_T5678.gen, formula = "Phase",
                                       p_adj_method = "holm", lib_cut = 1000,
                                       group = "Phase", struc_zero = TRUE, neg_lb = FALSE,
                                       tol = 1e-5, max_iter = 100, conserve = TRUE,
                                       alpha = 0.0001, global = TRUE)


#Results
res_out_ps_16S_succession_T5678 <- out_ps_16S_succession_T5678$res
res_out_ps_16S_succession_T5678.p = data.frame(out_ps_16S_succession_T5678$res$q_val)
res_out_ps_16S_succession_T5678.p.sig = data.frame(res_out_ps_16S_succession_T5678.p[res_out_ps_16S_succession_T5678.p$PhaseLate < 0.05,])
res_out_ps_16S_succession_T5678.p.sig$taxon_id = rownames(res_out_ps_16S_succession_T5678.p.sig)
colnames(res_out_ps_16S_succession_T5678.p.sig) = c("q_val", "taxon_id")
res_out_ps_16S_succession_T5678.p.sig = dplyr::arrange(res_out_ps_16S_succession_T5678.p.sig, desc(q_val))


#
df_lfc = data.frame(res_out_ps_16S_succession_T5678$lfc * res_out_ps_16S_succession_T5678$diff_abn, check.names = FALSE) %>% 
  tibble::rownames_to_column("taxon_id")
df_se = data.frame(res_out_ps_16S_succession_T5678$se * res_out_ps_16S_succession_T5678$diff_abn, check.names = FALSE) %>% 
  tibble::rownames_to_column("taxon_id")
colnames(df_se)[-1] = "SE"

df_fig_phase = df_lfc %>% 
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, PhaseLate, SE) %>%
  dplyr::filter(PhaseLate != 0) %>% 
  dplyr::arrange(desc(PhaseLate)) %>%
  dplyr::mutate(direct = ifelse(PhaseLate > 0, "Positive LFC", "Negative LFC"))
df_fig_phase$taxon_id = factor(df_fig_phase$taxon_id, levels = df_fig_phase$taxon_id)
df_fig_phase$direct = factor(df_fig_phase$direct, 
                           levels = c("Positive LFC", "Negative LFC"))
str(df_fig_phase)

#Get phylum for colouring
df_fig_phase.2 =df_fig_phase
rownames(df_fig_phase.2) = df_fig_phase.2$taxon_id
generalist = rownames(df_fig_phase.2)
df.test = subset_taxa(ps_16S_succession_T5678.gen, genus %in% generalist) %>%
  psmelt()
df.test = df.test[,c("phylum", "genus")]
colnames(df.test) = c("Phylum", "taxon_id")
#merge dataframe
df.ancom.16s.genus = merge(df_fig_phase, df.test, by = "taxon_id")

#Colorscheme with 13 colors
Purplerain <- c("#e0b62b","#f25a13", "#2b55e0", "#240785", "#9111ab", "#e33963","#580a82", "#f21395", "#7a202c", "#d1bb8a", "#e06f55", "#b78bd6", "#5674a8")
scales::show_col(Purplerain)

ggplot(data = df.ancom.16s.genus[abs(df.ancom.16s.genus$PhaseLate) > 0.1,], 
               aes(x = taxon_id, y = PhaseLate, fill = Phylum, color = Phylum)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = PhaseLate - SE, ymax = PhaseLate + SE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Significant differential abundant genera from Peak to Late phase") + 
  scale_fill_manual(values = Purplerain) +
  scale_color_manual(values = Purplerain) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))+
  coord_flip()


##########################################################
### ASV ANCOM with only structural zeros in late phase ###
##########################################################

#All ASVs
out_ps_16S_succession_T5678.asv <- ancombc(phyloseq = ps_16S_succession_T5678, formula = "Phase",
                                           p_adj_method = "holm", lib_cut = 1000,
                                           group = "Phase", struc_zero = TRUE, neg_lb = FALSE,
                                           tol = 1e-5, max_iter = 100, conserve = TRUE,
                                           alpha = 0.0001, global = TRUE)


#Results
res_out_ps_16S_succession_T5678.asv <- out_ps_16S_succession_T5678.asv$res

#P adjusted values
res_out_ps_16S_succession_T5678.p. = data.frame(out_ps_16S_succession_T5678.asv$res$q_val)
res_out_ps_16S_succession_T5678.p.sig = data.frame(res_out_ps_16S_succession_T5678.p[res_out_ps_16S_succession_T5678.p$PhaseLate < 0.05,])
res_out_ps_16S_succession_T5678.p.sig$taxon_id = rownames(res_out_ps_16S_succession_T5678.p.sig)
colnames(res_out_ps_16S_succession_T5678.p.sig) = c("q_val", "taxon_id")
res_out_ps_16S_succession_T5678.p.sig = dplyr::arrange(res_out_ps_16S_succession_T5678.p.sig, desc(q_val))

# Make df for plotting
df_lfc = data.frame(res_out_ps_16S_succession_T5678.asv$lfc * res_out_ps_16S_succession_T5678.asv$diff_abn, check.names = FALSE) %>% 
  tibble::rownames_to_column("taxon_id")
df_se = data.frame(res_out_ps_16S_succession_T5678.asv$se * res_out_ps_16S_succession_T5678.asv$diff_abn, check.names = FALSE) %>% 
  tibble::rownames_to_column("taxon_id")
colnames(df_se)[-1] = "SE"

df_fig_phase.asv = df_lfc %>% 
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, PhaseLate, SE) %>%
  dplyr::filter(PhaseLate != 0) %>% 
  dplyr::arrange(desc(PhaseLate)) %>%
  dplyr::mutate(direct = ifelse(PhaseLate > 0, "Positive LFC", "Negative LFC"))
df_fig_phase.asv$taxon_id = factor(df_fig_phase.asv$taxon_id, levels = df_fig_phase.asv$taxon_id)
df_fig_phase.asv$direct = factor(df_fig_phase.asv$direct, 
                             levels = c("Positive LFC", "Negative LFC"))

#get ASVs with structural zeros from late phase (the ones that disappear)
struc.zero = data.frame(out_ps_16S_succession_T5678.asv$zero_ind)
struc.zero = struc.zero %>% mutate(taxon_id = rownames(struc.zero))
struc.zero[struc.zero$structural_zero..Phase...Late. == TRUE,]
struc.zero[struc.zero$structural_zero..Phase...Peak. == TRUE,]

#Get taxonomy
tax.asv = data.frame(tax_table(ps_16S_succession_T5678))
tax.asv = tax.asv %>% mutate(taxon_id = rownames(tax.asv))

#merge tax data with structural zeros df
struc0late = left_join(struc.zero, tax.asv, "taxon_id")
struc0late = struc0late[struc0late$structural_zero..Phase...Late. == TRUE |struc0late$structural_zero..Phase...Peak. == TRUE,]

#Check what myxo that disappear
struc0late[struc0late$phylum == "Myxococcota" & struc0late$structural_zero..Phase...Late. == TRUE,]
str(struc0late)

#subset ancom results to only asvs with structural zeros in late 
  vectore = struc0late$taxon_id
  df_fig_phase.asv = df_fig_phase.asv[df_fig_phase.asv$taxon_id %in% vectore,]
  str(df_fig_phase.asv)
df.ancom.asv.struc0late = merge(df_fig_phase.asv, struc0late, by = "taxon_id")
str(df.ancom.asv.struc0late)

#Colorscheme with 13 colors for phylum plotting
Purplerain <- c("#e0b62b","#f25a13", "#2b55e0", "#240785", "#5674a8", "#9111ab", "#e33963","#580a82", "#f21395", "#7a202c", "#d1bb8a", "#e06f55", "#b78bd6", "#5674a8")
scales::show_col(Purplerain)
df.ancom.asv.struc0late$phylum = factor(df.ancom.asv.struc0late$phylum, levels = c("Proteobacteria", "Bacteroidota", "Bdellovibrionota",
                                                                                   "Planctomycetota", "Verrucomicrobiota","Myxococcota","Actinobacteriota","Firmicutes",
                                                                                   "Spirochaetota", "Fusobacteriota", "Campylobacterota", "Cyanobacteria","Deinococcota",
                                                                                   "Fibrobacterota"))
#Plot significant ASVs with structural zeros in late phase
ANCOM.struc0late.16S = ggplot(data = df.ancom.asv.struc0late[abs(df.ancom.asv.struc0late$PhaseLate) > 0.1,], 
       aes(x = taxon_id, y = PhaseLate, fill = phylum, color = phylum)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.1), alpha = 1, color = NA) +
  facet_wrap(.~phylum, ncol = 7, scales = "free_y")+
  geom_errorbar(aes(ymin = PhaseLate - SE, ymax = PhaseLate + SE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "\nLog fold change", 
       #title = "Significant differential abundant genera from Peak to Late phase"
       ) + 
  #scale_fill_manual(values = c("#5674a8","#b78bd6", "#e0b62b","#2b55e0","#9111ab","#f25a13","#240785","#580a82","#e06f55")) +
  scale_fill_manual(values = Purplerain) +
  geom_hline(yintercept = 0, color = "grey")+
  theme_bw(base_size = 12) + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold", color = "black"),
        #axis.text.y = element_text(color = "black", size = 8),
        strip.text = element_text(face = "bold", color = "black"),
        strip.background = element_rect(fill = "white"),
        axis.text.y = element_blank(),
        )+
  coord_flip()
ANCOM.struc0late.16S
ggsave(ANCOM.struc0late.16S, file = "Figures/ANCOM.struc0late.16S.tiff", width = 16, height = 6)




