library(metacoder)


#Heat tree
load("AmpliconAnalysis/16S/ps_16S.reduced.wTree.20062022.RData") # nasuh
ps_16S_filtered <- ps.ASV.reduced.16S 

ps_16S_filtered <- ps_16S_filtered %>%
  transform_sample_counts(function(x) x/sum(x))

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

##aggregate at genus
ps_16S_succession_T5678.gen = aggregate_taxa(ps_16S_succession_T5678, "genus")

#Get only tax to genus
tax_table(ps_16S_succession_T5678.gen) <- tax_table(ps_16S_succession_T5678.gen)[,1:6]

# Create metacoder environment
#TopASV50 = names(sort(taxa_sums(ps.TDA.NoTDA.day42), TRUE)[1:30])
#ps.day42.top50 = prune_taxa(TopASV50, ps.TDA.NoTDA.day42)

#Convert to metacoder object
metacoder <- parse_phyloseq(ps_16S_succession_T5678.gen)
#metacoder <- parse_phyloseq(physeq_100)
metacoder$data$tax_abund <- calc_taxon_abund(metacoder, data = "otu_table")
metacoder$data$tax_occ <- calc_n_samples(metacoder, "tax_abund", groups = "Phase")

#Compare groups
metacoder$data$diff_table <- compare_groups(metacoder, 
                                            data = "tax_abund", 
                                            cols = metacoder$data$sample_data$sample_id,
                                            groups = metacoder$data$sample_data$Phase)
#Adjust p-value
metacoder$data$diff_table$adjusted_p_value <- p.adjust(metacoder$data$diff_table$wilcox_p_value,
                                                       method = "fdr")

metacoder$data$diff_table$log2_median_ratio[metacoder$data$diff_table$adjusted_p_value > 0.05] <- 0

set.seed(41)
metacoder.plot42 <- heat_tree_matrix(metacoder,
                                   data = "diff_table",
                                   node_size = n_obs,
                                   node_label = taxon_names,
                                   node_color = log2_median_ratio,
                                   node_color_range = diverging_palette(),
                                   node_color_trans = "linear",
                                   node_color_interval = c(-5, 5),
                                   edge_color_interval = c(-5, 5),
                                   key_size = 0.7,
                                   node_label_size_range = c(0.025, 0.03),
                                   overlap_avoidance = 5,
                                   layout = "davidson-harel", # The primary layout algorithm
                                   initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                                   node_size_axis_label = "",
                                   node_color_axis_label = "Log2 ratio median proportions",
                                   output_file = "PhasePeakVsLate.pdf")
metacoder.plot42

