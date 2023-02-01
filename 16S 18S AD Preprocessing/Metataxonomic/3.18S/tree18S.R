## Phylogenetic tree and Unifrac distances ###
library("phyloseq")
library("ape")
library("ips")
library("seqinr")

# Phyloseq object to use (cleaned) 
ps_18S <- readRDS("AmpliconAnalysis/18S/phys.18s.rds")
ps_18S =  phyloseq::filter_taxa(ps_18S, function (x) {sum(x > 100) > 0}, prune=TRUE)

#Make DNAbin
seqs = as.DNAbin(refseq(ps_18S))

#alignen with mafft tool and putting on gap penalty of 1
seqs_align <- mafft(seqs, exec = "/usr/local/bin/mafft", 
                    options = c("--adjustdirection"), ep = 1) # ep = gap extension penalty.

View(as.character(seqs_align))
test = as.data.frame(as.character(seqs_align))


#visualize the alginment
image(seqs_align, cex.lab = 0.2)
table(nchar(as.character(refseq(ps_18S)))) # check length 

#collapsing the alignment into a format that we can make a fasta file of
seqs_alignX <- apply(as.character(seqs_align), 1, function(x) paste0(x, collapse = ""))
cat(file="asv_aligned.fasta", paste(paste0(">",names(seqs_alignX)), seqs_alignX, sep="\n"), sep="\n");

#correct ASV names
asv_aligned = seqinr::read.fasta("asv_aligned.fasta")
df.asv_aligned = data.frame(ASV=names(asv_aligned), Seqs=unlist(getSequence(asv_aligned, as.string=T)))
df.asv_aligned$ASV <- gsub("_R_","", as.character(df.asv_aligned$ASV))
cat(file="asv_aligned_correctnames.fasta", paste(paste0(">",df.asv_aligned$ASV), df.asv_aligned$Seqs, sep="\n"), sep="\n");

# Fast tree infers approximately-maximum-likelihood phylogenetic, GTR generalized time-reversible 
#runing Fasttree through the console by using the system2 function
system2("/Users/nnseh/FastTree", args = c("-gtr", "-nt", "-out", "asv_aligned.tree", "asv_aligned_correctnames.fasta"))

#importing the new tree 
treeNew <- read.tree("asv_aligned.tree")
plot(treeNew, cex = 0.5)
#Picking an outgroup automatically as the longest branch
pick_new_outgroup = function(tree.unrooted) {
  require("magrittr")
  require("data.table")
  require("ape") #ape::Ntip
  #tablify parts of tree that we need
  treeDT = 
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup = treeDT[which.max(length)]$id
  return(new.outgroup)
}

new_root = pick_new_outgroup(treeNew)

#rooting the new tree
treeNewR <- root(treeNew, outgroup = new_root, resolve.root = TRUE)

#plotting
plot(treeNewR, cex = 0.2)
#putting tree back in
ps.new = ps_18S
phy_tree(ps.new) <- phy_tree(treeNewR)
ps.new

#Save PS with tree
save(ps.new, file = "AmpliconAnalysis/18S/ps.18-asv.reduced.wTree.RData") 
  


# tree

# Phyloseq object to use (cleaned) 
ps_18S.new <- load("AmpliconAnalysis/18S/ps.18-asv.reduced.wTree.RData")
ps.new.RA = ps.new %>%
  phyloseq::transform_sample_counts(function(x) {x/sum(x)}*100 )
  
  
ps = filter_taxa(ps.new.RA, function(x) sum(x) > .5)
  #phyloseq::transform_sample_counts(function(x) {x/sum(x)}*100)

ps.class = phyloseq::tax_glom(ps, taxrank="class")

tax.df = data.frame(tax_table(ps.new.RA)) %>% replace_na(list(domain = "Uncharacterized",
                                                              phylum = "Uncharacterized", 
                                                              class = "Uncharacterized", 
                                                              order = "Uncharacterized",
                                                              family = "Uncharacterized", 
                                                              genus = "Uncharacterized",
                                                              species = "Uncharacterized",
                                                              strain = "Uncharacterized"))
tax_table(ps.new.RA) = as.matrix(tax.df)

ps.class = phyloseq::tax_glom(ps.new.RA, taxrank="class")



library(microbiome)
ps.class = microbiome::aggregate_taxa(ps.new.RA, "class")  

phyloseq::tax_table(ps.new.RA) = phyloseq::tax_table(ps.class)





plot_tree(ps.new.RA, color="timepoint", label.tips="order", size = 'abundance',
          nodelabf=nodeplotboot(), base.spacing=0.02, ladderize = "right")

tree = ggtree(ps.class, color="black", size=1.5, linetype=1, right=TRUE) + 
  #geom_text2(aes(subset=!isTip, label=label), hjust=-.2, size=4) +
  geom_tiplab(aes(label=class, color = class, 
                  fontface="bold.italic", size=4.5), hjust=-.1) +
  geom_tippoint(aes(color=family), size=6, alpha=1)+
  theme_tree(legend.position= "bottom") +
  #geom_text(aes(label=node), hjust=-.2) +
  scale_fill_manual(values = fam.color) +
  scale_color_manual(values = fam.color) +
  theme(legend.text = element_text(face = "italic", size = 8))+
  xlim(0, 0.8)
tree.top20



