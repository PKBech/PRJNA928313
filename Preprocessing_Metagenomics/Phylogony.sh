#### making phylogeny and metadata decorations

### Scan for hmm of RNA polymerase genes A and B

anvi-run-hmms -c CONTIGS_kaiju.db -H 00_HMMs/00_RNApol_genes_Universal_Markers/HMM_RNAP_ab/ -T 10

#### Make muscle aligment based on all RNAP A and B protein seq found in each MAG from CURATED_COLLECTION

anvi-get-sequences-for-hmm-hits -c CONTIGS_kaiju.db -p PROFILE.db -C BACTERIAL_CURATED_COLLECTION --hmm-sources HMM_RNAP_ab --return-best-hit --get-aa-sequences --align-with muscle --gene-names RNAP-b_all,RNAP-a_all --concatenate -o CONCATENATED_MARKER_PROTEINS.fa

#### Create tree in anvio based on alingment

anvi-gen-phylogenomic-tree -f CONCATENATED_MARKER_PROTEINS.fa -o tree_CONCATENATED_MARKER_PROTEINS

#### Make new directory with tree file and use --manual-mode


mkdir tree

mv tree_CONCATENATED_MARKER_PROTEINS tree/

cd tree/ 


### Make phyloseq tree clean by removing outliers in the tree, list the MAGs to remove and make a new collection and repeat above

anvi-import-collection -p PROFILE.db -c CONTIGS_kaiju.db BACTERIAL.txt -C BACTERIAL_CURATED_COLLECTION


anvi-interactive --manual-mode -t tree_CONCATENATED_MARKER_PROTEINS -p PROFILE.db -d ../metadata.txt

anvi-interactive --manual-mode -t tree_CONCATENATED_MARKER_PROTEINS -p PROFILE.db -d ../SUMMARY/bins_across_samples/detection.txt