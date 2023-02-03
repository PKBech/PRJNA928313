#### nblast ####

FASTQ_DIR='/home/perbec/PChem0002/Jyllinghavn_metagenomes2021/CONCOCT_Binned_DBs/SUMMARY_BACTERIAL/bin_by_bin/'

ls $FASTQ_DIR > temp
sample_list=$(cat temp)
mkdir new_fasta_headers_MAGs

for i in $sample_list
do
awk -v fname=$i '/^>/ {$0=$0 "_" fname}1' $FASTQ_DIR/$i/$i-contigs.fa > new_fasta_headers_MAGs/"$i"_new.fa
done

cat new_fasta_headers_MAGs/*fa > all_MAGs.fa

makeblastdb -in all_MAGs.fa -out all_MAGs -parse_seqids -dbtype nucl

makeblastdb -in final.contigs.fa -out final.contigs -parse_seqids -dbtype nucl

blastn -db all_MAGs -query ps_AD_filtered_succession_070622.fa -out AD_mapped_MAGs_070722.csv -outfmt 6 -num_threads 16
blastn -db final.contigs -query ps_AD_filtered_succession_070622.fa -out AD_mapped_final.contigs_070722.csv -outfmt 6 -num_threads 16
