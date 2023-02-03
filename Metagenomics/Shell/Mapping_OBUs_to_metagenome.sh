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

#Make nucl databases
#of MAGs
makeblastdb -in all_MAGs.fa -out all_MAGs -parse_seqids -dbtype nucl
#of the full co-assembly
makeblastdb -in final.contigs.fa -out final.contigs -parse_seqids -dbtype nucl

#Make blast
blastn -db all_MAGs -query ps_AD_filtered_OBU99_succession.fa -out AD_mapped_MAGs.csv -outfmt 6 -num_threads 16
blastn -db final.contigs -query ps_AD_filtered_OBU99_succession.fa -out AD_mapped_final.contigs.csv -outfmt 6 -num_threads 16
