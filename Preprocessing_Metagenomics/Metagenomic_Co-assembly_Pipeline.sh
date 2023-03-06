######################
#     QC with FastQC
######################

mkdir FastQC
cd FastQC

find *.fq.gz > list
sed 's/.fq.gz//g' list > list2
uniq list2 > sample_list
rm -f list*
sample_list=$(cat sample_list)
echo ${sample_list[@]}


#module load java/1.8.0  fastqc/0.11.8
for a in $sample_list
  do
    fastqc -o FastQC/ "$a".fq.gz -t 16
done  

##############################
#      Adapter removal
##############################

mkdir trimmed
cd trimmed

find *.fq.gz > list
sed 's/_..fq.gz//g' list > list2
uniq list2 > sample_list2
rm -f list*
sample_list=$(cat sample_list2)
echo ${sample_list2[@]}

for a in $sample_list2
  do
    AdapterRemoval --file1 "$a"_1.fq.gz --file2 "$a"_2.fq.gz \
    --basename "$a" --output1 trimmed/"$a"_filtered_1.fq.gz --output2 trimmed/"$a"_filtered_2.fq.gz \
    --adapter1 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
    --adapter2 GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
    --minlength 50 --threads 10 --gzip \
done

##############################
#      Remove duplicates
##############################

mkdir raw 
mkdir filtered

mv *_filtered* filtered
mv *.fq.gz raw

cd filtered

find *.fq.gz > temp
sed 's/_filtered_[1-2].fq.gz//g' temp > temp2
uniq temp2 > sample_list.txt
rm -f temp*
sample_list=$(cat sample_list.txt)


for a in $sample_list
do
cat "$a"_filtered_1.fq.gz| seqkit rmdup -j 20 -s -o "$a"_1_clean.fastq
cat "$a"_filtered_2.fq.gz| seqkit rmdup -j 20 -s -o "$a"_2_clean.fastq
done

###############################################
#      Re-pair fastq after rm duplicates 
###############################################
#Paired reads become disordered after dereplication and need to be re-paired again using the "repair.sh"

find *.fastq > temp
sed 's/_[1-2]_clean.fastq//g' temp > temp2
uniq temp2 > sample_list.txt
rm -f temp*
sample_list=$(cat sample_list.txt)

for a in $sample_list
do
repair.sh in="$a"_1_clean.fastq in2="$a"_2_clean.fastq out=Repaired/"$a"_1_repaired.fq out2=Repaired/"$a"_2_repaired.fq outsingle=Singletons/"$a"_singletons.fq 
rm "$a"_1_clean.fastq
rm "$a"_2_clean.fastq
done


####################################################
#		         		Filtering of Host DNA, using BWA	   		      	  
####################################################

#	Filter against:
#	Human       HG19
# Bryozoan GCA_914767715.1

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz -P host_ref/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/914/767/715/GCA_914767715.1_tzMemMemb1.1/GCA_914767715.1_tzMemMemb1.1_genomic.fna.gz -P host_ref/


mkdir Host_Removal
cd Host_removal

# 1) Build minimap Host_DB
module load samtools/1.9 bedtools/2.28.0 minimap2/2.6
minimap2 -d ref.mmi $REF_DIR/*

### Human Filtering with minimap
REF_DIR='/path/to/reference_genome'
FASTA='/path/to/fasta_files'
OUT='/path/to/directory_output'

cd $FASTA

cd $FASTA/
find *.fq > temp
sed 's/_[1-2]_repaired.fq//g' temp > temp2
uniq temp2 > sample_list.txt
rm -f temp*
sample_list=$(cat sample_list.txt)

#human filtering
Ref='human.fna'
for files in $sample_list
  do
    minimap2 -ax sr $REF_DIR/$Ref $FASTA/"$files"_1_repaired.fq $FASTA/"$files"_2_repaired.fq > $OUT/"$files"_mapped_and_unmapped.sam
      samtools view --threads 14 -bS $OUT/"$files"_mapped_and_unmapped.sam > $OUT/"$files"_mapped_and_unmapped.bam
        rm $OUT/"$files"_mapped_and_unmapped.sam
          samtools view --threads 14 -b -f 13 -F 1280 $OUT/"$files"_mapped_and_unmapped.bam > $OUT/"$files"_unmapped.bam
            rm $OUT/"$files"_mapped_and_unmapped.bam
              samtools sort --threads 14 -n $OUT/"$files"_unmapped.bam -o $OUT/"$files"_sorted.bam
                rm $OUT/"$files"_unmapped.bam
                  bedtools bamtofastq -i $OUT/"$files"_sorted.bam -fq $OUT/"$files"_nohuman_1.fastq -fq2 $OUT/"$files"_nohuman_2.fastq
                    pigz -p 20 $OUT/"$files"_nohuman_*.fastq                  
done


REF_DIR='/path/to/reference_genome'
FASTA='/path/to/fasta_files'
OUT='/path/to/directory_output'

#Bryozoan filtering
Ref='bryozoan.fna'
for files in $sample_list
  do
    minimap2 -ax sr $REF_DIR/$Ref $FASTA/"$files"_nohuman_1.fastq.gz $FASTA/"$files"_nohuman_2.fastq.gz > $OUT/"$files"_mapped_and_unmapped.sam
      samtools view --threads 16 -bS $OUT/"$files"_mapped_and_unmapped.sam > $OUT/"$files"_mapped_and_unmapped.bam
        rm $OUT/"$files"_mapped_and_unmapped.sam
          samtools view --threads 16 -b -f 13 -F 1280 $OUT/"$files"_mapped_and_unmapped.bam > $OUT/"$files"_unmapped.bam
            rm $OUT/"$files"_mapped_and_unmapped.bam
              samtools sort --threads 16 -n $OUT/"$files"_unmapped.bam -o $OUT/"$files"_sorted.bam
                rm $OUT/"$files"_unmapped.bam
                  bedtools bamtofastq -i $OUT/"$files"_sorted.bam -fq $OUT/"$files"_nohuman_nobryozoan_1.fastq -fq2 $OUT/"$files"_nohuman_nobryozoan_2.fastq
                    pigz -p 20 $OUT/"$files"_nohuman_nobryozoan_*.fastq                  
done

#############################################
#     Co-Assembly, using MegaHit  
 ############################################
 
module load megahit/1.1.1 anaconda2/4.4.0
FASTA_DIR='<path/to/FASTQ/>'
WORK_DIR='<path/to/WD/>'

FASTA=$(ls $FASTA_DIR/*.fq | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])')
# Check the FASTA list
echo $FASTA
# run megahit using preset meta-large since it is highly complex communities and the assembly will never finish if the preset are set to mega-sensitive
megahit -r $FASTA -o megahit_co-assembly -t 40 --preset meta-large --min-contig-len 1000

#############################################
#     Add taxonomy infomation with Kaiju 
# for supervised bining and refinement of bins
############################################

kaiju -t /path/to/kaiju/nodes.dmp -f /path/to/kaiju/nr_euk/kaiju_db_nr_euk.fmi -i gene_calls.fa -o gene_calls_nr.out -z 8 -v
addTaxonNames -t /path/to/kaiju/nodes.dmp -n /path/to/kaiju/names.dmp -i gene_calls_nr.out -o gene_calls_nr.names -r superkingdom,phylum,order,class,family,genus,species
anvi-import-taxonomy-for-genes -i gene_calls_nr.names -c CONTIGS.db -p kaiju --just-do-it




