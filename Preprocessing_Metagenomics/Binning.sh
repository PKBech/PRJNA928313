
############################################
#Building up contig database in Anvi'o
############################################

anvi-script-reformat-fasta megahit_co-assembly/final.contigs.fa -o megahit_co-assembly/final.contigs-fixed.fa -l 0 --simplify-names
anvi-gen-contigs-database -f megahit_co-assembly/final.contigs-fixed.fa -o CONTIGS.db -n 'CONTIG DB of Biofilm samples'


anvi-run-hmms -c CONTIGS.db --num-threads 16
anvi-run-ncbi-cogs -c CONTIGS.db --num-threads 16
anvi-run-scg-taxonomy -c CONTIGS.db -T 20

 
############################################
#           Mapping with bowtie2
############################################

#!/bin/sh
###Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=dtu_00031 -A dtu_00031
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N mapping_D00107
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=40
### Memory
#PBS -l mem=30gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 12 hours)
#PBS -l walltime=02:00:00:00

module load ngs tools bowtie2/2.4.2 samtools/1.14 bedtools/2.26.0 

#build an index for our contigs
bowtie2-build fixed.contigs.fa assembly

#Run mapping with bowtie2

FQ='path/to/fastq'
BAM='path/to/profiledbam'

find $FQ/*fq.gz > temp
sed 's/.fq.gz//g' temp > temp2
uniq temp2 > sample_list.txt
rm -f temp*
sample_list=$(cat sample_list.txt)

for b in $sample_list
      do
      bowtie2 -x assembly -q -1 contigs-fixed.fa $FQ/"$b"_nohuman_nobryozoan_1.fastq.gz $FQ/"$b"_nohuman_nobryozoan_2.fastq.gz
      --no-unal -p 40 -S $BAM/"$b"-raw.sam
        samtools view --threads 40 -b -o $BAM/"$b"_raw.bam $BAM/"$b"_raw.sam
        rm $BAM/"$b"-raw.sam #to much waste of space
        anvi-init-bam $BAM/"$b"_raw.bam -o $BAM/"$b"_out.bam
        rm $BAM/"$b"_raw.bam #to much waste of space
        anvi-profile -i $BAM/"$b"_out.bam -c CONTIGS.db -o "$b"_PROFILE.db -T 40 --cluster-contigs --profile-SCVs -M 1000 --write-buffer-size 1000
        rm $BAM/"$b"_out.bam
      done
      
      
############################################
#           Profiling in Anvi'o
############################################

### Merge profiles into a single profile database 
# please see http://merenlab.org/2016/06/22/anvio-tutorial-v2/ for further information

anvi-merge <profile_1.db> <profile_2.db> <profile_3.db> <profile_4.db> <profile_5.db> -o MERGED_Profile --enforce-hierarchical-clustering -c CONTIGS.db 

############################################
#      Automatic bining with CONCOCT
############################################

anvi-cluster-contigs -p PROFILE.db -c CONTIGS.db -C CONCOCT --driver concoct --just-do-it -T 40 --clusters 10
