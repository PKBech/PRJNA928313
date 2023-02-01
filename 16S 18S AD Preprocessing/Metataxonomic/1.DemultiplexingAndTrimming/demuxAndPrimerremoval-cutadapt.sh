### cutadapt demultiplexing #####
# 22032022


# demultiplex 
mkdir demux
cutadapt --no-indels -g file:barcodes.fasta -G file:barcodes.fasta -o demux/{name}_R1.fq.gz -p demux/{name}_R2.fq.gz Pool6_*_1.fq.gz Pool6_*_2.fq.gz -j 20


# Remove primers
mkdir primertrim
# Loop for removal of primers

for sample in $(cat samples)
do

    echo "On sample: $sample"
    
    cutadapt -g '^CCTACGGGNGGCWGCAG' -g '^GACTACHVGGGTATCTAATCC' -G '^CCTACGGGNGGCWGCAG' -G '^GACTACHVGGGTATCTAATCC' \
    -m 215 -M 285 -j 19 --discard-untrimmed \
    -o ${sample}_R1-trimmed.fq.gz -p ${sample}_R2-trimmed.fq.gz \
    demux/${sample}_R1.fq.gz demux/${sample}_R2.fq.gz \
    >> cutadapt_primer_trimming_stats.txt 2>&1

done

# Stats on cutadapt 

paste samples <(grep "passing" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")")

