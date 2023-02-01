#Read orientation for ADs

# Ensure prober read orientation
#module load python/v3.5.2 cutadapt/v2.6
samples=$(ls *_R1.fq.gz | sed 's/_R1.fq.gz//')
for sample in ${samples}
do
echo "On sample: $sample"
cutadapt -e 0.15 -g ^GCNTACNNNATNTACACNTCNGG -G ^NANGTCNCCNGTNCGGTA \
--discard-untrimmed \
-o ${sample}_1a_trimmed.fq.gz -p ${sample}_2a_trimmed.fq.gz \
${sample}_R1.fq.gz ${sample}_R2.fq.gz \
>> cutadapt_primer_trimming_stats.txt 2>&1
echo "On sample: $sample"
cutadapt -e 0.15 -g ^NANGTCNCCNGTNCGGTA -G ^GCNTACNNNATNTACACNTCNGG \
--discard-untrimmed \
-o ${sample}_1b_trimmed.fq.gz -p ${sample}_2b_trimmed.fq.gz \
${sample}_R1.fq.gz ${sample}_R2.fq.gz \
>> cutadapt_primer_trimming_stats.txt 2>&1

#Merge both files
cat ${sample}_1a_trimmed.fq.gz ${sample}_2b_trimmed.fq.gz > ${sample}_1_trimmed.fq.gz
cat ${sample}_2a_trimmed.fq.gz ${sample}_1b_trimmed.fq.gz > ${sample}_2_trimmed.fq.gz
rm ${sample}_1a_trimmed.fq.gz ${sample}_2b_trimmed.fq.gz ${sample}_2a_trimmed.fq.gz ${sample}_1b_trimmed.fq.gz
done
