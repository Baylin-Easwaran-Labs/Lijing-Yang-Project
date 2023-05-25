fastqc *gz

################
for i in  input1  chip1
do
trimmomatic PE -threads 4 \
${i}_1.fq.gz ${i}_2.fq.gz \
${i}_paired_clean_R1.fastq.gz \
${i}_unpair_clean_R1.fastq.gz \
${i}_paired_clean_R2.fastq.gz \
${i}_unpair_clean_R2.fastq.gz \
ILLUMINACLIP:/home/miniconda3/envs/chipseq/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10:1:true \
LEADING:3 TRAILING:3 \
SLIDINGWINDOW:4:20 MINLEN:50 TOPHRED33
done


###########


##########
fastqc *.fastq.gz

##############
for i in input1 chip1
do
bowtie2 -p 10 -x ./mm10/mm10 -1 /home/raw_data/${i}_paired_clean_R1.fastq.gz  -2  /home/raw_data/${i}_paired_clean_R2.fastq.gz  | samtools sort -O bam -@ 10 -o - > ${i}_align.bam
done


###############
for i in input1 chip1
do
samtools flagstat ${i}_align.bam > ${i}_align.bam.stat
echo${i}
cat ${i}_align.bam.stat
done


#############
for i in input1 chip1
do
samtools sort -n  ${i}_align.bam -o ${i}_align.sort.bam
samtools fixmate -m ${i}_align.sort.bam ${i}_align.fixmate.bam 
samtools sort  ${i}_align.fixmate.bam -o ${i}_align.positionsort.bam
samtools markdup -r ${i}_align.positionsort.bam ${i}_align.markdup.bam 
done



############################
for i in 1
do
macs2 callpeak -t chip${i}_align.markdup.bam -c input${i}_align.markdup.bam -f BAM -g mm -n chip${i} -B -q 0.01
echo${i}
done



##################
for i in 1
do
samtools index input${i}_align.markdup.bam input${i}_align.markdup.bam.bai
bamCoverage --normalizeUsing CPM -b input${i}_align.markdup.bam -o input${i}_align.markdup.bam.bw
done


for i in 1
do
samtools index chip${i}_align.markdup.bam chip${i}_align.markdup.bam.bai
bamCoverage --normalizeUsing CPM -b chip${i}_align.markdup.bam -o chip${i}_align.markdup.bam.bw
done

