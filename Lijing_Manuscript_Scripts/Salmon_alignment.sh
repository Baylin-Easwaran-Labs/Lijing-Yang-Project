#!/bin/sh

#  Salmon_alignment.sh
#
#
#  By hari easwaran 08/14/2020.


################################################
#Set options
################################################
#Set resource requirements
#$ -l mem_free=38.0G
#$ -l h_vmem=40.0G
#$ -l h_fsize=10G

#Name of the job
#$ -N Salmon_alignment

#Send email at end of job
#$ -m e
#$ -M heaswar2@jhmi.edu

#$ -e /Salmon/Analysis/Script_logs/Salmon_alignment_error.txt
#$ -o /Salmon/Analysis/Script_logs/Salmon_alignment_output.txt
################################################

################################################
#Salmon
#https://combine-lab.github.io/salmon/getting_started/
#https://salmon.readthedocs.io/en/latest/salmon.html
################################################
# Start script
date
echo "running Salmon_alignment.sh job now"
hostname

# Define path to fastq data.
#data_path=/RNA-Seq/X202SC20060371-Z01-F001/raw_data
#data_path=/RNA-Seq/Analysis/fastp_trimmed_data
data_path=/RNA-Seq/Analysis/trimGalore_trimmed_data

# Define path for output data.
#outputDir=/RNA-Seq/Analysis/SalmonAlignments/fastp_trimmed_Alignments
outputDir=/RNA-Seq/Analysis/SalmonAlignments/trimGalore_trimmed_Alignments
mkdir $outputDir

#for i in Y1
for i in Y1 Y10 Y11 Y12 Y13 Y14 Y15 Y16 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y9
do
fastq_file_R1=$data_path/$i/*1.fq.gz
fastq_file_R2=$data_path/$i/*2.fq.gz

sampleSpecificOutputDir=$outputDir/$i/
mkdir -p $sampleSpecificOutputDir

# Run Salmon alignment
## cd to drive containing salmon so that index path is available
cd /legacy/amber2/scratch/baylin/Hari/Programs/salmon/
salmon quant -i salmon_index_gencode.vM23 -l A -p 8 --gcBias -1 $fastq_file_R1 -2 $fastq_file_R2 -o $sampleSpecificOutputDir

done

echo "Done running all alignments"

date
