#!/bin/bash 

#SBATCH --partition=day 
#SBATCH --job-name=count_mut
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10GB 
#SBATCH --time=20:00:00 
#SBATCH --cpus-per-task=7

#module load BBMap 
#module load SAMtools 
#module load BEDTools 

name=${1} 
name2=${2} 
echo read1 ${name} 
echo read2 ${name2}

#echo DEDUPE
#clumpify.sh in1=${name} in2=${name2} out1=${name}.dedupe1 out2=${name2}.dedupe dedupe subs=0

#echo BBDUK 
#bbduk.sh in1=${name}.dedupe1 in2=${name2}.dedupe out=${name}.clean out2=${name2}.clean ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
#bbduk.sh in=${name}.clean in2=${name2}.clean out=${name}.trimmed out2=${name2}.trimmed ftl=10

#module load STAR 

#echo MAPPING 
#STAR --runMode alignReads \
#--genomeDir DMS_genome \
#--readFilesIn ${name}.trimmed ${name2}.trimmed \
#--outFileNamePrefix ${name} \
#--outSAMunmapped Within \
#--alignSJoverhangMin 999 \
#--alignIntronMax 1 \
#--alignIntronMin 999 \
#--runThreadN 16

#echo CONVERT TO BAM 
#samtools view -S -b ${name}Aligned.out.sam | samtools sort -o ${name}Aligned.out_sorted.bam 

#echo INDEX BAM 
#samtools index ${name}Aligned.out_sorted.bam 

#echo PILEUP 
#pileup.sh in=${name}Aligned.out.sam out=${name}_coverage.txt 

#echo ADD MD TAG 
#samtools calmd -b ${name}Aligned.out_sorted.bam short_long_MS2.fa > ${name}Aligned.out_sorted_MDtag.bam 

#echo INDEX BAM 
#samtools index ${name}Aligned.out_sorted_MDtag.bam 

module load Biopython 
module load Pysam 

echo COUNTING MUTATIONS 
python get_dms_PE.py ${name}Aligned.out_sorted_MDtag.bam

echo DONE
