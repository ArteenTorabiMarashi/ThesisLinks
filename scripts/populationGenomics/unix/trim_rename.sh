#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=64072M
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --output=trim.out



module load nixpkgs/16.09

module load trimmomatic/0.36


raw_dir=/home/arteen/projects/def-idworkin/arteen/rawReads/SociabilityAS_DNA

trim_dir=/home/arteen/projects/def-idworkin/arteen/SociabilityDNA/trimmed_reads

adapter=/home/arteen/projects/def-idworkin/arteen/SociabilityDNA/trimmed_reads/TruSeq3-PE-2.fa

files=(${raw_dir}/*_R1.fastq.gz)

for file in ${files[@]}
do
	name=${file}
	base=`basename ${name} _R1.fastq.gz`
	rename=$(awk -v X=$base -F',' '{ if ($1 == X) { print $2} }' /home/arteen/projects/def-idworkin/arteen/SociabilityDNA/trimmed_reads/renaming.csv)

	

	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 \
	${raw_dir}/${base}_R1.fastq.gz \
	${raw_dir}/${base}_R2.fastq.gz \
	${trim_dir}/${rename}_R1_PE.fastq.gz \
	${trim_dir}/${rename}_R1_SE.fastq.gz \
	${trim_dir}/${rename}_R2_PE.fastq.gz \
	${trim_dir}/${rename}_R2_SE.fastq.gz \
	ILLUMINACLIP:${adapter}:2:30:10 \
	LEADING:3 \
	TRAILING:3 \
	MAXINFO:40:0.5 \
	MINLEN:36
done
