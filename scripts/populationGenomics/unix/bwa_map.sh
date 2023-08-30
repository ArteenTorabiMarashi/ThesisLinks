#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=24072M
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output=bwa_map.out


module load StdEnv/2020
module load bwa/0.7.17


# Just so both genomics & RNA work are all the same version genome
genome=/home/arteen/projects/def-idworkin/arteen/SociabilityDNA/bwa_map/genome_index/all_ref.fa

trim_dir=/home/arteen/projects/def-idworkin/arteen/SociabilityDNA/trimmed_reads

map_dir=/home/arteen/projects/def-idworkin/arteen/SociabilityDNA/bwa_map/sam_files

files=(${trim_dir}/*_R1_PE.fastq.gz)

#For loop over every file
for file in ${files[@]} 
do
    name=${file}
    base=`basename ${name} _R1_PE.fastq.gz`
    bwa mem -t 16 -M ${genome} ${trim_dir}/${base}_R1_PE.fastq.gz ${trim_dir}/${base}_R2_PE.fastq.gz > ${map_dir}/${base}_bwa_PE.SAM
done
