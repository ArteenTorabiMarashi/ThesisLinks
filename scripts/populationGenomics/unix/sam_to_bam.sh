#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=34072M
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=sam_to_bam.out

module load StdEnv/2020
module load samtools/1.12

sam_dir=/home/arteen/projects/def-idworkin/arteen/SociabilityDNA/bwa_map/sam_files

bam_dir=/home/arteen/projects/def-idworkin/arteen/SociabilityDNA/bwa_map/bam_files

files=(${sam_dir}/*.SAM)
for file in ${files[@]}
do
	name=${file}
	base=`basename ${name} .SAM`
	
	samtools view -b -@8 ${sam_dir}/${base}.SAM | samtools sort -o ${bam_dir}/${base}.bam
done

