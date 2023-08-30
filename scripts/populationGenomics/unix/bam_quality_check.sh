#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=24072M
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --output=qualCheck_bam.out

module load StdEnv/2020
module load samtools/1.12

out_dir=/home/arteen/projects/def-idworkin/arteen/SociabilityDNA/bwa_map/bam_files/filter30

#temp storage 
in_dir=/home/arteen/projects/def-idworkin/arteen/SociabilityDNA/bwa_map/bam_files

files=(${in_dir}/*_bwa_PE.bam)

for file in ${files[@]} 
do 
	name=${file} 
	base=`basename ${name} _bwa_PE.bam`
	samtools view -b -q 30 -@ 8 ${in_dir}/${base}_bwa_PE.bam > ${out_dir}/${base}.bam
done 
