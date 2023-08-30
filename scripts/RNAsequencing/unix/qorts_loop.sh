#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=34072M
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=6
#SBATCH --output=qorts.out



trim_dir=~/projects/def-idworkin/arteen/SociabilityRNA/trimmed_reads/subset2

map_dir=~/projects/def-idworkin/arteen/SociabilityRNA/star/subset2

files=(${map_dir}/*Aligned.sortedByCoord.out.bam)


for file in ${files[@]}
do

name=${file}
sample_name=`basename ${name} Aligned.sortedByCoord.out.bam`

java -jar ~/projects/def-idworkin/arteen/SociabilityRNA/QC/qorts/QoRTs.jar QC \
        --generatePlots \
        --maxReadLength 101 \
        ${map_dir}/${sample_name}Aligned.sortedByCoord.out.bam \
        ~/projects/def-idworkin/arteen/SociabilityRNA/index/STAR_index/dmel-all-r6.38.gtf \
        ~/projects/def-idworkin/arteen/SociabilityRNA/QC/qorts/subset/${sample_name}

done


