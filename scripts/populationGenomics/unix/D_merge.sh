#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=15G
#SBATCH --time=03:00:00
#SBATCH --output=D_merge.out




in_dir=/home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/realignIndels/
out_dir=/home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/merge_run/merged_bam/

module load StdEnv/2020
module load samtools/1.15

samtools merge -@15 -o ${out_dir}/ER_D_merged.bam \
 ${in_dir}/ER_D_L1.bam ${in_dir}/ER_D_L2.bam \
 ${in_dir}/ER_D_L3.bam ${in_dir}/ER_D_L4.bam 
