#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=35G
#SBATCH --time=72:00:00
#SBATCH --output=mpileup.out



module load StdEnv/2020
module load samtools/1.15

samtools mpileup -Q 20 -q 20 -d 1200 \
-f /home/arteen/projects/def-idworkin/arteen/SociabilityDNA/scripts/gatk/formatted_genomes/dmel-all-chromosome-r6.38.fasta \
/home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/merge_run/merged_bam/*.bam \
-o /home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/merge_run/merged_mpileup/sociability_merged.mpileup
