#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=snp_calling.out


bash /home/arteen/projects/def-idworkin/arteen/SociabilityDNA/scripts/PoolSNP/PoolSNP-master/PoolSNP.sh \
mpileup=/home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/merge_run/merged_mpileup/sociability_merged.mpileup \
reference=/home/arteen/projects/def-idworkin/arteen/SociabilityDNA/scripts/gatk/formatted_genomes/dmel-all-chromosome-r6.38.fasta \
names=ANC,ER_C,ER_D,ER_U \
max-cov=0.98 \
min-cov=25 \
min-count=10 \
min-freq=0.01 \
miss-frac=0.2 \
jobs=30 \
BS=1 \
output=/home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/merge_run/snp_calling
