#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --output=detect_indels.out


module load python/2.7.18


python /home/arteen/projects/def-idworkin/arteen/SociabilityDNA/scripts/drosEU/DetectIndels.py \
--mpileup /home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/merge_run/merged_mpileup/sociability_merged.mpileup \
--minimum-count 20 \
--mask 5 \
| gzip > /home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/merge_run/filtering/InDel-positions_20.txt.gz