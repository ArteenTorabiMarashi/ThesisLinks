#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=10G
#SBATCH --time=03:00:00
#SBATCH --output=out_mpileup2sync.out


/home/arteen/projects/def-idworkin/arteen/SociabilityDNA/scripts/grenedalf/grenedalf-master/bin/grenedalf sync \
--pileup-path /home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/merge_run/merged_mpileup/sociability_merged.mpileup \
--out-dir /home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/merge_run/mpileup2sync
