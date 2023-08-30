#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=44072M
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=6
#SBATCH --output=bam_fastqc.out

module load StdEnv/2020

module load fastqc/0.11.9


cd /home/arteen/projects/def-idworkin/arteen/SociabilityDNA/bwa_map/bam_files/filter30

fastqc *.bam -o /home/arteen/projects/def-idworkin/arteen/SociabilityDNA/QC/bam_fastqc
