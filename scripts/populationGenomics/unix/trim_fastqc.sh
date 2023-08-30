#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=44072M
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=6
#SBATCH --output=trim_fastqc.out

module load StdEnv/2020

module load fastqc/0.11.9


cd /home/arteen/projects/def-idworkin/arteen/SociabilityDNA/trimmed_reads

fastqc *.fastq.gz -o /home/arteen/projects/def-idworkin/arteen/SociabilityDNA/QC/trim_fastqc
