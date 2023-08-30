#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=64072M
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=fastqc_DNA.out

module load StdEnv/2020

module load fastqc/0.11.9


cd /home/arteen/projects/def-idworkin/arteen/rawReads/SociabilityAS_DNA


fastqc *.fastq.gz -o /home/arteen/projects/def-idworkin/arteen/SociabilityDNA/QC/fastqc