#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=64072M
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=6
#SBATCH --output=bwa_index.out


module load StdEnv/2020
module load bwa/0.7.17


bwa index /home/arteen/projects/def-idworkin/arteen/SociabilityRNA/index/salmon/dmel-all-chromosome-r6.38.fasta
