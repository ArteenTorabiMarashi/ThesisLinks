#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=10G
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=10
#SBATCH --output=repeat_masker.out


module load StdEnv/2020
module load gcc/9.3.0
module load repeatmasker/4.1.1

RepeatMasker -pa 10 \
-lib /home/arteen/projects/def-idworkin/arteen/SociabilityDNA/scripts/gatk/formatted_genomes/dmel-transposons_idFixed-r6.38.fasta \
-gff /home/arteen/projects/def-idworkin/arteen/SociabilityDNA/scripts/gatk/formatted_genomes/dmel-all-chromosome-r6.38.fasta
