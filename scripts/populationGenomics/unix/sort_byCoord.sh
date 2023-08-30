#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --array=0-15
#SBATCH --mem-per-cpu=24072M
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=12
#SBATCH --output=sort_byCoord_%a.out

module load StdEnv/2020
module load samtools/1.15

in=/home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/fixmate
out=/home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/sorted_byCoord

declare -a forward=( ${in}/*.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .bam`

samtools sort -@ 12 -o ${out}/${base}_coordsorted.bam ${in}/${base}.bam
