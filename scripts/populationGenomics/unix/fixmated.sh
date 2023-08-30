#!/bin/bash


#SBATCH --account=def-idworkin
#SBATCH --array=0-15
#SBATCH --mem-per-cpu=24072M
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=12
#SBATCH --output=fixmate_%a.out

module load StdEnv/2020
module load samtools/1.15

in=/home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/sorted_byName
out=/home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/fixmate

declare -a forward=( ${in}/*.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .bam`

samtools fixmate -m -u -@ 12 -O bam ${in}/${base}.bam ${out}/${base}_fixmated.bam
