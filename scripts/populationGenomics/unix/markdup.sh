#!/bin/bash


#SBATCH --account=def-idworkin
#SBATCH --array=0-15
#SBATCH --mem-per-cpu=24072M
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=12
#SBATCH --output=markdup_%a.out

module load StdEnv/2020
module load samtools/1.15

in=/home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/sorted_byCoord
out=/home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/markdup

declare -a forward=( ${in}/*_namesorted_fixmated_coordsorted.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} _namesorted_fixmated_coordsorted.bam`

samtools markdup -l 150 -r -s -f stats.txt -d 2500 -@ 12 ${in}/${base}_namesorted_fixmated_coordsorted.bam ${out}/${base}.bam
