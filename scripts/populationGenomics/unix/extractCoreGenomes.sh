#!/bin/bash


#SBATCH --account=def-idworkin
#SBATCH --array=0-15
#SBATCH --mem-per-cpu=24072M
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=12
#SBATCH --output=extractCoreGenome_%a.out


module load samtools/1.15

in=/home/arteen/projects/def-idworkin/arteen/SociabilityDNA/bwa_map/bam_files/filter30
out=/home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/coreGenomeExtracted

declare -a forward=( ${in}/*.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .bam`

samtools view -h -@ 12 ${in}/${base}.bam | \
awk '{if ($3 == "2L" || $3 == "2R" || $3 == "3L" || $3 == "3R" || $3 == "4" || $3 == "X" || \
$2 == "SN:2L" || $2 == "SN:2R" || $2 == "SN:3L" || $2 == "SN:3R" || $2 == "SN:4" || $2 == "SN:X") {print $0}}' | \
samtools view -b -@ 12 -o ${out}/${base}.bam
