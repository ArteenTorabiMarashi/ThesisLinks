#!/bin/bash


#SBATCH --account=def-idworkin
#SBATCH --array=0-15
#SBATCH --mem-per-cpu=24072M
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=12
#SBATCH --output=gatk_mark_indels_%a.out

module load nixpkgs/16.09
module load gatk/3.8

#Path to input directory
in_dir=/home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/addReadGroups

#Path to output directory
out_dir=/home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/markIndels

genome=/home/arteen/projects/def-idworkin/arteen/SociabilityDNA/scripts/gatk/formatted_genomes/dmel-all-chromosome-r6.38.fasta

declare -a forward=( ${in_dir}/*.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} _RG.bam`

java -Xmx32g -jar "${EBROOTGATK}"/GenomeAnalysisTK.jar -I ${in_dir}/${base}_RG.bam \
-R ${genome} \
-T RealignerTargetCreator \
-o ${out_dir}/${base}.intervals
