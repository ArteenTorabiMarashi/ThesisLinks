#!/bin/bash


#SBATCH --account=def-idworkin
#SBATCH --array=0-15
#SBATCH --mem-per-cpu=24072M
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=12
#SBATCH --output=picard_RG_%a.out

module load StdEnv/2020
module load picard/2.26.3


### Fix the RGSM And RGPU calls here !!!! #######

in=/home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/markdup
out=/home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/addReadGroups

declare -a forward=( ${in}/*.bam )

R1=${forward[${SLURM_ARRAY_TASK_ID}]}

name=${R1}
base=`basename ${name} .bam`

java -jar -Xmx10g $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
	INPUT=${in}/${base}.bam \
	OUTPUT=${out}/${base}_RG.bam \
	SORT_ORDER=coordinate \
	RGID=library \
	RGLB=library \
	RGPL=illumina \
	RGSM=${base} \
	RGPU=library \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT
