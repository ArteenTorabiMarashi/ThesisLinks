#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=34072M
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=calcTin_check.out

module load python/3


virtualenv --no-download ~/scratch/arteen/mcmaster/rseqc
source ~/scratch/arteen/mcmaster/rseqc/bin/activate



cd ~/projects/def-idworkin/arteen/SociabilityRNA/QC/tin/

tin_dir=~/projects/def-idworkin/arteen/SociabilityRNA/QC/tin/

map_dir=~/projects/def-idworkin/arteen/SociabilityRNA/star/

files=(${map_dir}/*Aligned.sortedByCoord.out.bam)


for file in ${files[@]}
do
	name=${file}
	sample_name=`basename ${name} Aligned.sortedByCoord.out.bam`


	if [ ! -f "${tin_dir}/${sample_name}Aligned.sortedByCoord.out.tin.xls" ]; then
		echo "${sample_name}Aligned.sortedByCoord.out.tin.xls"
		echo "Does Not Exist, so running rseqc: "

		python3 ~/scratch/arteen/mcmaster/rseqc/bin/tin.py \
	-i ${map_dir} \
	-r ~/projects/def-idworkin/arteen/SociabilityRNA/index/dmel-all-r6.38_good.bed	
		
	
	else
		echo "Already Complete"
	fi

done