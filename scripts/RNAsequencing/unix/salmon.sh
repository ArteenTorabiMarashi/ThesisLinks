#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=64072M
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=salmon.out


module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3

module load salmon/1.4.0



index_dir=~/projects/def-idworkin/arteen/SociabilityRNA/index/salmon/salmon_index

trim_dir=~/projects/def-idworkin/arteen/SociabilityRNA/trimmed_reads

counts_dir=~/projects/def-idworkin/arteen/SociabilityRNA/salmon



files=(${trim_dir}/*_R1_PE.fastq.gz)


for file in ${files[@]}
do
	name=${file}
	sample_name=`basename ${name} _R1_PE.fastq.gz`
	

	if [ ! -f "${trim_dir}/${rename}_quant" ]; then
		echo "${rename}_quant"
		echo "Does Not Exist, so running salmon: "

			
		salmon quant -i ${index_dir} -l A \
			-1 ${trim_dir}/${sample_name}_R1_PE.fastq.gz \
			-2 ${trim_dir}/${sample_name}_R2_PE.fastq.gz \
			-p 6 --validateMappings --rangeFactorizationBins 4 \
			--seqBias --gcBias \
			-o ${counts_dir}/${sample_name}_quant

	
	else
		echo "Already Complete"
	fi

done
