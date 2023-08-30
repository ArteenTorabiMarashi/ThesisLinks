#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=44072M
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=star.out


module load  StdEnv/2020

module load star/2.7.9a

cd ~/projects/def-idworkin/arteen/SociabilityRNA/star

genome=~/projects/def-idworkin/arteen/SociabilityRNA/index/Star_index
trim_dir=~/projects/def-idworkin/arteen/SociabilityRNA/trimmed_reads
mapped_dir=~/projects/def-idworkin/arteen/SociabilityRNA/star

STAR --genomeLoad LoadAndExit \
	--genomeDir ${genome}


files=(${trim_dir}/*_R1_PE.fastq.gz)


for file in ${files[@]}
do
	name=${file}
	sample_name=`basename ${name} _R1_PE.fastq.gz`
	

	if [ ! -f "${mapped_dir}/${sample_name}Log.final.out" ]; then
		echo "${sample_name}Log.final.out"
		echo "Does Not Exist, so running star: "

			
		STAR --runThreadN 16 \
			--quantMode TranscriptomeSAM GeneCounts \
			--genomeDir ${genome} \
			--readFilesIn ${trim_dir}/${sample_name}_R1_PE.fastq.gz \
			${trim_dir}/${sample_name}_R2_PE.fastq.gz \
			--readFilesCommand zcat \
			--outFileNamePrefix ${mapped_dir}/${sample_name} \
			--outSAMtype BAM SortedByCoordinate

	
	else
		echo "Already Complete"
	fi

done

STAR --genomeLoad Remove \
	--genomeDir ${genome}
