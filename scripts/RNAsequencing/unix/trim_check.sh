#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=54072M
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --output=trim_check.out

module load nixpkgs/16.09

module load trimmomatic/0.36

raw_dir=~/scratch/arteen/mcmaster/rawData/SociabilityAS_RNA

trim_dir=~/projects/def-idworkin/arteen/SociabilityRNA/trimmed_reads

adapter=~/projects/def-idworkin/arteen/SociabilityRNA/trimmed_reads/TruSeq3-PE-2.fa

files=(${raw_dir}/*_R1.fastq.gz)


for file in ${files[@]}
do
	name=${file}
	base=`basename ${name} _R1.fastq.gz`
	rename=$(awk -v X=$base -F',' '{ if ($1 == X) { print $2} }' ~/projects/def-idworkin/arteen/SociabilityRNA/trimmed_reads/naming.csv)

	if [ ! -f "${trim_dir}/${rename}_R1_PE.fastq.gz" ]; then
		echo "${rename}_R1_PE.fastq.gz"
		echo "Does Not Exist, so trimming: "

		java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 \
  			${raw_dir}/${base}_R1.fastq.gz \
			${raw_dir}/${base}_R2.fastq.gz \
			${trim_dir}/${rename}_R1_PE.fastq.gz \
			${trim_dir}/${rename}_R1_SE.fastq.gz \
			${trim_dir}/${rename}_R2_PE.fastq.gz \
			${trim_dir}/${rename}_R2_SE.fastq.gz \
 			ILLUMINACLIP:${adapter}:2:30:10 LEADING:3 TRAILING:3 MAXINFO:20:0.2 MINLEN:36
	
	else
		echo "Already Complete"
	fi

done
