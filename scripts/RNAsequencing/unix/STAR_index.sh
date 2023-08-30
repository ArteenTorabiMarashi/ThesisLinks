#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=34072M
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=Star_index.out


module load  StdEnv/2020

module load star/2.7.9a

cd ~/projects/def-idworkin/arteen/SociabilityRNA/index/Star_index

STAR   --runMode genomeGenerate \
	--runThreadN 6 \
	--genomeDir ~/projects/def-idworkin/arteen/SociabilityRNA/index/Star_index \
	--genomeFastaFiles ~/projects/def-idworkin/arteen/SociabilityRNA/index/STAR_index/dmel-all-chromosome-r6.38.fasta \
	--sjdbGTFfile ~/projects/def-idworkin/arteen/SociabilityRNA/index/STAR_index/dmel-all-r6.38.gtf \
	--sjdbOverhang 99
