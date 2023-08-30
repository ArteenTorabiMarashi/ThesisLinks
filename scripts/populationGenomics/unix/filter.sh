#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --output=out_filter.out


module load python/2.7.18

python2.7 /home/arteen/projects/def-idworkin/arteen/SociabilityDNA/scripts/drosEU/FilterPosFromVCF.py \
--indel /home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/merge_run/filtering/InDel-positions_20.txt.gz \
--te /home/arteen/projects/def-idworkin/arteen/SociabilityDNA/scripts/gatk/formatted_genomes/dmel-all-chromosome-r6.38.fasta.out.gff \
--vcf /home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/merge_run/snp_calling/snp_calling.vcf \
> /home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/merge_run/filtering/SNPS_filtered.vcf
