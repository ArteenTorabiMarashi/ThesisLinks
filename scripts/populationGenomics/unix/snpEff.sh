#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH -A def-idworkin
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=out_spnEff.out



java -Xmx8g /home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/merge_run/snpEff/snpEff/snpEff.jar \
-ud -download Drosophila_melanogaster /home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/merge_run/snpEff/snpEff/CVD_highFst_lowCMH.bed > snp_eff.txt