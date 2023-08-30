#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=34072M
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --output=multiqc.out


cd ~/projects/def-idworkin/arteen/SociabilityRNA/QC/fastqc_trim


module load python/3.9

pip install multiqc

multiqc .
