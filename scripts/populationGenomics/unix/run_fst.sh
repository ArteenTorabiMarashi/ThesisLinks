#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=10G
#SBATCH --time=03:00:00
#SBATCH --output=out_fst.out


/home/arteen/projects/def-idworkin/arteen/SociabilityDNA/scripts/grenedalf/grenedalf-master/bin/grenedalf fst \
--window-type sliding \
--window-sliding-width 5000 \
--method unbiased-nei \
--pool-sizes 96 \
--sync-path /home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/merge_run/mpileup2sync/sync.sync \
--sample-name-list /home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/grenedalf/merged_run/sampleList.txt \
--out-dir /home/arteen/scratch/arteen/mcmaster/SociabilityDNA_outputs/merge_run/grenedalf/syncMpileup
