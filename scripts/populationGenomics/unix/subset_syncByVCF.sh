#!/bin/bash

#SBATCH --account=def-idworkin
#SBATCH --mem-per-cpu=2048M
#SBATCH --time=00:30:00
#SBATCH --output=sync_subsettedByVcf.out

## Either change your vcf & sync file names to: vcf.vcf & sync.sync or change the script to match your file names


## Split VCF by chromosome
sed -n '/2L/p' vcf.vcf > 2L.vcf
sed -n '/2R/p' vcf.vcf > 2R.vcf
sed -n '/3L/p' vcf.vcf > 3L.vcf
sed -n '/3R/p' vcf.vcf > 3R.vcf
sed -n '/X/p' vcf.vcf > X.vcf
sed '/2L/d' vcf.vcf | sed '/2R/d' | sed '/3L/d' | sed '/3R/d' | sed '/X/d' > 4.vcf

# Grab just the position IDs
awk -F'\t' '{print $2}' 2L.vcf > id_2L.txt
awk -F'\t' '{print $2}' 2R.vcf > id_2R.txt
awk -F'\t' '{print $2}' 3L.vcf > id_3L.txt
awk -F'\t' '{print $2}' 3R.vcf > id_3R.txt
awk -F'\t' '{print $2}' X.vcf > id_X.txt
awk -F'\t' '{print $2}' 4.vcf > id_4.txt



 ## Split Sync by chromosome
sed -n '/2L/p' sync.sync > 2L.sync
sed -n '/2R/p' sync.sync > 2R.sync
sed -n '/3L/p' sync.sync > 3L.sync
sed -n '/3R/p' sync.sync > 3R.sync
sed -n '/X/p' sync.sync > X.sync
sed '/2L/d' sync.sync | sed '/2R/d' | sed '/3L/d' | sed '/3R/d' | sed '/X/d' > 4.sync


# Run the awk line
awk -F'\t' 'NR==FNR{ids[$1]; next}$2 in ids' id_2L.txt 2L.sync > 2L_test.sync
awk -F'\t' 'NR==FNR{ids[$1]; next}$2 in ids' id_2R.txt 2R.sync > 2R_test.sync
awk -F'\t' 'NR==FNR{ids[$1]; next}$2 in ids' id_3L.txt 3L.sync > 3L_test.sync
awk -F'\t' 'NR==FNR{ids[$1]; next}$2 in ids' id_3R.txt 3R.sync > 3R_test.sync
awk -F'\t' 'NR==FNR{ids[$1]; next}$2 in ids' id_X.txt X.sync > X_test.sync
awk -F'\t' 'NR==FNR{ids[$1]; next}$2 in ids' id_4.txt 4.sync > 4_test.sync

# Stick em all back together
cat 2L_test.sync 2R_test.sync 3L_test.sync 3R_test.sync 4_test.sync X_test.sync > subsetted_sync.sync


# Clean out intermediate files
rm 2L.vcf 2R.vcf 3R.vcf 3L.vcf X.vcf 4.vcf
rm id_2L.txt id_2R.txt id_3L.txt id_3R.txt id_4.txt id_X.txt
rm 2L.sync 2R.sync 3L.sync 3R.sync 4.sync X.sync
rm 2L_test.sync 2R_test.sync 3L_test.sync 3R_test.sync 4_test.sync X_test.sync