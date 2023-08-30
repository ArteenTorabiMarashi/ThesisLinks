## Either change your vcf & bed file names to: vcf.vcf & bed.bed or change the script to match your file names

## Make sure the headers are off the VCF, cat it back on at the end

 ## Split bed by chromosome
sed -n '/2L/p' bed.bed > 2L.bed
sed -n '/2R/p' bed.bed > 2R.bed
sed -n '/3L/p' bed.bed > 3L.bed
sed -n '/3R/p' bed.bed > 3R.bed
sed -n '/X/p' bed.bed > X.bed
sed '/2L/d' bed.bed | sed '/2R/d' | sed '/3L/d' | sed '/3R/d' | sed '/X/d' > 4.bed

# Grab just the position IDs
awk -F'\t' '{print $2}' 2L.bed > id_2L.txt
awk -F'\t' '{print $2}' 2R.bed > id_2R.txt
awk -F'\t' '{print $2}' 3L.bed > id_3L.txt
awk -F'\t' '{print $2}' 3R.bed > id_3R.txt
awk -F'\t' '{print $2}' X.bed > id_X.txt
awk -F'\t' '{print $2}' 4.bed > id_4.txt

## Split VCF by chromosome
sed -n '/2L/p' vcf.vcf > 2L.vcf
sed -n '/2R/p' vcf.vcf > 2R.vcf
sed -n '/3L/p' vcf.vcf > 3L.vcf
sed -n '/3R/p' vcf.vcf > 3R.vcf
sed -n '/X/p' vcf.vcf > X.vcf
sed '/2L/d' vcf.vcf | sed '/2R/d' | sed '/3L/d' | sed '/3R/d' | sed '/X/d' > 4.vcf



# Run the awk line
awk -F'\t' 'NR==FNR{ids[$1]; next}$2 in ids' id_2L.txt 2L.vcf > 2L_test.vcf
awk -F'\t' 'NR==FNR{ids[$1]; next}$2 in ids' id_2R.txt 2R.vcf > 2R_test.vcf
awk -F'\t' 'NR==FNR{ids[$1]; next}$2 in ids' id_3L.txt 3L.vcf > 3L_test.vcf
awk -F'\t' 'NR==FNR{ids[$1]; next}$2 in ids' id_3R.txt 3R.vcf > 3R_test.vcf
awk -F'\t' 'NR==FNR{ids[$1]; next}$2 in ids' id_X.txt X.vcf > X_test.vcf
awk -F'\t' 'NR==FNR{ids[$1]; next}$2 in ids' id_4.txt 4.vcf > 4_test.vcf

# Stick em all back together
cat 2L_test.vcf 2R_test.vcf 3L_test.vcf 3R_test.vcf 4_test.vcf X_test.vcf > subsetted_vcf.vcf


# Clean out intermediate files
rm 2L.vcf 2R.vcf 3R.vcf 3L.vcf X.vcf 4.vcf
rm id_2L.txt id_2R.txt id_3L.txt id_3R.txt id_4.txt id_X.txt
rm 2L.bed 2R.bed 3L.bed 3R.bed 4.bed X.bed
rm 2L_test.vcf 2R_test.vcf 3L_test.vcf 3R_test.vcf 4_test.vcf X_test.vcf

## Dont forget to add your header back to the VCF!
