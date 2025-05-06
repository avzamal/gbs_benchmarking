find . -name "*_*_*_depth_5.vcf.gz" -exec cp -t all_vcf {} +


ls *.vcf.gz | awk 'BEGIN {FS="."}{print $1}' > filtering_M.txt

while IFS='' read -r ACC1 || [ -n "${ACC1}" ]; do
~/bin/bcftools-1.17/bcftools index "${ACC1}_depth_5_25_GT_ref_maf_order.vcf.gz"
comm -12 <(cat WGS_hard_ID_sorted.txt) <(sort "${ACC1}_ID.txt") | wc -l >> benchmarking_table_WGS_GBS_str.txt
wc -l "${ACC1}_ID.txt" >> benchmarking_table_calling_str.txt
~/bin/bedtools intersect -v -a "${ACC1}_depth_5_25_GT_ref_maf_order.vcf.gz" -b ~/DivSoy/Reference_genome/glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.bed > "${ACC1}_not_in_genes.txt"
wc -l "${ACC1}_not_in_genes.txt" >> benchmarking_table_not_in_genes_str.txt
~/bin/bedtools intersect -v -b "${ACC1}_depth_5_25_GT_ref_maf_order.vcf.gz" -a ~/DivSoy/Reference_genome/glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.bed > "${ACC1}_not_in_fragments.txt"
wc -l "${ACC1}_not_in_fragments.txt" >> benchmarking_table_not_in_fragments_str.txt
~/bin/bedtools intersect -v -a "${ACC1}_depth_5_25_GT_ref_maf_order.vcf.gz" -b ~/DivSoy/Reference_genome/glyma.Wm82.gnm4.ann1.T8TQ.repeatmasking.gff3 > "${ACC1}_not_in_repeats.txt"
wc -l "${ACC1}_not_in_repeats.txt" >> benchmarking_table_not_in_repeats_str.txt

done < filtering_M.txt

echo pair R_dosage > benchmarking_table_r_square_2.txt
echo pair NRD > benchmarking_table_NRD_2.txt
echo pair number TsTv > benchmarking_table_tstv_2.txt

for i in $(seq 1 63);
do
    for j in $(seq $((i+1))  63)
    do
        ACC1=$(awk -v var="$i" 'NR==var' filtering_M.txt)
        ACC2=$(awk -v var="$j" 'NR==var' filtering_M.txt)
        pair="${ACC1}_${ACC2}"
        ~/bin/bcftools-1.17/bcftools stats -s- "${ACC1}_depth_5_25_GT_ref_maf_order.vcf.gz" "${ACC2}_depth_5_25_GT_ref_maf_order.vcf.gz" > "${ACC1}_${ACC2}.txt"
        cat "${ACC1}_${ACC2}.txt" | grep -e "^GCsS" | LC_ALL=C awk -v pair=$pair '{x+=$11*$5; w+=$5; next} END{print pair,x/w}' >> benchmarking_table_r_square_2.txt
        cat "${ACC1}_${ACC2}.txt" | grep -e "^NRDs" | LC_ALL=C awk -v pair=$pair '{print pair,$3}' >> benchmarking_table_NRD_2.txt
        cat "${ACC1}_${ACC2}.txt" | grep -e "^TSTV" | tail -1 | LC_ALL=C awk -v pair=$pair '{print pair,$3+$4,$5}' >> benchmarking_table_tstv_2.txt
done
done