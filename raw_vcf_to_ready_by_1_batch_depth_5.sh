#!/bin/bash
#PBS -N vcf_to_ready
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -d .

LEFT_EDGE=$(expr $PBS_ARRAYID \* 1 + 1)
RIGHT_EDGE=$(expr $LEFT_EDGE + 1) 
COMMAND_1="${LEFT_EDGE}"p
sed -n $COMMAND_1 filtering_M.txt > "part_vcf_filtering_${PBS_ARRAYID}.txt"

while IFS='' read -r ITEM || [ -n "${ITEM}" ]; do

~/bin/bcftools-1.18/bcftools index "${ITEM}.vcf.gz"

~/bin/bcftools-1.18/bcftools view -m2 -M2 -v snps -Oz -o "${ITEM}_bi.vcf.gz" "${ITEM}.vcf.gz"

~/bin/bcftools-1.18/bcftools annotate --set-id '%CHROM\_%POS' -Oz -o "${ITEM}_bi_id.vcf.gz" "${ITEM}_bi.vcf.gz"

~/bin/bcftools-1.18/bcftools sort -Oz -o "${ITEM}_bi_id_sort.vcf.gz" "${ITEM}_bi_id.vcf.gz"

~/bin/bcftools-1.18/bcftools index "${ITEM}_bi_id_sort.vcf.gz"

~/bin/bcftools-1.18/bcftools view --include 'AVG(FORMAT/DP)>5' -Oz -o "${ITEM}_depth_5.vcf.gz" "${ITEM}_bi_id_sort.vcf.gz"

~/bin/bcftools-1.18/bcftools +/gss/home/a.zamalutdinov/bin/bcftools-1.18/plugins/setGT.so "${ITEM}_depth_5.vcf.gz" -- -t q -n . -i 'FMT/DP<3' > "${ITEM}_depth_5_GT.vcf"

~/bin/bcftools-1.18/bcftools norm --check-ref ws --fasta-ref /gss/home/a.zamalutdinov/DivSoy/Reference_genome/glyma.Wm82.gnm4.4PTR.genome_main.fna.gz -Oz -o "${ITEM}_depth_5_GT_ref.vcf.gz" "${ITEM}_depth_5_GT.vcf"

~/bin/bcftools-1.18/bcftools view -i'F_MISSING<0.4' --min-af 0.05 -Oz -o "${ITEM}_depth_5_GT_ref_maf.vcf.gz" "${ITEM}_depth_5_GT_ref.vcf.gz"

~/bin/bcftools-1.18/bcftools index "${ITEM}_depth_5_GT_ref_maf.vcf.gz"

~/bin/bcftools-1.18/bcftools view --regions-file contig_lengths.txt -S samples.txt "${ITEM}_depth_5_GT_ref_maf.vcf.gz" -Oz -o "${ITEM}_depth_5_25_GT_ref_maf_order.vcf.gz" --threads 1

~/bin/bcftools-1.18/bcftools index "${ITEM}_depth_5_25_GT_ref_maf_order.vcf.gz"

~/bin/bcftools-1.18/bcftools query -f '%ID\n' "${ITEM}_depth_5_25_GT_ref_maf_order.vcf.gz" > "${ITEM}_ID.txt"

~/bin/bcftools-1.17/bcftools stats -s- "${ITEM}_depth_5_25_GT_ref_maf_order.vcf.gz" > "${ITEM}_stats.txt"

done < "part_vcf_filtering_${PBS_ARRAYID}.txt"
