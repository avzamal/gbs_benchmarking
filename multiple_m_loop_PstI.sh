for i in $(seq 1 18);
do
MNUM=$(awk -v var="$i" 'NR==var' combinations.txt | awk '{print $1}')
REPEATNUM=$(awk -v var="$i" 'NR==var' combinations.txt | awk '{print $2}')
echo $MNUM
echo $REPEATNUM

mkdir "${MNUM}_${REPEATNUM}"

~/bin/bbmap/reformat.sh in1=PstI_R1.fastq.gz in2=PstI_R2.fastq.gz out1="${MNUM}_${REPEATNUM}/PstI_R1_subset.fastq.gz" out2="${MNUM}_${REPEATNUM}/PstI_R2_subset.fastq.gz" samplereadstarget=$MNUM

cd "${MNUM}_${REPEATNUM}"

cp ../PstI-stacks.tsv .

mkdir demultiplexed_PstI

/gss/home/a.zamalutdinov/bin/stacks-2.66/process_radtags -1 PstI_R1_subset.fastq.gz -2 PstI_R2_subset.fastq.gz -P -b PstI-stacks.tsv -o ./demultiplexed_PstI -q --renz-1 'pstI' --renz-2 'mspI' --threads 4 --retain-header

cd demultiplexed_PstI

cp ../../samples.txt .
cp ../../bbduk_filter_auto_for_pipelines.sh .
cp ../../bwa_to_bam_by_1_pipelines.sh .
chmod +x bbduk_filter_auto_for_pipelines.sh
chmod +x bwa_to_bam_by_1_pipelines.sh

./bbduk_filter_auto_for_pipelines.sh
./bwa_to_bam_by_1_pipelines.sh

~/bin/samtools-1.19/samtools mpileup -f /gss/home/a.zamalutdinov/DivSoy/Reference_genome/glyma.Wm82.gnm4.4PTR.genome_main.fna.gz 103.bam 10.bam 11.bam 25.bam 26.bam 36.bam 44.bam 45.bam 55.bam 60.bam 61.bam 6.bam 7.bam 8.bam 9.bam | java -jar ~/bin/VarScan.v2.4.6.jar mpileup2snp --min-coverage 2 --vcf-sample-list samples.txt --output-vcf 1 > "PstI_${MNUM}_${REPEATNUM}.vcf"

~/bin/bcftools-1.18/bcftools view -Oz -o "PstI_${MNUM}_${REPEATNUM}.vcf.gz" "PstI_${MNUM}_${REPEATNUM}.vcf"

~/bin/bcftools-1.18/bcftools index "PstI_${MNUM}_${REPEATNUM}.vcf.gz"

cd ../../

done