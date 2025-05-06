/*
* cat stacks.tsv | awk '{print $2}' > samples.txt
*/

params.reads = "/gss/home/a.zamalutdinov/DivSoy/Benchmarking/raw_GBS_reads/demultiplexed_PstI/*.{1,2}_paired_trim.fq.gz"
params.genome = "/gss/home/a.zamalutdinov/DivSoy/Reference_genome/glyma.Wm82.gnm4"
params.outdir = "/gss/home/a.zamalutdinov/DivSoy/Benchmarking/test_nextflow/results"
params.genome_fasta_gz = "/gss/home/a.zamalutdinov/DivSoy/Reference_genome/glyma.Wm82.gnm4.4PTR.genome_main.fna.gz"
params.genome_fasta = "/gss/home/a.zamalutdinov/DivSoy/Reference_genome/glyma.Wm82.gnm4.4PTR.genome_main.fna"
params.enzyme = "PstI"
params.vcfsamplelist = "/gss/home/a.zamalutdinov/DivSoy/Benchmarking/multiple_M/samples.txt"
params.popmap = "/gss/home/a.zamalutdinov/DivSoy/Benchmarking/18M/popmap.txt"
params.sample_map = "/gss/home/a.zamalutdinov/DivSoy/Benchmarking/18M/Map.sample_map"
params.intervals = "/gss/home/a.zamalutdinov/DivSoy/Benchmarking/18M/intervals.list"
params.sample_names = "/gss/home/a.zamalutdinov/DivSoy/Benchmarking/test_nextflow/samples_bcftools_correct.txt"

log.info """\
    Reads to bam pipeline
    ===================================
    reads         : ${params.reads}
    outdir        : ${params.outdir}
    genome_index  : ${params.genome}
    genome_fasta  : ${params.genome_fasta}
    enzyme        : ${params.enzyme}
    """
/*
process bbduk {
    memory '4 GB'
    cpus 1
    time '10m'
    executor = 'pbs'
    input: 
    tuple val(sample_id), path(reads)
    output:
    //path (*)
    tuple val(sample_id), path("*.{1,2}_paired_trim.fq.gz")
    script:
    """
    /gss/home/a.zamalutdinov/bin/bbmap/bbduk.sh in="$PWD/${reads[0]}" in2="$PWD/${reads[1]}" out="${sample_id}.1_paired_trim.fq.gz" \
    outm="${sample_id}.1_unpaired_trim.fq.gz" out2="${sample_id}.2_paired_trim.fq.gz" outm2="${sample_id}.2_unpaired_trim.fq.gz" \
    k=31 ref=artifacts,phix,adapters,lambda,pjet,mtst,kapa ordered cardinality qtrim=rl trimq=20 maq=25 tbo mink=11 ktrim=r \
    minlen=50 stats="${sample_id}.txt" t=1 -Xmx4g
    """
}
*/

process strobealign {
    memory { 25.GB * task.attempt }
    time { 48.h * task.attempt }

    errorStrategy { task.exitStatus in 134..143 ? 'retry' : 'terminate' }
    maxRetries 3
    //memory '1 MB'
    cpus 1
    //time '5m'
    //executor = 'pbs'
    publishDir "results", pattern: "*.bam" , mode: 'copy'
    publishDir "results", pattern: "*.sam" , mode: 'copy'
    input: 
    tuple val(sample_id), path(reads_trimmed)
    output:
    tuple val(sample_id), path("*.bam"), path("*.bai"), path("*.sam")
    script:
    """
    module load singularity/3.6.1
    singularity run ~/bin/strobealign_0.13.0.sif strobealign --use-index -t ${task.cpus} ${params.genome_fasta_gz} "${reads_trimmed[0]}" "${reads_trimmed[1]}" > "${sample_id}.sam"
    /gss/home/a.zamalutdinov/bin/samtools-1.14/samtools fixmate -O bam "${sample_id}.sam" "${sample_id}_fixmate.bam"
    /gss/home/a.zamalutdinov/bin/samtools-1.14/samtools sort -@ ${task.cpus} -O bam -o "${sample_id}.bam"  "${sample_id}_fixmate.bam"
    rm "${sample_id}_fixmate.bam"
    /gss/home/a.zamalutdinov/bin/samtools-1.14/samtools index "${sample_id}.bam"
    """
}
/*
process bwa_mem_2 {
    memory { 25.GB * task.attempt }
    time { 1.h * task.attempt }

    errorStrategy { task.exitStatus in 134..143 ? 'retry' : 'terminate' }
    maxRetries 3
    //memory '1 MB'
    cpus 1
    //time '5m'
    executor = 'pbs'
    publishDir "results", pattern: "*.bam" , mode: 'copy'
    input: 
    tuple val(sample_id), path(reads_trimmed)
    output:
    tuple val(sample_id), path("*.bam"), path("*.bai")
    script:
    """
    
    gzip -cd "${reads_trimmed[0]}" | head -10000 > F.fq
    gzip -cd "${reads_trimmed[1]}" | head -10000 > R.fq
    /gss/home/a.zamalutdinov/bin/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t ${task.cpus} ${params.genome_fasta_gz} F.fq R.fq > "${sample_id}.sam"
    /gss/home/a.zamalutdinov/bin/samtools-1.14/samtools fixmate -O bam "${sample_id}.sam" "${sample_id}_fixmate.bam"
    /gss/home/a.zamalutdinov/bin/samtools-1.14/samtools sort -@ ${task.cpus} -O bam -o "${sample_id}.bam"  "${sample_id}_fixmate.bam"
    rm "${sample_id}_fixmate.bam"
    /gss/home/a.zamalutdinov/bin/samtools-1.14/samtools index "${sample_id}.bam"
    """
}
*/
process deepvariant {
    //container = '/gss/home/a.zamalutdinov/bin/deepvariant_1.6.0.sif'
    //singularity.enabled = true
    
    memory { 19.GB * task.attempt }
    time { 24.h * task.attempt }
    cpus { 2 * task.attempt }
    //executor = 'pbs'
    publishDir "results/${params.enzyme}/str_DeepVariant", pattern: "*_renamed.DV.g.vcf.gz" , mode: 'copy'
    input: 
    tuple val(sample_id), path(bam), path(bai), path(sam)
    output:
    path("*_renamed.DV.g.vcf.gz")

    script:
    """
    module load singularity/3.6.1
    mkdir intermediate_results_dir_${sample_id}
    singularity exec -e /gss/home/a.zamalutdinov/bin/deepvariant_1.6.0.sif bash -c '/opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=${params.genome_fasta_gz} --reads=${bam} --sample_name ${sample_id} --intermediate_results_dir intermediate_results_dir_${sample_id} --output_vcf=${sample_id}.DV.vcf.gz --output_gvcf=${sample_id}.DV.g.vcf.gz --num_shards=${task.cpus}'
    echo ${sample_id} | awk '{print \$1}' > names.txt
    /gss/home/a.zamalutdinov/bin/bcftools-1.18/bcftools reheader -s names.txt "${sample_id}.DV.g.vcf.gz" > "${sample_id}_renamed.DV.g.vcf.gz"
    """

}

process bbcall {
    memory { 20.GB * task.attempt }
    time { 24.h * task.attempt }
    cpus { 8 * task.attempt }
    executor = 'pbs'
    publishDir "results/${params.enzyme}/str_BBcall", pattern: "*_str_BBcall.vcf.gz" , mode: 'copy'
    input: 
    tuple val(sample_id), path(bam), path(bai)
    output:
    path("*_str_BBcall.vcf.gz")
    script:
    """
    /gss/home/a.zamalutdinov/bin/bbmap/callvariants.sh -Xmx20000m in="${params.outdir}/103.sam","${params.outdir}/10.sam","${params.outdir}/11.sam","${params.outdir}/25.sam","${params.outdir}/26.sam","${params.outdir}/36.sam","${params.outdir}/44.sam","${params.outdir}/45.sam","${params.outdir}/55.sam","${params.outdir}/60.sam","${params.outdir}/6.sam","${params.outdir}/7.sam","${params.outdir}/8.sam","${params.outdir}/9.sam","${params.outdir}/61.sam" out="${params.enzyme}_str_BBcall.vcf.gz" ref=${params.genome_fasta_gz} multisample=t ploidy=2 threads=${task.cpus} maf=0.05 border=0
    """
}

process freebayes {
    memory { 20.GB * task.attempt }
    time { 48.h * task.attempt }
    cpus 1
    executor = 'pbs'
    publishDir "results/${params.enzyme}/str_freebayes", pattern: "*_str_freebayes.vcf.gz" , mode: 'copy'
    input: 
    tuple val(sample_id), path(bam), path(bai)
    output:
    path("*_str_freebayes.vcf.gz")
    script:
    """
    /gss/home/a.zamalutdinov/bin/bamaddrg-master/bamaddrg -b "${params.outdir}/103.bam" -s 103 -b "${params.outdir}/10.bam" -s 10 -b "${params.outdir}/11.bam" -s 11 -b "${params.outdir}/25.bam" -s 25 -b "${params.outdir}/26.bam" -s 26 -b "${params.outdir}/36.bam" -s 36 -b "${params.outdir}/44.bam" -s 44 -b "${params.outdir}/45.bam" -s 45 -b "${params.outdir}/55.bam" -s 55 -b "${params.outdir}/60.bam" -s 60 -b "${params.outdir}/61.bam" -s 61 -b "${params.outdir}/6.bam" -s 6 -b "${params.outdir}/7.bam" -s 7 -b "${params.outdir}/8.bam" -s 8 -b "${params.outdir}/9.bam" -s 9 | /gss/home/a.zamalutdinov/bin/freebayes-1.3.6 -f ${params.genome_fasta} -c > "${params.enzyme}_str_freebayes.vcf"
    /gss/home/a.zamalutdinov/bin/bcftools-1.18/bcftools view -Oz -o "${params.enzyme}_str_freebayes.vcf.gz" "${params.enzyme}_str_freebayes.vcf"
    """
}

process bcftools {
    memory { 20.GB * task.attempt }
    time { 24.h * task.attempt }
    cpus {4 * task.attempt }
    executor = 'pbs'
    publishDir "results/${params.enzyme}/str_bcftools", pattern: "*_str_bcftools.vcf.gz" , mode: 'copy'
    input: 
    tuple val(sample_id), path(bam), path(bai)
    output:
    path("*_str_bcftools.vcf.gz")
    script:
    """
    /gss/home/a.zamalutdinov/bin/bcftools-1.18/bcftools mpileup -Ou --threads ${task.cpus} -a AD,DP,ADF,ADR,SP,SCR -f ${params.genome_fasta_gz} ${params.outdir}/*.bam | /gss/home/a.zamalutdinov/bin/bcftools-1.18/bcftools call --threads ${task.cpus} -vmO z -o "${params.enzyme}_str_bcftools.vcf.gz"
    mv "${params.enzyme}_str_bcftools.vcf.gz" "${params.enzyme}_str_bcftools_before_fix.vcf.gz"
    ~/bin/bcftools-1.18/bcftools reheader -s "${params.sample_names}" "${params.enzyme}_str_bcftools_before_fix.vcf.gz" | ~/bin/bcftools-1.18/bcftools view -Oz -o "${params.enzyme}_str_bcftools.vcf.gz"
    """
}

process GATK {
    memory { 10.GB * task.attempt }
    time { 48.h * task.attempt }
    cpus 1
    //executor = 'pbs'
    publishDir "results/${params.enzyme}/str_GATK", pattern: "*.g.vcf.gz" , mode: 'copy'
    publishDir "results/${params.enzyme}/str_GATK", pattern: "*.g.vcf.gz.tbi" , mode: 'copy'
    input: 
    tuple val(sample_id), path(bam), path(bai), path(sam)
    output:
    tuple path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi")
    script:
    """
    module load singularity/3.6.1
    singularity exec -e /gss/home/a.zamalutdinov/bin/gatk_4.5.0.0.sif bash -c 'gatk AddOrReplaceReadGroups I=${sample_id}.bam O=${sample_id}.group.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20'
    /gss/home/a.zamalutdinov/bin/samtools-1.19/samtools index "${sample_id}.group.bam"
    singularity exec -e /gss/home/a.zamalutdinov/bin/gatk_4.5.0.0.sif bash -c 'gatk HaplotypeCaller -R ${params.genome_fasta_gz} -I ${sample_id}.group.bam -O ${sample_id}.g.vcf.gz -ERC GVCF --native-pair-hmm-threads ${task.cpus}'
    """
}
process varscan {
    memory { 10.GB * task.attempt }
    time { 24.h * task.attempt }
    cpus 1
    executor = 'pbs'
    publishDir "results/${params.enzyme}/str_VarScan", pattern: "*_str_VarScan.vcf.gz" , mode: 'copy'
    input: 
    tuple val(sample_id), path(bam), path(bai)
    output:
    path("*_str_VarScan.vcf.gz")
    script:
    """
    /gss/home/a.zamalutdinov/bin/samtools-1.19/samtools mpileup -f ${params.genome_fasta_gz} "${params.outdir}/103.bam" "${params.outdir}/10.bam" "${params.outdir}/11.bam" "${params.outdir}/25.bam" "${params.outdir}/26.bam" "${params.outdir}/36.bam" "${params.outdir}/44.bam" "${params.outdir}/45.bam" "${params.outdir}/55.bam" "${params.outdir}/60.bam" "${params.outdir}/61.bam" "${params.outdir}/6.bam" "${params.outdir}/7.bam" "${params.outdir}/8.bam" "${params.outdir}/9.bam" | java -jar /gss/home/a.zamalutdinov/bin/VarScan.v2.4.6.jar mpileup2snp --min-coverage 2 --vcf-sample-list "${params.vcfsamplelist}" --output-vcf 1 > "${params.enzyme}_str_VarScan.vcf"
    /gss/home/a.zamalutdinov/bin/bcftools-1.18/bcftools view -Oz -o "${params.enzyme}_str_VarScan.vcf.gz" "${params.enzyme}_str_VarScan.vcf"
    """
}

process stacks {
    memory { 20.GB * task.attempt }
    time { 24.h * task.attempt }
    cpus { 4 * task.attempt }
    //executor = 'pbs'
    publishDir "results/${params.enzyme}/str_stacks", pattern: "*_str_stacks.vcf.gz" , mode: 'copy'
    input: 
    tuple val(sample_id), path(bam), path(bai)
    output:
    path("*_str_stacks.vcf.gz")
    script:
    """
    /gss/home/a.zamalutdinov/bin/stacks-2.64/ref_map.pl --samples "${params.outdir}" --popmap ${params.popmap} -o ./ -T ${task.cpus}
    /gss/home/a.zamalutdinov/bin/stacks-2.64/populations -P ./ -O ./ --vcf
    /gss/home/a.zamalutdinov/bin/bcftools-1.18/bcftools sort -Oz -o "${params.enzyme}_str_stacks.vcf.gz" populations.snps.vcf
    """
}
process glnexus {
    memory { 20.GB * task.attempt }
    time { 4.h * task.attempt }
    cpus 8
    executor = 'pbs'
    publishDir "results/${params.enzyme}/str_DeepVariant", pattern: "*_str_DeepVariant.vcf.gz" , mode: 'copy'
    input: 
    path(reads)
    output:
    path("*_str_DeepVariant.vcf.gz")
    script:
    """
    module load singularity/3.6.1
    singularity exec -e /gss/home/a.zamalutdinov/bin/glnexus_v1.4.3.sif bash -c 'glnexus_cli --config DeepVariant -t ${task.cpus} -m 20 ${params.outdir}/${params.enzyme}/str_DeepVariant/*_renamed.DV.g.vcf.gz' | /gss/home/a.zamalutdinov/bin/bcftools-1.18/bcftools view -Oz --threads ${task.cpus} > "${params.enzyme}_str_DeepVariant.vcf.gz"
    """
}

process gatk_db {
    memory { 10.GB * task.attempt }
    time { 48.h * task.attempt }
    cpus 1
    executor = 'pbs'
    publishDir "results/${params.enzyme}/str_GATK", pattern: "*_str_GATK.vcf.gz" , mode: 'copy'
    input: 
    path(reads)
    output:
    path("*_str_GATK.vcf.gz")
    script:
    """
    module load singularity/3.6.1
    singularity run /gss/home/a.zamalutdinov/bin/gatk_4.5.0.0.sif bash -c 'gatk GenomicsDBImport --genomicsdb-workspace-path my_database --sample-name-map ${params.sample_map} --reader-threads ${task.cpus} -L ${params.intervals}'
    singularity run /gss/home/a.zamalutdinov/bin/gatk_4.5.0.0.sif bash -c 'gatk GenotypeGVCFs -R ${params.genome_fasta_gz} -V gendb://my_database -O ${params.enzyme}_str_GATK.vcf.gz'

    """
}

workflow {
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    strobealign(read_pairs_ch)
    //bwa_mem_2(read_pairs_ch)
    pool_ch = strobealign.out.collect()
    //pool_ch.view()
    deepvariant(strobealign.out)
    GATK(strobealign.out)
    bbcall(pool_ch)
    freebayes(pool_ch)
    bcftools(pool_ch)
    varscan(pool_ch)
    stacks(pool_ch)
    deepvariant_ch = deepvariant.out.collect()
    glnexus(deepvariant_ch)
    gatk_ch = GATK.out.collect()
    gatk_db(gatk_ch)
}



/*
 /Users/avzamal/Tools/BBMap_39.06/bbduk.sh in="/Users/avzamal/Downloads/test_nextflow/HindIII_R1_subset.fastq" in2="HindIII_R2_subset.fastq" out="Hind.1_paired_trim.fq.gz" \
    outm="Hind.1_unpaired_trim.fq.gz" out2="Hind.2_paired_trim.fq.gz" outm2="Hind.2_unpaired_trim.fq.gz" \
    k=31 ref=artifacts,phix,adapters,lambda,pjet,mtst,kapa ordered cardinality qtrim=rl trimq=20 maq=25 tbo mink=11 ktrim=r \
    minlen=50 stats="Hind.txt" t=1 -Xmx4g
singularity run ~/bin/strobealign_0.13.0.sif strobealign --use-index -t ${task.cpus} ${params.genome_fasta_gz} "${reads_trimmed[0]}" "${reads_trimmed[1]}" > "${sample_id}.sam"

/gss/home/a.zamalutdinov/bin/bwa/bwa mem -t 6 /gss/home/a.zamalutdinov/DivSoy/Reference_genome/glyma.Wm82.gnm4 "${ITEM}.1_paired_trim.fq.gz" "${ITEM}.2_paired_trim.fq.gz" > "${ITEM}.sam"

/gss/home/a.zamalutdinov/bin/samtools-1.14/samtools fixmate -O bam "${ITEM}.sam" "${ITEM}_fixmate.bam"
/gss/home/a.zamalutdinov/bin/samtools-1.14/samtools sort -@ 6 -O bam -o "${ITEM}.bam"  "${ITEM}_fixmate.bam"
/gss/home/a.zamalutdinov/bin/samtools-1.14/samtools index "${ITEM}.bam"

while IFS='' read -r ITEM || [ -n "${ITEM}" ]; do
mkdir "intermediate_results_dir_${ITEM}"
mkdir intermediate_results_dir_${sample_id}

COMMAND_2=$(echo "/opt/deepvariant/bin/run_deepvariant --model_type=WES --ref=/gss/home/a.zamalutdinov/Guar/Reference_genome/GCA_037177725.1_Cte_V1.0_genomic.fna.gz --reads=$ITEM.bam --intermediate_results_dir intermediate_results_dir_$ITEM --output_vcf=$ITEM.DV.vcf.gz --output_gvcf=$ITEM.DV.g.vcf.gz --num_shards=4")

singularity run /gss/home/a.zamalutdinov/bin/deepvariant_1.6.0.sif $COMMAND_2 &>"${ITEM}.txt"
done < samples.txt
singularity exec -e /gss/home/a.zamalutdinov/bin/deepvariant_1.6.0.sif bash -c '/opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=${params.genome}.fna.gz --reads=${bam} --intermediate_results_dir intermediate_results_dir_${sample_id} --output_vcf=${sample_id}.DV.vcf.gz --output_gvcf=${sample_id}.DV.g.vcf.gz --num_shards=1'
   
*/