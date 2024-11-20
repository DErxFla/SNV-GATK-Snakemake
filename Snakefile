# Snakefile

# Used directories
REF_GENOME = "Data/genome.fa"  # is already indexed
#genome_dict = "Data/genome.dict" # gatk CreateSequenceDictionary -R Data/genome.fa -O Data/genome.dict
known_sites = "Data/Homo_sapiens_assembly38.dbsnp138.vcf"
FastQ_Dir = "Data/samples"
fastqc_dir = "Data/fastqc"
aligned_reads = "Data/aligned_bam"
variants = "Data/variants"
SAMPLES = ["A", "B", "C"]

rule setup_directories:
    output:
        touch("Data/.setup_done")
    shell:
        """
        mkdir -p {fastqc_dir} {aligned_reads} {variants}
        touch Data/.setup_done
        """

rule all:
    input:
        "Data/.setup_done",
        expand(fastqc_dir + "/{sample}_fastqc.html", sample=SAMPLES),
        expand(fastqc_dir + "/{sample}_fastqc.zip", sample=SAMPLES),
        expand(aligned_reads + "/{sample}_sorted.bam", sample=SAMPLES),
        expand(aligned_reads + "/{sample}_marked_sorted.bam", sample=SAMPLES),
        expand(aligned_reads + "/{sample}_marked_sorted.bai", sample=SAMPLES),
        expand(aligned_reads + "/{sample}.recal.table", sample=SAMPLES),
        expand(aligned_reads + "/{sample}_recal_bqsr.bam", sample=SAMPLES)


rule fastqc:
    input:
        fastq = FastQ_Dir + "/{sample}.fastq"
    output:
        html = fastqc_dir + "/{sample}_fastqc.html",
        zip = fastqc_dir + "/{sample}_fastqc.zip"
    params:
        output_dir = fastqc_dir
    shell:
        """
        mkdir -p {params.output_dir}
        echo 'Run FastQC - QC'
        fastqc {input.fastq} -o {params.output_dir}
        """

# Reference genome should already be indexed (bwa index)
# Rule to map reads to the reference genome

rule map_reads:
    input:
        fastq = FastQ_Dir + "/{sample}.fastq",
        reference = REF_GENOME
    output:
        bam = aligned_reads + "/{sample}_aligned.bam"
    params:
        threads = 4
    log:
        aligned_reads + "/{sample}.map.log"
    shell:
        """
        bwa mem \
        -t {params.threads} \
        -R "@RG\\tID:{wildcards.sample}\\tPL:ILLUMINA\\tSM:{wildcards.sample}" \
        {input.reference} {input.fastq} 2> {log} | samtools view -bS - > {output.bam}
        """

# Rule to sort BAM files after alignment

rule sort_bam:
    input:
        unsorted_bam = aligned_reads + "/{sample}_aligned.bam"
    output:
        sorted_bam = aligned_reads + "/{sample}_sorted.bam"
    log:
        aligned_reads + "/{sample}.sort.log"
    shell:
        """
        samtools sort -o {output.sorted_bam} {input.unsorted_bam} 2> {log}
        """

# Rule to mark duplicates in sorted BAM files

rule mark_duplicates:
    input:
        sorted_bam = aligned_reads + "/{sample}_sorted.bam"
    output:
        marked_bam = aligned_reads + "/{sample}_marked_sorted.bam",
        index = aligned_reads + "/{sample}_marked_sorted.bai"
    params:
        memory = "4g"
    log:
        aligned_reads + "/{sample}.mark.log"
    shell:
        """
        gatk --java-options "-Xmx{params.memory}" MarkDuplicates \
        -I {input.sorted_bam} \
        -O {output.marked_bam} \
        -M Data/aligned_bam/{wildcards.sample}.metrics.txt \
        -CREATE_INDEX true \
        -VALIDATION_STRINGENCY LENIENT \
        -MAX_RECORDS_IN_RAM 1000000 2> {log}
        """
# base recalibration (optional)
rule base_recalibration:
    input:
        marked_bam = aligned_reads + "/{sample}_marked_sorted.bam",
        known = known_sites,
        reference = REF_GENOME
    output:
        recal = aligned_reads + "/{sample}.recal.table",
        recal_bam = aligned_reads + "/{sample}_recal_bqsr.bam"
    log:
        aligned_reads + "/{sample}.recal.log"
    shell:
        """
        set -e
        echo 'Running BaseRecalibrator for {sample}'
        gatk --java-options "-Xmx{params.memory}" BaseRecalibrator \
        -I {input.marked_bam} \
        -R {input.reference} \
        --known-sites {input.known} \
        -O {output.recal} 2>> {log}
        
        echo 'Applying BQSR for {sample}'
        gatk --java-options "-Xmx{params.memory}" ApplyBQSR \
        -I {input.marked_bam} \
        -R {input.reference} \
        --bqsr-recal-file {output.recal} \
        -O {output.recal_bam} 2>> {log}
        """

