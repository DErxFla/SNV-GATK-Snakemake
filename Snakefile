# Snakefile

# Used directories
REF_GENOME = "Data/genome.fa"  # is already indexed
FastQ_Dir = "Data/samples"
fastqc_dir = "Data/fastqc"
aligned_reads = "Data/aligned_bam"
SAMPLES = ["A", "B", "C"]

rule setup_directories:
    output:
        touch("Data/.setup_done")
    shell:
        """
        mkdir -p {fastqc_dir} {aligned_reads}
        touch Data/.setup_done
        """

rule all:
    input:
        "Data/.setup_done",
        expand(fastqc_dir + "/{sample}_fastqc.html", sample=SAMPLES),
        expand(fastqc_dir + "/{sample}_fastqc.zip", sample=SAMPLES),
        expand(aligned_reads + "/{sample}_sorted.bam", sample=SAMPLES),
        expand(aligned_reads + "/{sample}_marked_sorted.bam", sample=SAMPLES),
        expand(aligned_reads + "/{sample}_marked_sorted.bai", sample=SAMPLES)

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
        mkdir -p {aligned_reads}
        echo 'Run BWA - Mapping'
        bwa mem -t {params.threads} -R "@RG\\tID:{wildcards.sample}\\tPL:ILLUMINA\\tSM:{wildcards.sample}" \
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