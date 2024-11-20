# Snakefile

configfile: "config.yaml"

# Get variables from config file
# samples
samples = config["samples"]
# output directories 
FastQ_Dir = config["folders"]["FastQ_Dir"]
fastqc_dir = config["folders"]["fastqc_dir"]
aligned_reads = config["folders"]["aligned_reads"]
variants = config["folders"]["variants"]
# reference genome 38
reference_genome = config["input_files"]["reference_genome"]

# here is defined, what the final output should look like. 
# snakemake checks if all the outputs are in the folders they should be
# if something is missing --> error
rule all:
    input:
        "Data/.setup_done",
        expand(fastqc_dir + "/{sample}_fastqc.html", sample=samples),
        expand(fastqc_dir + "/{sample}_fastqc.zip", sample=samples),
        expand(aligned_reads + "/{sample}_sorted.bam", sample=samples),
        expand(aligned_reads + "/{sample}_marked_sorted.bam", sample=samples),
        expand(aligned_reads + "/{sample}_marked_sorted.bai", sample=samples),
        expand(aligned_reads + "/{sample}.g.vcf.gz", sample=samples),
        expand(aligned_reads + "/{sample}.g.vcf.gz.tbi", sample=samples),
        expand(variants + "/{sample}_snps.vcf.gz", sample=samples),
        expand(variants + "/{sample}_indels.vcf.gz", sample=samples)

# creates necessary directories for the outputs 
rule setup_directories:
    output:
        touch("Data/.setup_done")
    shell:
        """
        mkdir -p {fastqc_dir} {aligned_reads} {variants}
        touch Data/.setup_done
        """

# this rule runs FastQC on raw FastQ files to assess sequencing quality before mapping 
# (next step would be trimming if bad quality, here it wasn't the case)
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
        fastqc {input.fastq} -o {params.output_dir}
        """

# Maps reads to the reference genome using BWA and outputs an unsorted BAM.
rule map_reads:
    input:
        fastq = FastQ_Dir + "/{sample}.fastq",
        reference = reference_genome
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

# Sorts the BAM file to prepare it for duplicate marking and downstream analysis
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

# Marks duplicate reads in the sorted BAM file 
# CREATE_INDEX true --> index bam file for next step !
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
        -M {output.marked_bam}.metrics.txt \
        -CREATE_INDEX true \    
        -VALIDATION_STRINGENCY LENIENT \
        -MAX_RECORDS_IN_RAM 1000000 2> {log}
        """

# Calls variants using GATK's HaplotypeCaller and outputs a GVCF.
rule haplotype_caller:
    input:
        bam = aligned_reads + "/{sample}_marked_sorted.bam", 
        reference = reference_genome 
    output:
        gvcf = aligned_reads + "/{sample}.g.vcf.gz",  
        index = aligned_reads + "/{sample}.g.vcf.gz.tbi"  # index gvcf
    log:
        aligned_reads + "/{sample}_haplotypecaller.log"
    params:
        memory = "4g",
        ploidy = 2  # Adjust ploidy if working with non-diploid organisms !!
    shell:
        """
        gatk --java-options "-Xmx{params.memory}" HaplotypeCaller \
        -R {input.reference} \
        -I {input.bam} \
        -O {output.gvcf} \
        -ERC GVCF \
        --native-pair-hmm-threads 4 \
        --sample-ploidy {params.ploidy} \
        2> {log}
        """

# Select Variants

# Selects SNPs from the GVCF 
rule select_snps:
    input:
        gvcf = aligned_reads + "/{sample}.g.vcf.gz",
        index = aligned_reads + "/{sample}.g.vcf.gz.tbi",
        reference = reference_genome
    output:
        snps = variants + "/{sample}_snps.vcf.gz",
        snps_index = variants + "/{sample}_snps.vcf.gz.tbi"
    log:
        variants + "/{sample}_select_snps.log"
    params:
        memory = "4g"
    shell:
        """
        gatk --java-options "-Xmx{params.memory}" SelectVariants \
        -R {input.reference} \
        -V {input.gvcf} \
        --select-type-to-include SNP \
        -O {output.snps} \
        2> {log}
        """
# Selects INDELs from the GVCF
rule select_indels:
    input:
        gvcf = aligned_reads + "/{sample}.g.vcf.gz",
        index = aligned_reads + "/{sample}.g.vcf.gz.tbi",
        reference = reference_genome
    output:
        indels = variants + "/{sample}_indels.vcf.gz",
        indels_index = variants + "/{sample}_indels.vcf.gz.tbi"
    log:
        variants + "/{sample}_select_indels.log"
    params:
        memory="4g"
    shell:
        """
        gatk --java-options "-Xmx{params.memory}" SelectVariants \
        -R {input.reference} \
        -V {input.gvcf} \
        --select-type-to-include INDEL \
        -O {output.indels} \
        2> {log}
        """