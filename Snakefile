# Snakefile

configfile: "config.yaml"

import os
import glob

# List all .fastq files in the specified folder
fastq_folder = config["folders"]["fastq_folder"]
fastq_files = glob.glob(os.path.join(fastq_folder, "*.fastq"))
# Extract sample names from file names
# Assuming file names are like "sample1.fastq", "sample2.fastq", etc.
samples = [os.path.basename(f).replace(".fastq", "") for f in fastq_files]

# Extract sample names from file names
samples = [os.path.basename(f).replace(".fastq", "") for f in fastq_files]
# Allow focusing on a single sample dynamically via config for dag workflows
sample = config.get("sample", None)
if sample:
    samples = [sample]


# output directories 
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
        # "Data/.setup_done",
        expand(fastqc_dir + "/{sample}_fastqc.html", sample=samples),
        expand(fastqc_dir + "/{sample}_fastqc.zip", sample=samples),
        expand(variants + "/{sample}_snps.vcf", sample=samples),
        expand(variants + "/{sample}_indels.vcf", sample=samples),
        expand("dag/{sample}_workflow.png", sample=samples)

rule generate_dag:
    input:
        fastqc_dir + "/{sample}_fastqc.html"  # Any key output from the sample workflow
    output:
        dag = "dag/{sample}_workflow.png"  # Diagram output for the sample
    shell:
        """
        mkdir -p dag
        echo 'Generating DAG for sample {wildcards.sample}'
        snakemake --config sample="{wildcards.sample}" --dag | dot -Tpng > {output.dag}
        """

# this rule runs FastQC on raw FastQ files to assess sequencing quality before mapping 
# (next step would be trimming if bad quality, here it wasn't the case)
rule fastqc:
    input:
        fastq = fastq_folder + "/{sample}.fastq"
    output:
        html = fastqc_dir + "/{sample}_fastqc.html",
        zip = fastqc_dir + "/{sample}_fastqc.zip"
    params:
        output_dir = fastqc_dir
    shell:
        """
        echo 'Run FastQC - QC for {wildcards.sample}'
        fastqc {input.fastq} -o {params.output_dir}
        """

# Maps reads to the reference genome using BWA and outputs an unsorted BAM.
rule map_reads:
    input:
        fastq = fastq_folder + "/{sample}.fastq",
        reference = reference_genome
    output:
        sam = aligned_reads + "/{sample}_aligned.sam"
    log:
        aligned_reads + "/{sample}.map.log"
    shell:
        """
        echo 'Map reads to reference genome - {wildcards.sample}'
        bwa mem -t 6 \
        -R "@RG\\tID:{wildcards.sample}\\tPL:ILLUMINA\\tSM:{wildcards.sample}" \
        {input.reference} {input.fastq} > {output.sam} 2> {log}
        """

# Marks duplicate reads in the sorted BAM file 
rule mark_duplicates:
    input:
        unsorted_sam = aligned_reads + "/{sample}_aligned.sam"
    output:
        marked_bam = aligned_reads + "/{sample}_marked_sorted.bam"
    log:
        aligned_reads + "/{sample}.mark.log"
    shell:
        """
        echo 'Marking duplicates for {wildcards.sample}'
        gatk MarkDuplicatesSpark \
            -I {input.unsorted_sam} \
            -O {output.marked_bam} \
            2> {log}
        """

# Calls variants using GATK's HaplotypeCaller and outputs a GVCF.
rule haplotype_caller:
    input:
        bam = aligned_reads + "/{sample}_marked_sorted.bam", 
        reference = reference_genome 
    output:
        vcf = variants + "/{sample}.vcf"  
    log:
        variants + "/{sample}_haplotypecaller.log"
    shell:
        """
        echo 'Run HaplotypeCaller - {wildcards.sample}'
        gatk HaplotypeCaller -R {input.reference} -I {input.bam} -O {output.vcf}
        """

# Select Variants

# Selects SNPs from the VCF 
rule select_snps:
    input:
        vcf = variants + "/{sample}.vcf",
        reference = reference_genome
    output:
        snps = variants + "/{sample}_snps.vcf"
    log:
        variants + "/{sample}_select_snps.log"
    shell:
        """
        echo 'Selecting SNPs for {wildcards.sample}'
        gatk SelectVariants -R {input.reference} -V {input.vcf} --select-type SNP -O {output.snps}
        """

# Selects INDELs from the GVCF
rule select_indels:
    input:
        vcf = variants + "/{sample}.vcf",
        reference = reference_genome
    output:
        indels = variants + "/{sample}_indels.vcf"
    log:
        variants + "/{sample}_select_indels.log"
    shell:
        """
        gatk SelectVariants -R {input.reference} -V {input.vcf} --select-type INDEL -O {output.indels}
        """

