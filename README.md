# SNV-GATK-Snakemake

project in accordance....all names here

## Description
This project is a pipeline for variant calling using the GATK (Genomic Analysis Toolkit) and Snakemake
short explanation of what SNV is and why its important 

### Tools

## Workflow for variant calling using GATK and Snakemake
short explanation of each rule (bulletpoints
)
Data preparation: ​

setup_directories: creates output directories​

fastqc: runs quality control FASTQC over FASTQ files​

map_reads: maps reads to the reference genome using BWA​

sort_bam: sorts the unsorted BAM files using Samtools​

mark_duplicates: marks duplicate reads using GATK MarkDuplicates​

​

Variant calling: ​

haplotype_caller: calls variants using GATK Haplotype caller, outputs GVCF files​

select_snps: extracts SNPs using GATK SelectVariants ​

select_indels: extracts indels using GATK SelectVariants ​




### Visual Overview
![workflow_figure](rulegraph.svg)

## Future Improvements
base recalibration
indexing of reference genome 