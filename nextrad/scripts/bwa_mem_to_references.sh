#!/bin/bash
input_path=$1; host_genome=$2; symb_genome=$3; threads=$4
# Indexing host and symbiont genomes (if not already done)
#if [ ! -f $host_genome.bwt ]; then bwa index $host_genome; fi
#if [ ! -f $symb_genome.bwt ]; then bwa index $symb_genome; fi
# Create output directories
mkdir bwa_output_host bwa_output_symb
# Iterate over all FASTQ files and map to both host and symbiont genomes
for fastq_file in $input_path/*.fq.gz
do
    file_without_path=$( basename "$fastq_file" )
    sample_name=${file_without_path%.fq.gz}
    time=$(date +"%T")
    echo "$time - Mapping $sample_name against host and symbiont references..."
    bwa mem -t $threads -M $host_genome $fastq_file | \
        samtools sort -@ $threads -o bwa_output_host/$sample_name.bam
    bwa mem -t $threads -M $symb_genome $fastq_file | \
        samtools sort -@ $threads -o bwa_output_symb/$sample_name.bam
done
