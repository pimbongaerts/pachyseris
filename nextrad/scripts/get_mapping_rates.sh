#!/bin/bash
host_path=$1; sym_path=$2
# Output header
echo -ne "sample_name\thost_reads\tsymbiont_reads\tambiguous_reads\ttotal_reads"
# Iterate over all FASTQ files and map to both host and symbiont genomes
for host_bam_file in $host_path/*.bam
do
    # Set variables and empty temporary files
    file_without_path=$( basename "$host_bam_file" )
    symb_bam_file=${sym_path}/${file_without_path}
    sample_name=${file_without_path%_trimmed.bam}
    echo "" | tee temp_host.txt temp_symb.txt
    # Store lists of mapping reads for host and symbiont in temp files
    # Mapping defined as Q20 and not having these flags:
    #         0x4 (segment unmapped)
    #         0x100 (Secondary alignment)
    #         0x800 (supplementary alignment)
    samtools view -F0x904 -q20 $host_bam_file | cut -f1 | sort | uniq \
                  1> temp_host.txt
    samtools view -F0x904 -q20 $symb_bam_file | cut -f1 | sort | uniq \
                  1> temp_symb.txt
    # Determine uniquely mapping and ambiguous reads and output to STDOUT
    reads_totl="$(samtools view $host_bam_file | cut -f1 | sort | uniq | wc -l)"
    reads_host="$(comm -23 temp_host.txt temp_symb.txt | wc -l)"
    reads_symb="$(comm -13 temp_host.txt temp_symb.txt | wc -l)"
    reads_ambi="$(comm -12 temp_host.txt temp_symb.txt | wc -l)"
    echo -ne "$sample_name\t$reads_host\t$reads_symb\t$reads_ambi\t$reads_totl"
done