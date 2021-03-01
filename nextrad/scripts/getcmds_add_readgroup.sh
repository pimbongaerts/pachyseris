#!/bin/bash
input_folder=$1

# Variable with gatk path
export GATK_LOCAL_JAR=~/tools/gatk-4.0.1.2/gatk-package-4.0.1.2-local.jar

# Add compulsory Read Group (RG) header and sort/index bam
# (RG could have been supplied with -R parameter in bwa mem)
mkdir -p $input_folder\_final

for bam_file in $input_folder/*.bam
do
    file_without_path=$( basename "$bam_file" )
    sample_name=${file_without_path%_trimmed.bam}
    output_file=${input_folder}_final/${sample_name}.bam

    line="java -jar $GATK_LOCAL_JAR AddOrReplaceReadGroups"
    line="${line} -I $bam_file"
    line="${line} -O $output_file"
    line="${line} -RGID $sample_name"
    line="${line} -RGLB pachyseris"
    line="${line} -RGPL illumina"
    line="${line} -RGPU $sample_name"
    line="${line} -RGSM $sample_name"
    line="${line} -SORT_ORDER coordinate"
    line="${line} --CREATE_INDEX true"
    echo $line
done