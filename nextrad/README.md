# "Morphological stasis masks ecologically distinct coral species on tropical reefs"

TO BE COMPLETED: VCF files will be uploaded after review (to the `datafiles` folder).

## B - Reduced-representation sequencing (nextRAD)

Throughout the notebook there are references to the following environmental variables:

```shell
RAW_SEQ_PATH      # Folder with raw nextRAD sequence data
HOST_GENOME_FA    # Pachyseris genome FASTA (pspe_final_0.12.fasta)
SYMB_GENOME_FA    # Cladocopium genome fasta (sym_C1_genome.fasta)
GATK3_JAR         # JAR file for GATK3 (GenomeAnalysisTK.jar)
GATK4_JAR         # JAR file for GATK4 (gatk-package-4.0.1.2-local.jar)
PICARD_JAR        # JAR file for Picard (picard.jar)
MP_THREADS        # number of cores available for analyses (e.g. 32)
```

## B1 - Read filtering, mapping and variant calling

### Read filtering

[TrimGalore](https://github.com/FelixKrueger/TrimGalore) was used to trim Nextera adapters (with overlap of at least 5bp) and low-quality ends (below PHRED score of 20) , while discarding reads that become less than 30bp. After that, poor-performing samples (with less than 50MB of data) were removed.

```shell
$ find RAW_SEQ_PATH -name "*.fastq.gz" | parallel -j $MP_THREADS trim_galore -q 20 --length 30 \
                                                                             --stringency 5 \
                                                                             --nextera --phred33 {}
$ find . -name "*.fq.gz" -size -50M -delete; ls -1 *.fq.gz | wc -l
544 # number of remaining samples
```

### Mapping reads to genomes

We used [BWA-MEM](https://github.com/lh3/bwa) to map reads of all samples to the *Pachyseris speciosa* genome, and the *Cladocopium* ([Liu et al. 2018](https://www.nature.com/articles/s42003-018-0098-3)) genome to assess contamination using the [bwa_mem_to_reference.sh](scripts/bwa_mem_to_references.sh) script, and then summarized the mapping rates with the [get_mapping_rates.sh](scripts/get_mapping_rates.sh) script:

```shell
$ ./bwa_mem_to_references.sh $RAW_SEQ_PATH $HOST_GENOME_FA $SYMB_GENOME_FA $MP_THREADS
$ ./get_mapping_rates.sh bwa_output_host bwa_output_symb > pachy_mapping_rates_b1.txt
```

### Prepare BAM file for variant calling

[Samtools](http://www.htslib.org/) was used to index the *P. speciosa* genome fasta, and [CreateSequenceDictionary (Picard)](https://gatk.broadinstitute.org/hc/en-us/articles/360037422891-CreateSequenceDictionary-Picard-) to used to create a sequence dictionary:

```shell
$ samtools faidx $HOST_GENOME_FA
$ java -jar $PICARD_JAR CreateSequenceDictionary R=$HOST_GENOME_FA \
                                                 O=$GENOME_REF_PATH/pspe_final_0.12.dict
```

Using the [getcmds_add_readgroup.sh](scriopts/getcmds_add_readgroup.sh) script, we added compulsory Read Group (RG) headers with the [AddOrReplaceReadGroups](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-) command (script produces a list of commands to run in Picard/GATK; also removing `_trimmed` in output name) - executed using the [parallel](https://www.gnu.org/software/parallel/) tool:

```shell
$ parallel -j $MP_THREADS < <(getcmds_add_readgroup.sh bwa_output_host $GATK4_JAR)
```

### Variant calling with GATK

Create a file with a list of samples and run the [GATK UnifiedGenotyper](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php) for all samples at once (NOTE: not available in latest GATK version as replaced with HaplotypeCaller; however HaplotypeCaller too slow given the large number of samples):

```shell
$ find bwa_output_host_final -name \*.bam > bam_files.list
$ java -jar $GATK3_JAR -T UnifiedGenotyper \
                       -R $HOST_GENOME_FA \
                       -nt $MP_THREADS \
                       -nct 1 \
                       --genotype_likelihoods_model SNP \
                       -I bam_files.list \
                       -o pachy_b1_orig.vcf
```

### Variant QC filtering

Filtered variants based on GATK general recommendations for hard-filtering:

```shell 
$ java -jar $GATK3_JAR -T VariantFiltration \
                       -R $HOST_GENOME_FA \
                       -V pachy_b1_orig.vcf \
                       -o pachy_b1_temp1.vcf \
                       --filterExpression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || HaplotypeScore > 13.0" \
                       --filterName "std_gatk_hard_filter"
$ vcftools --vcf pachy_b1_temp1.vcf \
           --remove-filtered-all \
           --recode --recode-INFO-all --stdout \
           | gzip -c > pachy_b1_temp2.vcf.gz
$ grep "After" out.log
After filtering, kept 544 out of 544 Individuals
After filtering, kept 4839165 out of a possible 11065495 Sites
```

Reduce the dataset, by applying an additional hard-call filter on individual genotypes, using a minimum of 10x coverage and a genotyping quality of 30, and filtering for bi-allelelic SNPs with a minimum allele frequency (MAF) of 0.01. For now, retaining sites that are genotyped for at least 20% of samples (as to retain lineage-specific for downstream subsetting of the dataset):

```shell
$ vcftools --gzvcf pachy_b1_temp2.vcf.gz \
           --minDP 10 --maxDP 1000 --minQ 30 \
           --min-alleles 2 --max-alleles 2 \
           --maf 0.01 --max-missing 0.2 \
           --recode --recode-INFO-all --stdout > pachy_b1_temp3.vcf
$ grep "After" out.log
After filtering, kept 544 out of 544 Individuals
After filtering, kept 94439 out of a possible 4839165 Sites
```

Assess missing data (using [vcf_missing_data.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_missing_data.py) script) and remove very low-performing samples (<10% of data; note: will be filtered more stringently downstream):

```shell
$ vcf_missing_data.py pachy_b1_temp3.vcf > pachy_missing_data_b1.txt
$ vcftools --vcf pachy_b1_temp3.vcf \
           --remove <(tail +2 pachy_missing_data_b1.txt | awk '$5<10' | cut -f 1) \
           --mac 1 --recode --stdout | gzip -c > pachy_b1.vcf.gz
$ grep "After" out.log
After filtering, kept 516 out of 544 Individuals
After filtering, kept 94438 out of a possible 94439 Sites
```

## B2 - Preparation of overall dataset

### Eliminate low-performing samples and SNPs

Reduce dataset to SNPs that are genotyped for at least 80% of samples, and samples that are genotyped for at least 50% of SNPs: 

```shell
$ vcftools --gzvcf pachy_b1.vcf.gz \
           --max-missing 0.8 \
           --recode --stdout > pachy_b2_temp1.vcf
$ vcf_missing_data.py pachy_b2_temp1.vcf > pachy_missing_data_b2.txt
$ vcftools --vcf pachy_b2_temp1.vcf \
           --remove <(tail +2 pachy_missing_data_b2.txt | awk '$5<50' | cut -f 1) \
           --mac 1 --recode --stdout > pachy_b2_temp2.vcf
$ grep "After" out.log
After filtering, kept 501 out of 516 Individuals
After filtering, kept 8536 out of a possible 8536 Sites
```

### Assess genotyping accuracy of replicate samples

Assess "genotyping accuracy" observed for sperm (1x triplicate samples) versus regular sample replicates (7x duplicate samples) (all listed in [pachy_replicates_b2.txt](popfiles/pachy_replicates_b2.txt)) (using [vcf_clone_detect.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_clone_detect.py) script):

```shell
$ vcftools --vcf pachy_b2_temp2.vcf \
           --keep pachy_replicates_b2.txt \
           --recode --stdout > pachy_replicates_b2.vcf
$ vcf_clone_detect.py -v pachy_replicates_b2.vcf > pachy_reps_comparison_b2.txt
$ sed -n -e '/###3/,/###4/ p' pachy_reps_comparison_b2.txt | head -n 15
###3 - List of highest matches
99.84   0       [NA]    PSHSPERM41 vs PSHSPERM42        7592.5/7605     8385    7756
99.81   0.03    [NA]    PSHSPERM31 vs PSHSPERM42        7270.5/7284     8064    7756
99.8    0.01    [NA]    PSHSPERM31 vs PSHSPERM41        7552.5/7568     7890    8214
99.62   0.18    [NA]    PSRUDX0009 vs PSRUDH0009        7975.5/8006     8477    8065
99.6    0.02    [NA]    PSGUDH7699 vs PSGUDX7699        7697.0/7728     8048    8216
99.55   0.05    [NA]    PSCYDH8401 vs PSCYDX8401        7991.5/8028     8287    8277
99.5    0.05    [NA]    PSPGDX0957 vs PSPGDH0957        7741.0/7780     8514    7802
99.08   0.42    [NA]    PSCBDH8299 vs PSCBDX8299        7922.5/7996     8225    8307
99.02   0.06    [NA]    PSOSSH7184 vs PSOSSX7184        7118.5/7189     7641    8084
98.91   0.11    [NA]    PSGABX7659 vs PSGABH7659        7934.5/8022     8504    8054
98.0    -------------------- Potential threshold --------------------
93.64   5.27    [NA]    PSGUDH7699 vs PSHSPERM31        6725.5/7182     8043    7675
93.58   0.06    [NA]    PSHSPERM31 vs PSGUDX7699        6937.0/7413     7733    8216
93.53   0.05    [NA]    PSPGDX0957 vs PSHSPERM41        7573.5/8097     8419    8214
```

Retain the sample with the least missing data for each set of replicates:

```shell
$ vcftools --vcf pachy_b2_temp2.vcf \
           --remove <(sed -n -e '/###5/,// p' pachy_reps_comparison_b2.txt | grep -v "#") \
           --mac 1 --recode --stdout > pachy_b2_temp3.vcf
$ grep "After" out.log
After filtering, kept 492 out of 501 Individuals
After filtering, kept 8536 out of a possible 8536 Sites
```

### Detect and eliminate clones

Extract population assignment from sample name (using [popfile_from_vcf.py](https://github.com/pimbongaerts/radseq/blob/master/popfile_from_vcf.py) script) and assess potential clones (using [vcf_clone_detect.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_clone_detect.py) script):

```shell
$ popfile_from_vcf.py pachy_b2_temp3.vcf 3 5 > pachy_popfile_b2.txt
$ vcf_clone_detect.py -v pachy_b2_temp3.vcf \
                      -p pachy_popfile_b2.txt \
                      -o pachy_comparison_b2.txt > pachy_clone_detection_b2.txt
$ sed -n -e '/###2/,/###3/ p' pachy_clone_detection_b2.txt | head -n 15
###2 - Histogram (of pairwise genetic similarities)
 86       5 *****
 87    4774 *********************************************************************#
 88   44938 *********************************************************************#
 89   31556 *********************************************************************#
 90    3874 *********************************************************************#
 91   11786 *********************************************************************#
 92   15697 *********************************************************************#
 93    7525 *********************************************************************#
 94     572 *********************************************************************#
 95      22 **********************
 96       1 *
 97       4 ****
 98      15 ***************
 99      17 *****************
```

Lowest genotypic similarity between replicate samples was 98.91% (see above; `PSGABX7659` vs `PSGABH7659`), however when assessing histogram it seems safer to use a slightly lower cut-off to eliminate potential clones. Given that above 97% all matches are within-site comparisons, and below there start to be across-site comparisons, we refiltered using a 97.0% clone-threshold:

```shell
$ vcf_clone_detect.py -i pachy_comparison_b2.txt -t 97.0 > pachy_clone_detection_97_b2.txt
$ sed -n -e '/###3/,/##4/ p' pachy_clone_detection_97_b2.txt | grep -v "###4" | tail -n 12
98.37    0.05    [CYS]    PSCYSH8349 vs PSCYSH8357    5956.5/6055    7492    7099
97.98    0.39    [CYS]    PSCYSH8349 vs PSCYSH8345    6173.5/6301    7326    7511
97.93    0.05    [RDD]    PSRDDH0030 vs PSRDDH0036    6350.0/6484    6767    8253
97.92    0.01    [GDS]    PSGDSH8116 vs PSGDSH8115    6673.5/6815    8072    7279
97.2    0.72    [CBD]    PSCBDH8286 vs PSCBDH8290    5204.0/5354    5764    8126
97.0    -------------------- Manual threshold --------------------
96.65    0.55    [OKD-OSS]    PSOSSH7168 vs PSOKDH7049    4802.5/4969    5488    8017
95.4    1.25    [OSS]    PSOSSH7161 vs PSOSSH7177    6491.0/6804    7752    7588
95.25    0.15    [PGS-PVS]    PSPVSH1311 vs PSPGSH1059    7141.0/7497    8340    7693
95.23    0.02    [OSS]    PSOSSH7173 vs PSOSSH7177    7030.0/7382    8330    7588
95.18    0.05    [OSS]    PSOSSH7171 vs PSOSSH7185    7023.5/7379    7770    8145
```

Identify samples to be removed in non-clonal datasets (retaining sample with most data; except for the *P. rugosa* (and one misidentified *P. rugosa*: PSGMDH0861 (=> PRGMDH0861) which are retained at this stage):

```shell
$ sed -n -e '/###5/,// p' pachy_clone_detection_97_b2.txt \
      | grep -v "#" | grep -v "PR" \
      | grep -v "PSGMDH0861" > pachy_clones_to_remove_b2.txt
$ cat pachy_clones_to_remove_b2.txt <(echo -e "PRH0000001\nPRH0000002\nPSGMDH0861") \
      > pachy_clones_and_rugosa_to_remove_b2.txt
```

Eliminate duplicate genotypes (clones) and *P. rugosa* samples from dataset:

```shell
$ vcftools --vcf pachy_b2_temp3.vcf \
           --remove pachy_clones_and_rugosa_to_remove_b2.txt \
           --recode --stdout > pachy_b2.vcf
$ popfile_from_vcf.py pachy_b2.vcf 3 5 > pachy_popfile_b2.txt
$ grep "After" out.log
After filtering, kept 465 out of 492 Individuals
After filtering, kept 8536 out of a possible 8536 Sites
```

Dataset with *P. rugosa* but no duplicate genotypes (clones):

```shell
$ vcftools --vcf pachy_b2_temp3.vcf 
           --remove pachy_clones_to_remove_b2.txt 
           --recode --stdout > pachy_outgroup_b2.vcf
$ popfile_from_vcf.py pachy_outgroup_b2.vcf 3 5 > pachy_popfile_outgroup_b2.txt
$ grep "After" out.log
After filtering, kept 468 out of 492 Individuals
After filtering, kept 8536 out of a possible 8536 Sites
```

## B3 - Overall genetic structure

### B3a - Assess overall structure

Conduct a basic principal component analysis (PCA) on the overall dataset (excluding clones and *Pachyseris rugosa* samples) using the [adegenet](https://github.com/thibautjombart/adegenet/wiki) package in R (using RGB color to also visualize the 3rd PC; using [vcf2pca.R](scripts/vcf2pca.R]) script):

```shell
$ Rscript vcf2pca.R pachy_b2.vcf 5
```

Calculate genetic distances and construct a basic NJ for the overall dataset (excluding clones and *Pachyseris rugosa* samples) (using [vcf_gdmatrix.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_gdmatrix.py) and [gdmatrix2tree.py](https://github.com/pimbongaerts/radseq/blob/master/gdmatrix2tree.py) scripts):

```shell
$ vcf_gdmatrix.py pachy_b2.vcf pachy_popfile_b2.txt > pachy_gd_b3a.txt
$ gdmatrix2tree.py pachy_gd_b3a.txt pachy_gd_b3a.tre
```

Calculate genetic distances and construct a basic NJ for the overall dataset with *Pachyseris rugosa* samples (but no clones):

```shell
$ vcf_gdmatrix.py pachy_outgroup_b2.vcf pachy_popfile_outgroup_b2.txt > pachy_outgroup_gd_b3a.txt
$ gdmatrix2tree.py pachy_outgroup_gd_b3a.txt pachy_outgroup_gd_b3a.tre
```

### B3b - Snapclust

Extract five replicate datasets with random SNPs from main dataset, with each SNP separated by at least 2,500bp (to eliminate linked SNPs in close proximity) (using [vcf_single_snp.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_single_snp.py) script):

```shell
$ for i in {1..5}; do \
      vcf_single_snp.py pachy_b2.vcf -d 2500 > pachy_singlesnp$i\_b3b.vcf; \
  done
```

Evaluate number of clusters based on snapclust's AIC, BIC, and KIC (using [snapclust_chooseK.R](scripts/snapclust_chooseK.R) script) :

```shell
$ Rscript snapclust_chooseK.R pachy_singlesnp*.vcf
```

Snapclust was then ran for k=4-6 using the default recommended "Ward" algorithm to define initial group assignments,  with 50 iterations of the EM algorithm across the 5 replicates (using [snapclust.R](scripts/snapclust.R) script):

```shell
$ Rscript snapclust.R pachy_singlesnp*.vcf 4 6
```

### B3c - STRUCTURE

Run STRUCTURE [using the multi-threading wrapper [structure_mp](https://github.com/pimbongaerts/radseq/blob/master/structure_mp.py)]  and plot output (using [structure_mp_plot.py](https://github.com/pimbongaerts/radseq/blob/master/structure_mp_plot.py) script; note that the thinning threshold of 2,500bp is currently hard-coded in the script):

```shell
$ structure_mp.py pachy_b2.vcf pachy_popfile_b2.txt 6 20 20
Executing 20 parallel STRUCTURE runs for K = 6 ...10 reps DONE
Running CLUMPP on replicates for K = 6 ...
K = 6: MedMeaK 6 MaxMeaK 6 MedMedK 6 MaxMedK 6      pachy_b2.vcf
$ structure_mp_plot.py pachy_b2_1519698915 # (=structure_mp output name)
```

For each of the 20 CLUMPP-aligned replicates, check the maximum assignment value for each cluster, to assess whether there are runs with "ghost clusters":

```shell
$ for file in pachy_b2_1519698915/clumpp_K6_perm*.csv; do \
      echo -n "$file\t"; \
      for i in {3..9}; do \
          cut -d, -f$i $file | sort -nr | head -1 | xargs printf '%s\t'; \
      done; \
      echo ""; \
  done > pachy_b2_1519698915/structure_max_assignments_b3c.txt
```

Only retain the CLUMPP replicates that have at maximum assignment of at least 0.99 for every cluster (i.e. at least one sample has a full assignment to that cluster), copy those to a separate folder (`retained_runs`), and then produce a new summary using the mean of the assignment across replicates (using [csv_mean.py](https://github.com/pimbongaerts/radseq/blob/master/csv_mean.py) script; the `3` indicates that the data starts at the third column):

```shell
$ mkdir retained_runs
$ cat pachy_b2_1519698915/structure_max_assignments_b3c.txt \
      | awk '$2 > 0.99 && $3 > 0.99 && $4 > 0.99 && $5 > 0.99 && $6 > 0.99 && $7 > 0.99' \
      | cut -f1 \
      | while read file; \
        do \
            cp pachy_b2_1519698915/$file retained_runs/$file; \
        done
$ csv_mean.py retained_runs 3 > pachy_structure_k6_b3c.csv
```

Extract cluster assigments based on a 0.8 (lenient: `pachy_popfile_clusters_b3c.txt`) and 0.95 (stringent: `pachy_popfile_clusters_stringent_b3c.txt`) ancestry coefficient cut-off [using [popfile_from_clusters.py](https://github.com/pimbongaerts/radseq/blob/master/popfile_from_clusters.py) script] and rename cluster groups :

```shell
$ popfile_from_clusters.py pachy_structure_k6_b3c.csv  0.8 > pachy_popfile_clusters_b3c.txt
$ popfile_from_clusters.py pachy_structure_k6_b3c.csv  0.95 > pachy_popfile_clusters_stringent_b3c.txt
$ cat pachy_popfile_clusters_b3c.txt \
      | sed 's/CLUSTER_1/RED/g' | sed 's/CLUSTER_2/BL2/g' \
      | sed 's/CLUSTER_3/BLU/g' | sed 's/CLUSTER_4/ISR/g' \
      | sed 's/CLUSTER_5/GR2/g' | sed 's/CLUSTER_6/GRN/g' \
      > pachy_popfile_clusters_b3c.txt
$ cat pachy_popfile_clusters_stringent_b3c.txt \
      | sed 's/CLUSTER_1/RED/g' | sed 's/CLUSTER_2/BL2/g' \
      | sed 's/CLUSTER_3/BLU/g' | sed 's/CLUSTER_4/ISR/g' \
      | sed 's/CLUSTER_5/GR2/g' | sed 's/CLUSTER_6/GRN/g' \
      > pachy_popfile_clusters_stringent_b3c.txt
```

Quick summary of the assignments:

```shell
$ paste <(cut -f2 pachy_popfile_clusters_b3c.txt | sort | uniq -c | sort -rn) \
                <(cut -f2 pachy_popfile_clusters_stringent_b3c.txt | sort | uniq -c | sort -rn)
# Normal     Stringent
 146 GRN     136 GRN
 138 BLU     120 BLU
 107 RED     105 RED
  24 UNASSIGNED      55 UNASSIGNED
  20 ISR      20 ISR
  20 BL2      20 BL2
  10 GR2       9 GR2
$ grep -E "PSG|PSC" pachy_popfile_clusters_b3c.txt | grep UNASSIGNED | wc -l 
6 # Total number of "admixed" samples in GBR and WCS (normal)
$ grep -E "PSG|PSC" pachy_popfile_clusters_stringent_b3c.txt | grep UNASSIGNED | wc -l 
14 # Total number of "admixed" samples in GBR and WCS (stringent)
```

Check for correspondence with NJ tree by coloring samples according to STRUCTURE lenient (0.8) assignments (using [nexus_set_label_colors.py](https://github.com/pimbongaerts/radseq/blob/master/nexus_set_label_colors.py) script):

```shell
$ cat pachy_popfile_clusters_b3c.txt \
      | sed 's/GRN/#4daf4a/g' | sed 's/RED/#e41a1c/g' \
      | sed 's/BLU/#377eb8/g' | sed 's/ISR/#810f7c/g' \
      | sed 's/GR2/#238b45/g' | sed 's/BL2/#08519c/g' \
      | sed 's/UNASSIGNED/#000000/g' > pachy_popfile_cluster_colors_b3c.txt
$ nexus_set_label_colors.py pachy_gd_outgroup_b3a.tre \
                            pachy_popfile_cluster_colors_b3c.txt \
                            > pachy_gd_cluster_colors_b3c.tre
```

Combine the stringent assignments with the original mapping rates file, to be able to assess mapping rate differences between lineages:

```shell
$ join -a1 <(grep -v "sample_name" pachy_mapping_rates_b1.txt | sort ) \
           <(sort pachy_popfile_clusters_stringent_b3c.txt) \
           > pachy_mapping_rates_by_lineage_b3c.txt
```

## B4 - Cluster comparison

### B4a - Pairwise Fst between lineages and regions

Recreate a new VCF file with only stringently assigned samples (coancestry of >=0.95), grouped by cluster and region (excluding those with only 1 sample: `GRN_H` and `RED_P`), and refiltered for biallelic SNPs that are genotyped for at least half (0.5) of the samples in each region within each of the lineages (RED, GRN and BLU) (using [vcf_minrep_filter.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_minrep_filter.py) script). Note: FORMAT headers removed from vcf due to incompatibility issue with PyVCF.:

```shell
$ grep -v "UNA" pachy_popfile_clusters_stringent_b3c.txt \
       | awk '{print $1 "\t" $2 "_" substr($1,3,1)}' \
       | grep -vE "GRN_H|RED_P" > pachy_popfile_clusters_stringent_b4a.txt
$ vcftools --gzvcf pachy_b1.vcf.gz \
           --keep pachy_popfile_clusters_stringent_b4a.txt \
           --mac 1 --recode --stdout > pachy_clusters_b4a_temp1.vcf
$ grep -v "##FORMAT" pachy_clusters_b4a_temp1.vcf > pachy_clusters_b4a_temp2.vcf
$ vcf_minrep_filter.py pachy_clusters_b4a_temp2.vcf \
                       pachy_popfile_clusters_stringent_b4a.txt \
                       0.5 \
                       pachy_clusters_b4a.vcf
74344 out of 91999 SNPs failing 0.5 threshold
Lowest overall proportion genotyped of remaining SNPs: 0.5343137254901961
$ grep -v "#" pachy_clusters_b4a.vcf | wc -l
   17655 # SNPs remaining
$ grep -v "#" pachy_clusters_b4a.vcf | cut -f1 | sort | uniq | wc -l
    978 # Scaffolds with SNPs
```

Calculate Weir and Cockerham Fst between the clusters and regions (creating all pairwise combinations using the [itertools_combinations.py](https://github.com/pimbongaerts/radseq/blob/master/itertools_combinations.py) script):

```shell
# Create file with all popfile names (and exclude groups with small samplesize)
$ cut -f2 pachy_popfile_clusters_stringent_b4a.txt \
      | sort | uniq > pachy_popfile_clusternames_b4a.txt
# Create a separate popfile for each group
$ while read i; do \
      grep $i pachy_popfile_clusters_stringent_b4a.txt > $i.txt; \
  done <pachy_popfile_clusternames_b4a.txt
# Calculate mean Fst (using vcftools) for each pairwise comparison of groups
$ itertools_combinations.py pachy_popfile_clusternames_b4a.txt | wc -l
      45 # number of comparisons
$ itertools_combinations.py pachy_popfile_clusternames_b4a.txt \
      | awk '{system("vcftools --vcf pachy_clusters_b4a.vcf \
                               --weir-fst-pop "$1".txt \
                               --weir-fst-pop "$2".txt \
                               --out "$1"_vs_"$2)}'
```

Reformat outputted Fst values to CSV file so they can be imported into R:

```shell
$ grep "^Weir and Cockerham mean Fst" *.log \
      | awk -F "_vs_|.log|estimate: " '{ print $1 "," $2 "," $4 }' \
      > pachy_fsts_b4a.csv
```

### B4b - Genome-wide Fst between main lineages

#### Create new VCF focusing only on main lineages

Recreate a new VCF file with only stringently assigned samples (coancestry of >=0.95) belonging to the 4 main lineages, and focusing only on the `GBR`, `WCS` and `ISR` populations. Refiltering  for biallelic SNPs that are present in  at least 10 samples in each region within each lineage (so min. 70 indivs in total genotyped for each SNP) (using [vcf_minrep_filter_abs.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_minrep_filter_abs.py) script). Note: FORMAT headers removed from vcf due to incompatibility with PyVCF):

```shell
$ grep -v "UNA" pachy_popfile_clusters_stringent_b3c.txt \
      | grep -E "PSG|PSC|PSR" \
      | awk '{print $1 "\t" $2 "_" substr($1,3,1)}' \
      > pachy_popfile_clusters_b4b.txt
$ vcftools --gzvcf pachy_b1.vcf.gz \
           --keep pachy_popfile_clusters_b4b.txt \
           --mac 1 --recode --stdout \
           | grep -v "##FORMAT" > pachy_clusters_b4b_temp1.vcf
$ vcf_minrep_filter_abs.py pachy_clusters_b4b_temp1.vcf \
                           pachy_popfile_clusters_b4b.txt \
                           10 \
                           pachy_clusters_b4b_temp2.vcf
59371 out of 88658 SNPs failing 10.0 threshold
Lowest overall proportion genotyped of remaining SNPs: 0.29155313351498635
$ grep -v "#" pachy_clusters_b4b_temp2.vcf | wc -l
   29287 # SNPs remaining
$ grep -v "#" pachy_clusters_b4b_temp2.vcf | cut -f1 | sort | uniq | wc -l
    1128 # Scaffolds with SNPs
```

Quick summary of sample sizes for lineages split by region:

```shell
$ cut -f2 pachy_popfile_clusters_b4b.txt | sort | uniq -c
  47 BLU_C
  73 BLU_G
  29 GRN_C
  94 GRN_G
  20 ISR_R
  17 RED_C
  87 RED_G
```

#### Annotate VCF with gene regions and impacts

Build a new genome reference database in snpEff:

```shell
$ mkdir pspe_final_0.12 && cd "$_" # data directory with symlinks for snpeff
$ ln -s $GENOME_REF_PATH/pspe_0.12.maker_post_001.genes.gff genes.gff
$ ln -s $GENOME_REF_PATH/pspe_0.12.maker_post_001.proteins.fasta protein.fa
$ ln -s $GENOME_REF_PATH/pspe_final_0.12.fasta sequences.fa
$ echo "# Pachyseris speciosa genome, version 0.12\n
        data.dir = /home/pbongaerts/pachyseris/LFS/snpeff\npspe_final_0.12.genome : Pachyseris" \
      >> snpEff.config # add data directory to config file
$ snpEff build -gff3 -v pspe_final_0.12 # Create snpEff database
```

Annotate VCF file with snpEff:

```shell
$ snpEff pspe_final_0.12 pachy_clusters_b4b_temp2.vcf > pachy_annotated_b4b.vcf
```

Extract text tile with annotations from VCF:

```shell
$ grep -v "#" pachy_annotated_b4b.vcf | cut -f1,2,8 \
      | cut -d '|' -f1,2,3,4 | tr "|" "\\t" > pachy_annotations_b4b.txt
```

#### Calculate genome-wide Fst values

Calculate Weir and Cockerham Fst between the clusters and regions:

```shell
# Create file with all popfile names (and exclude groups with small samplesize)
$ cut -f2 pachy_popfile_clusters_b4b.txt | sort | uniq > pachy_popfile_clusternames_b4b.txt
# Create a separate popfile for each group
$ while read i; do \
      grep $i pachy_popfile_clusters_b4b.txt > $i.txt; \
  done < pachy_popfile_clusternames_b4b.txt
# Calculate mean Fst (using vcftools) for each pairwise comparison of groups
$ itertools_combinations.py pachy_popfile_clusternames_b4b.txt \
      | awk '{system("vcftools --vcf pachy_annotated_b4b.vcf \
                               --weir-fst-pop "$1".txt \
                               --weir-fst-pop "$2".txt \
                               --out "$1"_vs_"$2)}'
```

#### Create mapping reference to A. millepora pseudo chromosomes (for visualization only)

Use [RaGOO](https://github.com/malonge/RaGOO) to map P. speciosa scaffolds to A. millepora chromosomes (filtering fasta using [fasta_include.py](https://github.com/pimbongaerts/radseq/blob/master/fasta_include.py) and [fasta_exclude.py](https://github.com/pimbongaerts/radseq/blob/master/fasta_exclude.py) script, and extracting chromosome order with [vcf_ragoo_order.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_ragoo_order.py) script):

```shell
# Limit A. millepora reference genome to established chromosomes
$ fasta_include.py Amil.v2.01.chrs.fasta \
      <(grep "chr" Amil.v2.01.chrs.fasta | sed 's/>//') \
      > Amil_chroms_only.fasta
# Map P. speciosa draft to A. millepora chromosomes
$ ragoo.py pspe_final_0.12.fasta Amil_chroms_only.fasta
# Number of unmapped scaffolds
$ wc -l ragoo_output/orderings/Chr0_orderings.txt
     629 Chr0_orderings.txt
# Create ragoo fasta without Chr0 for dotplot:
$ fasta_exclude.py ragoo.fasta <(echo "Chr0_RaGOO") > ragoo_nochr0.fasta
# Create ordering index
$ vcf_ragoo_order.py BLU_C_vs_BLU_G.weir.fst ragoo_output/orderings | cut -f1-5 \
      > pachy_chrom_pos_order_b4b.txt
```

### B4c - Assessment of alternately fixed and outlier SNPs

#### Identify SNPs close to alternate fixation

Create hierarchical popfile with assignment of each sample to both lineage (`BLU`, `RED`, etc.) and region (`BLU_C`, `BLU_G`, etc.):

```shell
$ awk '{print $1 "\t" substr($2,1,3) "\t" substr($2,1,3) "_" substr($1,3,1)}' \
            pachy_popfile_clusters_stringent_b3c.txt > pachy_groups_b4b.txt
```

Identify SNPs that have are close to alternative fixation (allele frequency differential of at least 0.95) for each pairwise lineage comparison (using the `GBR` and `WCS` regions as replicates for RED, GRN and BLU linages) (using [vcf_afd_filter.py](https://github.com/pimbongaerts/radseq/blob/master/vcf_afd_filter.py) script):

```shell
$ vcf_afd_filter.py pachy_annotated_b4b.vcf \
      <(grep -E "PSG|PSC" pachy_groups_b4b.txt | grep -E "BLU|GRN") 0.95 \
      > pachy_fixed_BLUvsGRN.csv
$ vcf_afd_filter.py pachy_annotated_b4b.vcf \
      <(grep -E "PSG|PSC" pachy_groups_b4b.txt | grep -E "GRN|RED") 0.95 \
      > pachy_fixed_GRNvsRED.csv
$ vcf_afd_filter.py pachy_annotated_b4b.vcf \
      <(grep -E "PSG|PSC" pachy_groups_b4b.txt | grep -E "RED|BLU") 0.95 \
      > pachy_fixed_REDvsBLU.csv
$ vcf_afd_filter.py pachy_annotated_b4b.vcf \
      <(grep -E "PSG|PSC|PSR" pachy_groups_b4b.txt | grep -E "ISR|BLU") 0.95 \
      > pachy_fixed_ISRvsBLU.csv
$ vcf_afd_filter.py pachy_annotated_b4b.vcf \
      <(grep -E "PSG|PSC|PSR" pachy_groups_b4b.txt | grep -E "ISR|GRN") 0.95 \
      > pachy_fixed_ISRvsGRN.csv
$ vcf_afd_filter.py pachy_annotated_b4b.vcf \
      <(grep -E "PSG|PSC|PSR" pachy_groups_b4b.txt | grep -E "ISR|RED") 0.95 \
      > pachy_fixed_ISRvsRED.csv
```

Summarize number of SNPs close to alternate fixation:

```shell
$ for i in *.csv; do \
      echo $i; \
      grep "AFD_OUTLIER" $i | cut -d, -f3 | sort | uniq -c; \
  done | grep AFD_OUTLIER
  57 AFD_OUTLIER:BLU/GRN
  78 AFD_OUTLIER:RED/GRN
  62 AFD_OUTLIER:BLU/ISR
 108 AFD_OUTLIER:GRN/ISR
 207 AFD_OUTLIER:RED/ISR
  91 AFD_OUTLIER:BLU/RED
$ grep "AFD_OUTLIER" pachy_fixed_*.csv | cut -d: -f2 | cut -d, -f1,2 | sort | uniq | wc -l
    311 # Number of SNPs close to alternate fixation separating GRN/BLU/RED/ISR
$ grep "AFD_OUTLIER" pachy_fixed_BLUvsGRN.csv pachy_fixed_GRNvsRED.csv pachy_fixed_REDvsBLU.csv \
      | cut -d: -f2 | cut -d, -f1,2 | sort | uniq | wc -l
  146    # Number of SNPs close to alternate fixation separating GRN/BLU/RED
$ grep -v "#" pachy_annotated_b4b.vcf | wc -l
   29287    # Total number of SNPs
$ echo "(146/29287)*100" | bc -l
.49851469935466247800    # Percentage of SNPs that are alternatively fixed
```

#### Identify outlier SNPs with pcadapt

Create VCFs for pairwise comparisons between lineages:

```shell
# Create popfile for each pairwise comparison of lineages
$ itertools_combinations.py <(echo "BLU\nGRN\nRED\nISR") \
      | awk '{system("grep -E \"PSG|PSC|PSR\" pachy_groups_b4b.txt | grep -E \""$1"|"$2"\" > pachy_groups_"$1"_"$2".txt ")}'
# Create VCF for each pairwise comparison of lineages
$ itertools_combinations.py <(echo "BLU\nGRN\nRED\nISR") \
      | awk '{system("vcftools --vcf pachy_annotated_b4b.vcf --keep pachy_groups_"$1"_"$2".txt --recode --stdout > pachy_pcadapt_"$1"_vs_"$2".vcf ")}'
```

Identify outliers SNPs using pcadapt (Luu et al. 2016) for each pairwise comparison (using only the first PC and verifying that it separates only the lineages) based on q-values with an expected false discovery rate lower than 5% (using [pachyseris_pcadapt_between_lineages.R](scripts/pachyseris_pcadapt_between_lineages.R) script):

```shell
$ Rscript pachyseris_pcadapt_between_lineages.R
# Merge CHROM and POS from VCF with pcadapt output and filter outliers
$ for i in *_qvals.txt; do \
      paste -d" " <(grep -v "#" pachy_annotated_b4b.vcf | awk '{print $1 " " $2}') $i \
          | awk '$4 < 0.05' > pachy_outlier_$i; 
  done
$ wc -l pachy_outlier*.txt
     174 pachy_outlier_BLUvGRN_qvals.txt
    2671 pachy_outlier_BLUvISR_qvals.txt
     289 pachy_outlier_BLUvRED_qvals.txt
    2225 pachy_outlier_GRNvISR_qvals.txt
     170 pachy_outlier_GRNvRED_qvals.txt
    2738 pachy_outlier_REDvISR_qvals.txt
```

#### UniProt identification of all alternatively fixed and outlier SNPs

Obtain the UniProt ID and Name by uploading `uniprot_ids_b4c.txt` on the [UniProt](https://www.uniprot.org/uploadlists/) portal, and save the outcome as `tsv` file (with UniProt `Entry` as the first column, `Entry name` as second column, `Gene names` as the third column, and `Protein names` as fourth column).

```shell
$ head -2 uniprot_names.tab | cut -f1-4
Entry	Entry name	Gene names	Protein names
A0A075F932	SYT1_ANSCY	SYT1	Synaptotagmin-1 (Synaptotagmin I) (SytI)
$ join -1 2 -2 1 -o 1.1,1.2,2.3,2.4 -t $'\t' \
		<(sort -k2 pachy_gene2uniprot_b4c.txt) \
		<(sort uniprot_names.tab) \
    > pachy_gene2uniprotnames_b4c.txt
```

Create summary table ([pachy_divergent_snps.csv](outfiles/pachy_divergent_snps.csv); in `csv` format) of all alternatively fixed and outlier SNPs (using [pachyseris_summarize_divergent_snps.py](scripts/pachyseris_summarize_divergent_snps.py) script):

```shell
$ python3 pachyseris_summarize_divergent_snps.py
```

Create a subset table for only the alternatively fixed SNPs:

```shell
$ grep fix pachy_divergent_snps.csv > pachy_fixed_snps.csv
```

Output only the missense variants of the alternatively fixed SNPs in terms of the 3 main lineages:

```shell
$ grep -E 'fix_BLUvsGRN|fix_REDvsBLU|fix_GRNvsRED' pachy_divergent_snps.csv | grep missense | cut -d, -f3,8,10,11 | tr ',' '\t' | sed 's/fix_//g'  | sed 's/ISRvsRED//' | sed 's/ISRvsBLU//' | sed 's/ISRvsGRN//' | cut -d, -f1,2,3,8,10,11 | tr ',' '\t' | sed 's/fix_//g'  | sed 's/ISRvsRED//' | sed 's/ISRvsBLU//' | sed 's/ISRvsGRN//' # output manually formatted for readability
REDvsBLU;						pspe_0.1.m1.16807	#NA	NA
GRNvsRED;BLUvsGRN		pspe_0.1.m1.25527	#Rabep1		GTPase-binding effector protein 1 (Rabaptin-5) (Rabaptin-5alpha)
REDvsBLU;GRNvsRED		pspe_0.1.m1.27292	#Hecw2		ubiquitin-protein ligase HECW2 (EC 2.3.2.26) (HECT
GRNvsRED						pspe_0.1.m1.44145	#ZSWIM8		Zinc finger SWIM domain-containing protein 8
;REDvsBLU;BLUvsGRN	pspe_0.1.m1.48994	#Yap1 		Transcriptional coactivator YAP1 (Yes-associated protein 1) (Protein yorkie homolog) (Yes-associated protein YAP65 homolog)
REDvsBLU;						pspe_0.1.m1.51298	#dlx1a 		Homeobox protein Dlx1a (DLX-1) (Distal-less homeobox gene 1a)
REDvsBLU;GRNvsRED		pspe_0.1.m1.56819	#FGFR2		Fibroblast growth factor receptor 2 (FGFR-2) (EC 2.7.10.1) (PFR2)
REDvsBLU;BLUvsGRN		pspe_0.1.m1.59271	#Pappa		Pappalysin-1 (EC 3.4.24.79) (Insulin-like growth factor-dependent IGF-binding protein 4 protease) (IGF-dependent IGFBP-4 protease) (IGFBP-4ase) (Pregnancy-associated plasma protein A) (PAPP-A)
```

| Comparison     | Pachyseris gene   | UniProt | Gene   |
| -------------- | ----------------- | ------- | ------ |
| GRN vs BLU/RED | pspe_0.1.m1.25527 | O35551  | Rabep1 |
| BLU vs GRN/RED | pspe_0.1.m1.48994 | Q2EJA0  | Yap1   |
|                | pspe_0.1.m1.59271 | Q8R4K8  | Pappa  |
| RED vs GRN/BLU | pspe_0.1.m1.27292 | Q6I6G8  | Hecw2  |
|                | pspe_0.1.m1.56819 | Q91286  | Fgfr2  |
| GRN vs RED     | pspe_0.1.m1.44145 | A7E305  | Zswim8 |
| RED vs BLU     | pspe_0.1.m1.51298 | Q98875  | Dlx1a  |
|                | pspe_0.1.m1.16807 |         |        |

Output all the alternatively-fixed missense variants:

```shell
$ grep 'fix' pachy_divergent_snps.csv | grep missense | wc -l
24
$ grep 'fix' pachy_divergent_snps.csv | grep missense | cut -d, -f3,8,10,11 | tr ',' '\t'  # output manually formatted for readability
fix_ISRvsRED								pspe_0.1.m1.2147	#Gpr157	G-protein coupled receptor 157
fix_ISRvsRED								pspe_0.1.m1.19412	#Mfsd6 	Major facilitator superfamily domain-containing protein 6 (Macrophage MHC class I receptor 2)
fix_ISRvsGRN								pspe_0.1.m1.21556	#Prkn	E3 ubiquitin-protein ligase parkin (EC 2.3.2.31) (Parkin RBR E3 ubiquitin-protein ligase)
fix_ISRvsRED								pspe_0.1.m1.41119	#BRCA1	Breast cancer type 1 susceptibility protein homolog (EC 2.3.2.27) (RING-type E3 ubiquitin transferase BRCA1)
fix_ISRvsRED								pspe_0.1.m1.42452	#SACS	Sacsin (DnaJ homolog subfamily C member 29) (DNAJC29)
fix_GRNvsRED								pspe_0.1.m1.44145	#ZSWIM8	Zinc finger SWIM domain-containing protein 8
fix_ISRvsRED								pspe_0.1.m1.56705	#Pxdn Kiaa0230	Peroxidasin homolog (EC 1.11.1.7)
fix_ISRvsRED	pspe_0.1.m1.58064	#DICER1	Endoribonuclease Dicer (EC 3.1.26.3)
fix_ISRvsBLU;fix_ISRvsRED		pspe_0.1.m1.15762	#ALKBH7	Alpha-ketoglutarate-dependent dioxygenase alkB homolog 7
fix_REDvsBLU;fix_ISRvsRED		pspe_0.1.m1.16807	#NA	NA
fix_ISRvsBLU;fix_ISRvsRED		pspe_0.1.m1.33235	#TMEM143	Transmembrane protein 143
fix_ISRvsGRN;fix_ISRvsRED		pspe_0.1.m1.34598	#mybL Myb-like protein L
fix_ISRvsGRN;fix_ISRvsRED		pspe_0.1.m1.46418	#MYO10 Unconventional myosin-X (Unconventional myosin-10)
fix_REDvsBLU;fix_ISRvsRED		pspe_0.1.m1.51298	#dlx1a	Homeobox protein Dlx1a (DLX-1) (Distal-less homeobox gene 1a)
fix_ISRvsGRN;fix_ISRvsRED		pspe_0.1.m1.55671	#MYO15A	Unconventional myosin-XV (Unconventional myosin-15)
fix_ISRvsGRN;fix_ISRvsRED		pspe_0.1.m1.63692	#Ipo9 Importin-9 (Imp9) (Importin-9a) (Imp9a) (Importin-9b) (Imp9b) (Ran-binding protein 9) (RanBP9)
fix_ISRvsGRN;fix_ISRvsBLU;fix_ISRvsRED							pspe_0.1.m1.51059	#NA	NA
fix_REDvsBLU;fix_ISRvsRED;fix_GRNvsRED							pspe_0.1.m1.56819	#FGFR2	Fibroblast growth factor receptor 2 (FGFR-2) (EC 2.7.10.1) (PFR2)
fix_REDvsBLU;fix_ISRvsBLU;fix_BLUvsGRN							pspe_0.1.m1.59271	#Pappa	Pappalysin-1 (EC 3.4.24.79) (Insulin-like growth factor-dependent IGF-binding protein 4 protease) (IGF-dependent IGFBP-4 protease) (IGFBP-4ase) (Pregnancy-associated plasma protein A) (PAPP-A)
fix_REDvsBLU;fix_ISRvsRED;fix_GRNvsRED	pspe_0.1.m1.27292	Hecw2 #Kiaa1301 Nedl2	E3 ubiquitin-protein ligase HECW2 (EC 2.3.2.26) (HECT
fix_ISRvsBLU;fix_ISRvsRED;fix_GRNvsRED;fix_BLUvsGRN	pspe_0.1.m1.25527	#Rabep1 Rab5ep Rabpt5 Rabpt5a	Rab GTPase-binding effector protein 1 (Rabaptin-5) (Rabaptin-5alpha)
fix_ISRvsGRN;fix_REDvsBLU;fix_ISRvsRED;fix_BLUvsGRN	pspe_0.1.m1.48994	#Yap1 Yap Yap65	Transcriptional coactivator YAP1 (Yes-associated protein 1) (Protein yorkie homolog) (Yes-associated protein YAP65 homolog)
```

### B4e - GO-term enrichment analysis

#### GO-term enrichment analysis for alternatively fixed and outlier SNPs

Assess from the  [SwissProt blast](../genome/pspe.v0.12.blastSprot.outfmt6.gz) file how many of the *Pachyseris speciosa* genes have an annotation, and extract all matching UniProtID numbers:

```shell
$ cut -f1 pspe.v0.12.blastSprot.outfmt6 | sort | uniq | wc -l
   24693 # genes with SwissProt annotation
$ cut -d'|' -f2 pspe.v0.12.blastSprot.outfmt6 | sort | uniq > uniprot_ids_b4c.txt
```

Obtain the Gene Ontology (GO) IDs for each UniProt entry by uploading `uniprot_ids_b4c.txt` on the [UniProt](https://www.uniprot.org/uploadlists/) portal, and save the outcome as `tab`/`tsv` file (with UniProt `Entry` as the first column, `Entry name` as second column, `Gene Ontology IDs` as the third column, and `Gene ontology (GO` as fourth column). Then, extract just the SwissPort ID (first column) and GO (fourth column) columns. 

```shell
$ head -2 uniprot.tab | cut -f1-4
Entry    Entry name    Gene ontology IDs    Gene ontology (GO)
A0A075F932    SYT1_ANSCY    GO:0005509; GO:0005543; GO:0005544; GO:0005737; GO:0005886; GO:0007269; GO:0008021; GO:0016021; GO:0017158; GO:0030054; GO:0030672; GO:0042584; GO:0046883; GO:0048609; GO:0051592; GO:1903305; GO:1903861    cell junction [GO:0030054]; chromaffin granule membrane [GO:0042584]; cytoplasm [GO:0005737]; integral component of membrane [GO:0016021]; plasma membrane [GO:0005886]; synaptic vesicle [GO:0008021]; synaptic vesicle membrane [GO:0030672]; calcium ion binding [GO:0005509]; calcium-dependent phospholipid binding [GO:0005544]; phospholipid binding [GO:0005543]; multicellular organismal reproductive process [GO:0048609]; neurotransmitter secretion [GO:0007269]; positive regulation of dendrite extension [GO:1903861]; regulation of calcium ion-dependent exocytosis [GO:0017158]; regulation of hormone secretion [GO:0046883]; regulation of regulated secretory pathway [GO:1903305]; response to calcium ion [GO:0051592]
$ tail -n +2 uniprot.tab | cut -f1,4 > pachy_uniprot2go_b4c.txt
```

Create a file with the Pachyseris scaffold names and their corresponding SwissProt ID, and then one with the corresponding GOs (using [goterms_from_uniprot_blast.py](https://github.com/pimbongaerts/radseq/blob/master/goterms_from_uniprot_blast.py) script):

```shell
$ sed 's/|/./g' pspe.v0.12.blastSprot.outfmt6 \
      | cut -f1-2 | awk -F'.' '{print $1 "." $2 "." $3 "." $4 "\t" $6 }' \
      > pachy_gene2uniprot_b4c.txt
$ goterms_from_uniprot_blast.py pachy_gene2uniprot_b4c.txt \
                                pachy_uniprot2go_b4c.txt \
      > pachy_gene2go_b4c.txt
```

Create a list of GO IDs/terms and the genes assigned to that GO term using (using [create_gowinda_file.py](https://github.com/pimbongaerts/radseq/blob/master/create_gowinda_file.py) script), and create annoptation file from `.gff` to `.gtf` using script from Gowinda:

```shell
$ create_gowinda_file.py pachy_gene2go_b4c.txt > go2gene_ids_b4c.txt
$ python3 Gff2Gtf.py --input \
                     $GENOME_REF_PATH/pspe_0.12.maker_post_002.genes.gff3 \
      > pspe_0.12_genes.gtf
```

Create a file with all the SNPs and files with the subsets of different alternatively fixed SNPs:

```shell
# All SNPs in VCF
$ grep -v "#" pachy_annotated_b4a.vcf | cut -f1,2 > pachy_all_snps_b4c.txt
# Alternatively fixed SNPs
$ for i in pachy_fixed_*.csv; do \
		grep "AFD_OUTLIER" $i | awk -F',' '{print $1 "\t" $2}' > gowinda_$i;
  done
$ rename  's/csv/txt/' gowinda_pachy_fixed_*.csv
# Outlier SNPs (identified through pcadapt)
$ for i in pachy_outlier_*_qvals.txt; do \
		cut -f1,2 -d' ' $i | tr ' ' \\t > gowinda_$i;
	done
$ rename  's/_qvals//' gowinda_pachy_outlier_*.txt	
```

Run Gowinda for all three lineage comparisons:

```shell
$ for i in gowinda_pachy_*.txt; do \
		java -Xmx4g -jar ~/Gowinda-1.12.jar --snp-file pachy_all_snps_b4c.txt \
                                        --candidate-snp-file $i \
                                        --gene-set-file go2gene_ids_b4c.txt \
                                        --annotation-file pspe_0.12_genes.gtf \
                                        --simulations 1000000 \
                                        --min-significance 1 \
                                        --gene-definition updownstream1000 \
                                        --threads 32 \
                                        --mode gene \
                                        --min-genes 1 \
                                        --output-file output_$i;
	done
$ head output_gowinda* | cut -f1,5,6,7,8,9
==> output_gowinda_pachy_fixed_BLUvsGRN.txt <==
GO:0048368      0.4374086667    1       1       5       lateral mesoderm development
GO:0071149      0.4374086667    1       1       1       TEAD-2-YAP complex
GO:0071148      0.4374086667    1       1       3       TEAD-1-YAP complex
GO:0003015      0.4374086667    1       1       16      heart process
GO:2000737      0.4374086667    1       1       6       negative regulation of stem cell differentiation
GO:0072091      0.4374086667    1       1       6       regulation of stem cell proliferation
GO:0010837      0.4374086667    1       1       4       regulation of keratinocyte proliferation
GO:0061026      0.4374086667    1       1       2       cardiac muscle tissue regeneration
GO:0033148      0.4374086667    1       1       3       positive regulation of intracellular estrogen receptor signaling pathway
GO:0090263      0.4374086667    2       16      102     positive regulation of canonical Wnt signaling pathway

==> output_gowinda_pachy_fixed_GRNvsRED.txt <==
GO:0048557      0.3612798364    2       4       14      embryonic digestive tract morphogenesis
GO:0048264      0.3612798364    1       2       15      determination of ventral identity
GO:0001885      0.3612798364    1       2       10      endothelial cell development
GO:0035124      0.3612798364    1       2       6       embryonic caudal fin morphogenesis
GO:0034605      0.3612798364    2       9       61      cellular response to heat
GO:0030324      0.3612798364    2       11      74      lung development
GO:0046622      0.3612798364    1       1       4       positive regulation of organ growth
GO:0045629      0.3612798364    1       1       3       negative regulation of T-helper 2 cell differentiation
GO:0045627      0.3612798364    1       1       3       positive regulation of T-helper 1 cell differentiation
GO:0048513      0.3612798364    1       1       14      animal organ development

==> output_gowinda_pachy_fixed_ISRvsBLU.txt <==
GO:0060548      0.4032570000    2       7       52      negative regulation of cell death
GO:1990837      0.4783402000    1       1       12      sequence-specific double-stranded DNA binding
GO:0060317      0.4783402000    1       1       8       cardiac epithelial to mesenchymal transition
GO:0090403      0.4783402000    1       1       2       oxidative stress-induced premature senescence
GO:0043616      0.4783402000    1       1       7       keratinocyte proliferation
GO:0014068      0.5991239231    1       2       64      positive regulation of phosphatidylinositol 3-kinase signaling
GO:0032354      0.5991239231    1       1       2       response to follicle-stimulating hormone
GO:0007565      0.5991239231    1       1       23      female pregnancy
GO:1990089      0.5991239231    1       1       6       response to nerve growth factor
GO:1903818      0.5991239231    1       1       17      positive regulation of voltage-gated potassium channel activity

==> output_gowinda_pachy_fixed_ISRvsGRN.txt <==
GO:0060548      0.2738507000    2       7       52      negative regulation of cell death
GO:0005634      0.2738507000    15      536     5039    nucleus
GO:0030216      0.2738507000    2       5       23      keratinocyte differentiation
GO:0060317      0.2738507000    1       1       8       cardiac epithelial to mesenchymal transition
GO:0090403      0.2738507000    1       1       2       oxidative stress-induced premature senescence
GO:0043616      0.2738507000    1       1       7       keratinocyte proliferation
GO:0003774      0.2738507000    2       8       43      motor activity
GO:0016459      0.2738507000    2       7       34      myosin complex
GO:0048368      0.2738507000    1       1       5       lateral mesoderm development
GO:0071149      0.2738507000    1       1       1       TEAD-2-YAP complex

==> output_gowinda_pachy_fixed_ISRvsRED.txt <==
GO:0060449      0.3871710000    2       2       3       bud elongation involved in lung branching
GO:0045165      0.5121873071    2       3       49      cell fate commitment
GO:0060045      0.5121873071    2       3       24      positive regulation of cardiac muscle cell proliferation
GO:0043565      0.5121873071    5       34      437     sequence-specific DNA binding
GO:0042340      0.5121873071    2       2       14      keratan sulfate catabolic process
GO:0031090      0.5121873071    3       12      120     organelle membrane
GO:0048557      0.5121873071    2       4       14      embryonic digestive tract morphogenesis
GO:1902871      0.5121873071    1       1       1       positive regulation of amacrine cell differentiation
GO:0046533      0.5121873071    1       1       4       negative regulation of photoreceptor cell differentiation
GO:0070852      0.5121873071    2       5       23      cell body fiber

==> output_gowinda_pachy_fixed_REDvsBLU.txt <==
GO:0060449      0.0519570000    2       2       3       bud elongation involved in lung branching
GO:0060045      0.1180090000    2       3       24      positive regulation of cardiac muscle cell proliferation
GO:0048557      0.2329070000    2       4       14      embryonic digestive tract morphogenesis
GO:0008237      0.2473494286    3       18      78      metallopeptidase activity
GO:0047555      0.2473494286    2       5       15      3',5'-cyclic-GMP phosphodiesterase activity
GO:1902871      0.2473494286    1       1       1       positive regulation of amacrine cell differentiation
GO:0046533      0.2473494286    1       1       4       negative regulation of photoreceptor cell differentiation
GO:0048264      0.2599343214    1       2       15      determination of ventral identity
GO:0001885      0.2599343214    1       2       10      endothelial cell development
GO:0035124      0.2599343214    1       2       6       embryonic caudal fin morphogenesis

==> output_gowinda_pachy_outlier_BLUvGRN.txt <==
GO:0060548      0.0632200000    3       7       52      negative regulation of cell death
GO:0043202      0.0815120000    3       5       68      lysosomal lumen
GO:0042803      0.2885770000    9       99      1041    protein homodimerization activity
GO:0006027      0.3208397500    2       3       25      glycosaminoglycan catabolic process
GO:0070064      0.5927258000    2       5       13      proline-rich region binding
GO:0007288      0.6787820476    2       3       15      sperm axoneme assembly
GO:0050852      0.6787820476    2       7       95      T cell receptor signaling pathway
GO:0030216      0.6787820476    2       5       23      keratinocyte differentiation
GO:0060317      0.6787820476    1       1       8       cardiac epithelial to mesenchymal transition
GO:0043616      0.6787820476    1       1       7       keratinocyte proliferation

==> output_gowinda_pachy_outlier_BLUvISR.txt <==
GO:0009636      0.8450000000    7       8       72      response to toxic substance
GO:0006694      0.8823785750    6       7       56      steroid biosynthetic process
GO:0016514      0.8823785750    4       4       7       SWI/SNF complex
GO:0060548      0.8823785750    5       7       52      negative regulation of cell death
GO:0033588      0.8823785750    3       3       8       Elongator holoenzyme complex
GO:0019228      0.8823785750    6       7       70      neuronal action potential
GO:0030335      0.8823785750    12      22      260     positive regulation of cell migration
GO:0050807      0.8823785750    3       3       21      regulation of synapse organization
GO:0043005      0.8823785750    28      55      505     neuron projection
GO:0061512      0.8823785750    4       4       18      protein localization to cilium

==> output_gowinda_pachy_outlier_BLUvRED.txt <==
GO:0000792      0.1262650000    3       4       35      heterochromatin
GO:0005242      0.1262650000    2       2       8       inward rectifier potassium channel activity
GO:0060548      0.1547957500    3       7       52      negative regulation of cell death
GO:0099507      0.1547957500    2       2       9       ligand-gated ion channel activity involved in regulation of presynaptic membrane potential
GO:0060449      0.1847451111    2       2       3       bud elongation involved in lung branching
GO:0043202      0.1847451111    3       5       68      lysosomal lumen
GO:0071346      0.1847451111    3       5       54      cellular response to interferon-gamma
GO:0046654      0.1847451111    2       2       6       tetrahydrofolate biosynthetic process
GO:0046452      0.1847451111    2       2       4       dihydrofolate metabolic process
GO:0006029      0.2314382727    2       2       7       proteoglycan metabolic process

==> output_gowinda_pachy_outlier_GRNvISR.txt <==
GO:0070062      0.7575288824    45      90      810     extracellular exosome
GO:0016514      0.7575288824    4       4       7       SWI/SNF complex
GO:0033588      0.7575288824    3       3       8       Elongator holoenzyme complex
GO:0001228      0.7575288824    8       16      150     DNA-binding transcription activator activity, RNA polymerase II-specific
GO:0045652      0.7575288824    2       2       13      regulation of megakaryocyte differentiation
GO:0043005      0.7575288824    26      55      505     neuron projection
GO:0008022      0.7575288824    12      23      151     protein C-terminus binding
GO:0098609      0.7575288824    10      13      95      cell-cell adhesion
GO:0005635      0.7575288824    9       14      170     nuclear envelope
GO:0003774      0.7575288824    6       8       43      motor activity

==> output_gowinda_pachy_outlier_GRNvRED.txt <==
GO:0060484      0.5821490156    2       2       9       lung-associated mesenchyme development
GO:0030433      0.5821490156    3       10      64      ubiquitin-dependent ERAD pathway
GO:0044615      0.5821490156    2       2       8       nuclear pore nuclear basket
GO:0048557      0.5821490156    2       4       14      embryonic digestive tract morphogenesis
GO:1902871      0.5821490156    1       1       1       positive regulation of amacrine cell differentiation
GO:0046533      0.5821490156    1       1       4       negative regulation of photoreceptor cell differentiation
GO:0000792      0.5821490156    2       4       35      heterochromatin
GO:0016779      0.5821490156    2       3       10      nucleotidyltransferase activity
GO:0032447      0.5821490156    1       1       5       protein urmylation
GO:0034227      0.5821490156    1       1       5       tRNA thio-modification

==> output_gowinda_pachy_outlier_REDvISR.txt <==
GO:0043565      0.5020245000    20      34      437     sequence-specific DNA binding
GO:0002098      0.5020245000    4       4       13      tRNA wobble uridine modification
GO:0005249      0.5234854000    10      12      106     voltage-gated potassium channel activity
GO:0006813      0.5234854000    9       11      91      potassium ion transport
GO:0003774      0.5234854000    7       8       43      motor activity
GO:0032482      0.9032960719    7       9       68      Rab protein signal transduction
GO:0016514      0.9032960719    4       4       7       SWI/SNF complex
GO:0009880      0.9032960719    5       5       32      embryonic pattern specification
GO:0006694      0.9032960719    6       7       56      steroid biosynthetic process
GO:0006635      0.9032960719    5       6       38      fatty acid beta-oxidation
```

### B4e - Admixture

Subset all admixed samples (based on stringent assignment) and subset SNPs that are close to being alternatively fixed (AFD of >=0.95):

```shell
$ grep UNA pachy_popfile_clusters_stringent_b3c.txt | cut -f1 \
      > pachy_admixed_samples_b4d.txt
$ cat pachy_fixed_* | sort | uniq > pachy_fixed_all_b4d.txt
$ vcftools --gzvcf pachy_b1.vcf.gz \
           --keep pachy_admixed_samples_b4d.txt \
           --positions pachy_fixed_all_b4d.txt \
           --recode --stdout > pachy_admixed_b4f.vcf
$ grep "After" out.log
After filtering, kept 55 out of 516 Individuals
After filtering, kept 146 out of a possible 94438 Sites
```

Create separate assignment file for admixed samples (for visualization):

```shell
$ grep -Ff pachy_admixed_samples_b4d.txt \
                     pachy_structure_k6_b3c.csv \
      > pachy_structure_k6_admixed_b4d.csv
```

The seven “admixed” corals observed for Papua New Guinea had assignments to three or more clusters, likely reflecting clustering artefacts (e.g. due to admixture with an unsampled population/lineage), except for maybe two of the samples (PSPVDH1219 and PSPGSH1113). Seven samples from the Great Barrier Reef showed admixture of the “blue” lineage with the “dark green” cluster present only in Okinawa, which again likely represents a clustering artefact.

Assess potential F1 hybrid sample:

```shell
$ vcftools --vcf pachy_admixed_b4f.vcf \
           --indv PSGRBH8219 \
           --positions pachy_fixed_BLUvsGRN_b4c.txt \
           --recode --stdout \
      > pachy_8219_b4d.vcf
$ vcftools --vcf pachy_8219_b4d.vcf --hardy
$ paste <(grep -v "CHR" out.hwe| wc -l) \
        <(grep "0/0/0" out.hwe| wc -l) \
        <(grep "0/1/0" out.hwe| wc -l) \
        <(grep -E "1/0/0|0/0/1" out.hwe | wc -l)
# total            missing        heterozyg.   homozyg.
      57          18          35           4            
```

## B5 - Intra-lineage structuring

Split dataset in 3 clusters for further separate analyses (based on lenient admixture assignment; 0.8), and retain only SNPs genotyped for at least 80% of samples:

```shell
$ for i in BLU GRN RED; do \
      vcftools --gzvcf pachy_b1.vcf.gz \
                --keep <(grep $i pachy_popfile_clusters_b3c.txt) \
                --mac 1 --max-missing 0.8 --recode --stdout > pachy_$i\_b5a_temp1.vcf; 
  done
```

#### Overall dataset

For each dataset, filter out any sites that have an observed heterygosity of >0.5 (to filter out potential paralogs):

```shell
$ for i in BLU GRN RED; do \
        vcftools --vcf pachy_$i\_b5a_temp1.vcf \
                 --exclude-positions <(vcftools --vcf pachy_$i\_b5a_temp1.vcf \
                                                --hardy --stdout | tail +2 \
                     | awk -F'[\t/]' '($4 / ($3 + $4 +$5)) > 0.5 { print $1 "\t" $2 }') \
                 --recode --stdout > pachy_$i\_b5a.vcf; 
    done
$ for i in BLU GRN RED; do grep -v "#" pachy_$i\_b5a.vcf | wc -l; done
    6264
    6307
    6605
```

Create popfile for each dataset:

```shell
$ for i in BLU GRN RED; do \
      popfile_from_vcf.py pachy_$i\_b5a.vcf 3 5 > pachy_popfile_$i\_b5a.txt; \
  done
```

#### "Neutral" dataset

To identify a subset of SNPs that was more representative of neutral genomic diversity, we compiled a separate data set where we removed SNPS that were identified as outliers using pcadapt (Luu et al. 2016) without *a priori* population information and based on q-values with an expected false discovery rate lower than 10%:

```shell
# Pcadapt through the following script:
$ Rscript pachyseris_pcadapt.R
# Merge CHROM and POS from VCF with pcadapt output and filter outliers
$ for i in BLU GRN RED; do \
      paste -d" " <(grep -v "#" pachy_$i\_b5a.vcf | awk '{print $1 " " $2}') \
                  pachy_$i\_pcadapt_B5a.txt \
          | awk '$4 < 0.1' > pachy_$i\_outliers_b5a.txt; 
  done
$ for i in BLU GRN RED; do wc -l pachy_$i\_outliers_b5a.txt; done
      64 pachy_BLU_outliers_b5a.txt
     178 pachy_GRN_outliers_b5a.txt
      29 pachy_RED_outliers_b5a.txt
# Create neutral VCF datasets
$ for i in BLU GRN RED; do \
      vcftools --vcf pachy_$i\_b5a.vcf \
               --exclude-positions <(awk '{print $1 "\t" $2}' pachy_$i\_outliers_b5a.txt) \
               --recode --stdout \
          > pachy_$i\_neutral_b5a.vcf; \
  done
$ for i in BLU GRN RED; do grep -v "#" pachy_$i\_neutral_b5a.vcf | wc -l; done
    6200
    6129
    6576
```

#### Genetic structuring: PCA and DAPC

Conduct principal component analyses (PCA) and discriminant analyses of principal components (DAPC) for the three neutral and overall datasets using the [adegenet](https://github.com/thibautjombart/adegenet/wiki) package in R (using the [pachyseris_dapc_pca_b5b.R](scripts/pachyseris_dapc_pca_b5b.R) script):

```shell
$ Rscript pachyseris_dapc_pca_b5b.R
```

## B6 - Phylogenomic analyses

Select *Pachyseris* samples with the highest numbers of SNPs for each lineage:

```shell
$ vcf_missing_data.py pachy_outgroup_b2.vcf > pachy_outgroup_missing_data_b6.txt
$ join -a1 <(grep -v "INDIVIDUAL" pachy_outgroup_missing_data_b6.txt | sort ) \
           <(sort pachy_popfile_clusters_stringent_b3c.txt) \
           > pachy_outgroup_missing_data_by_lineage_b6.txt
$ for i in GR2 BL2 ISR; do \
		grep $i pachy_outgroup_missing_data_by_lineage_b6.txt  | sort -rn -k3 | head -n 6; \
	done > pachy_phylo_samples_b6.txt
$ for i in RED GRN BLU; do \
		grep PSG pachy_outgroup_missing_data_by_lineage_b6.txt | grep $i | sort -rn -k3 | head -n 6; \
	done >> pachy_phylo_samples_b6.txt
$ grep -E "PR|PSGMDH0861" pachy_outgroup_missing_data_by_lineage_b6.txt >> pachy_phylo_samples_b6.txt
```

 Gather *Pachyseris* samples for ipyrad analysis:

```shell
$ awk '{system("cp " $1 ".fastq.gz ../../pachy_phylo/orig_fq")}' pachy_phylo_samples_b6.txt
# we then manually included samples from Agaricia fragilis, Leptoseris (c.f.) glabra, and Pachyseris c.f. inatessa
```

Use [TrimGalore](https://github.com/FelixKrueger/TrimGalore) to trim Nextera adapters (with overlap of at least 5bp) and low-quality ends (below PHRED score of 20), while discarding reads that become less than 30bp, and never keeping more than the first 100bp (to trim sequences that were done in a 150bp run): 

```shell
$ mkdir trimmed_fq && cd "$_"
$ find ../orig_fq -name "*.fastq.gz" | parallel -j 32 trim_galore --hardtrim5 100 -q 20 --length 30 --stringency 5 --nextera --phred33 {}
$ rename -v 's/.100bp//' *.fq.gz
```

Create ipyrad params file, leave at default settings (apart from removing #8) and run [ipyrad](https://github.com/dereneaton/ipyrad):

```shell
$ ipyrad -n pachy_phylo
$ cut -f1 -d":" params-pachy_phylo.txt
------- ipyrad params file (v.0.9.62)-------------------------------------------
pachy_phylo                    ## [0] [assembly_name]
pachy_phylo/ipyrad 						 ## [1] [project_dir]
                               ## [2] [raw_fastq_path]
                               ## [3] [barcodes_path]
pachy_phylo/trimmed_fq/*.fq.gz ## [4] [sorted_fastq_path]
reference                      ## [5] [assembly_method]
pspe_final_0.12.fasta          ## [6] [reference_sequence]
rad                            ## [7] [datatype]
                               ## [8] [restriction_overhang]
5                              ## [9] [max_low_qual_bases]
33                             ## [10] [phred_Qscore_offset]
6                              ## [11] [mindepth_statistical]
6                              ## [12] [mindepth_majrule]
10000                          ## [13] [maxdepth]
0.85                           ## [14] [clust_threshold]
0                              ## [15] [max_barcode_mismatch]
2                              ## [16] [filter_adapters]
35                             ## [17] [filter_min_trim_len]
2                              ## [18] [max_alleles_consens]
0.05                           ## [19] [max_Ns_consens]
0.05                           ## [20] [max_Hs_consens]
4                              ## [21] [min_samples_locus]
0.2                            ## [22] [max_SNPs_locus]
8                              ## [23] [max_Indels_locus]
0.5                            ## [24] [max_shared_Hs_locus]
0, 0, 0, 0                     ## [25] [trim_reads]
0, 0, 0, 0                     ## [26] [trim_loci]
p, s, l, v                     ## [27] [output_formats]
                               ## [28] [pop_assign_file]
                               ## [29] [reference_as_filter]        
$ ipcluster start -n 32 --daemonize --profile="pachy_ipyrad"; sleep 60
$ ipyrad -p params-pachy_phylo.txt -s1234567 -c 32 --ipcluster pachy_ipyrad
#[...]
$ grep "//" pachy_phylo.loci| wc -l
98043 # number of loci
```

Run [tetrad](https://ipyrad.readthedocs.io/en/latest/API-analysis/cookbook-tetrad.html) analysis through interactive Python3 shell:

```python
>>> import ipyrad.analysis as ipa
>>> import ipyparallel as ipp
>>> converter = ipa.vcf_to_hdf5(name="pachy_tetrad", 
                                data="pachy_phylo.vcf", 
                                ld_block_size=2500)
>>> converter.run()
Indexing VCF to HDF5 database file
VCF: 713097 SNPs; 2344 scaffolds
[####################] 100% 0:07:01 | converting VCF to HDF5
HDF5: 713097 SNPs; 80108 linkage group
>>> ipyclient = ipp.Client(profile="pachy_ipyrad")
>>> print(len(ipyclient), 'cores')
>>> tet = ipa.tetrad(name="pachy_tetrad_ref", 
                     data="./analysis-vcf2hdf5/pachy_tetrad.snps.hdf5", 
                     nquartets=10e9, 
                     nboots=100)
loading snps array [48 taxa x 713097 snps]
max unlinked SNPs per quartet [nloci]: 80108
quartet sampler [full]: 194580 / 194580
>>> tet.run(auto=True)
```

Run RAXML using the`autoMRE` feature (automatic bootstrapping) and `-f a` option (ML search and bootstraps in the same run), and using the `GTRGAMMA` substition model:

```shell
$ raxmlHPC-PTHREADS-SSE3 -f a -T 32 -m GTRGAMMA -N autoMRE -x 12345 -p 54321 \
												 -n pachy_phylo_raxml -s pachy_phylo.phy
#[...]
Alignment has 1731480 distinct alignment patterns
#[...]
Overall execution time for full ML analysis: 676558.600422 secs or 187.932945 hours or 7.830539 days
$ head -n1 pachy_phylo.phy.reduced
48 10504278 # Total number of sites
```
