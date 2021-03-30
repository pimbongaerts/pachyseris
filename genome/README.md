# "Morphological stasis masks ecologically divergent coral species on tropical reefs"

## A - Genome assembly and annotation (PacBio)

We assembled and annotated a *de novo* reference genome using PacBio sequencing (~95x coverage) of a sperm sample from a *Pachyseris speciosa* colony from Orpheus Island on the Great Barrier Reef. 

## A1 - Genome assembly files

| file                                                         | checksum                         | description            |
| ------------------------------------------------------------ | -------------------------------- | ---------------------- |
| [pspe_final_0.12.fasta](pspe_final_0.12.fasta.gz)            | 16745e6d6ca50f29c7bd7755740ed88a | genome scaffolds       |
| [pspe_0.12.maker_post_002.genes.gff3](pspe_0.12.maker_post_002.genes.gff3.gz) | ba96f872693787865886f90a7782171b | gene models (GFF)      |
| [pspe_0.12.maker_post_002.proteins.fasta](pspe_0.12.maker_post_002.proteins.fasta.gz) | ff18770b4f293a9832f04441d541d587 | gene models (proteins) |
| [pspe_0.12.maker_post_002.transcripts.fasta](pspe_0.12.maker_post_002.transcripts.fasta.gz) | 115ed4bb5326ca96193be43b48f72b42 | transcripts            |
| [pspe.v0.12.blastSprot.outfmt6](pspe.v0.12.blastSprot.outfmt6.gz) |                                  |                        |

### A2 - Overview of DNA and RNA libraries

| Sample                              | Genomic DNA      |                     | RNA                 |
| ----------------------------------- | ---------------- | ------------------- | ------------------- |
| Sequencing platform                 | PacBio SMRT Cell | Illumina HiSeq 2500 | Illumina HiSeq 2500 |
| Number of Reads ( x 10<sup>6</sup>) | 8.03             | 203.6               | 238.7               |
| Read Length                         | 1 - 23 Kb*       | 2 x 250 bp          | 2 x 100 bp          |
| Insert size (bp)                    |                  | 450                 |                     |
| Total Bases ( x 10<sup>9</sup>)     | 84               | 101.8               | 47.8                |
| Genome  coverage†                   | 94.8             | 91.9                |                     |

\* = 5-95% range. The longest reads extend to 50 Kb. † = based on estimated genome size of 886.1 Mb

### A3 - **Estimation of genome size and heterozygosity rate**

| Program     | k-mer size | Genome size  (Mb) | heterozygosity  rate (%) |
| ----------- | ---------- | ----------------- | ------------------------ |
| sga.preqc   | 31         | 886.1             |                          |
| GenomeScope | 31         | 749.6             | 1.4                      |

### 