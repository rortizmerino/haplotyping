---
layout: default
title: Genotyping
nav_order: 4
has_children: true
permalink: /docs/genotyping
---

# Genotyping

Here you will find general info about genotyping methods. All these methods are based on Illumina DNA sequence reads and how they map, or align, into a particular genome reference. In these examples, the [Burrows-Wheeler Alignment (bwa) tool](https://bio-bwa.sourceforge.net/bwa.shtml) is used to get the alignments, [samtools](http://www.htslib.org/doc/samtools.html) is used to process them, and different methods within the [Genome Analysis Toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4) are used to process them. These methods are conviniently installed within the `vc` conda environment and their use is wrapped within the perl script [`varCaller_v4.pl`](https://github.com/rortizmerino/haplotyping/blob/main/scripts/varCaller_v4.pl). 

Once the `vc` conda environment is loaded

```bash
source activate vc
```

`varCaller_v4.pl` runs as follows:

```bash
varCaller_v4.pl <species> <reference> <path lib> <path run> <threads> <run>

species:   identifier for logfile, readgroups, and final joint vcf (e.g. Strain ID or CBS no)
reference: fasta file to create indexes and use for mapping
path lib:  path to fastq files
path run:  path to run directory
threads:   to be used for multithreaded steps
run:       0, print commands to log and exit; 1, run commands and print to log
```

One of the steps within `varCaller_v4.pl` runs `gatk DepthOfCoverage` which generates a table with the read coverage per genomic position. The script [`getCoveragePerNkb.pl`](https://github.com/rortizmerino/haplotyping/blob/main/scripts/getCoveragePerNkb.pl) takes those values into N kb non-overlapping windows to get more informative summarised values. `getCoveragePerNkb.pl` runs as follows:

```bash
getCoveragePerNkb.pl <reference> <path> <window> <version>

reference: fasta file to create indexes and use for mapping
path:      path to coverage files
window:    window size (e.g. 10000 for 10 kb)
version:   e.g. Scer_v0
```

Afterwards, the Single Nucleotide Variants (SNVs) must be filtered. In this case, [`fltrVCF.SNVs.pl`](https://github.com/rortizmerino/haplotyping/blob/main/scripts/fltrVCF.SNVs.pl) is used for that purpose by filtering out indels (2 or more nts), requiring absolute coverage >5 reads per variant called, variant coverage > 10% of the average and a Genome Quality GQ > 20. `fltrVCF.SNVs.pl` runs as follows:

```bash
fltrVCF.SNVs.pl <tag> <path>
tag:  to distinguish project, used in input files e.g. Scer_v0
path: to vcf file
```

A modified version can filter and print heterozygous positions only, this is [`fltrVCF.Hz.pl`](https://github.com/rortizmerino/haplotyping/blob/main/scripts/fltrVCF.Hz.pl) and runs just as `fltrVCF.SNVs.pl` above.

Output files from the codes above can be used to make predefined plots using R. Once the `r_env` conda environment is loaded

```bash
source activate r_env
```

Different Copy Number Variation (CNV) plots can be obtained depending on the type of reference genome used during the genotyping. In this project genome references used for CNV determination where mostly either two concatenated genomes, or the whole concatenated reference genomes for the Saccharomyces genus. [`plotAllSacchCov_bin_full.R`](https://github.com/rortizmerino/haplotyping/blob/main/scripts/plotAllSacchCov_bin_full.R) uses the whole concatenated reference genomes for the Saccharomyces genus, from there [`plotAllSacchCov_bin2.R`](https://github.com/rortizmerino/haplotyping/blob/main/scripts/plotAllSacchCov_bin2.R) can extract data relative to only 2 species (*S. cerevisiae* and *S. eubayanus* hardcoded), or the same information but only when 2 concatenated genomes where used as refererence [`plot2spCov_bin.R`](https://github.com/rortizmerino/haplotyping/blob/main/scripts/plot2spCov_bin.R). All of them run similarly:

```bash
Rscript --vanilla <script>.R <tag> <path> <winlen>
```

Another predefined way to analysed filtered SNVs consists of Allele Frequency plots. [`plotVCF.SNVs_n_cov_bin_32_i.R`](https://github.com/rortizmerino/haplotyping/blob/main/scripts/plotVCF.SNVs_n_cov_bin_32_i.R) produces these plots only for cases where *S. cerevisiae* and *S. eubayanus* were used, please be careful with hardcoded values. `plotVCF.SNVs_n_cov_bin_32_i.R` runs as follows:

```bash
Rscript --vanilla plotVCF.SNVs_n_cov_bin_32_i.R <tag> <path> <length>

example:  Rscript --vanilla plotVCF.SNVs_n_cov_bin_32_i.R Scen_Seub CPOPY0_Scen_Seub 1kb

optional: <first scaff> <last scaff> <lngstlen> <log|bin>
example:  Rscript --vanilla plotVCF.SNVs_n_cov_bin.R S1_np_ph_v0 CPOPY0_S1_np_ph_v0 1kb 1 100 1600000
```
