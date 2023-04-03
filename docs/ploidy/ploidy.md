---
layout: default
title: Ploidy
nav_order: 3
has_children: true
permalink: /docs/ploidy
---

# Genome complexity and ploidy determination

One of the first steps to start analysing a genome by its sequence consists of having an overview of its complexity before anything else. Here you will find general info about determining genome complexity and ploidy using Illumina DNA sequencing reads. 

## Genomescope

Genomescope is specially useful to generate an overview of complex polyploid genomes for which no reference exists yet. With a kmer-based genome profiling using Illumina DNA sequencing reads, it is able to estimate heterozygosity and genome size. These 2 parameters are required for most downstream analyses.

Yeast strains come in several ploidy levels. Before trying to estimate how complex a novel industrial isolate can be, it is useful to see the results of these methods using DNA sequencing read libraries obtained from strains with known ploidy. Some diploid(1N), triploid(2N), tetraploid(3N), and aneuploid datasets were therefore first analysed to generate standards. These datasets were obtained from a [2018 research article](https://doi.org/10.3389/fgene.2018.00123).

First, the datasets are downloaded as follows:

```bash
source activate sratk

# CBS7837 (2N)
fasterq-dump --split-files --skip-technical --verbose SRR3265396

# CBS2919 (3N)
fasterq-dump --split-files --skip-technical --verbose SRR3265389	

# CBS9564 (4N)
fasterq-dump --split-files --skip-technical --verbose SRR3265401

# YJM1098 (mostly 4N, mainly aneuploid)
fasterq-dump --split-files --skip-technical --verbose SRR3265463	
```

Once the files are downladed, they can be used with the [`gscoper.sh`](https://github.com/rortizmerino/haplotyping/blob/main/scripts/gscoper.sh) script to 1) run KMC to obtain a kmer distribution 2) run genomescope to get a profile and 3) run smudgeplot to get a ploidy profile. Below a sample command wich can be modified and used to run `gscoper.sh`

```bash
bash scripts/gscoper.sh -s <strain> -r ploidy/<strain> -p data/<strain>
```
