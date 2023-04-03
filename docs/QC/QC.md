---
layout: default
title: QC
nav_order: 5
has_children: true
permalink: /docs/QC
---

# Quality Control

Here you will find general info about Quality Control steps implemented for PacBio HiFi DNA sequencing reads. In general, these stepes can be considered as 1) diagnostic plots, 2) longest read subsetting, 3) random subsampling. 


LongQC was mainly used to generate diagnostic plots of read length distributions. To run it, a conda environment has first to be loaded, and then run as follows:

```bash
source activate LongQC

python <script location>/longQC.py sampleqc -x pb-sequel -o <output directory> --ncpu <number; normally 8> <input>.fastq
```

Warning: this software does not continue running if executed in the background with `Ctrl+Z` then `bg`, or even with `nohup`. Only method tested and working was with `screen`.

Filtlong was used to filter the longest reads in a library until they satisfied a coverage estimate. For this purpose, the genome size must first be calculated with genomescope or any other method. This calculation is the passed as an argument to run the program as follows:

```bash
<executable location>/filtlong --min_length 7000 --target_bases 960000000 <input>.fastq | gzip > <output>.fastq.gz
```

Examples to chose `--target_bases` are: 1) assuming 12Mb Scer ref genome (12 000 000 * 80 = 960 000 000), 2) assuming 27Mb genomescope (27 000 000 * 80 = 2 160 000 000), both of them meant to reach 80x coverage.

Rasusa was used similarly to Filtlong, to reach a certain coverage, altough it does a random subsampling to avoid throwing away small reads which are potentially informative. This method ended up being chosen over Filtlong and can be run as follows:

```bash
rasusa --genome-size 27m --coverage 80 --input <>.fastq --seed 123 --output <>_80x.fastq
#or 
rasusa --genome-size 27m --coverage 160 --input <>.fastq --seed 123 -v > <>_160x.fastq
```

Both `rasusa` and `filtlong` subsamplings were verified for their length distribution using `longQC`.
