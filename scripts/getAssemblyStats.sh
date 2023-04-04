#!/usr/bin/bash

SAMPLE=$1
SUFFIX=$2
METHOD=$3
INFILE="assembly/${METHOD}/${SAMPLE}_${SUFFIX}/${SAMPLE}.r2cat.${METHOD}.fasta"
OUTFILE="assembly/${METHOD}/${SAMPLE}_${SUFFIX}/${SAMPLE}.r2cat.${METHOD}.AssemblyStats.txt"

awk -v sample="$SAMPLE" 'BEGIN{FS=">| |=";print"#sample\tseq_name\tlength\tcov.\tcirc.\trepeat\tmult\talt_group\tgraph_path"}{if ($0 ~ />/ && $3 !~ /reverse/) printf"%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\n",sample,$2,$4,$8,$16,$14; if ($0 ~ />/ && $3 ~ /reverse/) printf"%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\n",sample,$2,$7,$11,$19,$17}' $INFILE > $OUTFILE
