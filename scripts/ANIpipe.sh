#!/usr/bin/bash

SAMPLE=$1
SUFFIX=$2
METHOD=$3
LENGTH=$4

source activate biopython
python scripts/fmtFasta.py -i assembly/${METHOD}/${SAMPLE}_${SUFFIX}/${SAMPLE}.r2cat.${METHOD}.fasta -o assembly/${METHOD}/${SAMPLE}_${SUFFIX}/
# fix table!!!
bash scripts/getAssemblyStats.sh ${SAMPLE} ${SUFFIX} ${METHOD}

source activate fastANI
python scripts/curateAssemblies.py -i assembly/${METHOD}/${SAMPLE}_${SUFFIX}/${SAMPLE}.r2cat.${METHOD}.AssemblyStats.txt -d assembly/${METHOD}/${SAMPLE}_${SUFFIX}

source activate r_env
Rscript --vanilla scripts/ANIhist.R assembly/${METHOD}/${SAMPLE}_${SUFFIX}/${SAMPLE}.r2cat.${METHOD}.ids_n_stats.txt assembly/${METHOD}/${SAMPLE}_${SUFFIX}/${SAMPLE}

source activate biopython
python scripts/fltrFasta.py -i assembly/${METHOD}/${SAMPLE}_${SUFFIX}/${SAMPLE}.r2cat.${METHOD}.fasta -o assembly/${METHOD}/${SAMPLE}_${SUFFIX}/ -l ${LENGTH}
mv assembly/${METHOD}/${SAMPLE}_${SUFFIX}/gt${LENGTH}kb.${SAMPLE}.r2cat.${METHOD}.fasta assembly/${METHOD}/${SAMPLE}_${SUFFIX}/${SAMPLE}.gt${LENGTH}kb.r2cat.${METHOD}.fasta

source activate dotplotter
cd assembly/${METHOD}/${SAMPLE}_${SUFFIX}
cp ${SAMPLE}.r2cat.${METHOD}.ids_n_stats.txt ${SAMPLE}.gt${LENGTH}kb.r2cat.${METHOD}.ids_n_stats.txt
perl ../../../scripts/compareAssemblies_ANI.pl ../../../references/genomes/Scer_noMt.fasta ${SAMPLE}.gt${LENGTH}kb.r2cat.${METHOD}.fasta

source activate r_env
cd ../../..
Rscript --vanilla scripts/Dotplot_n_cov.R assembly/${METHOD}/${SAMPLE}_${SUFFIX}/${SAMPLE}.gt${LENGTH}kb.r2cat.${METHOD}_vs_Scer_noMt.ANI.jpeg genotyping/v4/${SAMPLE}_Scen_Seho_coverage_per1kb_ScSe.jpeg ${SAMPLE}_gt${LENGTH}kb_${METHOD}_ScSe_dotplot_n_cov.jpeg
