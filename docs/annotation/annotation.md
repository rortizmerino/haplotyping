---
layout: default
title: Annotation
nav_order: 8
has_children: true
permalink: /docs/annotation
---

# Genome Annotation

Perhaps the best way to describe how the genome annotation was performed is to go trough a bash script wrapping the different funannotate functions used. WARNING funannotate is very complete but rather difficult to install so please have a look at the software section on it, does not hurt either to lookl at the funannotate documentation.

[FArun_07.sh](https://github.com/rortizmerino/haplotyping/blob/main/scripts/FArun_07.sh) has a real example, main idea is to load all software, run the cleaning step so assembled sequences are flagged for potential removal, and perform a soft masking as preparation. Those steps are somewhat straightforward so please just have a look at the main script.

Prediction:
This step uses the funannotate predict module and the main parameters used set the expected poidy to 3 (triploid), make an initial de novo run with genemark, then use augustus with S. cerevisiae reference proteins, and end up with a further refinement with BUSCO sets specific to Saccharomyces and Saccharomycetales.

```
funannotate predict -i ${sample}.masked.fa -o ${sample}${tag} --species ${sp_mod} --strain ${strain} --cpus $cpus --name ${tag} --cpus $cpus --ploidy 3 --protein_evidence Scer.uniprot.fasta --genemark_mode ES --augustus_species saccharomyces --busco_seed_species saccharomyces --busco_db saccharomycetales
```

Functional prediction:
Although this step can technichally be performed as a funannotate module, it is waaaay easier to do it separately. For this purpose interproscan has to be installed separately. This example only looks for PFAM and GO terms

```
interproscan.sh --applications Pfam -cpu 6 -f xml -goterms -i ${sample}${tag}/predict_results/${sp_mod}_${strain}.proteins.fa -o ${sample}${tag}.iprscan.xml
```

Functional annotation:
PFAM and GO terms have to then be merged with the gene prediction results as follows

```
funannotate annotate -i ${sample}${tag}/ --iprscan ${sample}${tag}.iprscan.xml --species ${sp_mod} --strain ${strain} --cpus $cpus --busco_db saccharomycetales
```

Afterwards, the comparison module can be used for 2 or more annotations. The following command shows how to use 2 parental annotations with a new one. WARNING: run_dnds is used for further validation steps and requrire extra attenntion during the installation so go have a look at the software section once more.

```
funannotate compare -i Scer_noMt.fmtd.gbk Saccharomyces_eubayanus_Seho.gbk Saccharomyces_pastorianus.gbk --out Scer_Seho_testfull --run_dnds estimate --cpus 6
```

