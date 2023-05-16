---
layout: default
title: Validation
nav_order: 9
has_children: true
permalink: /docs/validation
---

# Validation

As mentioned in the annotation section, funannotate compare can produce tables of synonymous substitution rates (dS). Main idea is that small number mean close evolutionary relations. The script [print_dS.py](https://github.com/rortizmerino/haplotyping/blob/main/scripts/print_dS.py) takes a funannotate compare output folder to print specific dS values qhich are then taken to [plot_dS.R](https://github.com/rortizmerino/haplotyping/blob/main/scripts/plot_dS.R) so they can visually assesed.

dS values are then useful to validate the annotations and assemblies by tracing the parental origin of genes, while ANI scores do the same for contigs/scaffolds. The script [modGff.pl](https://github.com/rortizmerino/haplotyping/blob/main/scripts/modGff.pl) merges these two indexes into the final annotation files so they can be used for validation and curation.

