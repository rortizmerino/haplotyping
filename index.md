---
title: Home
layout: home
---

Genome assembly and annotation of hybrid yeast strains using high-fidelity long sequencing reads
======
<br>
  Nowadays is commonly known that genome sequences are powerful tools to study and engineer living organisms. For instance, the yeast Saccharomyces cerevisiae, also known as brewer’s yeast, was the first eukaryotic organism whose genome was be fully sequenced and functionally annotated enabling a multitude of functional genomics efforts. As valuable as this knowledge has been, such reference genome mainly represents a laboratory strain artificially kept in a stable and simple configuration which is now known far from industrial strains with complex and plastic genomes. Here we used Pacific Biosciences high-fidelity (HiFi) reads, which are the most advanced genome sequencing technology available to date, to produce fully-resolved genome assemblies for 2 industrial yeast strains whose genomes are already known to be more complex than common laboratory strains. 
<br>
  The industrial S. pastorianus yeast strains used in this project are known to be interspecific hybrids, that is, they include genetic material from 2 different species, namely S. cerevisiae  and S. eubayanus. Furthermore, these strains present such genetic material in an unbalanced distribution and are referred as aneuploids, as opposed to the reference brewer’s yeast genome which is euploid and haploid by having 16 chromosomes in a single copy. Genomic sequences with either S. cerevisiae or S. eubayanus origin can be distinguished from each other given their distinct Average Nucleotide Identity (ANI). However, the different S. pastorianus strains have undergone different processes of genome evolution from which genomic regions from one parent, or the other, have been preferentially kept and independently modified by mutation and recombination events. These regions with the same parental origin cannot be as easily distinguished by their ANI alone. We therefore tried a “phasing” method in which HiFi long reads are classified by means of those slight differences and are then clustered into their different parental origins or haplotypes and analysed separately.
<br>
  As a result we have two alternative assemblies available per S. pastorianus strain, considering phasing or not. Both strategies have a high level of agreement with each other but differ in genomic regions representing copy number variation, and recombination events. These alternative assemblies have already been assigned to their corresponding parental origin and systematically named to facilitate their analysis. Furthermore, these assemblies have been annotated, meaning they include their encoded genes which were found by different comparative genomics strategies and include different representations of their cellular function (S. cerevisiae closest relative, functional domains, gene ontologies). These annotated assemblies have already been used as reference to analyse different evolved cell lines providing evidence of the genomic regions and haplotypes responsible for the phenotypes shown in parental and evolved industrial S. pastorianus strains.


----

[^1]: [It can take up to 10 minutes for changes to your site to publish after you push the changes to GitHub](https://docs.github.com/en/pages/setting-up-a-github-pages-site-with-jekyll/creating-a-github-pages-site-with-jekyll#creating-your-site).

[Just the Docs]: https://just-the-docs.github.io/just-the-docs/
[GitHub Pages]: https://docs.github.com/en/pages
[README]: https://github.com/just-the-docs/just-the-docs-template/blob/main/README.md
[Jekyll]: https://jekyllrb.com
[GitHub Pages / Actions workflow]: https://github.blog/changelog/2022-07-27-github-pages-custom-github-actions-workflows-beta/
[use this template]: https://github.com/just-the-docs/just-the-docs-template/generate
