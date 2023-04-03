---
layout: default
title: Phasing
nav_order: 7
has_children: true
permalink: /docs/Phasing
---

# Phasing

Nphase was chosen as method for phasing PacBio HiFi DNA sequencing reads given contact and helpful with the developer and its use in multiple research articles. Please have a look at the following links for firther details on the algorithm and use cases.

Available public versions lacked the option of not only extracting phased read clusters, but also the unphased reads theoretically representing those reads closest to the reference genome. Therefore, main scripts contained within the conda environment were modified to add this functionality. Modified scripts are [`nPhase.py`](https://github.com/rortizmerino/haplotyping/blob/main/scripts/nPhase.py), [`nPhaseCleaning.py`](https://github.com/rortizmerino/haplotyping/blob/main/scripts/nPhaseCleaning.py), [`nPhasePipeline_mod.py`](https://github.com/rortizmerino/haplotyping/blob/main/scripts/nPhasePipeline_mod.py) and [`nPhasePipelineFunctions.py`](https://github.com/rortizmerino/haplotyping/blob/main/scripts/nPhasePipelineFunctions.py) provided within this repository. It first has to run in default mode first as follows:

```bash
nphase pipeline --sampleName ${LIB}_${REFNAME} --reference $REFERENCE --output nPhaseOUT --longReads $LRLIBDIR/<>.fastq.gz --longReadPlatform pacbio --R1 $SRLIBDIR/<>_1.fastq.gz --R2 $SRLIBDIR/<>_2.fastq.gz --threads 8

date
echo "nphase finished, cleaning started"
nphase cleaning --sampleName ${LIB}_${REFNAME} --longReads $LRLIBDIR/<>.fastq.gz --resultFolder nPhaseOUT/${LIB}_${REFNAME}/
```

Warning: nPhase can be quite heavy on memory use, therefore, long read libraries were subsetted as detailed within the QC section.

Once the default phasing ran, then phased and unphased read clusters can be obtained as follows:

```bash
python <path>/nPhasePipeline_mod.py filter --sampleName ${LIB}_${REFNAME} --longReads $LRLIBDIR/<>.fastq.gz --clusteredTSV nPhaseOUT/${LIB}_${REFNAME}/Phased/*clusterReadNames.tsv --output nPhaseOUT/${LIB}_${REFNAME}/Phased/FastQ/

# for "cleaned" reads"
python <path>/nPhasePipeline_mod.py filter --sampleName ${LIB}_${REFNAME} --longReads $LRLIBDIR/<>.fastq.gz --clusteredTSV nPhaseOUT/${LIB}_${REFNAME}/Phased/Cleaned/${LIB}_${REFNAME}/${LIB}_${REFNAME}_cleaned.clusterReadNames.tsv --output nPhaseOUT/${LIB}_${REFNAME}/Phased/Cleaned/${LIB}_${REFNAME}/cleanFastQ/
```

Warning: all provided nPhase scripts must be located within the same <path>
