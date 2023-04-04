---
layout: default
title: Assembly
nav_order: 6
has_children: true
permalink: /docs/Assembly
---

# Assembly

Here you will find general info about genome assembly methods implemented with long DNA sequencing reads (either PacBio HiFi or Oxford Nanopore Technologies). There are multiple different genome assemblers, `canu` and `flye` were used here with slightly different purposes. In general terms, `canu` was found useful for full read sets (up to 800x coverage) whereas `flye` gave better results for smaller datasets such as those obtained from phased read clusters (around 40x coverage). While setting up the assembly runs is rather straightforward, some post-processing steps are shown here to facilitate their analysis.

Canu has a heavy memory use, therefore assembly runs with full read sets (640-800x coverage) were launched in FAT nodes with a minimum of 240 Gb RAM. A sample command looks like:

```bash
canu -pacbio-hifi ${LIBDIR}/${LIB}_hifi.fastq genomeSize=28m -d $RUNDIR -p ${LIB} useGrid=false maxThreads=16
```
Warning: `genomeSize` here was estimated using genomescope as described in the QC section. `maxThreads` is generally recommended to be kept between 8-16. If you consider yourself as an hpc pro, set `useGrid=true` and configure it to use your favourite scheduler pretty sure admins will love you.

Flye is also quite straightforward to run, however, has very different performance with small and large datasets. Considering the example above, 640x coverage of a 28 Mb genome required 720 Gb RAM. A sample command looks like:

```bash
flye --pacbio-hifi ${LIBDIR}/${LIB}_hifi_640x.fastq --genome-size 28m --out-dir $RUNDIR --threads 16 --keep-haplotypes
```
Warning: See warnings on `--genome-size` and `--threads` above. `--keep-haplotypes` is specific to flye but its effect has not been tested thoroughly, is more for peace of mind in this case.

Flye was also used iteratively to separately assemble unphased, and phased read clusters obtained by nPhase as described in the Phasing section. A sample bash snippet is shown below:

```bash 
# set variables
LIB=<>
RUNDIR=${LIB}_002
LIBDIR=nPhaseOUT/<>_Scer/Phased/FastQ

# set run directory and go there
if [ ! -d "$RUNDIR" ]; then
 mkdir -p $RUNDIR
fi
cd $RUNDIR

# START
date
# gscope kmer genome size ~27.9 
flye --pacbio-hifi ${LIBDIR}/${LIB}_Scer_unphased.fasta --genome-size 28m --out-dir $RUNDIR/unphased --threads 16 --keep-haplotypes
date

FILES=($(ls ${LIBDIR}/*fasta))
for FILE in "${FILES[@]}"; do
 FILENAME=$(echo $FILE | awk -F / '{ printf "%s",$NF }')
 OUT=$(echo $FILENAME | awk -F _ '{ printf "%s_%s_%s_%s",$3,$4,$5,$6 }')

 if [ ! -d "$RUNDIR/$OUT" ]; then
  mkdir -p $RUNDIR/$OUT
 fi

 TMP=""
 if [[ $FILENAME =~ merged ]]; then
  TMP=$(echo $FILENAME | awk -F _ '{ printf "%s_%s",$6,$7 }')
 elif [[ $FILENAME =~ cluster ]]; then
  TMP=$(echo $FILENAME | awk -F _ '{ printf "%s",$6 }')
 fi
 TMP2=$(echo $TMP | sed s/".fasta"/".fastq.gz"/ )
 FQ=($(ls ${LIBDIR}/*${TMP2}))

 # phased
 # Sizes according to Scen (CEN.PK ONT) reference genome
 if [[ $OUT =~ Scer_01 ]]; then 
  echo "Fasta = $FILE"
  echo "Fastq = $FQ"
  echo -e "Assembling $FQ into $RUNDIR/$OUT\n"
  flye --pacbio-hifi $FQ --genome-size 0.21m --out-dir $RUNDIR/$OUT --threads 16
 fi
 
 ...
 
 if [[ $OUT =~ Scer_16 ]]; then 
  echo "Fasta = $FILE"
  echo "Fastq = $FQ"
  echo -e "Assembling $FQ into $RUNDIR/$OUT\n"
  flye --pacbio-hifi $FQ --genome-size 0.96m --out-dir $RUNDIR/$OUT --threads 16
 fi
done
 
```
Warning: Please note this snippet was made to work on an nPhase output folder structure already processed to print unphased reads as described in the Phasing section. It also has hardcoded references and parameters specific to *S. cerevisiae* (Scer) as its genome was used as reference. Phased read clusters are obtained given each reference chromosome, only chromosomes 1 and 16 are shown here and their corresponding sizes are used to set `--genome-size` parameters accordingly.

Next step is ordering the scaffolds/contigs according to their synteny and parental origin. As strains used here were already known to be *S. cerevisiae* x *S. eubayanus* hybrids, a concatenated reference containing both genomes was used. A quick and dirty option for this syntenic ordering is to use r2cat. Pull requests are welcome in case a command-line option is available. Until then please visit this link <https://bibiserv.cebitec.uni-bielefeld.de/cgcat?id=cgcat_r2cat> to download and run the stand alone tool. Next steps work regardless of this syntenic ordering but actually doing it makes life easier.


