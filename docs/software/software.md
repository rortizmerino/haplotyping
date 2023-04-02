---
layout: default
title: Software
nav_order: 2
has_children: true
has_toc: true
permalink: /docs/software
---

# Software

Here you will find general info about Software 

Most software used in this project was installed and run trough conda environments. Please check the conda documentation in their source page. Below are listed environment names and tools installed within them. Tools and environments are then used, mentioned and linked within the other sections in this repository.

## Conda

The specific conda version used here can be obtained with the commands below

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh
bash Miniconda3-py38_4.12.0-Linux-x86_64.sh
conda init bash
```
Some useful commands are 

```bash
conda info --envs
conda list
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

## sra-toolkit

This toolkit can be used to download public datasets from the National Centre for Biotechnology Information (NCBI), specifically, from the Sequence Read Archive (SRA)

```bash
conda create --name sratk
source activate sratk
conda install -c bioconda sra-tools
```

## longqc

LongQC is a tool for the data quality control of the PacBio and ONT long reads. It was mainly used to make diagnostic plots altough it has many more functionalities. Check out https://github.com/yfukasawa/LongQC most instructions below were taken from there, except for a few workarounds

```bash
conda create --name LongQC
source activate LongQC
	
conda install python=3.9
conda install numpy
conda install scipy
conda install matplotlib
conda install scikit-learn
conda install pandas
conda install jinja2
conda install h5py
	
conda install -c bioconda pysam
conda install -c bioconda edlib
conda install python-edlib -y

# in case they are missing
# sudo apt install make
# sudo apt-get install libz-dev
	
git clone https://github.com/yfukasawa/LongQC.git
cd LongQC/minimap2-coverage && make

# if the following error is found 
# ImportError: libcrypto.so.1.1: cannot open shared object file: No such file or directory
# workaround https://github.com/merenlab/anvio/issues/1479

# some quick and dirty tricks, sorry if you cannot do the same but this library can be annoying
cp /data/raul/software/miniconda3/envs/vc/lib/libcrypto.so.1.1 /data/raul/software/miniconda3/envs/LongQC/lib/.
cd /data/raul/software/miniconda3/envs/LongQC/lib/
rm libcrypto.so
ln -s libcrypto.so.1.1 libcrypto.so
ln -s libcrypto.so.1.1 libcrypto.so.1.0.0
```

## filtlong

Filtlong is a tool for filtering long reads by quality and length. It can take a set of long reads and produce a smaller, arguably better subset. Check out https://github.com/rrwick/Filtlong

```bash
git clone https://github.com/rrwick/Filtlong.git
cd Filtlong
make -j
#export PATH=$PATH:/scratch/tmp/rauort/Filtlong/bin
```

## rasusa

This tool takes its name from random sub-sampling. As opossed to filtlong, it better preserves the read length ditribution as it does not only take the X-longest reads but takes a random sample which satisfies the condition of a certain amount of bases (i.e. target coverage)

```bash
conda create -n rasusa
conda activate rasusa
conda install rasusa
```

## vc

These tools were used for genotyping within a couple scripts described in the corresponding section. They involve read mapping, variant calling, and coverage calculation.

```bash
conda create -n vc
conda activate vc
conda install -c bioconda java-jdk=8
conda install -c bioconda bwa
conda install -c bioconda samtools

wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip
unzip gatk-4.3.0.0.zip
```

## R

R here is mainly used for plotting. Installing it within conda has the advantage of avoiding the need of setting up alternative/custom package locations. Please checkout https://docs.anaconda.com/anaconda/user-guide/tasks/using-r-language/

```bash
conda create --name r_env
conda activate r_env
	
conda install r-base
# contains some of the so-called tidyverse tools
conda install r-essentials
	
# to be executed within R (literally just type R once the environment is loaded)
install.packages("BiocManager")
library(BiocManager)
# used for Copy Number Variation (CNV) plots 
BiocManager::install("DNAcopy")
#library(DNAcopy)

# used to make compound plots
install.packages("jpeg")
#library(jpeg)	
```

## dotplotter

This environment just needs base R and mummer to generate dotplots. Here they are maintained separately from the other R environment and used within a Perl wrapper since it seems they do not communicate well with other tools/languages.

```bash
# quick and dirty mode
conda create --name dotplotter --clone r_env

# otherwise
conda create --name dotplotter
conda activate dotplotter
conda install r-base
	
conda install -c bioconda perl-bioperl-core
conda install -c bioconda mummer

# bad communication example
# conda install -c conda-forge biopython 
# seems to mess with something in R 
# conda error while loading shared libraries: libreadline.so.6
```

## fastani

Fast implementation of an Average Nucleotide Identity calculation used by NCBI

```bash
conda create --name fastANI
conda activate fastANI
conda install -c bioconda fastani

conda install python=3.8
conda install -c conda-forge biopython
```

## circos

Original Perl-based implementation of circos, please check http://circos.ca/ for documentation. Includes mummer as it is used within wrapper perl scripts to get alignments and coordinates to display in the plots. 

```bash
conda create --name circos
conda activate circos
	
conda install -c bioconda circos
conda install -c bioconda perl-bioperl
conda install -c bioconda perl-bioperl-core
conda install -c bioconda mummer
```

## biopython

Every now and then a clean biopython install makes life easier.

```bash
conda create --name biopython
conda activate biopython
	
conda install python=3.8
conda install -c conda-forge biopython
```

## genomescope2

Estimate genome heterozygosity, repeat content, and size from sequencing reads using a kmer-based statistical approach. Check https://www-nature-com.tudelft.idm.oclc.org/articles/s41467-020-14998-3 and https://github.com/tbenavi1/genomescope2.0

```bash
conda create --name genomescope2
conda activate genomescope2
	
conda install genomescope2
```

## smudgeplot

This one is meant to be part of genomescope2 but it is actually not (yet?). Can be used with either jellyfish or KMC3 to do the kmer calculation.

```bash
conda create --name smudgeplot
conda activate smudgeplot
	
conda install -c bioconda smudgeplot
conda install -c bioconda jellyfish
	
mkdir KMC3
cd KMC3
wget https://github.com/refresh-bio/KMC/releases/download/v3.2.1/KMC3.2.1.linux.tar.gz
tar zxvf KMC3.2.1.linux.tar.gz
./KMC3/bin/kmc
```

## nPhase

Please checkout https://github.com/OmarOakheart/nPhase for detailed info. This particular version does not use an insane amount of memory as newer versions do. It was also modified to print out unphased reads, this is a personal modification so please check the phasing section for more details.

```bash
conda create -n nPhaseMod python=3.8
conda activate nPhaseMod
conda install -c oakheart nphase=1.1.10
# fix plotnine issue
conda install matplotlib=3.5.3
```

## iprscan

Interproscan can be used to obtain a plethora of functional terms from biological sequences. Here is mainly used to complement funannotate, please check https://interproscan-docs.readthedocs.io/en/latest/index.html for more details

```bash
mkdir /scratch/rortizmerino/software/interproscan
cd /scratch/rortizmerino/software/interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.60-92.0/interproscan-5.60-92.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.60-92.0/interproscan-5.60-92.0-64-bit.tar.gz.md5
md5sum -c interproscan-5.60-92.0-64-bit.tar.gz.md5
tar -pxvzf interproscan-5.60-92.0-*-bit.tar.gz

# to setup databases within a compute cluster running SLURM (such as delftblue)
srun --job-name="ipr_setup" --partition=compute --time=00:30:00 --ntasks=2 --cpus-per-task=1 --mem-per-cpu=10GB --pty bash
python3 setup.py interproscan.properties

# delftblue specific, otherwise install openjdk v11 somehow else
module load openjdk/11.0.12_7
```

## funannotate

Warning: this is one of the most complicated installations ever, the steps below worked here but there is absolutelly no warranty they will work elsewhere, please check https://funannotate.readthedocs.io/en/latest/ and good luck!

```bash
# define a large enough folder to hold the packages
conda config --add pkgs_dirs /scratch/rortizmerino/conda/pkgs
conda config --show
	
conda create --name funannotate --prefix /scratch/rortizmerino/conda/funannotate "python>=3.6,<3.9"
conda activate /scratch/rortizmerino/conda/funannotate
conda install funannotate
	
# reveal missing dependencies, might have to run a few times troughout installation
funannotate check --show-versions

# missing local::lib not installed, install with cpanm local::lib
conda install -c bioconda perl-local-lib

funannotate setup -d /scratch/rortizmerino/software/funannotate_db
funannotate setup --busco_db fungi -d /scratch/rortizmerino/software/funannotate_db
	
export FUNANNOTATE_DB=/scratch/rortizmerino/software/funannotate_db
	
# get GeneMark-ES/ET/EP ver 4.71_lic , LINUX 64 kernel 3.10 - 5 !!!!!!!
http://topaz.gatech.edu/GeneMark/license_download.cgi
	
cd /scratch/rortizmerino/software/
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_BySRT/gmes_linux_64_4.tar.gz
gunzip gm_key_64.gz
cp gm_key_64 ~/.gm_key
tar zxvf gmes_linux_64_4
cd gmes_linux_64_4/
bash check_install.bash
perl change_path_in_perl_scripts.pl /scratch/rortizmerino/conda/funannotate/bin/perl 
	
# add before running job !!!!!!!!!!!!!!!!!!!!
export GENEMARK_PATH=/scratch/rortizmerino/software/gmes_linux_64_4
export PATH="/scratch/rortizmerino/software/gmes_linux_64_4:$PATH"
export PATH="/scratch/rortizmerino/software/gmes_linux_64_4/gmes_petap.pl:$PATH"

# install augustus (<v3.3.3) separately!!!!
wget https://github.com/Gaius-Augustus/Augustus/releases/download/v3.3.3/augustus-3.3.3.tar.gz
tar zxf augustus-3.3.3.tar.gz
cd Augustus
make augustus
	
export PATH=/scratch/rortizmerino/software/augustus-3.3.3/bin:/scratch/rortizmerino/software/augustus-3.3.3/scripts:$PATH
export AUGUSTUS_CONFIG_PATH=/scratch/rortizmerino/software/augustus-3.3.3/config/
# remove augustus within conda 
conda remove augustus -n funannotate --force
	
funannotate setup -b saccharomycetales

ete3 upgrade-external-tools paml
# target directory? [/home/rortiz/.etetoolkit/]: y
```

## mitohifi_env

Assemble and annotate mitochondrial chromosome(s) from PacBio HiFi data. Checkout https://github.com/marcelauliano/MitoHiFi

```bash
git clone https://github.com/RemiAllio/MitoFinder.git
cd MitoFinder
./install.sh

export PATH="/data/raul/software/MitoFinder:$PATH"

cd ..
git clone https://github.com/marcelauliano/MitoHiFi.git
cd MitoHiFi

conda env create -n mitohifi_env -f mitohifi_env.yml 
source activate mitohifi_env

export PATH="/data/raul/software/MitoFinder:$PATH"
export PATH="/tudelft.net/staff-bulk/tnw/BT/imb/sequence_data/bin/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin:$PATH"
export PATH="/data/raul/software/MitoHiFi:$PATH"
```
