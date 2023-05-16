#!/bin/sh

#SBATCH --job-name="FArun_07"
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

#SBATCH --partition=compute		# compute | gpu | memory
#SBATCH --time=24:00:00			# max 120h
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G		# 6 * 20 = 120G

#SBATCH --mail-type=ALL

# load modules, environments, binaries, ...
module load 2022r2 openjdk/11.0.12_7

module load miniconda3/4.12.0
conda activate /scratch/rortizmerino/conda/funannotate

export PATH="/scratch/rortizmerino/software/augustus-3.3.3/bin:/scratch/rortizmerino/software/augustus-3.3.3/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH=/scratch/rortizmerino/software/augustus-3.3.3/config/
export FUNANNOTATE_DB=/scratch/rortizmerino/software/funannotate_db
export GENEMARK_PATH=/scratch/rortizmerino/software/gmes_linux_64_4
export PATH="/scratch/rortizmerino/software/gmes_linux_64_4:$PATH"
export PATH="/scratch/rortizmerino/software/gmes_linux_64_4/gmes_petap.pl:$PATH"

RUNDIR=/scratch/rortizmerino/annotation/FA_07

# set run directory and go there
if [ ! -d "$RUNDIR" ]; then
 mkdir -p $RUNDIR
fi
cd $RUNDIR

cp /scratch/rortizmerino/genomes/LM07.r2cat.canu.fasta LM07.fa
cp /scratch/rortizmerino/genomes/Scer.uniprot.fasta .

# START
date
 sample="LM07"
species="Saccharomyces pastorianus"
 sp_mod="${species// /_}"
 strain="LM07"
    tag="LM07full"
  fasta="LM07"
   cpus=6

date
echo "cleaning started"
cmd="funannotate clean -i ${fasta}.fa -o ${sample}.clean.fa"
echo "CMD: " $cmd
$cmd
echo "cleaning finished"

date
echo "masking started"
cmd="funannotate mask -i ${fasta}.fa -o ${sample}.masked.fa --cpus $cpus"
echo "CMD: " $cmd
$cmd
echo "masking finished"

awk '{if ($0 ~ /^>/) {print $1} else {print} }' ${sample}.masked.fa > tmp.fa
#awk '{if ($0 ~ /^>/) {gsub("LM07_002_Scer_",""); gsub("mergedCluster_contig_",""); gsub("unphased_contig_",""); print $1} else {print} }' /scratch/rortizmerino/annotation/FA_05/LM07.masked.fa > tmp.fa
mv tmp.fa ${sample}.masked.fa

date
echo "prediction started"
cmd="funannotate predict -i ${sample}.masked.fa -o ${sample}${tag} --species ${sp_mod} --strain ${strain} --cpus $cpus --name ${tag} --cpus $cpus --ploidy 3 --protein_evidence Scer.uniprot.fasta --genemark_mode ES --augustus_species saccharomyces --busco_seed_species saccharomyces --busco_db saccharomycetales"
echo "CMD: " $cmd
$cmd
echo "prediction finished"

date
echo "interpro started"
cmd="/scratch/rortizmerino/software/interproscan/interproscan-5.60-92.0/interproscan.sh --applications Pfam -cpu 6 -f xml -goterms -i ${sample}${tag}/predict_results/${sp_mod}_${strain}.proteins.fa -o ${sample}${tag}.iprscan.xml"
echo "CMD: " $cmd
$cmd
echo "interpro finished"

date
echo "annotation started"
cmd="funannotate annotate -i ${sample}${tag}/ --iprscan ${sample}${tag}.iprscan.xml --species ${sp_mod} --strain ${strain} --cpus $cpus --busco_db saccharomycetales"
echo "CMD: " $cmd
$cmd
echo "annotation finished"

date
echo "bye"
# END

date
echo "doei!"
