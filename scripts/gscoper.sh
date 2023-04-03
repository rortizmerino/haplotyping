
Help()
{
 # Display help
 echo "Usage: gscoper.sh -s <strain> -r <rundir> -p <path>"
 echo
 echo " s: strain ID"
 echo " r: directory in which to run the analysis"
 echo " p: path to input data (fastq files)"
 echo " h: print this help and exit"
 echo
}

while getopts s:r:p:h: option
do
 case "${option}" in
  s) STRAIN=${OPTARG};;
  r) RUNDIR=${OPTARG};;
  p) PATH2DATA=${OPTARG};;
  h) Help
     exit;;
  *) echo "ERROR: invalid option"
     Help
     exit;;
 esac
done

if [ -z $1 ]
then 
 echo "ERROR: arguments needed"
 Help
 exit
fi

######
date
echo "running gscoper.sh -s ${STRAIN} -r ${RUNDIR} -p ${PATH2DATA}"

KMCpath="/home/nfs/rortizmerino/seqdata/raul/haplotyping/software/KMC3/bin"

if [ ! -d "${RUNDIR}" ]; then
 mkdir -p ${RUNDIR}
fi
cd ${RUNDIR}
pwd

# create a directory for temporary files
mkdir tmp logs

# create a file with both the raw read files
ls ${PATH2DATA}/*.fastq > FILES

# run kmc
${KMCpath}/kmc -k21 -t16 -m64 -ci1 -cs10000 @FILES kmer_counts tmp
${KMCpath}/kmc_tools transform kmer_counts histogram kmer_k21.hist -cx10000
${KMCpath}/kmc_dump -ci15 -cx3000 kmer_counts kmer_k21.dump

source activate genomescope2
genomescope2 -i kmer_k21.hist -k 21 -p 2 -o . -n ${STRAIN}_genomescope
conda deactivate

source activate smudgeplot
L=$(smudgeplot.py cutoff kmer_k21.hist L)
U=$(smudgeplot.py cutoff kmer_k21.hist U)
echo
echo "-> L:$L U:$U"
echo
smudgeplot.py hetkmers -o kmer_pairs < kmer_k21.dump
smudgeplot.py plot -o ${STRAIN} -t "${STRAIN}" kmer_pairs_coverages.tsv
conda deactivate

# for tetraploids looking like diploids only
source activate genomescope2
genomescope2 -i kmer_k21.hist -k 21 -p 4 -o . -l 10 -n ${STRAIN}_genomescope_p4
conda deactivate

source activate smudgeplot
smudgeplot.py plot -o ${STRAIN}_L10_n17 -n 17 -t "${STRAIN}" kmer_pairs_coverages.tsv
